import os
import pathlib
import tempfile
from abc import ABC, abstractmethod
from typing import Callable

import pandas as pd
from loguru import logger
from openmm.app import PDBxFile
from pdbfixer import PDBFixer
from tqdm.auto import tqdm

from loopbuilder.convert import extract_segment_from_mmcif, join_segments
from loopbuilder.segment import Segment, SegmentModel
from loopbuilder.score import Scorer, Filter
from loopbuilder.typing import StrPath


class Builder(ABC):
    """Base class for loop builders

    For concrete builders, it should normaly be sufficient to implement the
    `build_segment` method. This should contain the logic for constructing
    a single `SegmentModel` for a given `Segment`.

    A default implementation of the higher order `build` method is provided,
    This takes care of organizing the building process for multiple segments
    and models (making calls to `build_segment`), also including filtering and scoring.

    How missing segments are found in the first place is defined in `find_segments`.
    Full control over these segments can alternatively be achieved through manipulation
    of the `segments` attribute, which is a list of `Segment` objects. The `build` method
    will call `find_segments` only if no segments are present.
    """

    def __init__(
        self,
        structure_file: StrPath,
        output_directory: StrPath,
        working_directory: StrPath | None = None,
        scorers: list[Scorer | Callable[[SegmentModel], None]] | None = None,
        filters: list[Filter | Callable[[SegmentModel], bool]] | None = None,
        progress_bar: bool = True,
    ):
        """Initialize the builder

        Args:
            structure_file: Path to the input structure file
            output_directory: Path to the output directory. Successful models will be saved here.

        Keyword args:
            working_directory: Path to a working directory. All trial models and intermediate files
                will be saved here. If `None`, uses a temporary working directory.
            scorers: List of scorers to use on generated segment models
            filters: List of filters to use on generated segment models
            progress_bar: Show a progress bar for the model building process
        """
        if working_directory is not None:
            working_directory = pathlib.Path(working_directory).resolve()
        self.working_directory = working_directory
        self.structure_file = pathlib.Path(structure_file).resolve()
        self.output_directory = pathlib.Path(output_directory).resolve()
        self.scorers = scorers or []
        self.filters = filters or []
        self.segments: list[Segment] = []
        self.progress_bar = progress_bar

    def find_segments(self) -> None:
        """Find missing segments in the input structure

        NOTE: Uses PDBFixer at the moment but this part could also be
            implemented for example in Biopython if needed.

        NOTE: Ignores missing terminal residues
        """

        self.segments = []
        fixer = PDBFixer(filename=str(self.structure_file))
        fixer.findMissingResidues()

        chains = list(fixer.topology.chains())
        non_terminal = {}
        for key, value in fixer.missingResidues.items():
            chain_index, missing_index = key
            if missing_index == 0:
                continue

            if missing_index == len(chains[chain_index]):
                continue

            non_terminal[key] = value

        for i, ((chain_index, residue_start_index), residue_names) in enumerate(non_terminal.items(), 1):
            chain_residues = list(chains[chain_index].residues())
            segment = Segment(
                identifier=f"loop_{i}",
                chain_index=chain_index,
                chain_name=chains[chain_index].id,
                residue_start_index=residue_start_index,
                residue_start_seqid=int(chain_residues[residue_start_index - 1].id) + 1,
                residue_index_offset=int(chain_residues[0].id),
                residue_names=residue_names,
                parent_structure_file=self.structure_file,
                models=[],
            )
            self.segments.append(segment)

    def build(self, n: int = 1, *, max_tries: int | None = None) -> list[Segment] | None:
        """Build n models for the input structure

        Keyword args:
            n: Number of models (passing the `filters`) to build
            max_tries: Maximum number of models to build before aborting
                even if `n` has not been reached. If `None`, defaults to `n * 10`.
            scorers: List of scorers to use
            filters: List of filters to use

        Returns:
            List of `Segment` objects with the models built
        """

        if self.working_directory is not None:
            working_directory = self.working_directory
            working_directory.mkdir(parents=True, exist_ok=True)
        else:
            temp_dir = tempfile.TemporaryDirectory()
            working_directory = pathlib.Path(temp_dir.name)

        # NOTE: Make sure the working directory is writable for other processes
        os.chmod(working_directory, 0o777)

        if n < 1:
            logger.warning("Building 0 models. Exiting.")
            return

        if max_tries is None:
            max_tries = n * 10

        if max_tries < 1:
            logger.warning("Using 0 trial models. Exiting.")
            return

        logger.info(f"Building n={n} models for {self.structure_file} (max_tries={max_tries})")
        self.output_directory.mkdir(parents=True, exist_ok=True)
        logger.info(f"Saving output to {self.output_directory}")
        logger.info(f"Using working directory {working_directory}")

        scorer_str = "\n".join([f"  {scorer!r}" for scorer in self.scorers])
        if scorer_str:
            scorer_str = ":\n" + scorer_str
        logger.info(f"Using {len(self.scorers)} scorer(s){scorer_str}")

        filter_str = "\n".join([f"  {filter!r}" for filter in self.filters])
        if filter_str:
            filter_str = ":\n" + filter_str
        logger.info(f"Using {len(self.filters)} filter(s){filter_str}")

        if not self.segments:
            logger.info("Looking for segments")
            self.find_segments()

        if not self.segments:
            logger.warning("No segments found")
            return self.segments

        logger.info(f"Found {len(self.segments)} segments")

        if self.progress_bar:
            it_segments = tqdm(self.segments, desc="Building segments", unit="segment")
        else:
            it_segments = self.segments

        for segment in it_segments:
            logger.info(f"Building models for segment {segment.identifier}")
            n_success = 0
            n_tries = 0

            while True:
                if n_success >= n:
                    logger.success(f"Reached the target number of models (success_rate={n_success / n_tries:.2%})")
                    break

                if n_tries >= max_tries:
                    logger.warning(f"Reached the maximum number of tries (success_rate={n_success / n_tries:.2%})")
                    break

                n_tries += 1
                segment_model = self.build_segment(segment, trial_id=str(n_tries), working_directory=working_directory)

                segment_start = segment.residue_start_seqid
                segment_end = segment_start + len(segment) - 1
                model_structure_file = segment_model.structure_file.with_stem(f"{segment_model.structure_file.stem}_segment")
                extract_segment_from_mmcif(
                    segment_model.structure_file,
                    model_structure_file,
                    residue_indices={segment_start, segment_end},
                    chain_id=segment.chain_name,
                )
                segment_model.structure_file = model_structure_file

                for scorer in self.scorers:
                    scorer(segment_model)

                logger.info(f"Scored trial model {n_tries} for segment {segment.identifier}: {segment_model.scores}")

                for filter_ in self.filters:
                    if not filter_(segment_model):
                        logger.info(f"Trial model {n_tries} for segment {segment.identifier} failed filter {filter_}")
                        break
                else:
                    n_success += 1
                    model_structure_file = (
                        self.output_directory / f"{segment.parent_structure_file.stem}_{segment.identifier}_{n_success}.cif"
                    )
                    segment_model.structure_file.rename(model_structure_file)
                    segment_model.structure_file = model_structure_file
                    segment_model.index = n_success
                    segment.models.append(segment_model)
                    logger.success(f"Built model {n_success} for segment {segment.identifier}")

            if segment.models:
                model_structure_files = [m.structure_file for m in segment.models]
                joined_model_structure_file = self.output_directory / f"{segment.parent_structure_file.stem}_{segment.identifier}.cif"
                join_segments(
                    model_structure_files,
                    joined_model_structure_file,
                )
                for model in segment.models:
                    model.structure_file = joined_model_structure_file
                for file in model_structure_files:
                    file.unlink(missing_ok=True)


        segment_df = pd.DataFrame(
            self.segments,
            columns=[
                "identifier",
                "chain_index",
                "chain_name",
                "residue_start_index",
                "residue_start_seqid",
                "residue_index_offset",
                "residue_names",
                "parent_structure_file",
            ],
        )
        segment_df.to_csv(self.output_directory / "segments.csv", index=False)

        models = [m for segment in self.segments for m in segment.models]
        model_df = pd.DataFrame(
            models,
        )
        model_df.to_csv(self.output_directory / "models.csv", index=False)

        return self.segments

    @abstractmethod
    def build_segment(self, segment: Segment, trial_id: str, working_directory: pathlib.Path) -> SegmentModel:
        """Build a model for a segment

        Args:
            segment: `Segment` object for which to build a model
            trial_id: Trial ID for the model, used to generate a unique (temporary) output filename
            working_directory: Path to a working directory. All trial models and intermediate files
                will be saved here.

        Returns:
            A `SegmentModel` object
        """


class PDBFixerBuilder(Builder):
    def build_segment(self, segment: Segment, trial_id: str, working_directory: pathlib.Path) -> SegmentModel:
        """Build a segment model using PDBFixer

        Args:
            segment: `Segment` object for which to build a model
            trial_id: Trial ID for the model, used to generate a unique (temporary) output filename
            working_directory: Path to a working directory. All trial models and intermediate files
                will be saved here.

        Returns:
            A `SegmentModel` object
        """

        fixer = PDBFixer(str(segment.parent_structure_file))
        fixer.missingResidues = {(segment.chain_index, segment.residue_start_index): segment.residue_names}

        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        model_structure_file = (
            working_directory / f"{segment.parent_structure_file.stem}_{segment.identifier}_{trial_id}.cif"
        )
        with open(model_structure_file, "w") as fp:
            # NOTE: Write out full structure here for scoring. Segments will be split off later.
            PDBxFile.writeFile(fixer.topology, fixer.positions, file=fp, keepIds=True)

        return SegmentModel(identifier=segment.identifier, structure_file=model_structure_file, scores={})
