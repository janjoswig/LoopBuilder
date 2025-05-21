import pathlib
import tempfile
from abc import ABC, abstractmethod
from loguru import logger
from openmm.app import PDBxFile
from pdbfixer import PDBFixer
from tqdm.auto import tqdm

from loopbuilder.segment import Segment, SegmentModel
from loopbuilder.score import Scorer, Filter
from loopbuilder.typing import StrPath


class Builder(ABC):
    """Base class for loop builders"""

    def __init__(
        self,
        structure: StrPath,
        output_directory: StrPath,
        working_directory: StrPath | None = None,
        scorers: list[Scorer] | None = None,
        filters: list[Filter] | None = None,
    ):
        """Initialize the builder

        Args:
            structure: Path to the input structure file
            output_directory: Path to the output directory. Successful models will be saved here.

        Keyword args:
            working_directory: Path to a working directory. All trial models and intermediate files
                will be saved here. If `None`, uses a temporary working directory.
            scorers: List of scorers to use (optional)
            filters: List of filters to use (optional)
        """
        if working_directory is not None:
            working_directory = pathlib.Path(working_directory).resolve()
        self.working_directory = working_directory
        self.structure = pathlib.Path(structure).resolve()
        self.output_directory = pathlib.Path(output_directory).resolve()
        self.scorers = scorers or []
        self.filters = filters or []
        self.segments: list[Segment] = []

    def find_segments(self) -> None:
        """Find missing segments in the input structure

        NOTE: Uses PDBFixer at the moment but this part can also be
            implemented for example in Biopython if needed.
        """

        self.segments = []
        fixer = PDBFixer(filename=str(self.structure))
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

        for i, ((chain_index, residue_start_index), residue_names) in enumerate(non_terminal.items()):
            segment = Segment(
                identifier=f"loop_{i}",
                chain_index=chain_index,
                chain_name=chains[chain_index].id,
                residue_start_index=residue_start_index,
                residue_index_offset=int(next(chains[chain_index].residues()).id),
                residue_names=residue_names,
                parent_structure_file=self.structure,
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

        if n < 1:
            logger.warning("Building 0 models. Exiting.")
            return

        if max_tries is None:
            max_tries = n * 10

        if max_tries < 1:
            logger.warning("Using 0 trial models. Exiting.")
            return

        logger.info(f"Building n={n} models for {self.structure} (max_tries={max_tries})")
        self.output_directory.mkdir(parents=True, exist_ok=True)
        logger.info(f"Saving output to {self.output_directory}")

        scorer_str = "\n".join([f"  {scorer!r}" for scorer in self.scorers])
        logger.info(f"Using {len(self.scorers)} scorers\n{scorer_str}")

        filter_str = "\n".join([f"  {filter!r}" for filter in self.filters])
        logger.info(f"Using {len(self.filters)} filters\n{filter_str}")

        if not self.segments:
            logger.info("Looking for segments")
            self.find_segments()

        if not self.segments:
            logger.warning("No segments found")
            return self.segments

        logger.info(f"Found {len(self.segments)} segments")

        for segment in tqdm(self.segments):
            logger.info(f"Building models for segment {segment.identifier}")
            n_success = 0
            n_tries = 0

            if n_success >= n:
                logger.success(f"Reached the target number of models (success_rate={n_success / n_tries:.2%})")
                break

            if n_tries >= max_tries:
                logger.warning(f"Reached the maximum number of tries (success_rate={n_success / n_tries:.2%})")
                break

            segment_model = self.build_segment(segment, trial_id=str(n_tries), working_directory=working_directory)
            n_tries += 1

            for scorer in self.scorers:
                scorer.score(segment_model)

            for filter in self.filters:
                if not filter.filter(segment_model):
                    continue

            n_success += 1
            segment.models.append(segment_model)
            logger.success(f"Built model {n_success} for segment {segment.identifier}")

        return self.segments

    @abstractmethod
    def build_segment(self, segment: Segment, trial_id: str, working_directory: pathlib.Path) -> SegmentModel:
        """Build a model for a segment

        Args:
            segment: `Segment` object to build a model for
            trial_id: Trial ID for the model, used to generate a unique output filename
            working_directory: Path to a working directory. All trial models and intermediate files
                will be saved here.

        Returns:
            A `SegmentModel` object
        """


class PDBFixerBuilder(Builder):
    def build_segment(self, segment: Segment, trial_id: str, working_directory: pathlib.Path) -> SegmentModel:
        """Build a segment model using PDBFixer

        Args:
            segment: `Segment` object to build a model for
            trial_id: Trial ID for the model, used to generate a unique output filename

        Returns:
            A `SegmentModel` object
        """

        fixer = PDBFixer(str(segment.parent_structure_file))
        fixer.findMissingResidues = {(segment.chain_index, segment.residue_start_index): segment.residue_names}

        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        model_structure_file = (
            working_directory
            / f"{segment.parent_structure_file}_segment_{segment.chain_index}_{segment.residue_start_index}_{trial_id}.cif"
        )
        with open(model_structure_file, "w") as fp:
            PDBxFile.writeFile(fixer.topology, fixer.positions, file=fp, keepIds=True)

        return SegmentModel(identifier=segment.identifier, structure_file=model_structure_file, scores={})
