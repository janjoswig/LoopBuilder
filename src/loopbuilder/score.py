import shlex
import subprocess
from abc import ABC, abstractmethod
from typing import Any

from loopbuilder.segment import SegmentModel
from loopbuilder.typing import StrPath


class SegmentModelEvaluatorBase:
    """Method-providing common base for `Scorer` and `Filter` classes
    
    On init, the `parameters` dictionary will be populated with default parameters
    from (in that order):
        1. `SegmentModelEvaluatorBase.__default_parameters`
        2. `<Subclass>.__default_parameters`
        3. User-provided keyword arguments
    """

    __default_parameters = {"cleanup": True}

    def __init__(self, identifier: str | None = None, **kwargs: Any):
        if identifier is None:
            identifier = self.__class__.__name__
        self.identifier = identifier

        subclass_default_parameters = getattr(self, f"_{self.__class__.__name__}__default_parameters", {})
        self.parameters = {**SegmentModelEvaluatorBase.__default_parameters, **subclass_default_parameters, **kwargs}

    def __repr__(self) -> str:
        param_str = ", ".join(f"{k}={v}" for k, v in self.parameters.items())
        if param_str:
            param_str += ", "
        return f"{self.identifier}({param_str})"


class Scorer(ABC, SegmentModelEvaluatorBase):
    """Base class for loop scorers"""

    @abstractmethod
    def score(self, model: SegmentModel) -> None:
        """Score a model

        Compute arbitrary quantities for a model and store them in the model's `scores` dictionary.
        """

    def __call__(self, model: SegmentModel) -> None:
        self.score(model)


class Filter(ABC, SegmentModelEvaluatorBase):
    """Base class for loop filters"""

    @abstractmethod
    def filter(self, model: SegmentModel) -> bool:
        """Filter a model

        Look at arbitrary quantities for a model. Return `True` if the model passes the filter
        and `False` otherwise.
        """

    def __call__(self, model: SegmentModel) -> bool:
        return self.filter(model)


class MolProbityScorer(Scorer):
    """MolProbity-based loop scoring"""

    __default_parameters = {
        "executable": "molprobity.molprobity",
    }

    def score(self, model: SegmentModel) -> None:
        """Score a model using MolProbity

        This method runs the MolProbity executable on the model's structure file and parses the output.
        The scores are stored in the model's `scores` dictionary. Currently read scores are:

          * `"clashscore"`: Clash score (0 is best)
          * `"ramachandran_outliers"`: Ramachandran outlier rate
          * `"rotamer_outliers"`: Rotamer outlier rate
          * `"cbeta_deviations"`: Absolute number of C-beta deviations
          * `"rms_bonds"`: Root-mean-square deviation of bond lengths from optimal values
          * `"rms_angles"`: Root-mean-square deviation of bond angles from optimal values
          * `"molprobity_score"`: MolProbity total score (0 is best)

        If an error occurs while running MolProbity, the error message is stored in the model's `scores` dictionary
        under the key `"molprobity_error"`.
        """

        exe = self.parameters["executable"]

        # NOTE: This is a bit of a hack to optionally run MolProbity through a Docker container
        # NOTE: MolProbity writes its main output file (molprobity.out) to the current working directory
        docker_image = self.parameters.get("docker_image")
        if docker_image is not None:
            structure_dir = str(model.structure_file.parent).replace('\\', '/')
            cmd = f"docker run -v {structure_dir}:/data -w /data {docker_image} {exe} {model.structure_file.name}"
            cwd = None
        else:
            cmd = f"{exe} {model.structure_file.name}"
            cwd = model.structure_file.parent

        result = subprocess.run(shlex.split(cmd), capture_output=True, encoding="utf-8", cwd=cwd)

        if result.returncode != 0:
            model.scores["molprobity_error"] = result.stderr
            return

        output_file = model.structure_file.with_name("molprobity.out")
        scores = self.parse_output(output_file)
        model.scores.update(scores)
        
        if self.parameters.get("cleanup", False):
            output_file.unlink(missing_ok=True)
            output_file.with_name("molprobity_probe.txt").unlink(missing_ok=True)
            output_file.with_name("molprobity_coot.py").unlink(missing_ok=True)

    def parse_output(self, output_file: StrPath) -> dict[str, Any]:
        """Parse the output of MolProbity and extract scores"""

        score_key_map = {
            # Key in molprobity.out, used key, converter
            "Clashscore": ("clashscore", float),
            "Ramachandran outliers": ("ramachandran_outliers", lambda x: round(float(x.strip(" %")) / 100, 4)),
            "Rotamer outliers": ("rotamer_outliers", lambda x: round(float(x.strip(" %")) / 100, 4)),
            "C-beta deviations": ("cbeta_deviations", float),
            "RMS(bonds)": ("rms_bonds", float),
            "RMS(angles)": ("rms_angles", float),
            "MolProbity score": ("molprobity_score", float),
        }

        scores = {}
        with open(output_file, "r") as fp:
            # Find summary section
            for line in fp:
                if not line.startswith("==="):
                    continue

                try:
                    section = line.split()[1]
                except IndexError:
                    continue

                if section == "Summary":
                    break

            for line in fp:
                try:
                    key, value = line.split("=")
                except ValueError:
                    continue

                try:
                    key_, convert = score_key_map[key.strip()]
                except KeyError:
                    continue

                scores[key_] = convert(value.strip())
        return scores
