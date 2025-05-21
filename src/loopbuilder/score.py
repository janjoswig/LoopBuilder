from abc import ABC, abstractmethod
from typing import Any

from loopbuilder.segment import SegmentModel


class SegmentModelEvaluatorBase:
    __default_parameters = {}

    def __init__(self, identifier: str | None = None, **kwargs: Any):
        if identifier is None:
            identifier = self.__class__.__name__
        self.identifier = identifier
        self.parameters = {**SegmentModelEvaluatorBase.__default_parameters, **self.__default_parameters, **kwargs}

    def __repr__(self) -> str:
        param_str = ", ".join(f"{k}={v}" for k, v in self.parameters.items())
        if param_str:
            param_str += ", "
        return f"{self.__class__.__name__}({self.identifier}{param_str})"


class Scorer(ABC, SegmentModelEvaluatorBase):
    """Base class for loop scorers"""

    @abstractmethod
    def score(self, model: SegmentModel) -> None:
        """Score a model

        Compute arbitrary quantities for a model and store them in the model's `scores` dictionary.
        """


class Filter(ABC, SegmentModelEvaluatorBase):
    """Base class for loop filters"""

    @abstractmethod
    def filter(self, model: SegmentModel) -> bool:
        """Filter a model

        Look at arbitrary quantities for a model. Return `True` if the model passes the filter
        and `False` otherwise.
        """
