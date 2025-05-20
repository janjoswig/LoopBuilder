import pathlib
from abc import ABC, abstractmethod
from loguru import logger

from loopbuilder.typing import StrPath


class Builder(ABC):
    """Base class for loop builders"""

    def __init__(self, structure: StrPath, output_directory: StrPath):
        self.structure = pathlib.Path(structure)
        self.output_directory = pathlib.Path(output_directory)

    @abstractmethod
    def build(self, n: int = 1) -> None:
        """Build n models for the input structure"""
