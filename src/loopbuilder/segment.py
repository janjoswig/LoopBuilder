import pathlib
from dataclasses import dataclass
from typing import Any


@dataclass
class SegmentModel:
    identifier: str
    structure_file: pathlib.Path
    scores: dict[str, Any]
    index: int = 0


@dataclass
class Segment:
    identifier: str
    chain_index: int
    chain_name: str
    residue_start_index: int
    residue_start_seqid: int
    residue_index_offset: int
    residue_names: list[str]
    parent_structure_file: pathlib.Path
    models: list[SegmentModel]

    def __len__(self) -> int:
        return len(self.residue_names)
