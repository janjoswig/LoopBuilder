import gemmi

from loopbuilder.typing import StrPath


def extract_segment_from_mmcif(
    input_file: StrPath,
    output_file: StrPath,
    *,
    residue_indices: set[int],
    chain_id: str,
) -> None:
    """Extract a subset of consecutive residue entries from one into another CIF file

    Args:
        input_file: Path to the input CIF file
        output_file: Path to the output CIF file

    Keyword args:
        residue_indices: Set of residue indices to keep in the output file.
            It is sufficient to specify the first and last residue index.
        chain_id: Chain ID to filter residues from
    """

    doc = gemmi.cif.read(str(input_file))
    block = doc[0]

    table = block.find("_atom_site.", ['label_asym_id', 'label_seq_id'])

    # NOTE: Cannot iterate over table and delete rows at the same time
    keep_rows = []
    for row in table:
        if (row[0] == chain_id) and (int(row[1]) in residue_indices):
            keep_rows.append(row.row_index)

    start, end = min(keep_rows), max(keep_rows)
    if start > 0:
        del table[:start]
    if end < len(table) - 1:
        # NOTE: Deleting rows shifts row indices, so subtract `start` from `end`
        del table[end - start + 1:]

    doc.write_file(str(output_file))


def join_segments(
    input_files: list[StrPath],
    output_file: StrPath,
) -> None:
    """Join multiple CIF files into a single CIF file with multiple models

    NOTE: Uses a naive concatenation of the text files for now as how to do this
        in `gemmi` with proper parsing is somewhat obscure. If this becomes a
        problem, a `gemmi`-based solution should be revisited.

    Args:
        input_files: List of paths to the input CIF files.
        output_file: Path to the output CIF file.

    Keyword args:
        chain_id: Chain ID to filter residues from
    """

    with open(str(output_file), "w") as fpo:
        with open(str(input_files[0]), "r") as fpi:
            for line in fpi:
                fpo.write(line)

            # NOTE: Assumes last line is last ATOM entry of last model
            prev_model_count = int(line.split()[-1])

        for input_file in input_files[1:]:
            with open(str(input_file), "r") as fpi:
                for line in fpi:
                    if line.startswith("ATOM"):
                        model_count = int(line.split()[-1])
                        fpo.write(" ".join(line.split()[:-1]) + f" {model_count + prev_model_count}\n")
                prev_model_count = model_count + prev_model_count