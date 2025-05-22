import pathlib
import shutil
import tempfile

import mdtraj

from loopbuilder.convert import join_segments


SEGMENT_MODEL_1 = """
data_cell
_cell.length_a 93.4570
_cell.length_b 93.4570
_cell.length_c 166.6790
_cell.angle_alpha 90.0000
_cell.angle_beta 90.0000
_cell.angle_gamma 90.0000

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x_esd
_atom_site.Cartn_y_esd
_atom_site.Cartn_z_esd
_atom_site.occupancy_esd
_atom_site.B_iso_or_equiv_esd
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 3250 N N . GLU B ? 600 . -7.9148 21.0821 35.5826 0.0 0.0 ? ? ? ? ? . 600 GLU B N 1
ATOM 3251 C CA . GLU B ? 600 . -9.4031 21.1704 36.1498 0.0 0.0 ? ? ? ? ? . 600 GLU B CA 1
ATOM 3252 C CB . GLU B ? 600 . -8.8760 21.6280 37.4816 0.0 0.0 ? ? ? ? ? . 600 GLU B CB 1
ATOM 3253 C CG . GLU B ? 600 . -10.0507 22.0078 38.3734 0.0 0.0 ? ? ? ? ? . 600 GLU B CG 1
ATOM 3254 C CD . GLU B ? 600 . -10.0608 23.2786 39.2084 0.0 0.0 ? ? ? ? ? . 600 GLU B CD 1
ATOM 3255 O OE1 . GLU B ? 600 . -9.6763 24.3452 38.6682 0.0 0.0 ? ? ? ? ? . 600 GLU B OE1 1
ATOM 3256 O OE2 . GLU B ? 600 . -10.4644 23.1338 40.3786 0.0 0.0 ? ? ? ? ? . 600 GLU B OE2 1
ATOM 3257 C C . GLU B ? 600 . -10.5877 20.0467 36.3904 0.0 0.0 ? ? ? ? ? . 600 GLU B C 1
ATOM 3258 O O . GLU B ? 600 . -10.4480 18.9385 35.8789 0.0 0.0 ? ? ? ? ? . 600 GLU B O 1
"""

SEGMENT_MODEL_2 = """
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x_esd
_atom_site.Cartn_y_esd
_atom_site.Cartn_z_esd
_atom_site.occupancy_esd
_atom_site.B_iso_or_equiv_esd
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 3250 N N . GLU B ? 600 . -7.8856 21.1628 35.6874 0.0 0.0 ? ? ? ? ? . 600 GLU B N 1
ATOM 3251 C CA . GLU B ? 600 . -9.2370 21.0021 36.5602 0.0 0.0 ? ? ? ? ? . 600 GLU B CA 1
ATOM 3252 C CB . GLU B ? 600 . -8.9416 21.4748 38.0327 0.0 0.0 ? ? ? ? ? . 600 GLU B CB 1
ATOM 3253 C CG . GLU B ? 600 . -9.2919 22.9694 38.4026 0.0 0.0 ? ? ? ? ? . 600 GLU B CG 1
ATOM 3254 C CD . GLU B ? 600 . -9.0269 24.3883 37.7703 0.0 0.0 ? ? ? ? ? . 600 GLU B CD 1
ATOM 3255 O OE1 . GLU B ? 600 . -8.7756 25.3357 36.9418 0.0 0.0 ? ? ? ? ? . 600 GLU B OE1 1
ATOM 3256 O OE2 . GLU B ? 600 . -8.9361 23.5737 36.8202 0.0 0.0 ? ? ? ? ? . 600 GLU B OE2 1
ATOM 3257 C C . GLU B ? 600 . -10.8806 21.0456 36.2740 0.0 0.0 ? ? ? ? ? . 600 GLU B C 1
ATOM 3258 O O . GLU B ? 600 . -10.9994 21.7169 35.2300 0.0 0.0 ? ? ? ? ? . 600 GLU B O 1
"""


def test_join_segments():
    tmpdir = tempfile.mkdtemp(dir=".")
    tmpdir = pathlib.Path(tmpdir)
    input_file_1 = tmpdir / "segment_model_1.cif"
    input_file_2 = tmpdir / "segment_model_2.cif"
    output_file = tmpdir / "joined_segments.cif"

    with open(input_file_1, "w") as f:
        f.write(SEGMENT_MODEL_1)

    with open(input_file_2, "w") as f:
        f.write(SEGMENT_MODEL_2)

    join_segments([input_file_1, input_file_2], output_file)

    with open(output_file, "r") as f:
        content = f.read()

    assert SEGMENT_MODEL_1 in content

    traj = mdtraj.load(str(output_file))

    assert traj.n_frames == 2
    assert traj.n_atoms == 9

    shutil.rmtree(tmpdir, ignore_errors=True)