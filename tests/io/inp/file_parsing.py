# tests/io/inp/test_file_parsing.py
"""INP Parser file parsing tests."""

import os
import pytest
import pygcmc


def test_read_gcmc_inp():
    """Test reading the actual gcmc.inp file"""
    # Get the path to the test data directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    gcmc_inp_path = os.path.join(current_dir, "..", "..", "data", "gcmc.inp")
    
    # Parse the file
    param = pygcmc.io.INPParser.parse_file(gcmc_inp_path)

    # Test file paths
    assert len(param.file_info.par_files) == 3
    assert param.file_info.par_files == ["ffnonbonded.itp", "silcs.itp", "nbfix.itp"]
    assert len(param.file_info.fragment_top_files) == 9
    assert param.file_info.fragment_top_files == [
        "mol/benx.itp", "mol/prpx.itp", "mol/dmee.itp", 
        "mol/meoh.itp", "mol/form.itp", "mol/imia.itp",
        "mol/acey.itp", "mol/mamy.itp", "mol/sol.itp"
    ]
    assert param.file_info.atomtype_file == "atomtypes.atp"
    assert param.file_info.monomer_dir == "mol"
    assert param.file_info.topology_file == "test.top"
    assert param.file_info.input_pdb_file == "test.pdb"
    assert param.file_info.output_top_file == "test.1.gc.0.top"
    assert param.file_info.output_pdb_file == "test.1.gc.0.pdb"
    assert param.file_info.map_prefix == "test_silcs.gc.1.0"

    # Test space parameters
    assert param.space_info.grid_spacing == pytest.approx(1.0)
    assert param.space_info.box_size[0] == pytest.approx(36.736)
    assert param.space_info.box_size[1] == pytest.approx(40.850)
    assert param.space_info.box_size[2] == pytest.approx(49.379)
    assert param.space_info.cutoff == pytest.approx(12.0)
    assert param.space_info.gc_center[0] == pytest.approx(33.368)
    assert param.space_info.gc_center[1] == pytest.approx(35.425)
    assert param.space_info.gc_center[2] == pytest.approx(39.690)
    assert param.space_info.sys_center[0] == pytest.approx(33.368)
    assert param.space_info.sys_center[1] == pytest.approx(35.425)
    assert param.space_info.sys_center[2] == pytest.approx(39.690)

    # Test fragment parameters
    assert len(param.file_info.fragment_names) == 9
    assert param.file_info.fragment_names == [
        "benx", "prpx", "dmee", "meoh", "form", 
        "imia", "acey", "mamy", "sol"
    ]
    assert len(param.fragment_info.conc_list) == 9
    assert param.fragment_info.conc_list == pytest.approx([
        0.25, 0.25, 0.25, 0.25, 0.25, 
        0.25, 0.25, 0.25, 55.00
    ])
    assert len(param.fragment_info.muex_list) == 9
    assert param.fragment_info.muex_list == pytest.approx([
        -0.79, 1.96, -1.79, -5.36, -10.92,
        -14.18, -97.31, -68.49, -5.60
    ])

    # Test MC parameters
    assert param.mc_info.print_freq == 1000
    assert param.mc_info.mc_steps == 10000

    # Test bias parameters
    assert not param.bias_info.use_cavity_bias
    assert not param.bias_info.use_conf_bias

    # Test basic info parameters
    assert param.basic_info.init_cycle
    assert not param.basic_info.conserve_fragments
    assert not param.file_info.generate_maps


def test_inp_parser_file(tmp_path):
    """Test parsing INP parameters from a file."""
    # Create a temporary input file
    inp_content = """par:ffnonbonded.itp
par:silcs.itp
par:nbfix.itp
fragitp:mol/benx.itp
fragitp:mol/prpx.itp
atomtypes:atomtypes.atp
monomerdir:mol
top:test.top
pdb:test.pdb
protitp:test.top
grid_dx:   1.000
box_size:  36.736   40.850   49.379
cutoff:  12.000
gc_center:  33.368   35.425   39.690
sys_center:  33.368   35.425   39.690
fragname:     benx   prpx   dmee   meoh   form   imia   acey   mamy    sol 
fragconc:     0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.25  55.00 
fragmuex:    -0.79   1.96  -1.79  -5.36 -10.92 -14.18 -97.31 -68.49  -5.60 
nprint:1000
conserve_frags:no
map_generation:no
op_top:test.1.gc.0.top
op_pdb:test.1.gc.0.pdb
initcycle:yes
map_filename_prefix:test_silcs.gc.1.0
mcsteps:10000
use_cavity_bias:no
use_conf_bias:no"""

    # Write test input file
    inp_file = tmp_path / "test.inp"
    inp_file.write_text(inp_content)

    # Parse the file
    param = pygcmc.io.INPParser.parse_file(str(inp_file))

    # Test file info parameters
    assert len(param.file_info.par_files) == 3
    assert param.file_info.par_files[0] == "ffnonbonded.itp"
    assert len(param.file_info.fragment_top_files) == 2
    assert param.file_info.fragment_top_files[0] == "mol/benx.itp"
    assert param.file_info.atomtype_file == "atomtypes.atp"
    assert param.file_info.monomer_dir == "mol"
    assert param.file_info.topology_file == "test.top"
    assert param.file_info.input_pdb_file == "test.pdb"
    assert param.file_info.output_top_file == "test.1.gc.0.top"
    assert param.file_info.output_pdb_file == "test.1.gc.0.pdb"
    assert param.file_info.map_prefix == "test_silcs.gc.1.0"
    assert not param.file_info.generate_maps

    # Test space info parameters
    assert param.space_info.grid_spacing == pytest.approx(1.0)
    assert param.space_info.box_size[0] == pytest.approx(36.736)
    assert param.space_info.box_size[1] == pytest.approx(40.850)
    assert param.space_info.box_size[2] == pytest.approx(49.379)
    assert param.space_info.cutoff == pytest.approx(12.0)
    assert param.space_info.gc_center[0] == pytest.approx(33.368)
    assert param.space_info.gc_center[1] == pytest.approx(35.425)
    assert param.space_info.gc_center[2] == pytest.approx(39.690)
    assert param.space_info.sys_center[0] == pytest.approx(33.368)
    assert param.space_info.sys_center[1] == pytest.approx(35.425)
    assert param.space_info.sys_center[2] == pytest.approx(39.690)

    # Test fragment info parameters
    assert len(param.file_info.fragment_names) == 9
    assert param.file_info.fragment_names[0] == "benx"
    assert len(param.fragment_info.conc_list) == 9
    assert param.fragment_info.conc_list[0] == pytest.approx(0.25)
    assert param.fragment_info.conc_list[-1] == pytest.approx(55.00)
    assert len(param.fragment_info.muex_list) == 9
    assert param.fragment_info.muex_list[0] == pytest.approx(-0.79)
    assert param.fragment_info.muex_list[-1] == pytest.approx(-5.60)

    # Test MC info parameters
    assert param.mc_info.print_freq == 1000
    assert param.mc_info.mc_steps == 10000

    # Test bias info parameters
    assert not param.bias_info.use_cavity_bias
    assert not param.bias_info.use_conf_bias

    # Test basic info parameters
    assert param.basic_info.init_cycle
    assert not param.basic_info.conserve_fragments