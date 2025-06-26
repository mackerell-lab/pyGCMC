# tests/io/inp/test_string_parsing.py
"""INP Parser string parsing tests."""

import pytest
import pygcmc


def test_inp_parser_string():
    """Test parsing INP parameters from a string."""
    # Test parsing from string
    inp_content = """top:test.top
pdb:test.pdb
grid_dx:1.0
box_size:10.0 10.0 10.0
cutoff:12.0
fragname:benx prpx
fragconc:0.25 0.25
fragmuex:-0.79 1.96
mcsteps:1000"""

    param = pygcmc.io.INPParser.parse_string(inp_content)

    assert param.file_info.topology_file == "test.top"
    assert param.file_info.input_pdb_file == "test.pdb"
    assert param.space_info.grid_spacing == pytest.approx(1.0)
    assert param.space_info.cutoff == pytest.approx(12.0)
    assert len(param.file_info.fragment_names) == 2
    assert len(param.fragment_info.conc_list) == 2
    assert len(param.fragment_info.muex_list) == 2
    assert param.mc_info.mc_steps == 1000