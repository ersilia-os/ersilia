import os
import sys

from ersilia.io.input import GenericInputAdapter

root = os.path.abspath(os.path.dirname(__file__))
inputs_path = os.path.abspath(os.path.join(root, "..", "inputs"))
sys.path.append(inputs_path)
from ..inputs.compound_list import smiles_list as compound_list_input
from ..inputs.compound_lists import smiles_lists as compound_lists_input
from ..inputs.compound_single import smiles as compound_single_input
from ..inputs.compound_singles import smiles as compound_singles_input


def test_compound_single():
    adapter = GenericInputAdapter(input_type="compound", input_shape="single")
    inp_csv = os.path.abspath(os.path.join(inputs_path, "compound_single.csv"))
    inp_json = os.path.abspath(os.path.join(inputs_path, "compound_single.json"))
    inp_py = compound_single_input
    d_csv = [d for d in adapter.adapt_one_by_one(inp_csv)]
    d_json = [d for d in adapter.adapt_one_by_one(inp_json)]
    d_py = [d for d in adapter.adapt_one_by_one(inp_py)]
    assert d_csv == d_json == d_py


def test_compound_singles():
    adapter = GenericInputAdapter(input_type="compound", input_shape="single")
    inp_csv = os.path.abspath(os.path.join(inputs_path, "compound_singles.csv"))
    inp_json = os.path.abspath(os.path.join(inputs_path, "compound_singles.json"))
    inp_py = compound_singles_input
    d_csv = [d for d in adapter.adapt_one_by_one(inp_csv)]
    d_json = [d for d in adapter.adapt_one_by_one(inp_json)]
    d_py = [d for d in adapter.adapt_one_by_one(inp_py)]
    assert d_csv == d_json == d_py
