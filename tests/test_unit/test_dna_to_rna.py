#write the script to test the dna_to_rna function 

import pytest 
import logging
from GeneBankProject.ToolKit.modules.dna_to_rna import convert_dna_to_rna

def integer_value():
    with pytest.raises(TypeError):
        convert_dna_to_rna(10)

def non_base_value():
    with pytest.raises(ValueError):
        convert_dna_to_rna(agctagctx)

def no_value():
    dna = ""
    expected = convert_dna_to_rna(dna)
    assert dna == expected

def expected_string():
    rna = "acguacguacguacgaucgau"
    expected = convert_dna_to_rna(dna)
    assert rna == expected

def mixed_case():
    assert convert_dna_to_rna("AgCt") == "agcu"

def test_long_sequence():
    dna = "ATGC" * 1000
    rna = convert_dna_to_rna(dna)
    assert len(rna) == len(dna)
    assert "t" not in rna
    assert "u" in rna

@pytest.mark.parametrize("dna,expected", [
    ("A", "a"),
    ("T", "u"),
    ("C", "c"),
    ("G", "g"),
    ("a", "a"),
    ("t", "u"),
    ("c", "c"),
    ("g", "g"),
])
def test_single_base_conversion(dna, expected):
    assert convert_dna_to_rna(dna) == expected