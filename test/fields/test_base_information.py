
import os
import pytest
import sys
from ersilia.hub.content.base_information import BaseInformation

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.github/scripts")))

from airtableops import ReadmeMetadata, AirtableMetadata

@pytest.fixture
def valid_data():
    data = BaseInformation()
    data.title = "Molecular weight calculator"
    data.description = "The model is simply an implementation of the function Descriptors.MolWt of the chemoinformatics package RDKIT. \
    It takes as input a small molecule (SMILES) and calculates its molecular weight in g/mol."
    data.identifier = "eos3b5e"
    data.slug = "molecular-wight"
    data.input = ["Compound"]
    data.input_shape = "Single"
    data.task = "Regression"
    data.source = "Online"
    data.output = "Other value"
    data.output_type = ["Float"]
    data.output_shape = "Single"
    data.interpretation = "This is a regression model."
    data.computational_performance_one = 10.11
    data.computational_performance_ten = 14.12
    data.computational_performance_hund = 50.45
    data.publication = "https://publication.url"
    data.source_code = "https://source.code"
    data.contributor = "sample-contributor"
    data.license = "MIT"
    return data

@pytest.fixture
def readme_metadata(valid_data):
    return ReadmeMetadata(model_id="eos3b5e")

def test_write_information(readme_metadata, valid_data):
    generated_text = readme_metadata.write_information(valid_data, readme_path=None) 
    
    assert "## Identifiers" in generated_text
    assert "## Characteristics" in generated_text
    assert "* [Publication]" in generated_text
    assert "* [Source Code]" in generated_text
    assert "## Citation" in generated_text

def test_write_information_to_file(readme_metadata, valid_data):
    test_file_path = "test_readme.md"
    
    readme_metadata.write_information(valid_data, readme_path=test_file_path)
    
    assert os.path.exists(test_file_path)
    
    with open(test_file_path, "r") as file:
        content = file.read()
        assert "## Identifiers" in content
    
    os.remove(test_file_path)