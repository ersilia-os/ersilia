# Inputs for Testing

This folder contains various input files used for testing different components of the Ersilia project. Each test may require specific inputs in different formats, and these are organized to maintain clarity and separation.

## Folder Structure

- **Test-Specific Folders**: Each test has its own dedicated folder, named after the test. Inside these folders, you'll find input files specific to that test.
  - For example, the folder `test_inputs/` contains inputs for the `test_inputs.py` test.

## Adding New Test Inputs

When adding new tests, create a folder named after the test and place all relevant input files inside. Make sure to follow the format conventions mentioned above.

## Input Formats

The input files are available in different formats to accommodate various needs during testing:

* `.csv`: Input in tabular format.
* `.json`: Input in JSON format.
* `.py`: These files are used to mimick Python variables, i.e. useful when we use Ersilia in the Python API.

## Chemistry

Some tests such as test_inputs require chemical structures of molecules (drugs), which are following (in SMILES format):

```
CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O # artemisin
C1=CN=CC=C1C(=O)NN # isoniazid
CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O # tenofovir
CC(=O)OC1=CC=CC=C1C(=O)O # aspirin
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O # ibuprofen
CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C # remdesivir
COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5 # cephalotaxin
```