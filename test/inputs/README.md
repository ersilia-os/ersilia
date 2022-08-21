# Inputs for testing

In this folder we outline different kinds of inputs for testing.

## Input formats

* `.csv`: Input in tabular format.
* `.json`: Input in JSON format.
* `.py`: These files are used to mimick Python variables, i.e. useful when we use Ersilia in the Python API.

## Chemistry

The molecules (drugs) considered in this files are the following (in SMILES format):

```
CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O # artemisin
C1=CN=CC=C1C(=O)NN # isoniazid
CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O # tenofovir
CC(=O)OC1=CC=CC=C1C(=O)O # aspirin
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O # ibuprofen
CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C # remdesivir
COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5 # cephalotaxin
```