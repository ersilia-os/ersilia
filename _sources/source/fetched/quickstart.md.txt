# Quick start

Ersilia ML models are stored in GitHub and Open Science Framework. Once `ersilia` package is installed, you can fetch each model individually.

Each model has a unique identifier in the format `eos0abc`.

## Command Line Interface

First, download your model of interest. It will be stored locally at `~/eos`.

```bash
ersilia fetch eos0abc
```

You can check the catalog of models available in your computer.

```bash
ersilia catalog --local
```

For more information, check the model card.
```bash
ersilia card eos0abc
```

Serve your model. A URL will be displayed, together with the APIs available for the model.
```bash
ersilia serve eos0abc
```

To run the model, use your API of choice. For example, `predict` for Caffeine.
```bash
ersilia api eos0abc predict "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
```

Don't forget to close the model when you are finished.
```bash
ersilia close eos0abc
```

If you don't want to use this model anymore, remove it from your computer.
```bash
ersilia delete eos0abc
```

## Python Package

Ersilia can be also run as a Python package.

```python
from ersilia import ErsiliaModel

mdl = ErsiliaModel('eos0abc') # fetches model if not available locally
mdl.serve()
mdl.predict("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
mdl.close()
```
