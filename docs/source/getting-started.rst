Getting started
===============

Install Ersilia from PyPI. Docker should be installed and running, since models
are fetched from DockerHub by default.

.. code-block:: bash

   # install the Ersilia CLI and Python package
   pip install ersilia

   # check that the installation worked
   ersilia --help

Find a model
------------

Every model in the Ersilia Model Hub has an identifier such as ``eos4e40``.
Browse the models and their identifiers in the
`Ersilia Model Hub catalog <https://catalog.ersilia.io>`_, or list them from
the terminal:

.. code-block:: bash

   # list all models available in the Ersilia Model Hub
   ersilia catalog --hub

Command line
------------

Fetch and serve a model, run it on a single-column CSV of SMILES, and close it
when you're done:

.. code-block:: bash

   # download the model (antibiotic activity prediction, Stokes et al. 2020)
   ersilia fetch eos4e40

   # start the model server in this terminal
   ersilia serve eos4e40

   # generate 5 example inputs (SMILES) for the model
   ersilia example -n 5 -o input.csv

   # run predictions and save them to a CSV file
   ersilia run -i input.csv -o output.csv

   # stop the model server
   ersilia close

See :doc:`cli` for all commands and options.

Python
------

.. code-block:: python

   from ersilia.api import Model

   # create a handle for the model
   model = Model("eos4e40")

   # download the model
   model.fetch()

   # serve the model inside the block; it is closed automatically at the end
   with model:
       # run predictions on a list of SMILES; returns a pandas DataFrame
       df = model.run(["CCO", "c1ccccc1"])

See :doc:`python-api` for the full API. For installation details and user
guides, see the `Ersilia Book <https://ersilia.gitbook.io/ersilia-book/>`_.
