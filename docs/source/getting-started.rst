Getting started
===============

Install Ersilia from PyPI. Docker should be installed and running, since models
are fetched from DockerHub by default.

.. code-block:: bash

   pip install ersilia
   ersilia --help

Command line
------------

Fetch and serve a model, run it on a single-column CSV of SMILES, and close it
when you're done:

.. code-block:: bash

   ersilia fetch eos3b5e
   ersilia serve eos3b5e
   ersilia example -n 5 -o input.csv
   ersilia run -i input.csv -o output.csv
   ersilia close

See :doc:`cli` for all commands and options.

Python
------

.. code-block:: python

   from ersilia.api import Model

   model = Model("eos3b5e")
   model.fetch()
   with model:
       df = model.run(["CCO", "c1ccccc1"])

See :doc:`python-api` for the full API. For installation details and user
guides, see the `Ersilia Book <https://ersilia.gitbook.io/ersilia-book/>`_.
