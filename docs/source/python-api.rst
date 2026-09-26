Python API
==========

The public Python API consists of two classes, :class:`~ersilia.api.Model` and
:class:`~ersilia.api.Catalog`. They mirror the CLI commands of the same name.

.. code-block:: python

   from ersilia.api import Model

   model = Model("eos4e40")
   model.fetch()
   with model:
       df = model.run(["CCO", "c1ccccc1"])

Model
-----

.. autoclass:: ersilia.api.Model
   :members: fetch, serve, run, close, info, example, delete, is_fetched

Catalog
-------

.. autoclass:: ersilia.api.Catalog
   :members: catalog
