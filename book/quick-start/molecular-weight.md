# Molecular weight

If you explore the [Ersilia Model Hub](https://airtable.com/shrUcrUnd7jB9ChZV), you will find a list of computational assets that you can use, for instance, to calculate or predict properties for your molecules of interest. A simple example of this would be the calculation of **molecular weight**. Let's use this as a toy example to quickly test the tool.

First, you have to **fetch** the model from our repository. The molecular weight calculator is named `molecular-weight`:

```bash
# download and install
ersilia fetch molecular-weight
```

Then, you can **serve** this tool so that it becomes available as an API:

```bash
# start a service
ersilia serve molecular-weight
```

Our **API** of interest is `calculate`. We can use it to calculate the molecular weight of Aspirin.

```bash
# calculate molecular weight of Aspirin using SMILES
ersilia api calculate -i "CC(=O)OC1=CC=CC=C1C(=O)O"
```

Now that we are done, we can **close** the service.

```
# close service
ersilia close
```

This is it! Of course, in real life you don't need a repository like the Ersilia Model Hub to calculate molecular weight. This is overkill - [RDKit ](https://www.rdkit.org)or any other chemoinformatics tool has built-in functions for this. But we hope that this example illustrates how Ersilia can be used. For a more detailed and realistic example, please refer to the **Antibiotic activity prediction** tutorial in the next page.
