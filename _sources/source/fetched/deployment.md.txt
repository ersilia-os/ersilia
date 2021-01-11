# Model deployment

Ersilia is powered by the open-source [BentoML](http://bentoml.ai) library. BentoML is an outstanding tool that makes it seamlessly easy to bring ML models to productions.
Accordingly, most functionalities available from the Ersilia CLI are directly inherited by the BentoML command line tool.

The following is a limited list of deployment functionalities inherited directly from BentoML. For more details, please visit
the [BentoML docs](https://docs.bentoml.org/en/latest/index.html).

## Get an overview of your local repository

By default, BentoML stores all model files in the `~/bentoml` directory.
This directory contains all the code, configs and files required for deployment.

```bash
ersilia list
```

For a prettier view of the list in your browser, you can try the Yatai service following [these instructions](https://docs.bentoml.org/en/latest/quickstart.html#save-prediction-service-for-distribution).

## Command-line interface

```bash
ersilia run eos0aaa --input
```

## Web based UI

To start a REST API server locally, simply use the `serve` command:

```bash
ersilia serve eos0aaa
```

The `eos0aaa` model is served at `localhost:5000`. You can go to [http://localhost:5000](http://example.com) and use a simple UI.

## Check BentoML!

BentoML many more functionalities and all of those can be applied to Ersilia models. Visit the [BentoML docs](https://docs.bentoml.org/en/latest/quickstart.html) to know more.
