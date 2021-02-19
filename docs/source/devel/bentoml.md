# Using Bento to pack and serve models

Core concepts are very clearly explained in the following [documentation](https://docs.bentoml.org/en/latest/concepts.html)

All credit goes to [this](https://towardsdatascience.com/pytorch-bentoml-heroku-the-simple-stack-65196aad181c) article.

BentoML is a tool that helps you build production-ready API endpoints for your ML models. All the ML frameworks are supported like TF, Keras, Pytorch, Scikit-learn and more. It is still relatively new and under active development, though the current APIs are stable enough for production use. Since its a new tool I’ll tell you more.

## Core Concepts

So the main idea of BentoML is to help Data Science teams ship prediction services faster while making it easy to test and deploy models following the best practices from DevOps. It's almost like a wrapper for your model which creates an API endpoint and making it easier to deploy. The following are the important concepts you use to create any prediction service using BentoML.

### Bento Service

This is the base class for building our prediction service and all the services we write inherited all its properties from this class. we define the properties specific to the API like what type of data the endpoint expects, what the dependencies are, how the model handles the data we get from the endpoint and more. In short, all the information on how to create the endpoint and pack the model is inferred from the attributes of this class.

### Model Artifacts

Artifacts are the trained models which are packed using bento. BentoML supports different ML frameworks and these are to be handled in their own ways. The model artifact handles serialization and deserialization automatically according to the artifact chosen corresponding to the ML framework used. For a complete list of artifacts supported check out the doc
For a full list of model artifacts, check out the [documentation](https://docs.bentoml.org/en/latest/api/artifacts.html).

### Input Adapters (former Handlers)

For a full list of input adapters, see [documentation](https://docs.bentoml.org/en/latest/api/adapters.html#).
These specify the type of data the model is expecting. I can be JSON, Pandas dataframe, Tensorflow tensor, images etc. The complete list of handlers is given in the doc.
