---
description: High-level information about the packaging of models using Ersilia Pack
---

# Model packaging

To package models in Ersilia, there are two essential components, namely the [Ersilia model template](../model-contribution/model-template.md), and the [Ersilia Pack](https://github.com/ersilia-os/ersilia-pack) tool.  Ersilia Pack creates a FastAPI app based on the information available in a given model repository, based on the template.

**Ersilia Pack** is a lightweight, plug-and-play framework to serve AI models from the [Ersilia Model Hub](https://ersilia.io) via an API.

It allows developers to:

* Package models with their dependencies
* Run them as web services
* Submit data and retrieve predictions
* Monitor jobs and system performance
* Access detailed model metadata

[Detailed documentation](https://github.com/ersilia-os/ersilia-pack) on Ersilia Pack is available in the corresponding GitHub repository. Here, we simply provide a high-level overview of the tool.

### Typical Workflow

1. **Lint a model repository** to check it has the correct structure
2. **Package the model** into a self-contained bundle
3. **Serve the model via API** locally or in the cloud
4. **Interact with the model** through documented endpoints
