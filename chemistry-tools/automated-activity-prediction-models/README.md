---
description: >-
  We are developing AutoML tools for chemistry data to facilitate adoption of
  AI/ML
---

# Automated activity prediction models

Quick baseline modeling of chemistry data can be done with [LazyQSAR](https://github.com/ersilia-os/lazy-qsar), our fast modelling tool. LazyQSAR produces light-weight models for binary classification and regression tasks.&#x20;

{% content-ref url="light-weight-automl-with-lazyqsar.md" %}
[light-weight-automl-with-lazyqsar.md](light-weight-automl-with-lazyqsar.md)
{% endcontent-ref %}

Our flagship AutoML tool for chemistry is [ZairaChem](https://github.com/ersilia-os/zaira-chem). This Python library offers robust ensemble-based modeling capabilities applicable to a wide range of modeling scenarios. At the moment, ZairaChem is focused on binary classification and regression tasks.

{% content-ref url="accurate-automl-with-zairachem.md" %}
[accurate-automl-with-zairachem.md](accurate-automl-with-zairachem.md)
{% endcontent-ref %}

In addition, we have developed a model distillation pipeline named [Olinda](https://github.com/ersilia-os/olinda) aimed at producing light, interoperable models in [ONNX](https://onnx.ai/) format.

{% content-ref url="model-distillation-with-olinda.md" %}
[model-distillation-with-olinda.md](model-distillation-with-olinda.md)
{% endcontent-ref %}
