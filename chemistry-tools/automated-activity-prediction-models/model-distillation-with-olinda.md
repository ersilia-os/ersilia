---
description: Olinda is a model distillation tool for chemistry data
---

# Model distillation with Olinda

[Olinda](https://github.com/ersilia-os/olinda) is a generic cheminformatics model distillation library. It can automatically distill models from Pytorch, Tensorflow, ONNX and models fetched directly from the [Ersilia Model Hub](https://ersilia.io/model-hub).

## Getting started

Olinda is available on PyPi and can be installed using pip.

```python
from olinda import distill

student_model = distill(your_model)
```

### How the distillation works?

The distillation function first downloads a reference SMILES dataset if it is not already present. It then generates featurized inputs using the reference SMILES dataset for training the student model. Next it uses the provided teacher model to generate input-output pairs. The input-output pairs together with the featurized inputs constitute the training dataset for the student model. Finally a suitable architecture for the student model is selected using heuristics and the selected model is trained using the training dataset.

1. Generate reference SMILES dataset
2. Generate a training dataset using the given teacher model
3. Search a suitable architecture for the student model
4. Train student model

During the distillation process, helpful messages and progress bars are printed to keep the user informed. In the case of a crash or process interruption the distillation process can be resumed automatically. It caches all the intermediate results in a local directory (`xdg_home() / olinda`).

### Distillation customization

The distillation API is very flexible and covers a wide varietry of use cases. User can easily customize the distillation behavior by passing parameters to the `distill` function.

```python
def distill(
    model: Any,
    featurizer: Optional[Featurizer],
    working_dir: Path,
    clean: bool = False,
    tuner: ModelTuner = AutoKerasTuner(),
    reference_smiles_dm: Optional[ReferenceSmilesDM] = None,
    featurized_smiles_dm: Optional[FeaturizedSmilesDM] = None,
    generic_output_dm: Optional[GenericOutputDM] = None,
) -> pl.LightningModule:
    """Distill models.

    Args:
        model (Any): Teacher Model.
        featurizer (Optional[Featurizer]): Featurizer to use.
        working_dir (Path): Path to model workspace directory.
        clean (bool): Clean workspace before starting.
        tuner (ModelTuner): Tuner to use for selecting and optimizing student model.
        reference_smiles_dm (Optional[ReferenceSmilesDM]): Reference SMILES datamodules.
        featurized_smiles_dm (Optional[FeaturizedSmilesDM]): Reference Featurized SMILES datamodules.
        generic_output_dm (Optional[GenericOutputDM]): Precalculated training dataset for student model.

    Returns:
        pl.LightningModule: Student Model.
    """
```

