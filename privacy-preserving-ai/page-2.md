---
description: >-
  This page describes our ChemXor library, a tool to build privacy-preserving
  machine learning models for drug discovery. We believe that encrypted assets
  can foster collaboration
---

# Encryption of AI/ML models

## Privacy preserving AI/ML for drug discovery

We are developing an open-source **privacy-preserving machine learning platform** for drug discovery. It has been widely argued that artificial intelligence and machine learning (AI/ML) can transform the pharmaceutical industry. However, AI/ML techniques are bound to the availability of training datasets, oftentimes restricted by intellectual property (IP) liabilities. As a result, a wealth of proprietary experimental screening results remains inaccessible to researchers and impossible to share without compromising the IP of the companies.

The current project offers a solution to this problem. We propose that **sensitive experimental results can be shared securely in the form of AI/ML models**, which retain the essential properties of the dataset but do not display the identity of the screened compounds. Sharing encrypted AI/ML tools instead of datasets may enable new forms of collaboration between pharma, biotech and academia, and offers a new means of contribution to Open Science.

## Fully homomorphic encryption for AI/ML models

Fully homomorphic encryption (FHE) allows **computation on encrypted data** without leaking any information about the encrypted data. More succinctly:

$$
FHE(a·b) = FHE(a)·FHE(b)
$$

The result of the computation can only be decrypted by the party that holds the decryption key.

### The current state of FHE

Fully homomorphic encryption is still a nascent area of research in the field of cryptography compared to other established cryptographic techniques. It was theorised back in the 70s but the first practical breakthrough happened in 2009 with the seminal thesis of Craig Gentry on this topic. Since then, various new schemes have been proposed to improve usability and performance. Some of the major schemes used in practice are discussed below.

| Scheme | Description                                                                                                                                                                                                            |
| ------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| BFV    | Brakerski/Fan-Vercauteren scheme (BFV) is a Ring-Learning with Errors (RLWE)-based cryptosystem. It supports arithmetic over integers.                                                                                 |
| BGV    | Brakerski-Gentry-Vaikuntanathan (BGV) scheme is also a Ring-Learning with Errors (RLWE)-based cryptosystem and offers capabilities similar to the BFV scheme.                                                          |
| CKKS   | Cheon-Kim-Kim-Song (CKKS) scheme supports addition and multiplication over real numbers but yields only approximate results. It is best suited for application in machine learning.                                    |
| TFHE   | THFE is a fast fully homomorphic encryption scheme over the torus. It is based on GSW and its ring variants. It supports arithmetic over integers. More information can be found [here](https://tfhe.github.io/tfhe/). |

### **Application of FHE in AI/ML**

CryptoNets (2016) was the first paper to show that FHE can be successfully used to encrypt a machine learning model to perform encrypted evaluation and training. The model was hand-tuned and used Microsoft's SEAL library to implement FHE functions.

However, the adoption of FHE in AI/ML applications is still low despite its enormous potential in enabling privacy-preserving ML-as-a-Service systems with strong theoretical security guarantees. FHE still suffers from performance issues which makes it practically infeasible to use with large machine learning models. We'll discuss the challenges further in the following sections.

#### **The selection of encryption parameters is not trivial**

Selecting encryption parameters depends on the computation being performed. So, it takes some trial and error in each case. There are some projects trying to solve this problem by using a compiler approach. Look [here](https://github.com/microsoft/EVA) for more details.

#### **The time complexity of computation scales poorly with input size**

The current generation of FHE libraries suffers from severe performance issues. As the input size increases, the evaluation time quickly becomes infeasibly large. This limits the size of input matrices to an ML model.

#### **Poor integration of current FHE libraries with popular ML frameworks**

FHE libraries are not well integrated with the rest of the machine learning ecosystem. For example, TenSeal tensors are not interoperable with Pytorch tensors.

#### **Poor support for hardware accelerator backends in FHE libraries to speed up the computation**

None of the major FHE libraries implements a CUDA backend for computation. So, GPUs cannot be used to speed up computations.

#### **Poor community support**

FHE community is still small which results in poor documentation and limited worked-out examples.

### Our work

We have created a Python library ([ChemXor](https://github.com/ersilia-os/chemxor)). It provides a set of pre-tuned model architectures for evaluating FHE(Fully homomorphic encryption) encrypted inputs. These models can be trained as normal Pytorch models. It also provides convenient functions to quickly query and host these models as a service with strong privacy guarantees for the end-user. It is built on top of TenSEAL and Pytorch.

**Encryption**

The encryption context for input data is based on the third-party TenSeal library. TenSeal is currently using the CKKS encryption schema but it can be adapted to incorporate other encryption schemas as well. Computational performance is extremely sensitive to CKKS parameters. For the built-in model architectures available in ChemXor, we already provide manually tuned CKKS parameters. As a result, ChemXor provides a straightforward API to perform this otherwise laborious FHE step.

**Partitioned models**

FHE inputs also suffer from fixed multiplication depth. After a certain number of multiplication operations, the noise in the input grows too large. This limits the number of layers that a neural network can have. To overcome this problem, ChemXor encrypted models are partitioned. After a certain number of multiplications, the output is sent back to the user. The user decrypts the output, recovers the plain text and encrypts it again to send back to the model to continue execution. ChemXor provides functions to do all of this automatically.



#### **Getting started with ChemXor**

ChemXor is available on PyPi and can be installed using pip.

```bash
pip install chemxor
```

#### Tutorials

[Model Training](https://github.com/ersilia-os/chemxor/blob/main/notebooks/train\_olindanet\_models.ipynb)

[Model Evaluation and Serving](https://github.com/ersilia-os/chemxor/blob/main/notebooks/fhe\_olindanet\_models.ipynb)

#### Model selection and training

At the moment, one can choose from 3 pre-tuned models.

* `OlindaNetZero` : Slimmest model with one convolution and 3 linear layers
* `OlindaNet`: Model with two convolutions and 4 linear layers
* `OlindaOneNet`: Model with four convolutions and 4 linear layers

These models accept a `32 x 32` input and can be configured to produce a signle or multiple outputs.&#x20;

```python
from chemxor.models import OlindaNetZero, OlindaNetOne, OlindaNet

# model for regression
model = OlindaNetZero(output = 1)
```

The model is a normal Pytorch Lightning module which is compatible with Pytorch `NN` module.

#### Dataset Preparation

ChemXor provides two generic Pytorch Lightning Datamodules (Regression, Classification) that can be used to train and evaluate the models. These Datamodules expects raw data as CSV files with two columns (target, SMILES).&#x20;

```python
from chemxor.data import OlindaCDataModule, OlindaRDataModule

dm_regression = OlindaRDataModule(csv_path="path/to/csv")

# Use the threshold value to automatically create categorical 
# classes from the target column of the CSV
dm_classification = OlindaCDataModule(csv_path="path/to/csv", threshold=[0.5])

```

The DataModules will take care of converting the `smiles` input to `32 x 32` images. &#x20;

#### Model training

It is recommended to use a Pytorch Lightning trainer to train the models. Although a normal Pytorch training loop can also be used.&#x20;

```python
import pytorch_lightning as pl

# Save the best 3 checkpoints based on validation loss
checkpoint_callback = pl.callbacks.ModelCheckpoint(
        dirpath="path/to/save/checkpoints",
        save_top_k=3,
        monitor="VAL_Loss",
    )
trainer = pl.Trainer(callbacks=[checkpoint_callback], accelerator="auto")
trainer.fit(model=model, datamodule=data_module)
```

#### **FHE models**

After training, the models can be wrapped using their specific FHE wrappers to process FHE inputs. FHE wrappers will take care of Tenseal context parameters and keys management.

```python
from chemxor.models import OlindaNetZero, OlindaNetOne, OlindaNet
from chemxor.models import FHEOlindaNetZero, FHEOlindaNetOne, FHEOlindaNet

model = OlindaNetZero(output = 1)
model.load("path/to/checkpoint")
fhe_model = FHEOlindaNetZero(model=model)
```

#### **FHE inputs evaluation**

The Datamodules can generate Pytorch dataloaders that produce encrypted inputs for the model.

```python
from chemxor.data import OlindaCDataModule, OlindaRDataModule

dm_regression = OlindaRDataModule(csv_path="path/to/csv")
dm_regression.setup("test")
enc_data_loader = dm_classification.enc_dataloader(context=fhe_model.context)
enc_sample = next(iter(enc_data_loader))
```

Also, the FHE models are partitioned to control multiplicative depth. So, the forward function is modified to accept a step parameter. For testing, The FHE model can be evaluated locally as follows:

```python
from chemxor.utils import process_fhe_input

output = enc_sample[0]
for step in fhe_model.steps:
    output = fhe_model(output, step)
    dec_out = output.decrypt()
    output = process_fhe_input(
                    dec_out,
                    fhe_model.pre_process[step],
                    fhe_model.context
                )

# final decryted output
decrypted_output = output.decrypt()
```

This process can automated using a utility function provided by ChemXor

```python
from chemxor.utils import evaluate_fhe_model

decrypted_output = evaluate_fhe_model(fhe_model, enc_sample[0])
```

#### Serve models

FHE Models can be served in the form of a Flask app as follows:

```python
from chemxor.service import PartitionNetServer

fhe_model_server = PartitionNetServer(fhe_model)

if __name__ == "__main__":
    fhe_model_server.run()
```

ChemXor's Pre defined Models can also be served using the CLI

```bash
chemxor serve olida|olinda_zero|olinda_one 
```

#### Query models

We can then query models with this simple command:

```bash
chemxor query [MODEL_URL] [INPUT_SMILES_STRING]
```

### **Future work**

It might be possible to offload the encrypted model evaluation to the client with the help of re-encryption proxy schemes. It will eliminate the need for hosting models.
