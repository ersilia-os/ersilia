# Fully homomorphic encryption for AI/ML models

### FHE crash course

FHE allows computation on encrypted data without leaking any information about the encrypted data. More succinctly\


FHE(a \* b) = FHE(a) \* FHE(b)



The result of the computation can only be decrypted by the party that holds the decryption key.

### The current state of FHE

Fully homomorphic encryption is still a nascent area of research in the field of cryptography compared to other established cryptographic techniques. It was theorised back in the 70s but the first practical breakthrough happened in 2009 with the seminal thesis of Craig Gentry on this topic. Since then, various new schemes have been proposed to improve usability and performance. Some of the major schemes used in practice are discussed below.



**BFV**

Brakerski/Fan-Vercauteren scheme (BFV) is a Ring-Learning with Errors (RLWE)-based cryptosystem. It supports arithmetic over integers.

**BGV**

Brakerski-Gentry-Vaikuntanathan (BGV) scheme is also a Ring-Learning with Errors (RLWE)-based cryptosystem and offers capabilities similar to the BFV scheme.

**CKKS**

Cheon-Kim-Kim-Song (CKKS) scheme supports addition and multiplication over real numbers but yields only approximate results. It is best suited for application in machine learning.

**TFHE**

THFE is a fast fully homomorphic encryption scheme over the torus. It is based on GSW and its ring variants. It supports arithmetic over integers. More information can be found [here](https://tfhe.github.io/tfhe/).\
\


**Application of FHE in AI/ML**

CryptoNets (2016) was the first paper to show that FHE can be successfully used to encrypt a machine learning model to perform encrypted evaluation and training. The model was hand-tuned and used Microsoft's SEAL library to implement FHE functions.

However, the adoption of FHE in AI/ML applications is still low despite its enormous potential in enabling privacy-preserving ML-as-a-Service systems with strong theoretical security guarantees. FHE still suffers from performance issues which makes it practically infeasible to use with large machine learning models. We'll discuss the challenges further in the following sections.

**The selection of encryption parameters is not trivial**

Selecting encryption parameters depends on the computation being performed. So, it takes some trial and error in each case. There are some projects trying to solve this problem by using a compiler approach. Look [here](https://github.com/microsoft/EVA) for more details.\


**The time complexity of computation scales poorly with input size**

The current generation of FHE libraries suffers from severe performance issues. As the input size increases, the evaluation time quickly becomes infeasibly large. This limits the size of input matrices to an ML model.\


**Poor integration of current FHE libraries with popular ML frameworks**

FHE libraries are not well integrated with the rest of the machine learning ecosystem. For example, TenSeal tensors are not interoperable with Pytorch tensors.\


**Poor support for hardware accelerator backends in FHE libraries to speed up the computation**

None of the major FHE libraries implements a CUDA backend for computation. So, GPUs cannot be used to speed up computations.\


**Poor community support**

FHE community is still small which results in poor documentation and limited worked-out examples.\


### Our Work

We have created a Python library(ChemXor) for training and evaluating PyTorch models on FHE(Fully homomorphic encryption) encrypted inputs with no code changes to the original model. It also provides convenient functions to quickly query and host these models as a service with strong privacy guarantees for the end-user. It is built on top of TenSEAL, Pytorch and ONNX.

#### **Getting started with chemxor**

Chemxor is available on PyPi and can be installed using pip.

```bash
pip install chemxor
```

**Encryption context**

We first need to create an encryption context to begin encrypting models and inputs. The TenSeal library is used to create encryption contexts. In the example below, we are using the CKKS encryption scheme.&#x20;

```python
import tenseal as ts

# Create a tenseal context
context = ts.context(
    ts.SCHEME_TYPE.CKKS,
    poly_modulus_degree=8192,
    coeff_mod_bit_sizes=[60, 40, 40, 60]
)
context.global_scale = pow(2, 40)
context.generate_galois_keys()
```

There are other encryption schemes available with their own strengths and weaknesses. Choosing parameters for encryption schemes is not trivial and requires trial and error. More detailed documentation for available encryption schemes and their parameters is available [here](https://github.com/Microsoft/SEAL#examples) and [here](https://github.com/OpenMined/TenSEAL/tree/main/tutorials). There are other projects ([EVA](https://github.com/Microsoft/EVA) compiler for CKKS) that are trying to automate the selection of parameters for specific encryption schemes. However, this is currently out of scope for ChemXor.

**Encrypted Datasets**

ChemXor provides functions to easily convert your Pytorch datasets to Encrypted datasets.

```python
from torch.utils.data import DataLoader
from chemxor.data_modules.enc_dataset import EncDataset

# Use the context that we created earlier
enc_pytorch_dataset = EncDataset(context, pytorch_dataset)

# The encrypted datasets can also be used to create dataloaders
DataLoader(enc_pytorch_dataset, batch_size=None)
```

`EncDataset` class is a wrapper that modifies that **`__getitem__`** method of the `Dataset` class from Pytorch. It encrypts the items using the provided `context` before returning the items.

ChemXor also provides `EncConvDataset` class, a variant to `EncDataset` class for inputs that undergo convolution operations.

```python
from torch.utils.data import DataLoader
from chemxor.data_modules.enc_conv_dataset import EncConvDataset

# Use the context that we created earlier
enc_pytorch_dataset = EncConvDataset(context, pytorch_dataset, kernel_size, stride)

# The encrypted datasets can also be used to create dataloaders
DataLoader(enc_osm_train, batch_size=None)
```

It uses image-to-column encoding of inputs to speed up computation. More details on this topic can be found [here](https://github.com/OpenMined/TenSEAL/blob/main/tutorials/Tutorial%204%20-%20Encrypted%20Convolution%20on%20MNIST.ipynb).

> `EncDataset` and `EncConvDataset` does not encrypt the data on disk. Items are encryted lazily on the fly as needed.

**Encrypted models**

ChemXor can automatically convert Pytorch models to models that can be evaluated on encrypted inputs. However, evaluating any arbitrary converted model on encrypted inputs can take an infeasibly long time. This is a major limitation of FHE at the moment.

```python
import tenseal as ts
from chemxor.crypt import FHECryptor
from chemxor.model.cryptic_sage import CrypticSage

# Create a tenseal context
context = ts.context(
    ts.SCHEME_TYPE.BFV,
    poly_modulus_degree=4096,
    plain_modulus=1032193
)

# Initialize the FHECryptor with tenseal context
fhe_cryptor = FHECryptor(context)

# Use any Pytorch Model
model = CrypticSage()

# Convert model using the FHECryptor
converted_model = fhe_cryptor.convert_model(model, dummy_input)

# Converted model is still a Pytorch lightning module
# So use it as usual for evaluating encrypted inputs
enc_output = converted_model(enc_input)
```

ChemXor first converts the Pytorch model to an ONNX model. This ONNX model is then used to create an equivalent function chain that can process encrypted inputs. The resulting converted model is still a Pytorch model with a modified forward function. We are still working on supporting all the operations in the ONNX spec. But, some of the operations might not be available at the time of release.

It is also possible to manually wrap an existing Pytorch model class to make it compatible with encrypted inputs. This is the recommended approach for now as the automatic conversion is not mature yet. There are several models with their encrypted wrappers in ChemXor that can be used as examples.

```python
# Pytorch lightning model
# Adapted from https://github.dev/OpenMined/TenSEAL/blob/6516f215a0171fd9ad70f60f2f9b3d0c83d0d7c4/tutorials/Tutorial%204%20-%20Encrypted%20Convolution%20on%20MNIST.ipynb
class ConvNet(pl.LightningModule):
    """Cryptic Sage."""

    def __init__(self: "ConvNet", hidden: int = 64, output: int = 10) -> None:
        """Init."""
        super().__init__()
        self.hidden = hidden
        self.output = output
        self.conv1 = nn.Conv2d(1, 4, kernel_size=7, padding=0, stride=3)
        self.fc1 = nn.Linear(256, hidden)
        self.fc2 = nn.Linear(hidden, output)

    def forward(self: "ConvNet", x: Any) -> Any:
        """Forward function.

        Args:
            x (Any): model input

        Returns:
            Any: model output
        """
        x = self.conv1(x)
        # the model uses the square activation function
        x = x * x
        # flattening while keeping the batch axis
        x = x.view(-1, 256)
        x = self.fc1(x)
        x = x * x
        x = self.fc2(x)
        return x

# Encrypted wrapper
# Adapted from https://github.dev/OpenMined/TenSEAL/blob/6516f215a0171fd9ad70f60f2f9b3d0c83d0d7c4/tutorials/Tutorial%204%20-%20Encrypted%20Convolution%20on%20MNIST.ipynb
class EncryptedConvNet(pl.LightningModule):
    """Encrypted ConvNet."""

    def __init__(self: "EncryptedConvNet", model: ConvNet) -> None:
        """Init."""
        super().__init__()

        self.conv1_weight = model.conv1.weight.data.view(
            model.conv1.out_channels,
            model.conv1.kernel_size[0],
            model.conv1.kernel_size[1],
        ).tolist()
        self.conv1_bias = model.conv1.bias.data.tolist()

        self.fc1_weight = model.fc1.weight.T.data.tolist()
        self.fc1_bias = model.fc1.bias.data.tolist()

        self.fc2_weight = model.fc2.weight.T.data.tolist()
        self.fc2_bias = model.fc2.bias.data.tolist()

    def forward(self: "EncryptedConvNet", x: Any, windows_nb: int) -> Any:
        """Forward function.

        Args:
            x (Any): model input
            windows_nb (int): window size.

        Returns:
            Any: model output
        """
        # conv layer
        enc_channels = []
        for kernel, bias in zip(self.conv1_weight, self.conv1_bias):
            y = x.conv2d_im2col(kernel, windows_nb) + bias
            enc_channels.append(y)
        # pack all channels into a single flattened vector
        enc_x = ts.CKKSVector.pack_vectors(enc_channels)
        # square activation
        enc_x.square_()
        # fc1 layer
        enc_x = enc_x.mm(self.fc1_weight) + self.fc1_bias
        # square activation
        enc_x.square_()
        # fc2 layer
        enc_x = enc_x.mm(self.fc2_weight) + self.fc2_bias
        return enc_x
```

A few things to note here:

* We converted Pytorch tensors to a list in the encrypted wrapper. This is required as Pytorch tensors are not compatible with TenSeal encrypted tensors.
* We are not using the standard ReLU activation. CKKS encryption scheme cannot evaluate non-linear piecewise functions. So, either alternative activation functions can be used or polynomial approximations of non-linear activation functions can be used.

#### Serve models

```python
from chemxor.service import create_model_server

# `create_model_server` returns a flask app
flask_app = create_model_server(model, dummy_input)

if __name__ == "__main__":
    flask_app.run()
```

#### Query models

```bash
chemxor query -i [input file path] [model url]
```

#### Distilled models

To overcome the performance limitations of FHE, we used ChemXor to create simpler distilled models from larger complex models.  accept inputs as molecules encoded as 32 x 32 images and predict the properties of these molecules. These models are hosted using AWS lambdas and are available for public use.\


**Future Work**

It might be possible to offload the encrypted model evaluation to the client with the help of Re-encryption proxy schemes. It will eliminate the need for hosting models.
