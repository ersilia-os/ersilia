# Windows Installation of Ersilia 
### Prerequisites 
**Linux operating system** 

For now, Ersilia does not work natively on Windows operating systems. To use Ersilia on a Windows system, you have to install a Linux operating system. 

There are two ways you can use a Linux operating system in Windows:

- You can use [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install) to install a Linux distribution. (Ubuntu is the default for Windows Subsystem for Linux). 
You can also use [Visual studio code](https://code.visualstudio.com/docs/remote/wsl) with Windows Subsystem for Linux.

**Note: Windows Subsystem for Linux is not available on Windows 8.1 and previous versions.**

- Or you can download a Virtual Machine for Windows hosts ([Oracle Virtualbox](https://www.virtualbox.org/wiki/Downloads) or [VMWare Workstation Player](https://www.vmware.com/products/workstation-player.html)) and download a Linux distribution (for example [Ubuntu](https://ubuntu.com/download/desktop)) to enable you use the Linux operating system in Windows.

**Conda**

Install [Conda for Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) to be able to work in a conda environment to install Ersilia. 

**Python**

You will need [Python](https://docs.python-guide.org/starting/install3/linux/) to be able to install Ersilia from Python Package Index (PyPI)

**Git**

Install the version control tool [Git](https://git-scm.com/download/linux) to be able to clone and install Ersilia if you decide not to install it as a package from Python Package Index. 

### Installation of Ersilia

After you are done with the prerequisites, proceed to install Ersilia.

1. Create and activate Ersilia environment
```
conda create -n ersilia python=3.7
conda activate ersilia
```
2. Install Ersilia
```
#from pypi
pip install ersilia

#or latest from github
pip install git+https://github.com/ersilia-os/ersilia.git
```
3. To test if Ersilia is properly installed, type:
```
ersilia --help
```

