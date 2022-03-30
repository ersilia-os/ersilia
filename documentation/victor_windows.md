
INSTALLATION STEPS FOR ERSILIA USING WINDOWS

The Ersilia model cannot be installed on windows. It can only be installed in a linux based operating system. One of the best ways to do this is to install Ubuntu on windows with the aid of a virtual box. They steps below shows how to install Ubuntu on windows which will enable us install the Ersilia model.

Step 1:
Download and install the [Oracle Virtualbox](https://www.virtualbox.org/wiki/Downloads).

Step 2:
Download a Linux distribution for example [Ubuntu](https://ubuntu.com/download/desktop) to enable you use the Linux operating system in Windows. Don’t install yet.

Step 3:
Start the virtualbox and click on the ‘New’ symbol.
 







Step 4:
Give the name field Ubuntu and click on Next
 

Step 5:
Allocate RAM to the OS. I recommend 2gig RAM.  Your PC should have a RAM size  of atleast 6 gig so you could allocate 2gig to the OS.

 

Step 5:
Create a virtual disk. This will serve as the hard disk of the virtual Linux system. It is where the virtual system will store its files.

 










Step 6:
Choose a hard disk file type. The VDI type is given as a default.
 











Step 7:
Choose the method for creating the hard disk. The 'dynamically allocated' is given as a default.
 








Step 8:
Choose the hard disk size you would like to allocate to the OS. I recommend allocating a disk size of at least 20gig.
 









Step 9:
Double tap the Ubuntu to boot it. You can also boot it by tapping it once and clicking start (The green arrow at the top right).
 









Step 10:
Click on Install Ubuntu .
 















Step 11
Select the language of your keyboard and click continue
 















Step 12
Select  normal installation and click continue
 














Step 13:
Select ‘Erase disk and install Ubuntu’. Don’t worry it won’t delete anything on your windows operating system. You are using a virtual disk space of 20 gig that we created in the previous steps. It won’t impact the real operating system.
 













Step 14:
Click continue
 



Step 15:
Select your country and click continue.
 





Step 16:
Enter your details. Try to use a password that you can remember. Afterwards click continue.
 


You are almost done. It will take about 5 to 10 minutes to complete the installation.
 




Once the installation finishes restart the virtual system.
 


Click on your profile to enter your details and open the Ubuntu OS. There you have it, you now have Ubuntu running on your windows.

 




Step 17:      INSTALL DEPENDENCIES

Dependencies are software packages or programs that’s essential for a particular code or command to work. You will need to install the following dependencies in your Ubuntu OS so that our Ersilia model could work.

**Anaconda**

Install [Anaconda for Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) to be able to work in a conda environment to install Ersilia. 

**Python**

You will need [Python](https://docs.python-guide.org/starting/install3/linux/) to be able to install Ersilia from Python Package Index (PyPI)

**Git**

Install the version control tool [Git](https://git-scm.com/download/linux) to be able to clone and install Ersilia if you decide not to install it as a package from Python Package Index. 

**Docker**

Install [Docker]( https://runnable.com/docker/install-docker-on-linux) on Ubuntu. Docker helps all AI/ML programs to run efficiently.

### Installation of Ersilia

After you are done installing the dependencies, proceed to install Ersilia. Navigate to the Ubuntu terminal and run the following commands:
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
Either way works.

3. To test if Ersilia is properly installed, type:
```
ersilia --help
```





