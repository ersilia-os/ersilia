---
description: >-
  A few interesting learning resources related to drug discovery, machine
  learning, and global health
---

# External resources

## Join our Slack community

There is an enormous amount of resources to learn about **drug discovery**, **artificial intellicence & machine learning**, and **global health**.

We continuously screen the literature, social media and the web in search of relevant content related these areas of research. Feel free to join our Slack community (message us at [hello@ersilia.io](mailto:hello@ersilia.io) and we will send you an invitation link). In particular, you may be interested in the **#literature** and **#technology** and **#media** channels.

## Drug discovery materials

### Wellcome Centre Anti-Infectives Research (WCAIR)

The WCAIR in Dundee offers a set of free [online courses](https://wcair.dundee.ac.uk/training/training-resources/) to support drug discovery research. We highly recommend the following:

* [Medicinal Chemistry](https://wcair.dundee.ac.uk/training/training-resources/medchem/)
* [The Basics of DMPK](https://wcair.dundee.ac.uk/training/training-resources/basicsofdmpk/)
* [In Conversation](https://wcair.dundee.ac.uk/training/training-resources/in-conversation-with/)

### Literature reviews

To get started with AI/ML applied to drug discovery, you can read these two review articles published in [_Nature Reviews Drug Discovery_](https://www.nature.com/nrd/):

* [Schneider et al. (2020)](https://www.nature.com/articles/s41573-019-0050-3). "Rethinking drug design in the artificial intelligence era"
* [Vamathevan et al. (2019)](https://www.nature.com/articles/s41573-019-0024-5). "Applications of machine learning in drug discovery and development"

## Tools to get started in programming

### **Python**&#x20;

Python is a popular beginner-friendly general-purpose programming language. It can be used for web development, automation, data analysis and machine learning. With an easy-to-read syntax, Python has easily become popular amongst scientists, accountants and programmers.

### **Google Colab**&#x20;

Google Colaboratory (Colab) is a Jupyter notebook that allows users to write and execute Python code for free in Google cloud. At Ersilia, we develop our lessons in Colab because it does not require any installation on the student's local server and it uses the computational capacity of Google cloud, removing any requirements for the student’s computers.

Colab's free plan allows for 12 hours of runtime and 16GB of storage. These features are more than sufficient for the purposes of training and education. In addition, Colab can allocate GPU computing time which helps in machine learning by processing large amounts of data in a short time.

#### **How to access Colab**

The only requisite to access Colab is to have a Google Account. If you don’t have an account in Google, you can [sign up](https://accounts.google.com/) to create one. Next, head to [Colab](https://colab.research.google.com/) and install the extension for your google account (free and only required the first time that Colab is used). Colab notebooks are stored on your Google Drive. To create a new notebook, go to the desired folder in your Drive, right click and select More > Cola&#x62;**.**

#### **Colab Sessions**&#x20;

**Running Colab Notebooks**&#x20;

Colab is an online Jupyter Notebook. Jupyter notebooks are the most widely interactive platform for computing. They are composed of individual cells, each one containing a piece of code that can be executed independently from the others. In Colab, you can click on the left of each cell the “play” button to execute the code in this cell.

**Importing files to Colab**&#x20;

Many times we need to import files, such as an excel table, into our Colab notebook. There are several ways to do so, here we explain in detail two of them:

**By mounting Google Drive**

You can access files in your Google Drive by mounting your drive in Colab’s runtime virtual machine

* Import drive and use the .mount function

```python
from google.colab import drive
drive.mount('/content/drive')
```

When you run the above command, an authentication page will come up seeking your permission to allow Colab access files in your drive. Once you have approved, you should get the output that reads, “Mounted at /content/drive”

* Once the drive is mounted, the folder icon on the left panel of the Colab page will open your drive hierarchy. Locate the file of interest in drive and copy the file path (right click → copy path).
* Load the file in Colab using the copied path. For example, to import a .csv file into a Pandas Dataframe, you would use the following code:Use the code below to read the file (assuming it is an Excel file in .csv format) using read\_csv() and the copied file path pasted in the parentheses

```python
import pandas as pd
df = pd.read_csv("/content/drive/MyDrive/data/excel_file.csv")
```

If the file is imported successfully, there will be no error message.

**By importing files from Github**

* Copy the URL of the .csv file hosted on Github. The URL should link to the file in raw format.&#x20;
* Paste the link and import the file.

```python
url = 'copied_csv_file_in_raw_form_Github_link'
df = pd.read_csv(url)
```

If the file import is successful, there will be no error message.

