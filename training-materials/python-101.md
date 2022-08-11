# Python 101

The goal of this course is to give an overview of the fundamentals of the Python programming language and encourage students to explore the possibilities of applying basic computational pipelines in their daily work

**Objectives:**

* To teach biomedical scientists the basics of Python&#x20;
* To introduce Google Colab as a cloud-based solution for coding&#x20;
* To show how Python and its libraries are used in data science and machine learning

### Content of Python 101 course

**Session 1: Python programming language**&#x20;

The basics of Python as a programming language will be covered in this session. Topics include, data types (strings, numbers, booleans, lists, sets, dictionaries), variables, string formatting, conditionals, loops and functions

**Session 2: Pandas**&#x20;

Pandas is a very convenient package specifically built to work with data formats such as .xlxs and .csv. In this session we dive into Pandas and its functions that will enable us to work with large datasets.

**Session 3: Matplotlib**&#x20;

Matplotlib is a Python library to represent information from data in the form of charts, graphs and plots. In this session, we will go over some of the most widely used plots.

**Session 4: Supervised Machine Learning**&#x20;

Supervised machine learning uses well-labeled data to train an algorithm (i.e the output is known). The algorithm learns if it is doing right by comparing the predicted output vs the real output. This session covers the steps to train and evaluate a model using Scikitlearn. Scikitlearn is a Python library that consists of classification, regression and clustering algorithms for machine learning.

### How to import Python Packages in Colab&#x20;

Python has libraries that offer standard solutions to problems that programmers encounter regularly. These libraries are broad and offer a wide range of facilities to make programming easier and faster for the programmer. There are libraries for mathematics, visualisations, data analysis and manipulation, machine learning and many more.

The Python packages that will be used in the courses include, Numpy, Pandas, Matplotlib and Seaborn. To use these packages, you have to install them first, and then import them into your runtime on Colab

Once you install the packages using pip (Python’s package installer), you do not need to run the command again.

**Numpy**&#x20;

Numpy is a library that offers useful mathematical functions and computing tools with an easy to understand syntax for programmers from any background or with any experience level to benefit from. It is already available in Colab and doesn’t need to be installed.

```python
import numpy as np
```

**Pandas**&#x20;

Pandas is a powerful and user friendly library used to study, explore and exercise control over data. It doesn’t come pre-installed in Colab so it needs to be installed first.

```python
!pip install pandas 

#you have to import the library to use it
import pandas as pd  #pd is a common abbreviation for pandas
```

**Matplotlib**&#x20;

Matplotlib is a comprehensive library used to create interactive visualisations, create publication quality plots, export to many file formats all in Python. It has a module called Pyplot that makes plotting a graph easy and fun. Matplotlib has to be installed as well before you import matplotlib and the pyplot module.

```python
!pip install matplotlib

import matplotlib.pyplot as plt  #plt is a common alias for matplotlib.pyplot
```

**Seaborn**&#x20;

Seaborn is a package built on Matplotlib that produces nicer plots. Its structure is based off Matplotlib, but you should now call the plot functions off Seaborn (sns) instead of Matplotlib (plt). You have to install the Seaborn package as well.

```python
!pip install seaborn

import seaborn as sns
```

### **How to access course materials**&#x20;

The materials for this course are located in a publicly accessible folder on Google Drive under the folders: Data Exercises Exercises\_solved Theory It is convenient to have the Courses folder in your drive.&#x20;

Below are the steps to copy the Courses folder to your drive.&#x20;

* Access the [Courses folder](https://drive.google.com/drive/folders/1m7MdRpC9Ic7ch8NlFhS1C4fVolooCR4k) on Google drive. It will appear in your "Shared with me" folder automatically.&#x20;
* In the Courses folder, right click on the Python 101 folder.&#x20;
* In the menu that opens, click on Download to download the folder to your local machine&#x20;
* Once the folder has been downloaded, go to “My Drive”&#x20;
* Click on the icon that reads “New” and click on “File Upload”&#x20;
* In the dialog box that opens, navigate to the just downloaded file’s location and click on the folder. Then click on “Upload”.&#x20;

The folder should upload into your drive and, you should have the Python 101 course folder with its files in your drive.

### Summary of the sessions&#x20;

**Session 1: Python Programming Language**&#x20;

* Comments: to add text inside a code cell, you must comment it using a hash (#) for a single line or three quotes ("""...""") for docstrings. Comments will not be run as code, and they are really helpful for following the code.
*   Data types:

    *   String: text is passed as a string (str) in python, always inside single ('...') or double quotes ("...") Numbers

        * Integers (int): whole numbers (positive or negative)
        * Floating Point Numbers (float): real numbers with decimal positions


    * Booleans: a boolean (bool) value is either a True or False&#x20;
    * Lists: a list is an ordered collection of items which can include different data types inside. Each item in a list is recognized by its position. Keep in mind that lists start at position 0 in python, not 1. A list is defined by square brackets \[1,2,2,3]&#x20;
    * Sets: A set is an unordered collection of unique items. It is defined by curly brackets {1,2,3}&#x20;
    * Dictionaries: A dictionary is an unordered conatiner of Key:Value pairs. Keys must be unique, and can map to one or more values. A dictionary is defined by curly brackets {Key:Value}, for example {"Name":"Gemma", "Age":30}

    **Lists and dictionaries are mutable, i.e the elements inside can change, but sets are immutable.**
* Variables: a variable stores data values, and its name must be UNIQUE. Variables are simply assigned by the “=” operator, and can be reassigned at any moment.&#x20;
* Conditionals: The usual conditions in mathematics like <, >, ==, <=, >= are used in Python. “If… else” statements are also used as conditionals in Python. The way conditions work is, if a condition is true, a block of code should be executed.&#x20;
* Loops: loops are used to go over a piece of code multiple times. There are for loops and while loops.&#x20;
* List comprehensions: list comprehensions are shorter lines of code used to create lists based on an existing list.&#x20;
* Functions: writing functions is a way to create reusable code. After writing the code, you only need to ‘call’ it when you want to use it and not write it again.

**Session 2: Pandas**&#x20;

* Series: the basic object of a pandas dataframe is a Series (a column of indexed elements). The Series method is called off pandas (pd) by using pd.Series&#x20;
* Dataframes: a dataframe is a collection of Series (columns) identified by an index.&#x20;
* Missing values: often times, data contains missing values. Missing values are converted to None or NaN by Pandas. These must be handled with before proceeding to further data analysis. The two most common options for dealing with missing values are:&#x20;
  * Eliminating rows or columns with missing values&#x20;
  * Imputing the null values
* Working with large dataframes: to work with large dataframes, you have to import the dataset.&#x20;
* Exporting data: Data can be exported as csv, excel (xlsx) or html files. SQL integration is also possible but it depends on the SQL API

**Session 3: Matplotlib & Seaborn**&#x20;

* Plots and subplots: most of the plots are under Matplotlib’s pyplot module.and the pyplot module is imported using its popular alias ‘plt’.
* Categorical plots: these plots are used to visualise categorical data (data that can be classified into groups). They include boxplots, barplots and violin plots&#x20;
* Distribution plots: these plots are used to display how data is distributed or spread. They include lineplots, scatterplots and histograms.&#x20;
* Seaborn: If you have time, explore the seaborn library, a package built on matplotlib that produces nicer plots. Its structure is based off Matplotlib, but you should now call the plot functions off Seaborn (with the alias sns) instead of Matplotlib (plt)
* Geographical maps: this is a very specific function to plot data in a geographical map. You can read more about [Choropleth maps](https://plotly.com/python/choropleth-maps/#reference) and the plotly library in its documentation.

**Session 4: Supervised Machine Learning.**&#x20;

* Steps to train a model:&#x20;
  1. Data collection and processing.&#x20;
  2. Division of training data in Train and Test sets.&#x20;
  3. Use the train set to train the model.&#x20;
  4. Predict an output for the test set and compare the predicted vs real results.&#x20;
  5. Improve the model until we are satisfied with the performance on the test set&#x20;
* Types of supervised ML models&#x20;
  * Classification: Classification problems are characterized by having categorical output (i.e: active, inactive), so the model tries to predict to which class the input belongs. It can include several classes, is not limited to a binary classification.&#x20;
  * Regression: Regression problems are characterized by continuous variables, where the model tries to predict the exact value of the input (i.e: IC50 of a specific compound)
* Evaluation of supervised ML models&#x20;
  * Classification metrics: We obtain the following data to evaluate the model:&#x20;
    * Y\_pred: predictions on the test set&#x20;
    * Y\_real: real outcome of the test set&#x20;
  * Accuracy: number of correct predictions divided by the total number of predictions (TP/(lenY)). For example, if we have predicted correctly 5 out of 10 data points --> Accuracy = 50%
  * Precision: identification of only real positives (with a 100% precision, a model does not classify any negatives as positive) --> TP/(TP+FP)&#x20;
  * Recall: identification of all positives (with a 100% recall, no positive is classified as negative, but some negatives might be included in the positives) --> TP / (TP+FN)&#x20;
  * Confusion matrix: plots the real vs the predicted values in a table, to easily obtain the FP, TP, TN, FN values.
* Regression metrics: In a regression task, we obtain as error the difference between the predicted value and the real value (i.e: IC50real=0.1, IC50pred = 0.5 --> error of 0.4).&#x20;
  * Mean Absolute Error: mean of the absolute values of errors.&#x20;
  * Mean Squared Error: mean of the square error. By squaring, larger errors are contributing more and therefore the model punishes them.&#x20;
  * Root Mean Squared Error: root of the mean of square error to simplify interpretation (by using MSE, we also square the units which makes them difficult to interpret).&#x20;
  * R-square: coefficient of determination, the amount of variance explained by the model (from 0 to 1, the closer to 1, the better our model is)
* Sklearn Package: Sklearn package is a popular Python package for machine learning. It contains algorithms for machine learning techniques like classification, regression and clustering.
