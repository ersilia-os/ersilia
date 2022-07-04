# OPEN SOURCE LICENSES

**MODEL INCORPORATION GUIDELINES**

**Documentation guidelines**

Each model in the Ersilia Model Hub is contained within an individual repository. The README file of the repository should cover the essential information about it, using the provided [template](https://ersilia-os/eos-template).

In summary:

**The template to follow when adding a new model**

Information on a model is required to understand how it works and to get a backstory on it. Please follow the guide to know the information to provide and how to structure it. The example model used here is the Avalon model.

1. Model title: the title of the model will be the heading of the README file. The format to write a heading in Github is “# Avalon fingerprints” and it will look as so:

![](../.gitbook/assets/0)

1. Model identifiers: this is the section to include information that identifies the model. It is a subheading and will be written in the format “## Model identifiers”. The fields needed include,

* Slug: It is a one-word reference to the model
* Ersilia ID: a generated unique model identifier
* Tags: a few words that help to further identify the model.

The result should look as so:![](../.gitbook/assets/1)

1. Model description: this section includes information that describes the model. It is a subheading and has the following fields

* Input: an example is SMILES
* Output: {unit and description of output)
* Model type: (Regression or Classification)
* Training set: (number of compounds and link to the training data)
* Mode of training: (is it pretrained? that is, were the checkpoints downloaded and used to train the model? or is it retrained? that is, was it trained from scratch with updated data)

![](../.gitbook/assets/2)

1. Source code: this is the fourth subheading and contains information on the new model's source publication and source code

* Cite the source publication.
* Code: include the link to the source code
* Checkpoints: include the link to the checkpoints used if the model is a pretrained model

It will look like this:

![](../.gitbook/assets/3)

1. License and copyright notice: In this section, state the licences used which are the GPL v3 license used by Ersilia and the license used by the source code, if any exists. Use this guide to license new models in Ersilia's model hub. The result should look like this:

![](../.gitbook/assets/4)

1. History: here, state the date when the model was downloaded and incorporated into Ersilia. This is how it should look’

![](../.gitbook/assets/5)

**How to deal with the different Open Source Licences when incorporating code developed by third parties in the Ersilia Model Hub.**

In the open-source community, while open-source software is available to all, there is still the need for licences to guide users on the rights and permissions in the use of the software. At Ersilia, we licence code developed by Ersilia under GPLv3.

This article is a guide to licence third-party code to be incorporated into the Ersilia Model Hub.

**Terms:**

* Licence: according to [Merriam Webster](https://www.merriam-webster.com/dictionary/licence), is a grant by the holder of a copyright or patent to another of any of the rights embodied in the copyright or patent short of an assignment of all rights
* Open Source (OS) Licence: according to[ Wikipedia](https://en.wikipedia.org/wiki/Open-source\_license), is a type of licence for computer software and other products that allows the source code, blueprint, or design to be used, modified and/or shared under defined terms and conditions
* Source Code: The [Linux Information Project](http://www.linfo.org/source\_code.html) defines source code as the version of software as it is originally written (i.e., typed into a computer) by a human in plain text (i.e., human-readable alphanumeric characters).
* Copyright: according to [Merriam Webster](https://www.merriam-webster.com/dictionary/copyright), is the exclusive legal right to reproduce, publish, sell, or distribute the matter and form of something (such as a literary, musical, or artistic work)

**Types of Open Source Licences**

OS licences fall under two types: Permissive licence and Copyleft licence.

* Permissive licences give users fewer restrictions when using source code. Users can take the permissive-licensed software, make it their own through changes or additions, and distribute that modified program with only a handful of conditions.[\[1\]](https://fossa.com/blog/what-do-open-source-licenses-even-mean/) Some popular examples of permissive licences include MIT, Apache, BSD etc.
* Copyleft licences are more restrictive than permissive licences. Generally, they require that any derivative work of the copyleft-licensed software be released under the same licence as the original.[\[2\]](https://fossa.com/blog/what-do-open-source-licenses-even-mean/) Some examples of Copyleft licences are AGPL, GPL, LGPL, Mozilla etc.

**Summary of some famous examples of the two types of OS licences**

**Permissive Licences**

* **MIT licence** [\[3\]](https://opensource.org/licenses/MIT) lets users use source code without restrictions and limitations, even for commercial purposes, on the condition that the original copyright and licence notice are included in all copies and all substantial copies of the source code. However, the MIT licence is compatible with most licences like the GPLs.

Ersilia guidelines:

* Keep the original licence
* Add a GPLv3 licence.
* **Apache licence** [\[4\]](https://opensource.org/licenses/Apache-2.0) lets users, without restriction, use source code with its licence as long as the original copyright and licence notices are included in the copies, modified files carry the notice that they have been modified and include a NOTICE file, if there is one, in your copy. The Apache licence is compatible with most licences but it is not compatible with GPLv2.

Ersilia guidelines:

* Keep the original licence
* Indicate changes made in modified files
* Include the NOTICE file, if there is one
* Add the GPLv3 licence when incorporating the model to Ersilia Model Hub
* **BSD 3-Clause licence** [\[5\]](https://opensource.org/licenses/BSD-3-Clause) is short for Berkeley Software Distribution 3-Clause licence. It allows users to use source code with its licence without restrictions and limitations as long as the original licence is kept and the source code’s author and contributors are not used to promote the distributed derivative works.

**Copyleft Licenses**

* **GNU Public Licence (GPL) version 3** [\[6\]](https://opensource.org/licenses/GPL-3.0) requires that derivative works of the source code with its licence keep the licence, modified files should be marked as changed and the software should not be privatised or become proprietary

Ersilia guidelines:

* The GPLv3 should be maintained
* Modified files should carry the notice that they have been modified
* **GNU Affero General Public Licence (AGPL)** [\[7\]](https://opensource.org/licenses/AGPL-3.0) lets users distribute source code and even charge for the software. It requires that derivative works maintain the AGPL and modified files carry the notice of being modified.

Ersilia guidelines:

* The AGPL should be maintained
* Modified files should carry the notice that they have been modified
* **GNU Lesser General Public Licence (LGPL)** [\[8\]](https://opensource.org/licenses/LGPL-3.0) permits users to copy and distribute verbatim copies of the licensed source code but users should not change the source code.

Ersilia guidelines:

* The LGPL should be maintained
* Modified files should carry the notice that they have been modified
* **Mozilla Public Licence** [\[9\]](https://opensource.org/licenses/MPL-2.0) requires that derivative works should maintain the licence. It permits distribution, modification, and even using the software for commercial purposes.

Ersilia guidelines:

* Keep the Mozilla Public Licence

**How to license a third-party model in Ersilia**

1. Check the licence in the source code to know its permissions
2. In the model’s repository on Ersilia’s Github page, include the source code’s licence in the folder that has the source code of the model first.
3. Then include the GPL v3 licence in the main folder of the model by following these steps:
   1. In the main repository, click on create “Add file” then “Create New file” ![](../.gitbook/assets/6)
   2. Type in “LICENSE” as the name of the file

![](../.gitbook/assets/7)

*
  1. On the right side of the page, you should see “Choose a license template”. Click on it.
  2. Select “GNU General Public License v3.0”![](../.gitbook/assets/8)
  3. Click on “Review and submit”

![](../.gitbook/assets/9)

*
  1. Click on “Commit new file” at the bottom of the page that will be loaded![](../.gitbook/assets/10)

1. On the model’s README, include the License and copyright notice of both the GPL v3 licence and the source code’s licence.

If the licence you are dealing with is not described here please refer to the sources used to write these guidelines and contact the Ersilia Team ([hello@ersilia.io](mailto:hello@ersilia.io)):

1. Mahak Bandi. All About Open Source Licenses. _FOSSA_ [https://fossa.com/blog/what-do-open-source-licenses-even-mean/](https://fossa.com/blog/what-do-open-source-licenses-even-mean/)
2. Mahak Bandi. All About Open Source Licenses. _FOSSA_ [https://fossa.com/blog/what-do-open-source-licenses-even-mean/](https://fossa.com/blog/what-do-open-source-licenses-even-mean/)
3. The MIT License _Open Source Initiative_ [https://opensource.org/licenses/MIT](https://opensource.org/licenses/MIT)
4. Apache License, Version 2.0. _Open Source Initiative_ [https://opensource.org/licenses/Apache-2.0](https://opensource.org/licenses/Apache-2.0)
5. The 3-Clause BSD License _Open Source Initiative_ [https://opensource.org/licenses/BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause)
6. GNU General Public License _Open Source Initiative_ [https://opensource.org/licenses/GPL-3.0](https://opensource.org/licenses/GPL-3.0)
7. GNU Affero General Public License _Open Source Initiative_ [https://opensource.org/licenses/AGPL-3.0](https://opensource.org/licenses/AGPL-3.0)
8. GNU Lesser General Public License _Open Source Initiative_ [https://opensource.org/licenses/LGPL-3.0](https://opensource.org/licenses/LGPL-3.0)
9. Mozilla Public License _Open Source Initiative_ [https://opensource.org/licenses/MPL-2.0](https://opensource.org/licenses/MPL-2.0)
