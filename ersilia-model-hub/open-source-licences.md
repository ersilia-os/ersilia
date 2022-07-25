---
description: >-
  How to deal with the different Open Source Licences when incorporating code
  developed by third parties in the Ersilia Model Hub.
---

# Open Source Licences

In the open-source community, while open-source software is available to all, there is still the need for licences to guide users on the rights and permissions in the use of the software. At Ersilia, we release our code using a GPLv3 Licence.

This article is a guide to licence third-party code to be incorporated into the Ersilia Model Hub.

### **Terms**

* Licence: grant by the holder of a copyright or patent to another of any of the rights embodied in the copyright or patent short of an assignment of all rights ([Merriam Webster](https://www.merriam-webster.com/dictionary/licence))
* Open Source (OS) Licence: type of licence for computer software and other products that allows the source code, blueprint, or design to be  freely used, modified, and shared ([Open Source Initiative](https://opensource.org/licenses))
* Source Code: version of software as it is originally written (i.e., typed into a computer) by a human in plain text (i.e., human-readable alphanumeric characters) ([Linux Information Project](http://www.linfo.org/source\_code.html)).
* Copyright: exclusive legal right to reproduce, publish, sell, or distribute the matter and form of something ([Merriam Webster](https://www.merriam-webster.com/dictionary/copyright)).

### **Types of Open Source Licences**

OS licences fall under two types: Permissive licence and Copyleft licence.

* Permissive licences give users fewer restrictions when using source code. Users can take the permissive-licensed software, make it their own through changes or additions, and distribute that modified program with only a handful of conditions.[\[1\]](https://fossa.com/blog/what-do-open-source-licenses-even-mean/) Some popular examples of permissive licences include MIT, Apache, BSD etc.
* Copyleft licences are more restrictive than permissive licences. Generally, they require that any derivative work of the copyleft-licensed software be released under the same licence as the original.[\[2\]](https://fossa.com/blog/what-do-open-source-licenses-even-mean/) Some examples of Copyleft licences are AGPL, GPL, LGPL, Mozilla etc.

### **Permissive Licences**

* **MIT licence** [\[3\]](https://opensource.org/licenses/MIT) allows use of the source code without restrictions and limitations, even for commercial purposes, on the condition that the original copyright and licence notice are included in all copies and all substantial copies of the source code. The MIT licence is compatible with most licences like the GPLs. _When incorporating code MIT-licenced, please keep the original copyright and licence notice on the folder where the source code is and reference it in the README file._
* **Apache licence** [\[4\]](https://opensource.org/licenses/Apache-2.0) allows use of the source code without restrictions as long as the original copyright and licence notices are included in the copies, modified files carry the notice that they have been modified and include a NOTICE file, if there is one, in your copy. The Apache licence is compatible with most licences but it is not compatible with GPLv2. _When incorporating code MIT-licenced, please keep the original copyright and licence notice on the folder where the source code is and reference it in the README file. If modifications have been done to the original code, create a NOTICE file explaining such modifications._
* **BSD 3-Clause licence** [\[5\]](https://opensource.org/licenses/BSD-3-Clause) is short for Berkeley Software Distribution 3-Clause licence. It allows users to use source code with its licence without restrictions and limitations as long as the original licence is kept and the source code’s author and contributors are not used to promote the distributed derivative works. _When incorporating code BSD-3-licenced, please keep the original copyright and licence notice on the folder where the source code is and reference it in the README file._

### **Copyleft Licenses**

* **GNU Public Licence (GPL) version 3** [\[6\]](https://opensource.org/licenses/GPL-3.0) requires that derivative works of the source code with its licence keep the licence, modified files should be marked as changed and the software should not be privatised or become proprietary. This is the preferred licence for Ersilia's software, since it ensures code will remain open to all. _If incorporating code already licensed under GPLv3, simply add the LICENSE.md file on the root of the repository and create  NOTICE.md file for any modifications._
* **GNU Affero General Public Licence (AGPL)** [\[7\]](https://opensource.org/licenses/AGPL-3.0) lets users distribute source code and even charge for the software. It requires that derivative works maintain the AGPL and modified files carry the notice of being modified. _When incorporating code AGPL-licenced, please keep the original copyright and licence notice on the folder where the source code is and reference it in the README file. Also create a NOTICE.md file in the relevant folder for any modifications._
* **GNU Lesser General Public Licence (LGPL)** [\[8\]](https://opensource.org/licenses/LGPL-3.0) permits users to copy and distribute verbatim copies of the licensed source code but users should not change the source code. _When incorporating code LGPL-licenced, please keep the original copyright and licence notice on the folder where the source code is and reference it in the README file._
* **Mozilla Public Licence** [\[9\]](https://opensource.org/licenses/MPL-2.0) requires that derivative works should maintain the licence. It permits distribution, modification, and even using the software for commercial purposes, but the original code and its modifications must maintain the MLP licence. _When incorporating MPL-licenced code, maintain the licence notice for the specific folders where source code is._

If the licence you are dealing with is not described here please refer to the sources used to write these guidelines and contact the Ersilia Team ([hello@ersilia.io](mailto:hello@ersilia.io)):

1. Mahak Bandi. All About Open Source Licenses. _FOSSA_ [https://fossa.com/blog/what-do-open-source-licenses-even-mean/](https://fossa.com/blog/what-do-open-source-licenses-even-mean/)
2. Mahak Bandi. All About Open Source Licenses. _FOSSA_ [https://fossa.com/blog/what-do-open-source-licenses-even-mean/](https://fossa.com/blog/what-do-open-source-licenses-even-mean/)
3. [The MIT License](https://opensource.org/licenses/MIT). _Open Source Initiative._ Accessed June 2022
4. Apache License, Version 2.0. _Open Source Initiative_ [https://opensource.org/licenses/Apache-2.0](https://opensource.org/licenses/Apache-2.0)
5. The 3-Clause BSD License _Open Source Initiative_ [https://opensource.org/licenses/BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause)
6. GNU General Public License _Open Source Initiative_ [https://opensource.org/licenses/GPL-3.0](https://opensource.org/licenses/GPL-3.0)
7. GNU Affero General Public License _Open Source Initiative_ [https://opensource.org/licenses/AGPL-3.0](https://opensource.org/licenses/AGPL-3.0)
8. GNU Lesser General Public License _Open Source Initiative_ [https://opensource.org/licenses/LGPL-3.0](https://opensource.org/licenses/LGPL-3.0)
9. Mozilla Public License _Open Source Initiative_ [https://opensource.org/licenses/MPL-2.0](https://opensource.org/licenses/MPL-2.0)
