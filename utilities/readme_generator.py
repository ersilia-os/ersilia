from __future__ import print_function
import sys
import json

class cfile(file):
    def __init__(self, name, mode = 'r'):
        self = file.__init__(self, name, mode)

    def h1(self, string):
        self.writelines("# " + string + '\n')
        return None

    def h2(self, string):
        self.writelines("## " + string + '\n')
        return None

    def writeLine(self, string):
        self.writelines(string + '\n')
        return None

    def beginTable(self):
        self.writeLine("| | |")
        self.writeLine("|-|-|")
        return None

    def tableRow(self, parentKey, key, value, schema = True):
        if schema:
            if key in parentKey.keys():
                self.writelines("| " + key +  " | " + value + " | " + '\n')
        else:
            self.writelines("| " + key +  " | " + value + " | " + '\n')
        return None


def generateReadMe(data):
    fid = cfile('README.md', 'w')
    fid.h1(data["meta"]["name"])
    fid.writeLine("This repository hosts the contributor source files for the " +  data["meta"]["name"] + " model. ModelHub" +
                   " integrates these files into an engine and controlled runtime environment. A unified API allows" +
                   " for out-of-the-box reproducible implementations of published models. For more information, please" +
                  " visit [www.modelhub.ai](http://modelhub.ai/) or contact us [info@modelhub.ai](mailto:info@modelhub.ai).")
    #
    fid.h2("meta")
    meta = data["meta"]
    fid.beginTable()
    fid.tableRow(data, "id", data["id"])
    fid.tableRow(meta, "application_area", meta["application_area"])
    fid.tableRow(meta, "task", meta["task"])
    fid.tableRow(meta, "task_extended", meta["task_extended"])
    fid.tableRow(meta, "data_type", meta["data_type"])
    fid.tableRow(meta, "data_source", meta["data_source"])
    #
    if "publication" in data.keys():
        fid.h2("publication")
        pub = data["publication"]
        fid.beginTable()
        fid.tableRow(pub, "title", pub["title"])
        fid.tableRow(pub, "source", pub["source"])
        fid.tableRow(pub, "url", pub["url"])
        fid.tableRow(pub, "year", str(pub["year"]))
        fid.tableRow(pub, "authors", pub["authors"])
        fid.tableRow(pub, "abstract", pub["abstract"])
        fid.tableRow(pub, "google_scholar", pub["google_scholar"])
        fid.tableRow(pub, "bibtex", pub["bibtex"])
    #
    fid.h2("model")
    model = data["model"]
    fid.beginTable()
    fid.tableRow(model, "description", model["description"])
    fid.tableRow(model, "provenance", model["provenance"])
    fid.tableRow(model, "architecture", model["architecture"])
    fid.tableRow(model, "learning_type", model["learning_type"])
    fid.tableRow(model, "format", model["format"])
    fid.tableRow(model, "I/O", "model I/O can be viewed [here](contrib_src/model/config.json)", schema = False)
    fid.tableRow(model, "license", "model license can be viewed [here](contrib_src/license/model)", schema = False)

    fid.h2("run")
    fid.writeLine("To run this model and view others in the collection, view the instructions on [ModelHub](http://app.modelhub.ai/).")

    fid.h2("contribute")
    fid.writeLine("To contribute models, visit the [ModelHub docs](https://modelhub.readthedocs.io/en/latest/).")

    fid.close()
    print ("readme.md generated.")


if __name__ == "__main__":
    try:
        with open(sys.argv[1]) as f:
            data = json.load(f)
        generateReadMe(data)
    except Exception as e:
        print (e)
