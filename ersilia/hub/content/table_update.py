import os
import tempfile

from .catalog import CatalogTable

counter_path = os.path.join(tempfile.mkdtemp(), "model_counter")
# start of script - read or initialise counter
try:
    with open(counter_path, "r") as count_in:
        counter = int(count_in.read())

except FileNotFoundError:
    counter = 0


class table(object):
    def __init__(self, catalog):
        self.catalog = catalog
        self.page_size = 15
        self.counter = counter

    def initialise(self):
        table_no = 0
        data = self.catalog.data
        page_size = self.page_size
        R = data[table_no * page_size : table_no * page_size + page_size]
        # sets counter to next table i.e. table 1
        with open(counter_path, "w") as count_in:
            count_in.write(str(1))
        return CatalogTable(data=R, columns=self.catalog.columns)

    def next_table(self):

        table_no = counter
        data = self.catalog.data
        page_size = self.page_size
        R = data[table_no * page_size : table_no * page_size + page_size]
        if table_no * page_size + page_size > len(data) < len(data) / page_size:
            print("Last Table! ersilia catalog --previous to go back")
        with open(counter_path, "w") as count_out:
            count_out.write(str(counter + 1))
        # if table finished, freezes the counter to the last table_no and returns
        if len(data) / page_size <= counter:
            with open(counter_path, "w") as count_out:
                count_out.write(str(counter))
            return "Table Ended! Please, ersilia catalog --previous to go back"
        return CatalogTable(data=R, columns=self.catalog.columns)

    def prev_table(self):
        table_no = counter
        data = self.catalog.data
        page_size = self.page_size

        table_no = table_no - 2

        if table_no == 0:
            print("Table 0. ersilia catalog --next to go to next table")
        R = data[table_no * page_size : table_no * page_size + page_size]
        with open(counter_path, "w") as count_out:
            count_out.write(str(counter - 1))
        # if already table 0, frezes the counter to 1 and return
        if counter < 1:
            with open(counter_path, "w") as count_out:
                count_out.write(str(1))
                return "Table ended! ersilia catalog --next to go next "
        return CatalogTable(data=R, columns=self.catalog.columns)
