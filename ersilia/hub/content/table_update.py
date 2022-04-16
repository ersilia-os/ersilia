from .catalog import CatalogTable




class table(object):
    def __init__(self, catalog, table_no):
        self.catalog = catalog
        self.table_no = table_no
        self.page_size = 5
    
    def initialise(self):
        table_no = 0
        data = self.catalog.data
        page_size = self.page_size
        R = data[table_no * page_size : table_no * page_size + page_size]
        return CatalogTable(data=R, columns=self.catalog.columns)
    
    def next_table(self):
        table_no = self.table_no
        data = self.catalog.data
        page_size = self.page_size
        R = data[table_no * page_size : table_no * page_size + page_size]
        if table_no * page_size + page_size > len(data):
            print("Last Table! ersilia catalog --previous to go back")

        return CatalogTable(data = R, columns = self.catalog.columns)

    def prev_table(self):
        table_no = self.table_no
        data = self.catalog.data
        page_size = self.page_size

        table_no = table_no - 2

        if table_no * page_size + page_size > len(data):
            table_no = table_no - 2
        if table_no < 1:
            table_no = 0
            print('Table 1. ersilia catalog --next to go to next table')
        R = data[table_no * page_size : table_no * page_size + page_size]
        return CatalogTable(data = R, columns = self.catalog.columns)

