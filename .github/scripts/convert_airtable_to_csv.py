import os
import csv
import re
from airtable import Airtable

current_directory = os.path.dirname(os.path.abspath(__file__))

data_directory = os.path.join(current_directory, '../../ersilia/hub/content/data')
os.makedirs(data_directory, exist_ok=True)
file_path= os.path.abspath(os.path.join(data_directory, 'models.tsv'))

def convert_airtable_to_tsv(airtable_api_key, airtable_base_id,airtable_table_id,file_path):
 
    base = Airtable(airtable_base_id, airtable_table_id, airtable_api_key)
    records = base.get_all()
    with open(file_path, 'w') as f:
        print("Writing to {0}".format(file_path))
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(records[0]["fields"].keys())
        for record in records:
            clean_values= []
            for value in record["fields"].values():
                clean_value= re.sub(r'[\r\n]', '', str(value))
                clean_values.append(clean_value)
            writer.writerow(record["fields"].values())
    
if __name__ == "__main__":

    print("Getting environmental variables")
    airtable_api_key = os.environ.get('AIRTABLE_API_KEY')
    airtable_base_id = os.environ.get('AIRTABLE_BASE_ID')
    airtable_table_id = os.environ.get('AIRTABLE_TABLE_NAME')

    print("Converting AirTable base to CSV file")
    convert_airtable_to_tsv(airtable_api_key, airtable_base_id,airtable_table_id, file_path)
