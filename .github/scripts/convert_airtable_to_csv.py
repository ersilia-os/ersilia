import os
from airtable import Airtable
import csv

current_directory = os.path.dirname(os.path.abspath(__file__))

data_directory = os.path.join(current_directory, '../../ersilia/hub/content/data')
os.makedirs(data_directory, exist_ok=True)
file_path= os.path.abspath(os.path.join(data_directory, 'models.csv'))

def convert_airtable_to_csv(airtable_api_key, airtable_base_id,airtable_table_id,file_path):
 
    base = Airtable(airtable_base_id, airtable_table_id, airtable_api_key)
    records = base.get_all()
    with open(file_path, 'w') as f:
        print("Writing to {0}".format(file_path))
        writer = csv.writer(f)
        writer.writerow(records[0]["fields"].keys())
        for record in records:
            writer.writerow(record["fields"].values())
    
if __name__ == "__main__":

    print("Getting environmental variables")
    airtable_api_key = os.environ.get('AIRTABLE_API_KEY')
    airtable_base_id = os.environ.get('AIRTABLE_BASE_ID')
    airtable_table_id = os.environ.get('AIRTABLE_TABLE_NAME')

    print("Converting AirTable base to CSV file")
    convert_airtable_to_csv(airtable_api_key, airtable_base_id,airtable_table_id, file_path)
