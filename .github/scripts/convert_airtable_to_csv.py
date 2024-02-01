import os
from airtable import Airtable
import csv

airtable_api_key = os.environ.get('AIRTABLE_API_KEY')
airtable_base_id = os.environ.get('AIRTABLE_BASE_ID')
airtable_table_id = os.environ.get('AIRTABLE_TABLE_NAME')

def convert_airtable_to_csv(airtable_api_key, airtable_base_id,airtable_table_id):

    base = Airtable(airtable_base_id, airtable_table_id, airtable_api_key)
    records = base.get_all()

    script_directory = os.path.dirname(os.path.abspath(__file__))
    full_output_path = os.path.join(script_directory, "../../hub/content/data/models_from_airtable.csv")
    os.makedirs(os.path.dirname(full_output_path), exist_ok=True)

    csv_file = os.path.join(full_output_path)

    with open(csv_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(records[0]["fields"].keys())
        for record in records:
            writer.writerow(record["fields"].values())
    
if __name__ == "__main__":

    airtable_api_key = os.environ.get('AIRTABLE_API_KEY')
    airtable_base_id = os.environ.get('AIRTABLE_BASE_ID')
    airtable_table_id = os.environ.get('AIRTABLE_TABLE_NAME')

    convert_airtable_to_csv(airtable_api_key, airtable_base_id,airtable_table_id)