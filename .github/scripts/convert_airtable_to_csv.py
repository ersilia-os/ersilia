import sys
import os
from airtable import Airtable
import csv


def convert_airtable_to_csv(airtable_api_key, airtable_base_id,airtable_table_id,output_path):

    base = Airtable(airtable_base_id, airtable_table_id, airtable_api_key)
    records = base.get_all()

    script_directory = os.path.dirname(os.path.abspath(__file__))
    full_output_path = os.path.join(script_directory, "../../hub/content/data/")
    os.makedirs(os.path.dirname(full_output_path), exist_ok=True)

    csv_file = os.path.join(full_output_path, output_path)

    with open(csv_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(records[0]["fields"].keys())
        for record in records:
            writer.writerow(record["fields"].values())
    
if __name__ == "__main__":
    airtable_api_key = sys.argv[1]
    airtable_base_id= sys.argv[2]
    airtable_table_id= sys.argv[3]
    output_path = sys.argv[4]

    convert_airtable_to_csv(airtable_api_key, airtable_base_id,airtable_table_id,output_path)