from genericpath import exists
import json
import csv
import time
from ..cli import echo
from ..hub.delete.delete import ModelFullDeleter

# in days :
model_usage_lim   = 30 
model_cleanup_lim = 7

def seconds_to_days(s):
    return s / (24 * 3600) 

def main():
    ts_dict = {}
    if not exists('last_cleaned.json'):
        ts_dict['timestamp'] = str(time.time())
        with open('last_cleaned.json', 'w') as outfile:
            json.dump(ts_dict, outfile)
    else:
        current_ts = time.time()
        with open('last_cleaned.json') as json_file:
            ts_dict = json.load(json_file)
            ts = float(ts_dict['timestamp'])
            if (seconds_to_days(current_ts - ts))>model_cleanup_lim: 
                with open("fetched_models.txt") as infile:
                    fetched_models = dict(csv.reader(infile))
                for model_id in fetched_models: 
                    if (seconds_to_days(current_ts-float(fetched_models[model_id])))>model_usage_lim:
                        del_unused_model(model_id)

def del_unused_model(model_id):
        # model_id = ModelBase(model).model_id
        md = ModelFullDeleter()
        if md.needs_delete(model_id):
            echo("Deleting model {0}".format(model_id))
            md.delete(model_id)
            echo(
                ":collision: Model {0} deleted successfully!".format(model_id),
                fg="green",
            )

        # Remove deleted model's entry from fetched_models.txt
        with open("fetched_models.txt") as infile:
            models = dict(csv.reader(infile))
        infile.close()
        del models[model_id]
        with open('fetched_models.txt', 'w') as f:
            for key, values in models.items():
                f.write(f"{key},{values}\n")
                      