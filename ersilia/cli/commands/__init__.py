import click
from bentoml.cli.click_utils import BentoMLCommandGroup
from ... import __version__
from ... import logger
from ..echo import Silencer
from pathlib import Path
from ...utils import cron
import json
import csv
import time


# in days :
model_usage_lim   = 30 
model_cleanup_lim = 7

def seconds_to_days(s):
    return s / (24 * 3600)

@click.group(cls=BentoMLCommandGroup)
@click.version_option(version=__version__)
@click.option(
    "-v",
    "--verbose",
    default=False,
    is_flag=True,
    help="Show logging on terminal when running commands.",
)
@click.option(
    "-s",
    "--silent",
    default=False,
    is_flag=True,
    help="Do not echo any progress message.",
)
def ersilia_cli(verbose, silent):
    """
    Ersilia CLI
    """
    if verbose:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)
    silencer = Silencer()
    silencer.speak()  # To reset default
    if silent:
        silencer.silence()
    
    ts_dict = {}
    if not Path('last_cleaned.json').exists():
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
                for m_name in fetched_models: 
                    if (seconds_to_days(current_ts-float(fetched_models[m_name])))>model_usage_lim:
                        cron.del_unused_model(model_id=m_name)



