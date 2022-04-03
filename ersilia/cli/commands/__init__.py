import click
from bentoml.cli.click_utils import BentoMLCommandGroup
from ... import __version__
from ... import logger
from ..echo import Silencer

import json
import csv
import time

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
    if not exists('last_cleaned.json'):
        ts_dict['timestamp'] = str(time.time())
        with open('last_cleaned.json', 'w') as outfile:
            json.dump(ts_dict, outfile)
    else:
        current_ts = time.time()
        with open('last_cleaned.json') as json_file:
            ts_dict = json.load(json_file)
            ts = float(ts_dict['timestamp'])
            if(current_ts - ts) > 604800:
                with open("fetched_models.txt") as infile:
                    fetched_models = dict(csv.reader(infile))
                for m_name in fetched_models: 
                    if( current_ts - float(fetched_models[m_name])) > 2592000:
                        del_cmd = 'ersilia delete ' + m_name
                        test_op = os.system(del_cmd)
