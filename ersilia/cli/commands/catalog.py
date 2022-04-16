import click



from . import ersilia_cli
from ...hub.content.catalog import ModelCatalog
from ...hub.content.search import ModelSearcher
from ...hub.content.table_update import table

# define storage file path based on script path (__file__)
import os
counter_path = os.path.join(os.path.dirname(__file__), 'my_counter')
# start of script - read or initialise counter 
try:
    with open(counter_path, 'r') as count_in:
        counter = int(count_in.read())
             
except FileNotFoundError:
    counter = 0


def catalog_cmd():
    """Creates catalog command"""
    # Example usage: ersilia catalog
    @ersilia_cli.command(help="List a catalog of models")
    @click.option(
        "-l",
        "--local",
        is_flag=True,
        default=False,
        help="Show catalog of models available in the local computer",
    )
    
    @click.option("-t", "--text", default = None, type=click.STRING, help ="Shows the  model related to input keyword")
    @click.option("-m", "--mode", default = None, type=click.STRING, help = "Shows the  model trained via input mode")
    @click.option("-n", "--next",is_flag=True, default=False,  help = "Shows the next table")
    @click.option("-p", "--previous",is_flag=True, default=False, help = "Shows previous table")
    def catalog(local=False, search=None, text=None , mode=None, next = False, previous = False):

        mc = ModelCatalog()
        if not local:
            if not (next and previous):
                catalog = mc.hub()
                catalog = table(catalog, 0).initialise()
                with open(counter_path, 'w') as count_in:
                       count_in.write(str(1))
            if next:
                catalog = mc.hub()
                catalog = table(catalog, counter).next_table()
                with open(counter_path, 'w') as count_out:
                     count_out.write(str(counter+1))
            
            if previous:
                catalog = mc.hub()
                catalog = table(catalog, counter).prev_table()
                with open(counter_path, 'w') as count_out:
                     count_out.write(str(counter-1))
                if counter < 1:
                    with open(counter_path, 'w') as count_out:
                        count_out.write(str(1))
                     


        else:
            catalog = mc.local()
        
        if text:
            catalog = ModelSearcher(catalog).search(text)
        if mode:
            catalog = ModelSearcher(catalog).search(mode)
        click.echo(catalog)


        