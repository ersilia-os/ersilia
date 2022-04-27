import click



from . import ersilia_cli
from ...hub.content.catalog import ModelCatalog
from ...hub.content.search import ModelSearcher
from ...hub.content.table_update import table


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
            if not (next or previous) :
                catalog = mc.hub()
                catalog = table(catalog).initialise()
                
            if next:
                catalog = mc.hub()
                catalog = table(catalog).next_table()
                
            
            if previous:
                catalog = mc.hub()
                catalog = table(catalog).prev_table()
                
        else:
            catalog = mc.local()
        
        if text:
            catalog = ModelSearcher(catalog).search(text)
        if mode:
            catalog = ModelSearcher(catalog).search(mode)
        click.echo(catalog)


        