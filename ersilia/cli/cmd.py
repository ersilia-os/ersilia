

class Command(object):

    def __init__(self):
        pass

    def api(self):
        from .commands.api import api_cmd
        api_cmd()

    def auth(self):
        from .commands.auth import auth_cmd
        auth_cmd()

    def card(self):
        from .commands.card import card_cmd
        card_cmd()

    def catalog(self):
        from .commands.catalog import catalog_cmd
        catalog_cmd()

    def close(self):
        from .commands.close import close_cmd
        close_cmd()

    def conda(self):
        from .commands.conda import conda_cmd
        conda_cmd()

    def delete(self):
        from .commands.delete import delete_cmd
        delete_cmd()

    def deploy(self):
        from .commands.deploy import deploy_cmd
        deploy_cmd()

    def dockerize(self):
        from .commands.dockerize import dockerize_cmd
        dockerize_cmd()

    def fetch(self):
        from .commands.fetch import fetch_cmd
        fetch_cmd()

    def serve(self):
        from .commands.serve import serve_cmd
        serve_cmd()

    def setup(self):
        from .commands.setup import setup_cmd
        setup_cmd()

    def store(self):
        from .commands.store import store_cmd
        store_cmd()
