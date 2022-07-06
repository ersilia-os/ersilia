import importlib


class Command(object):
    def __init__(self):
        pass

    def api(self):
        m = importlib.import_module("ersilia.cli.commands.api")
        m.api_cmd()

    def auth(self):
        m = importlib.import_module("ersilia.cli.commands.auth")
        m.auth_cmd()

    def card(self):
        m = importlib.import_module("ersilia.cli.commands.card")
        m.card_cmd()

    def catalog(self):
        m = importlib.import_module("ersilia.cli.commands.catalog")
        m.catalog_cmd()

    def clear(self):
        m = importlib.import_module("ersilia.cli.commands.clear")
        m.clear_cmd()

    def close(self):
        m = importlib.import_module("ersilia.cli.commands.close")
        m.close_cmd()

    def delete(self):
        m = importlib.import_module("ersilia.cli.commands.delete")
        m.delete_cmd()

    def example(self):
        m = importlib.import_module("ersilia.cli.commands.example")
        m.example_cmd()

    def fetch(self):
        m = importlib.import_module("ersilia.cli.commands.fetch")
        m.fetch_cmd()

    def publish(self):
        m = importlib.import_module("ersilia.cli.commands.publish")
        m.publish_cmd()

    def serve(self):
        m = importlib.import_module("ersilia.cli.commands.serve")
        m.serve_cmd()

    def setup(self):
        m = importlib.import_module("ersilia.cli.commands.setup")
        m.setup_cmd()
