import importlib


class Command(object):
    """
    Command class to dynamically import and execute CLI commands.
    """

    def __init__(self):
        pass

    def catalog(self):
        """
        Display the catalog.
        """
        m = importlib.import_module("ersilia.cli.commands.catalog")
        m.catalog_cmd()

    def uninstall(self):
        """
        Uninstall the application.
        """
        m = importlib.import_module("ersilia.cli.commands.uninstall")
        m.uninstall_cmd()

    def close(self):
        """
        Close the application.
        """
        m = importlib.import_module("ersilia.cli.commands.close")
        m.close_cmd()

    def delete(self):
        """
        Delete the application.
        """
        m = importlib.import_module("ersilia.cli.commands.delete")
        m.delete_cmd()

    def example(self):
        """
        Show an example.
        """
        m = importlib.import_module("ersilia.cli.commands.example")
        m.example_cmd()

    def info(self):
        """
        Display information.
        """
        m = importlib.import_module("ersilia.cli.commands.info")
        m.info_cmd()

    def fetch(self):
        """
        Fetch data.
        """
        m = importlib.import_module("ersilia.cli.commands.fetch")
        m.fetch_cmd()

    def publish(self):
        """
        Publish data.
        """
        m = importlib.import_module("ersilia.cli.commands.publish")
        m.publish_cmd()

    def run(self):
        """
        Execute the command.
        """
        m = importlib.import_module("ersilia.cli.commands.run")
        m.run_cmd()

    def stop(self):
        """
        Stop the command.
        """
        m = importlib.import_module("ersilia.cli.commands.stop")
        m.stop_cmd()

    def restart(self):
        """
        Restart the command.
        """
        m = importlib.import_module("ersilia.cli.commands.restart")
        m.restart_cmd()

    def serve(self):
        """
        Serve the application.
        """
        m = importlib.import_module("ersilia.cli.commands.serve")
        m.serve_cmd()

    def setup(self):
        """
        Set up the application.
        """
        m = importlib.import_module("ersilia.cli.commands.setup")
        m.setup_cmd()

    def test(self):
        """
        Test the application.
        """
        m = importlib.import_module("ersilia.cli.commands.test")
        m.test_cmd()

    def dump(self):
        """
        Dump precalculation.
        """
        m = importlib.import_module("ersilia.cli.commands.dump")
        m.dump_cmd()
