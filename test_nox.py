import os
from test.playground.runner import NoxRunner

BASE_PATH = "test/playground"

config_path = os.path.join(os.getcwd(), BASE_PATH, "config.yml")
noxfile = os.path.join(os.getcwd(), BASE_PATH, "noxfile.py")
if __name__ == "__main__":
    runner = NoxRunner(config_path=config_path, noxfile=noxfile)
    runner.setup()
    runner.test_from_github()
    runner.execute_all()
