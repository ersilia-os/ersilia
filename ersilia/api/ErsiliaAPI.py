from commands import run

class ErsiliaAPI:
    def __init__(self, model_id):
        self.model_id = model_id

    def run(self, input, batch_size):
        run.run(input, batch_size)
