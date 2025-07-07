from .commands import run, serve

class ErsiliaAPI:
    def __init__(self, model_id):
        self.model_id = model_id
        serve.serve(model=model_id,
                    port=None,
                    track=False,
                    tracking_use_case= "local",
                    enable_local_cache=True,
                    local_cache_only=False,
                    cloud_cache_only=False,
                    cache_only=False,
                    max_cache_memory_frac=None,
                    )

    def run(self, input, batch_size):
        run.run(input, batch_size)
