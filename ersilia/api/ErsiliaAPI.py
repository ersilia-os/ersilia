from .commands import run, serve, close, info

class ErsiliaAPI:
    def __init__(self, model_id):
        self.model_id = model_id
        serve_dict = serve.serve(model=model_id,
                    port=None,
                    track=False,
                    tracking_use_case= "local",
                    enable_local_cache=True,
                    local_cache_only=False,
                    cloud_cache_only=False,
                    cache_only=False,
                    max_cache_memory_frac=None,
                    )
        print(serve_dict)

    def run(self, input, batch_size):
        print(run.run(self.model_id, input, batch_size))

    def close(self):
        close.close(self.model_id)
        print(f"Model {self.model_id} closed.")
    
    def info(self):
        info.info(self.model_id)
