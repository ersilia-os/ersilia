from .commands import run, serve, close, info, example, fetch, delete

class ErsiliaAPI:
    def __init__(self, model_id):
        self.model_id = model_id
        print(fetch.fetch(model=model_id,
        overwrite=True,
        from_dir=None,
        from_github=False,
        from_dockerhub=True,
        version="latest", 
        from_s3=False,
        from_hosted=False,
        with_fastapi=True,
        with_bentoml=None,
        hosted_url=None,
        ))
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
        print(close.close(self.model_id))
    
    def info(self):
        info.info(self.model_id)

    def example(self, n_samples, random, deterministic):
        print(example.example(self.model_id, n_samples, random, deterministic))

    def delete(self):
        print(delete.delete(self.model_id))

