from .commands import close, delete, example, fetch, info, run, serve


class ErsiliaAPI:
    def __init__(self, model_id):
        self.model_id = model_id
        
    def fetch(self):
        fetch.fetch(
            model=self.model_id,
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
            verbose=False,
        )

    def serve(self):
        serve.serve(
            self.model_id,
            port=None,
            track=False,
            tracking_use_case="local",
            enable_local_cache=True,
            local_cache_only=False,
            cloud_cache_only=False,
            cache_only=False,
            max_cache_memory_frac=None,
            verbose=False,
        )

    def run(self, input, output, batch_size):
        print(run.run(self.model_id, input, output, batch_size))

    def close(self):
        close.close(self.model_id)

    def info(self):
        info.info(self.model_id)

    def example(self, file_name, simple, random, n_samples, deterministic):
        print(
            example.example(
                self.model_id, file_name, simple, random, n_samples, deterministic
            )
        )

    def delete(self):
        delete.delete(self.model_id, verbose=False)

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()