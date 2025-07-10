from .commands import close, delete, example, fetch, info, run, serve


class ErsiliaAPI:
    """
    Python API wrapper for interacting with Ersilia Model Hub models.

    This class provides a programmatic interface to run, serve, fetch, and manage
    machine learning models from the Ersilia Model Hub. It wraps the existing CLI-based 
    functionality in a clean Pythonic interface and supports use in context managers.

    Parameters
    ----------
    model_id : str
        The unique identifier of the model to be managed and executed.

    Methods
    -------
    fetch():
        Downloads the specified model and its dependencies from DockerHub.
        
    serve():
        Serves the specified model locally to prepare for inference.

    run(input, output, batch_size):
        Runs the model on the given input and writes predictions to the output file.

    close():
        Terminates the model server and cleans up associated resources.

    info():
        Prints metadata and technical details about the specified model.

    example(file_name, simple, random, n_samples, deterministic):
        Generates example input files for the model.

    delete():
        Deletes the specified model and its local artifacts.

    __enter__():
        Context manager entry — automatically serves the model.

    __exit__(exc_type, exc_value, traceback):
        Context manager exit — automatically closes the served model.

    Usage Example
    -------------
    >>> from ersilia.api.create_api import ErsiliaAPI
    >>> molecular_weight = ErsiliaAPI("eos3b5e")
    >>> molecular_weight.fetch()
    >>> molecular_weight.serve()
    >>> with molecular_weight as model:
    >>>     model.info() 
    """
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

    def __enter__(self):
        self.serve()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
