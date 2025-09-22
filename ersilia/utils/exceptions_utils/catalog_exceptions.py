from .exceptions import ErsiliaError

# ruff: noqa: D101, D102


class CatalogErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running catalog command"
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)


class LocalCatalogDestBundleModelsMismatch(ErsiliaError):
    def __init__(self, dest_model_ids, bundle_model_ids):
        self.message = (
            "The models in the $EOS/dest folder do not match the models in the "
            "$EOS/bundles folder. Models in both folders must be the same, "
            "otherwise, there may be inconsistencies with the catalog and the "
            "management of Ersilia.\n"
            "It is possible that this issue is caused by an error at fetching time, "
            "which caused a model to be in the dest folder but not in the bundles folder.\n"
            "To solve this issue, please delete manually the problematic models."
        )
        self.hints = (
            f"Dest models: {dest_model_ids}, Bundle models: {bundle_model_ids}\n"
            f"Models in dest but not in bundles: {sorted(set(dest_model_ids) - set(bundle_model_ids))}\n"
            f"Models in bundles but not in dest: {sorted(set(bundle_model_ids) - set(dest_model_ids))}\n"
            "To solve this issue, please delete manually the problematic models."
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class LocalCatalogIdleModelsInConda(ErsiliaError):
    def __init__(self, idle_model_ids):
        self.message = (
            "There are idle models in the local catalog that have a conda environment "
            "but are not registered in the catalog. This may cause inconsistencies "
            "with the catalog and the management of Ersilia.\n"
        )
        self.hints = (
            f"Idle models with conda env: {idle_model_ids}\n"
            "To solve this issue, please delete manually the problematic models."
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class LocalCatalogIdleModelsInDocker(ErsiliaError):
    def __init__(self, idle_model_ids):
        self.message = (
            "There are idle models in the local catalog that have a docker image "
            "but are not registered in the catalog. This may cause inconsistencies "
            "with the catalog and the management of Ersilia.\n"
        )
        self.hints = (
            f"Idle models with docker image: {idle_model_ids}\n"
            "To solve this issue, please delete manually the problematic models."
        )
        ErsiliaError.__init__(self, self.message, self.hints)
