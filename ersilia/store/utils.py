class OutputSource():
    LOCAL_ONLY = "local-only"
    CLOUD_ONLY = "cloud-only"
    # CLOUD_FIRST = "cloud-first"
    ALL = [
        LOCAL_ONLY,
        CLOUD_ONLY,
        #CLOUD_FIRST
        ]

    @classmethod
    def is_local(cls, option):
        return option == cls.LOCAL_ONLY

    @classmethod
    def is_cloud(cls, option):
        return option in (
            cls.CLOUD_ONLY,
            #cls.CLOUD_FIRST
            )

from pydantic import BaseModel
class InferenceStoreApiPayload(BaseModel):
    model: str
    inputs: list[str] = [] # validation error if inputs is None e.g. if inchi->smiles fails