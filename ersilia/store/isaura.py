import pandas as pd

from ..core.session import Session
from ..default import DEFAULT_PRIVATE_BUCKET, DEFAULT_PUBLIC_BUCKET, ERSILIA_BUCKET
from ..utils.logging import logger as log

try:
    from isaura.helpers import logger
    from isaura.manage import (
        IsauraChecker,
        IsauraCopy,
        IsauraInspect,
        IsauraReader,
        IsauraWriter,
    )
except ImportError:
    IsauraChecker = None
    IsauraReader = None
    IsauraInspect = None
    logger = None

# ruff: noqa: D101, D102


class IsauraStore:
    def __init__(self):
        self.checker = None
        self.reader = None
        self.check_dict = {}

        if not self.is_installed():
            return

        logger.set_verbosity(log.verbosity)

        session = Session(config_json=None)
        model_id = session.current_model_id()
        current_status = session.current_store_status()
        self.access = current_status[2]
        nn = current_status[-1]
        bucket = self.resolve_default_bucket(self.access)
        version = self.latest_version(model_id, bucket)
        self.model_id = model_id
        self.version = version

        self.reader = IsauraReader(
            bucket=bucket,
            model_id=model_id,
            model_version=version,
            approximate=nn,
            input_csv=None,
        )
        self.writer = IsauraWriter(
            input_csv=None,
            model_id=self.model_id,
            model_version=self.version,
            access=self.access,
            bucket=ERSILIA_BUCKET,
        )
        self.copier = IsauraCopy(
            model_id=self.model_id, model_version=self.version, bucket=ERSILIA_BUCKET
        )

    @staticmethod
    def is_installed():
        return (
            IsauraChecker is not None
            and IsauraReader is not None
            and IsauraInspect is not None
            and IsauraWriter is not None
            and IsauraCopy is not None
            and logger is not None
        )

    def check(self, vs):
        cks = IsauraChecker(
            bucket=DEFAULT_PUBLIC_BUCKET[0],
            model_id=self.model_id,
            model_version=self.version,
        ).seen_many(vs)
        missed = [k for k, v in cks.items() if not v[0]]
        if missed:
            ckss = IsauraChecker(
                bucket=DEFAULT_PRIVATE_BUCKET[0],
                model_id=self.model_id,
                model_version=self.version,
            ).seen_many(missed)
            for k, v in ckss.items():
                cks[k] = v

        self.check_dict = cks
        cks = {k: v[0] for k, v in cks.items()}
        return cks

    def write(self, df):
        self.writer.write(df=df)
        self.copier.copy()

    def read(self, vs):
        if not self.reader:
            return pd.DataFrame(columns=["input"])

        priv_inputs = [
            inp
            for inp, v in self.check_dict.items()
            if v[-1] == DEFAULT_PRIVATE_BUCKET[0] and v[0]
        ]
        pub_inputs = [
            inp
            for inp, v in self.check_dict.items()
            if v[-1] == DEFAULT_PUBLIC_BUCKET[0] and v[0]
        ]

        frames = []
        if priv_inputs:
            self.reader.bucket = DEFAULT_PRIVATE_BUCKET[0]
            priv_df = self.reader.read(df=pd.DataFrame(priv_inputs, columns=["input"]))
            frames.append(priv_df)

        if pub_inputs:
            self.reader.bucket = DEFAULT_PUBLIC_BUCKET[0]
            pub_df = self.reader.read(df=pd.DataFrame(pub_inputs, columns=["input"]))
            frames.append(pub_df)

        if not frames:
            return pd.DataFrame(columns=["input"])

        df = pd.concat(frames, ignore_index=True)

        if "input" in df.columns:
            df["input"] = pd.Categorical(df["input"], categories=vs, ordered=True)
            df = df.sort_values("input").reset_index(drop=True)
        return df

    def get_bucket_records(self, bucket):
        if not self.is_installed():
            return []
        insp = IsauraInspect(model_id="_", model_version="_", cloud=False)
        return insp.inspect_models(bucket, prefix_filter="")

    def latest_version(self, model, bucket):
        if not self.is_installed() or model is None:
            return None

        records = self.get_bucket_records(bucket)
        versions = []

        for r in records:
            name = r["model"]
            if not name.startswith(model + "/"):
                continue
            parts = name.split("/")
            if len(parts) < 2:
                continue
            v = parts[1]
            if len(v) > 1 and v[0] == "v" and v[1:].isdigit():
                versions.append(int(v[1:]))

        return f"v{max(versions)}" if versions else "v1"

    def resolve_default_bucket(self, access):
        if access is None or access == "public":
            return DEFAULT_PUBLIC_BUCKET[0]
        return DEFAULT_PRIVATE_BUCKET[0]
