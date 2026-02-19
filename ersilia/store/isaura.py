import time
from datetime import datetime, timezone

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
    IsauraWriter = None
    IsauraCopy = None
    logger = None

# ruff: noqa: D101, D102


class IsauraStore:
    def __init__(self):
        self.checker = None
        self.reader = None
        self.check_dict = {}

        session = Session(config_json=None)
        current_status = session.current_store_status()
        read_store, write_store = current_status[0], current_status[1]
        installed = self.is_installed()
        if not installed or (not read_store and not write_store):
            return

        logger.set_verbosity(log.verbosity)
        model_id = session.current_model_id()
        self.access = current_status[2]
        nn = current_status[3]
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
            model_id=self.model_id,
            model_version=self.version,
            bucket=ERSILIA_BUCKET,
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
        t0 = time.perf_counter()
        ts = datetime.now(timezone.utc).isoformat(timespec="milliseconds")
        log.debug(
            f"IsauraStore.check start ts={ts} n={len(vs)} model_id={self.model_id} version={self.version}"
        )
        try:
            t = time.perf_counter()
            cks = IsauraChecker(
                bucket=DEFAULT_PUBLIC_BUCKET[0],
                model_id=self.model_id,
                model_version=self.version,
            ).seen_many(vs)
            log.debug(
                f"IsauraStore.check seen_many public n={len(vs)} dt={(time.perf_counter()-t):.6f}s"
            )

            t = time.perf_counter()
            missed = [k for k, v in cks.items() if not v[0]]
            log.debug(
                f"IsauraStore.check missed n={len(missed)} dt={(time.perf_counter()-t):.6f}s"
            )

            if missed:
                t = time.perf_counter()
                ckss = IsauraChecker(
                    bucket=DEFAULT_PRIVATE_BUCKET[0],
                    model_id=self.model_id,
                    model_version=self.version,
                ).seen_many(missed)
                log.debug(
                    f"IsauraStore.check seen_many private n={len(missed)} dt={(time.perf_counter()-t):.6f}s"
                )

                t = time.perf_counter()
                cks.update(ckss)
                log.debug(f"IsauraStore.check merge dt={(time.perf_counter()-t):.6f}s")

            self.check_dict = cks
            out = {k: v[0] for k, v in cks.items()}
            log.debug(
                f"IsauraStore.check done dt_total={(time.perf_counter()-t0):.6f}s"
            )
            return out
        except Exception as e:
            log.error(
                f"IsauraStore.check failed dt_total={(time.perf_counter()-t0):.6f}s err={e}"
            )
            self.check_dict = {}
            return {}

    def write(self, df):
        t0 = time.perf_counter()
        n = len(df) if hasattr(df, "__len__") else None
        log.debug(
            f"IsauraStore.write start rows={n} model_id={getattr(self, 'model_id', None)} version={getattr(self, 'version', None)}"
        )

        t = time.perf_counter()
        self.writer.write(df=df)
        log.debug(f"IsauraStore.write writer_write dt={(time.perf_counter() - t):.6f}s")

        t = time.perf_counter()
        self.copier.copy()
        log.debug(f"IsauraStore.write copier_copy dt={(time.perf_counter() - t):.6f}s")

        log.debug(f"IsauraStore.write done dt_total={(time.perf_counter() - t0):.6f}s")

    def read(self, vs):
        if not self.reader:
            log.debug("IsauraStore.read no_reader returning_empty_df")
            return pd.DataFrame(columns=["input"])

        t0 = time.perf_counter()
        log.debug(
            f"IsauraStore.read start n={len(vs) if vs is not None else None} model_id={self.model_id} version={self.version}"
        )

        try:
            t = time.perf_counter()
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
            log.debug(
                f"IsauraStore.read split priv={len(priv_inputs)} pub={len(pub_inputs)} "
                f"check_dict={len(self.check_dict)} dt={(time.perf_counter() - t):.6f}s"
            )

            frames = []

            if priv_inputs:
                t = time.perf_counter()
                self.reader.bucket = DEFAULT_PRIVATE_BUCKET[0]
                log.debug(
                    f"IsauraStore.read set_bucket bucket={self.reader.bucket} dt={(time.perf_counter() - t):.6f}s"
                )

                t = time.perf_counter()
                priv_df = self.reader.read(
                    df=pd.DataFrame(priv_inputs, columns=["input"])
                )
                log.debug(
                    f"IsauraStore.read reader_read bucket={self.reader.bucket} rows={len(priv_df)} dt={(time.perf_counter() - t):.6f}s"
                )
                frames.append(priv_df)

            if pub_inputs:
                t = time.perf_counter()
                self.reader.bucket = DEFAULT_PUBLIC_BUCKET[0]
                log.debug(
                    f"IsauraStore.read set_bucket bucket={self.reader.bucket} dt={(time.perf_counter() - t):.6f}s"
                )

                t = time.perf_counter()
                pub_df = self.reader.read(
                    df=pd.DataFrame(pub_inputs, columns=["input"])
                )
                log.debug(
                    f"IsauraStore.read reader_read bucket={self.reader.bucket} rows={len(pub_df)} dt={(time.perf_counter() - t):.6f}s"
                )
                frames.append(pub_df)

            t = time.perf_counter()
            if not frames:
                log.debug(
                    f"IsauraStore.read no_frames dt={(time.perf_counter() - t):.6f}s dt_total={(time.perf_counter() - t0):.6f}s"
                )
                return pd.DataFrame(columns=["input"])
            log.debug(
                f"IsauraStore.read frames_ok frames={len(frames)} dt={(time.perf_counter() - t):.6f}s"
            )

            t = time.perf_counter()
            df = pd.concat(frames, ignore_index=True)
            log.debug(
                f"IsauraStore.read concat rows={len(df)} cols={len(df.columns)} dt={(time.perf_counter() - t):.6f}s"
            )

            t = time.perf_counter()
            if "input" in df.columns:
                order = {}
                for i, v in enumerate(vs):
                    if v not in order:
                        order[v] = i
                df["_order"] = df["input"].map(order)
                df["_order"] = df["_order"].fillna(len(order)).astype(int)
                df = (
                    df.sort_values("_order")
                    .drop(columns=["_order"])
                    .reset_index(drop=True)
                )
            log.debug(
                f"IsauraStore.read reorder_sort dt={(time.perf_counter() - t):.6f}s"
            )

            log.debug(
                f"IsauraStore.read done rows={len(df)} dt_total={(time.perf_counter() - t0):.6f}s"
            )
            return df
        except Exception as e:
            log.error(
                f"IsauraStore.read failed model_id={self.model_id} version={self.version} "
                f"dt_total={(time.perf_counter() - t0):.6f}s err={e}"
            )
            return pd.DataFrame(columns=["input"])

    def get_bucket_records(self, bucket):
        t0 = time.perf_counter()
        if not self.is_installed():
            log.debug(
                f"IsauraStore.get_bucket_records not_installed bucket={bucket} dt={(time.perf_counter() - t0):.6f}s"
            )
            return []
        insp = IsauraInspect(model_id="_", model_version="_", cloud=False)
        t1 = time.perf_counter()
        log.debug(
            f"IsauraStore.get_bucket_records inspect_init dt={(t1 - t0):.6f}s bucket={bucket}"
        )
        recs = insp.inspect_models(bucket, prefix_filter="")
        t2 = time.perf_counter()
        log.debug(
            f"IsauraStore.get_bucket_records inspect_models n={len(recs) if recs is not None else None} dt={(t2 - t1):.6f}s"
        )
        return recs

    def latest_version(self, model, bucket):
        t0 = time.perf_counter()
        if not self.is_installed() or model is None:
            log.debug(
                f"IsauraStore.latest_version skip installed={self.is_installed()} model={model} bucket={bucket} dt={(time.perf_counter() - t0):.6f}s"
            )
            return None

        records = self.get_bucket_records(bucket)
        t1 = time.perf_counter()
        log.debug(
            f"IsauraStore.latest_version got_records n={len(records)} dt={(t1 - t0):.6f}s model={model} bucket={bucket}"
        )

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

        out = f"v{max(versions)}" if versions else "v1"
        log.debug(
            f"IsauraStore.latest_version resolved version={out} candidates={len(versions)} dt_total={(time.perf_counter() - t0):.6f}s"
        )
        return out

    def resolve_default_bucket(self, access):
        bucket = (
            DEFAULT_PUBLIC_BUCKET[0]
            if access is None or access == "public"
            else DEFAULT_PRIVATE_BUCKET[0]
        )
        log.debug(f"IsauraStore.resolve_default_bucket access={access} bucket={bucket}")
        return bucket
