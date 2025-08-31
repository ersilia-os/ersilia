import datetime
import os
import re

import validators

try:
    from validators import ValidationFailure
except ImportError:
    from validators import ValidationError as ValidationFailure

from ...utils.identifiers.model import ModelIdentifier

# ruff: noqa


_SEMVER_REGEX = re.compile(
    r"""^
    (0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)    
    (?:-([0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*))?   
    (?:\+([0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*))?   
    $""",
    re.VERBOSE,
)


class BaseInformationValidator:
    @staticmethod
    def _read_default_fields(field):
        root = os.path.dirname(os.path.abspath(__file__))
        filename = field.lower().replace(" ", "_")
        file_path = os.path.join(root, "metadata", filename + ".txt")
        with open(file_path, "r") as f:
            valid_field = f.read().split("\n")
        return valid_field

    @staticmethod
    def is_valid_url(u: str) -> bool:
        r = validators.url(u)
        return False if isinstance(r, ValidationFailure) else bool(r)

    @staticmethod
    def is_semver(version: str, allow_v_prefix: bool = True) -> bool:
        if not isinstance(version, str):
            return False
        s = version.strip()
        if allow_v_prefix and s[:1] in ("v", "V"):
            s = s[1:]
        if not allow_v_prefix and s != version.strip():
            return False
        return _SEMVER_REGEX.match(s) is not None

    @staticmethod
    def is_numeric(x) -> bool:
        try:
            float(x)
            return True
        except Exception:
            return False

    @staticmethod
    def to_numeric(x):
        x = str(x)
        return float(x) if "." in x else int(float(x))

    @staticmethod
    def to_list(x):
        if isinstance(x, (list, tuple)):
            return list(x)
        s = str(x).strip().strip("'\"")
        if (s.startswith("[") and s.endswith("]")) or (
            s.startswith("(") and s.endswith(")")
        ):
            items = [i.strip() for i in s[1:-1].split(",")]
            return items
        return [s]

    @classmethod
    def validate_identifier(cls, v) -> bool:
        return ModelIdentifier().is_valid(v)

    @classmethod
    def validate_slug(cls, s) -> bool:
        return s.lower() == s and 5 <= len(s) <= 60

    @classmethod
    def validate_status(cls, s) -> bool:
        return s in cls._read_default_fields("Status")

    @classmethod
    def validate_title(cls, t) -> bool:
        return 10 <= len(t) <= 300

    @classmethod
    def validate_description(cls, d) -> bool:
        return bool(d) and len(d) >= 200

    @classmethod
    def validate_mode(cls, m) -> bool:
        return m is None or m in cls._read_default_fields("Mode")

    @classmethod
    def validate_source(cls, s) -> bool:
        return s in cls._read_default_fields("Source")

    @classmethod
    def validate_source_type(cls, st) -> bool:
        return st in cls._read_default_fields("Source Type")

    @classmethod
    def validate_input(cls, inp) -> bool:
        return all(i in cls._read_default_fields("Input") for i in cls.to_list(inp))

    @classmethod
    def validate_input_shape(cls, shp) -> bool:
        return shp is None or shp in cls._read_default_fields("Input Shape")

    @classmethod
    def validate_input_dimension(cls, d) -> bool:
        return isinstance(d, int) and d >= 1

    @classmethod
    def validate_task(cls, t) -> bool:
        return isinstance(t, str) and t in cls._read_default_fields("Task")

    @classmethod
    def validate_subtask(cls, st) -> bool:
        return isinstance(st, str) and st in cls._read_default_fields("Subtask")

    @classmethod
    def validate_biomedical_area(cls, ba) -> bool:
        return all(
            v in cls._read_default_fields("Biomedical Area") for v in cls.to_list(ba)
        )

    @classmethod
    def validate_target_organism(cls, to) -> bool:
        return all(
            v in cls._read_default_fields("Target Organism") for v in cls.to_list(to)
        )

    @classmethod
    def validate_output(cls, o) -> bool:
        return all(v in cls._read_default_fields("Output") for v in cls.to_list(o))

    @classmethod
    def validate_output_type(cls, ot) -> bool:
        return ot is None or all(
            v in cls._read_default_fields("Output Type") for v in cls.to_list(ot)
        )

    @classmethod
    def validate_output_shape(cls, osz) -> bool:
        return osz is None or osz in cls._read_default_fields("Output Shape")

    @classmethod
    def validate_output_dimension(cls, d) -> bool:
        return isinstance(d, int) and d >= 1

    @classmethod
    def validate_output_consistency(cls, c) -> bool:
        return c in cls._read_default_fields("Output Consistency")

    @classmethod
    def validate_interpretation(cls, it) -> bool:
        return 10 <= len(it) <= 300

    @classmethod
    def validate_tag(cls, tg) -> bool:
        return all(v in cls._read_default_fields("Tag") for v in cls.to_list(tg))

    @classmethod
    def validate_publication(cls, url) -> bool:
        return not url or cls.is_valid_url(url)

    @classmethod
    def validate_publication_type(cls, pt) -> bool:
        return pt in cls._read_default_fields("Publication Type")

    @classmethod
    def validate_publication_year(cls, y) -> bool:
        current = datetime.date.today().year
        return isinstance(y, int) and 1900 <= y <= current

    @classmethod
    def validate_source_code(cls, url) -> bool:
        return not url or cls.is_valid_url(url)

    @classmethod
    def validate_license(cls, lic) -> bool:
        return lic in cls._read_default_fields("License")

    @classmethod
    def validate_contributor(cls, c) -> bool:
        return not c or isinstance(c, str)

    @classmethod
    def validate_incorporation_date(cls, d) -> bool:
        if not d:
            return True
        try:
            iso = datetime.datetime.fromisoformat(d).date().isoformat()
            return d == iso
        except (ValueError, TypeError):
            return False

    @classmethod
    def validate_dockerhub(cls, url) -> bool:
        return not url or url.startswith("https://hub.docker.com/r/ersiliaos/")

    @classmethod
    def validate_docker_architecture(cls, arch) -> bool:
        return all(
            v in cls._read_default_fields("Docker Architecture")
            for v in cls.to_list(arch or [])
        )

    @classmethod
    def validate_s3(cls, url) -> bool:
        return not url or url.startswith(
            "https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/"
        )

    @classmethod
    def validate_model_size(cls, sz) -> bool:
        return sz is None or cls.is_numeric(sz)

    @classmethod
    def validate_environment_size(cls, sz) -> bool:
        return sz is None or cls.is_numeric(sz)

    @classmethod
    def validate_image_size(cls, sz) -> bool:
        return sz is None or cls.is_numeric(sz)

    @classmethod
    def validate_computational_performance(cls, v) -> bool:
        return v is None or cls.is_numeric(v)

    @classmethod
    def validate_deployment(cls, dep) -> bool:
        return all(
            v in cls._read_default_fields("Deployment") for v in cls.to_list(dep or [])
        )

    @classmethod
    def validate_both_identifiers(cls, ident, slug) -> bool:
        return bool(ident and slug)
