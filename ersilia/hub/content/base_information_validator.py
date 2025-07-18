import datetime

import validators

try:
    from validators import ValidationFailure
except ImportError:
    from validators import ValidationError as ValidationFailure

from ...utils.identifiers.model import ModelIdentifier

# ruff: noqa: D101, D102


class BaseInformationValidator:
    @staticmethod
    def is_valid_url(u: str) -> bool:
        r = validators.url(u)
        return False if isinstance(r, ValidationFailure) else bool(r)

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
        return s in cls.read_list("Status")

    @classmethod
    def validate_title(cls, t) -> bool:
        return 10 <= len(t) <= 300

    @classmethod
    def validate_description(cls, d) -> bool:
        return bool(d) and len(d) >= 200

    @classmethod
    def validate_mode(cls, m) -> bool:
        return m is None or m in cls.read_list("Mode")

    @classmethod
    def validate_source(cls, s) -> bool:
        return s in cls.read_list("Source")

    @classmethod
    def validate_source_type(cls, st) -> bool:
        return st in cls.read_list("Source Type")

    @classmethod
    def validate_input(cls, inp) -> bool:
        return all(i in cls.read_list("Input") for i in cls.to_list(inp))

    @classmethod
    def validate_input_shape(cls, shp) -> bool:
        return shp is None or shp in cls.read_list("Input Shape")

    @classmethod
    def validate_input_dimension(cls, d) -> bool:
        return isinstance(d, int) and d >= 1

    @classmethod
    def validate_task(cls, t) -> bool:
        return isinstance(t, str) and t in cls.read_list("Task")

    @classmethod
    def validate_subtask(cls, st) -> bool:
        return isinstance(st, str) and st in cls.read_list("Subtask")

    @classmethod
    def validate_biomedical_area(cls, ba) -> bool:
        return all(v in cls.read_list("Biomedical Area") for v in cls.to_list(ba))

    @classmethod
    def validate_target_organism(cls, to) -> bool:
        return all(v in cls.read_list("Target Organism") for v in cls.to_list(to))

    @classmethod
    def validate_output(cls, o) -> bool:
        return all(v in cls.read_list("Output") for v in cls.to_list(o))

    @classmethod
    def validate_output_type(cls, ot) -> bool:
        return ot is None or all(
            v in cls.read_list("Output Type") for v in cls.to_list(ot)
        )

    @classmethod
    def validate_output_shape(cls, osz) -> bool:
        return osz is None or osz in cls.read_list("Output Shape")

    @classmethod
    def validate_output_dimension(cls, d) -> bool:
        return isinstance(d, int) and d >= 1

    @classmethod
    def validate_output_consistency(cls, c) -> bool:
        return c in cls.read_list("Output Consistency")

    @classmethod
    def validate_interpretation(cls, it) -> bool:
        return 10 <= len(it) <= 300

    @classmethod
    def validate_tag(cls, tg) -> bool:
        return all(v in cls.read_list("Tag") for v in cls.to_list(tg))

    @classmethod
    def validate_publication(cls, url) -> bool:
        return not url or cls.is_valid_url(url)

    @classmethod
    def validate_publication_type(cls, pt) -> bool:
        return pt in cls.read_list("Publication Type")

    @classmethod
    def validate_publication_year(cls, y) -> bool:
        current = datetime.date.today().year
        return isinstance(y, int) and 1900 <= y <= current

    @classmethod
    def validate_source_code(cls, url) -> bool:
        return not url or cls.is_valid_url(url)

    @classmethod
    def validate_license(cls, lic) -> bool:
        return lic in cls.read_list("License")

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
            v in cls.read_list("Docker Architecture") for v in cls.to_list(arch or [])
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
        return all(v in cls.read_list("Deployment") for v in cls.to_list(dep or []))

    @classmethod
    def validate_both_identifiers(cls, ident, slug) -> bool:
        return bool(ident and slug)
