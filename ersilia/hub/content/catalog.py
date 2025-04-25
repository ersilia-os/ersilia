"""See available models in the Ersilia Model Hub"""

import csv
import json
import os

from ... import ErsiliaBase
from ...db.hubdata.interfaces import JsonModelsInterface
from ...default import MODEL_SOURCE_FILE, TableConstants
from ...tools.bentoml.exceptions import BentoMLException
from ...utils.identifiers.model import ModelIdentifier
from ...utils.terminal import run_command
from .card import ModelCard

try:
    import webbrowser
except ModuleNotFoundError:
    webbrowser = None

try:
    from github import Github
except ModuleNotFoundError:
    Github = None


class CatalogTable(object):
    """
    Class to handle the catalog table of models.

    This class provides methods to convert the model catalog into various formats,
    such as list of dictionaries, JSON, and table format. It also provides methods
    to write the catalog table to a file.

    Parameters
    ----------
    data : list
        The data of the catalog table.
    columns : list
        The columns of the catalog table.
    """

    def __init__(self, data, columns):
        self.data = data
        self.columns = columns

    def as_list_of_dicts(self) -> list:
        """
        Convert the catalog table to a list of dictionaries.

        Returns
        -------
        list
            The catalog table as a list of dictionaries.
        """
        R = []
        for r in self.data:
            d = {}
            for i, c in enumerate(self.columns):
                d[c] = r[i]
            R += [d]
        return R

    def as_json(self) -> str:
        """
        Convert the catalog table to JSON format.

        Returns
        -------
        str
            The catalog table in JSON format.
        """
        R = self.as_list_of_dicts()
        return json.dumps(R, indent=4)

    def generate_separator_line(self, left, middle, right, horizontal, widths):
        """
        Generates a separator line for the table based on the given borders and column widths.

        Parameters
        ----------
        left : str
            The character to use for the left border of the line.
        middle : str
            The character to use between columns (as separators).
        right : str
            The character to use for the right border of the line.
        horizontal : str
            The character used to draw the horizontal border.
        widths : list of int
            A list of column widths to determine how much horizontal space each column takes.

        Returns
        -------
        str
            The formatted separator line as a string.
        """
        return left + middle.join(horizontal * (width + 2) for width in widths) + right

    def as_table(self) -> str:
        """
        Convert the catalog table to table format.

        Returns
        -------
        str
            The catalog table in table format.
        """
        column_widths = [
            max(
                len(str(item)) if item is not None else 0
                for item in [col] + [row[i] for row in self.data]
            )
            for i, col in enumerate(self.columns)
        ]

        table_constants = TableConstants

        row_format = table_constants.COLUMN_SEPARATOR.join(
            f"{{:<{width}}}" for width in column_widths
        )

        table = (
            self.generate_separator_line(
                table_constants.TOP_LEFT,
                table_constants.TOP_MIDDLE,
                table_constants.TOP_RIGHT,
                table_constants.HORIZONTAL,
                column_widths,
            )
            + "\n"
        )
        table += (
            table_constants.VERTICAL
            + table_constants.CELL_PADDING
            + row_format.format(*self.columns)
            + table_constants.CELL_PADDING
            + table_constants.VERTICAL
            + "\n"
        )
        table += (
            self.generate_separator_line(
                table_constants.MIDDLE_LEFT,
                table_constants.MIDDLE_MIDDLE,
                table_constants.MIDDLE_RIGHT,
                table_constants.HORIZONTAL,
                column_widths,
            )
            + "\n"
        )

        for index, row in enumerate(self.data):
            row = [str(item) if item is not None else "" for item in row]
            table += (
                table_constants.VERTICAL
                + table_constants.CELL_PADDING
                + row_format.format(*row)
                + table_constants.CELL_PADDING
                + table_constants.VERTICAL
                + "\n"
            )

            if index < len(self.data) - 1:
                table += (
                    self.generate_separator_line(
                        table_constants.MIDDLE_LEFT,
                        table_constants.MIDDLE_MIDDLE,
                        table_constants.MIDDLE_RIGHT,
                        table_constants.HORIZONTAL,
                        column_widths,
                    )
                    + "\n"
                )

        table += self.generate_separator_line(
            table_constants.BOTTOM_LEFT,
            table_constants.BOTTOM_MIDDLE,
            table_constants.BOTTOM_RIGHT,
            table_constants.HORIZONTAL,
            column_widths,
        )

        return table

    def write(self, file_name: str):
        """
        Write the catalog table to a file.

        Parameters
        ----------
        file_name : str
            The name of the file to write the catalog table to.
        """
        with open(file_name, "w") as f:
            if file_name.endswith(".tsv"):
                delimiter = "\t"
            else:
                delimiter = ","  # Default to CSV always
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(self.columns)
            for r in self.data:
                writer.writerow(r)

    def __str__(self):
        return self.as_json()

    def __repr__(self):
        return self.__str__()


class ModelCatalog(ErsiliaBase):
    """
    Class to handle the model catalog.

    This class provides methods to manage the model catalog, including adding, updating,
    and retrieving models.

    Attributes
    ----------
    LESS_FIELDS : list
        List of fields with less information.
    MORE_FIELDS : list
        List of fields with more information.
    """

    LESS_FIELDS = ["Identifier", "Slug"]
    MORE_FIELDS = LESS_FIELDS + [
        "Title",
        "Task",
        "Output",
        "Output Shape",
    ]

    def __init__(self, config_json=None, less=True):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.mi = ModelIdentifier()
        self.less = less

    def _is_eos(self, s):
        if self.mi.is_valid(s):
            if not self.mi.is_test(s):
                return True
        return False

    def _get_item(self, card, item):
        if "card" in card:
            card = card["card"]
        if item.lower() in card:
            return card[item.lower()]
        elif item.capitalize() in card:
            return card[item.capitalize()]
        elif item.title() in card:
            return card[item.title()]
        else:
            return None

    def _get_title(self, card):
        return self._get_item(card, "title")

    def _get_slug(self, card):
        return self._get_item(card, "slug")

    def _get_status(self, card):
        return self._get_item(card, "status")

    def _get_input(self, card):
        return self._get_item(card, "input")[0]

    def _get_output(self, card):
        return self._get_item(card, "output")[0]

    def _get_model_source(self, card):
        model_id = self._get_item(card, "identifier")
        model_source_file = os.path.join(self._model_path(model_id), MODEL_SOURCE_FILE)
        if os.path.exists(model_source_file):
            with open(model_source_file) as f:
                return f.read().rstrip()
        else:
            return None

    def _get_service_class(self, card):
        if "service_class" in card:
            return card["service_class"]
        if "Service_class" in card:
            return card["Service_class"]
        return None

    def airtable(self):
        """
        List models available in AirTable Ersilia Model Hub base.
        """
        if webbrowser:
            webbrowser.open("https://airtable.com/shrUcrUnd7jB9ChZV")  # TODO Hardcoded

    def _get_catalog(self, columns: list, model_cards: list):
        R = []
        columns = ["Index"] + columns

        idx = 0
        for card in model_cards:
            status = self._get_status(card)
            slug = self._get_slug(card)
            if status == "In Progress":
                continue
            if "test" in slug:
                continue
            idx += 1
            r = [idx]
            for field in columns[1:]:
                if field == "Model Source":
                    r += [self._get_model_source(card)]
                else:
                    r += [self._get_item(card, field)]
            # Sort R by Identifier

            R += [r]
        R = sorted(R, key=lambda x: x[0])
        return CatalogTable(data=R, columns=columns)

    def hub(self):
        """List models available in Ersilia model hub from the S3 JSON"""
        ji = JsonModelsInterface()
        models = ji.items_all()
        columns = self.LESS_FIELDS if self.less else self.MORE_FIELDS
        table = self._get_catalog(columns, models)
        return table

    def local(self) -> CatalogTable:
        """
        List models metadata from the local repository.

        Returns
        -------
        CatalogTable
            The catalog table containing the models available locally.
        """
        mc = ModelCard()
        columns = self.LESS_FIELDS if self.less else self.MORE_FIELDS + ["Model Source"]
        cards = []
        for model_id in os.listdir(self._bundles_dir):
            if not self._is_eos(model_id):
                continue
            card = mc.get(model_id)
            if not card:
                continue
            cards += [card]
        table = self._get_catalog(columns, cards)
        return table

    def bentoml(self) -> CatalogTable:
        """
        List models available as BentoServices.

        Returns
        -------
        CatalogTable
            The catalog table containing the models available as BentoServices.
        """
        try:
            stdout, stderr, returncode = run_command(["bentoml", "list"], quiet=True)
            if returncode != 0:
                raise BentoMLException(f"BentoML list failed: {stderr}")

            # Process stdout to build CatalogTable
            output_lines = stdout.split("\n")
            if not output_lines or len(output_lines) == 1:
                return CatalogTable(data=[], columns=[])  # Return empty table

            # Extract columns and values
            columns = ["BENTO_SERVICE", "AGE", "APIS", "ARTIFACTS"]
            header = output_lines[0]
            values = output_lines[1:]

            # Parse table data
            cut_idxs = [header.find(col) for col in columns]
            R = []
            for row in values:
                r = []
                for i, idx in enumerate(zip(cut_idxs, cut_idxs[1:] + [None])):
                    r.append(
                        row[idx[0] : idx[1]].rstrip()
                        if idx[1]
                        else row[idx[0] :].rstrip()
                    )
                R.append([r[0].split(":")[0]] + r)

            return CatalogTable(data=R, columns=["Identifier"] + columns)

        except BentoMLException:
            raise
        except Exception as e:
            self.logger.error(f"Unexpected error: {str(e)}")
            raise BentoMLException(f"Failed to fetch BentoML models: {str(e)}")
