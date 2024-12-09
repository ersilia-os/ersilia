"""See available models in the Ersilia Model Hub"""

import subprocess
import requests
import shutil
import os
import json
import csv
from .card import ModelCard
from ... import ErsiliaBase
from ...utils.identifiers.model import ModelIdentifier
from ...db.hubdata.interfaces import JsonModelsInterface
from ...default import BENTOML_PATH, MODEL_SOURCE_FILE
from ...default import TableConstants
from ... import logger

try:
    import webbrowser
except ModuleNotFoundError as err:
    webbrowser = None

try:
    from github import Github
except ModuleNotFoundError as err:
    Github = None


class CatalogTable(object):
    def __init__(self, data, columns):
        self.data = data
        self.columns = columns

    def as_list_of_dicts(self):
        R = []
        for r in self.data:
            d = {}
            for i, c in enumerate(self.columns):
                d[c] = r[i]
            R += [d]
        return R

    def as_json(self):
        R = self.as_list_of_dicts()
        return json.dumps(R, indent=4)

    def generate_separator_line(self, left, middle, right, horizontal, widths):
        """
        Generates a separator line for the table based on the given borders and column widths.

        The line starts with a 'left' border, followed by repeated 'horizontal'
        sections for each column's width, joined by 'middle' separators, and ends
        with a 'right' border.

        Args:
            left (str): The character to use for the left border of the line.
            middle (str): The character to use between columns (as separators).
            right (str): The character to use for the right border of the line.
            horizontal (str): The character used to draw the horizontal border.
            widths (list[int]): A list of column widths to determine how much
                                horizontal space each column takes.

        Returns:
            str: The formatted separator line as a string.
        """
        return left + middle.join(horizontal * (width + 2) for width in widths) + right

    def as_table(self):
        """
        Returns the catalog data in table format. The method calculates the
        column widths dynamically by determining the maximum width of each column,
        based on the data and column headers.

        A row format string is then constructed using the column widths,
        specifying that each cell is left-aligned and padded with a column seperator sperating the rows.

        The method starts by constructing the top border using the
        'generate_separator_line' helper. It then adds the headers,
        formatted to fit the column widths, followed by a separator line
        also created by the helper function.

        A 'for' loop iterates over the data rows, adding each row to the
        table with borders and padding. After each row, a separator line is inserted.
        Finally, the bottom border is added using the helper function, completing the table.

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

    def write(self, file_name):
        with open(file_name, "w") as f:
            if file_name.endswith(".csv"):
                delimiter = ","
            elif file_name.endswith(".tsv"):
                delimiter = "\t"
            else:
                return None
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(self.columns)
            for r in self.data:
                writer.writerow(r)

    def __str__(self):
        return self.as_json()

    def __repr__(self):
        return self.__str__()


class ModelCatalog(ErsiliaBase):
    def __init__(self, config_json=None, only_identifier=True):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.mi = ModelIdentifier()
        self.only_identifier = only_identifier

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

    def _get_model_source(self, model_id):
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
        """List models available in AirTable Ersilia Model Hub base"""
        if webbrowser:
            webbrowser.open("https://airtable.com/shrUcrUnd7jB9ChZV") #TODO Hardcoded


    def hub(self):
        """List models available in Ersilia model hub from the S3 JSON"""
        mc = ModelCard()
        ji = JsonModelsInterface()
        models = ji.items_all()
        R = []
        for model in models:
            status = self._get_status(model)
            if status == "In Progress":
                continue
            
            identifier = self._get_item(model, "identifier")
            if self.only_identifier:
                R += [[identifier]]
            else:
                slug = self._get_slug(model)
                title = self._get_title(model)
                R += [[identifier, slug, title]]

        columns = ["Identifier"] if self.only_identifier else ["Identifier", "Slug", "Title"]
        return CatalogTable(R, columns=columns)
                
    def local(self):
        """List models available locally"""
        mc = ModelCard()
        R = []
        logger.debug("Looking for models in {0}".format(self._bundles_dir))
        if self.only_identifier:
            R = []
            for model_id in os.listdir(self._bundles_dir):
                if not self._is_eos(model_id):
                    continue
                R += [[model_id]]
            columns = ["Identifier"]
        else:
            for model_id in os.listdir(self._bundles_dir):
                if not self._is_eos(model_id):
                    continue
                card = mc.get(model_id)
                slug = self._get_slug(card)
                title = self._get_title(card)
                status = self._get_status(card)
                inputs = self._get_input(card)
                output = self._get_output(card)
                model_source = self._get_model_source(model_id)
                service_class = self._get_service_class(card)
                R += [
                    [
                        model_id,
                        slug,
                        title,
                        status,
                        inputs,
                        output,
                        model_source,
                        service_class,
                    ]
                ]
            columns = [
                "Identifier",
                "Slug",
                "Title",
                "Status",
                "Input",
                "Output",
                "Model Source",
                "Service Class",
            ]
        logger.info("Found {0} models".format(len(R)))
        if len(R) == 0:
            return CatalogTable(data=[], columns=columns)
        return CatalogTable(data=R, columns=columns)

    def bentoml(self):
        """List models available as BentoServices"""
        try:
            result = subprocess.run(
                ["bentoml", "list"], stdout=subprocess.PIPE, env=os.environ, timeout=10
            )
        except Exception as e:
            shutil.rmtree(BENTOML_PATH)
            return None
        result = [r for r in result.stdout.decode("utf-8").split("\n") if r]
        if len(result) == 1:
            return
        columns = ["BENTO_SERVICE", "AGE", "APIS", "ARTIFACTS"]
        header = result[0]
        values = result[1:]
        cut_idxs = []
        for col in columns:
            cut_idxs += [header.find(col)]
        R = []
        for row in values:
            r = []
            for i, idx in enumerate(zip(cut_idxs, cut_idxs[1:] + [None])):
                r += [row[idx[0] : idx[1]].rstrip()]
            R += [[r[0].split(":")[0]] + r]
        columns = ["Identifier"] + columns
        return CatalogTable(data=R, columns=columns)
