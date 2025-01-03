import csv
import os

# TODO: For now, only explicitly tabular results are returned. We could, in principle, output any other result

MAX_LOG_COLUMNS = 10


class TabularResultLogger(object):
    """
    Class to log tabular results.

    This class handles logging of results from tabular files such as CSV and TSV.

    Parameters
    ----------
    None
    """

    def __init__(self):
        pass

    def tabulate(self, result, identifier=None, model_id=None):
        """
        Tabulate the results from a file.

        Parameters
        ----------
        result : str
            Path to the result file.
        identifier : str, optional
            Identifier to include in the results.
        model_id : str, optional
            Model ID to include in the results.

        Returns
        -------
        list or None
            List of tabulated results or None if the file is not tabular.
        """
        if self._is_tabular_file(result):
            if result.endswith(".h5"):
                return False  # TODO include HDF5 compatibility
            delimiter = self._get_delimiter(result)
            with open(result, "r") as f:
                reader = csv.reader(f, delimiter=delimiter)
                h = []
                if identifier:
                    h += ["identifier"]
                if model_id:
                    h += ["model_id"]
                h += next(reader)[:MAX_LOG_COLUMNS]
                R = [h]
                for r in reader:
                    s = []
                    if identifier:
                        s += [identifier]
                    if model_id:
                        s += [model_id]
                    s += r[:MAX_LOG_COLUMNS]
                    R += [s]
                return R
        return None

    def _is_tabular_file(self, s):
        if type(s) is str:
            if not os.path.exists(s):
                return False
            if s.endswith(".csv"):
                return True
            if s.endswith(".tsv"):
                return True
            if s.endswith(".h5"):
                return True
        else:
            return False

    def _get_delimiter(self, s):
        if s.endswith(".csv"):
            return ","
        if s.endswith(".tsv"):
            return "\t"
        return None
