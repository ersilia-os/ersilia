"""Search for models"""

import re

import numpy as np

from .catalog import CatalogTable


class ModelSearcher(object):
    """
    Class to handle model search operations.

    This class provides methods to search for models in the catalog table
    using text and mode of training.

    Parameters
    ----------
    catalog : CatalogTable
        The catalog table containing the models.
    """

    """This class is used for searching through the catalog table

    Attributes:
    catalog : list and description of model on the hub or in local

    """

    def __init__(self, catalog):
        self.catalog = catalog

    def levenshtein_ratio_and_distance(self, s: str, t: str) -> float:
        """
        Algorithm to fuzzy match strings using Levenshtein distance.

        Parameters
        ----------
        s : str
            The first string.
        t : str
            The second string.

        Returns
        -------
        float
            The Levenshtein distance ratio between the two strings.
        """
        # Initialize matrix of zeros
        s = str(s)
        t = str(t)
        t = t.lower()
        rows = len(s) + 1
        cols = len(t) + 1
        distance = np.zeros((rows, cols), dtype=int)

        # Populate matrix of zeros with the indeces of each character of both strings
        for i in range(1, rows):
            for k in range(1, cols):
                distance[i][0] = i
                distance[0][k] = k

        # Iterate over the matrix to compute the cost of deletions,insertions and/or substitutions
        for col in range(1, cols):
            for row in range(1, rows):
                if s[row - 1] == t[col - 1]:
                    cost = 0  # If the characters are the same in the two strings in a given position [i,j] then the cost is 0
                else:
                    cost = 2

                distance[row][col] = min(
                    distance[row - 1][col] + 1,  # Cost of deletions
                    distance[row][col - 1] + 1,  # Cost of insertions
                    distance[row - 1][col - 1] + cost,
                )  # Cost of substitutions

        # Computation of the Levenshtein Distance Ratio
        Ratio = ((len(s) + len(t)) - distance[row][col]) / (len(s) + len(t))
        return Ratio

    def search_text(self, s: str) -> CatalogTable:
        """
        Search using text and return the closest matching string.

        Parameters
        ----------
        s : str
            The text to search for.

        Returns
        -------
        CatalogTable
            The catalog table containing the closest matching models.
        """
        idxs = set()
        s = s.lower()
        ratio_list = []
        data = self.catalog.data
        for i, r in enumerate(data):
            string_ratio = []
            ratio = self.levenshtein_ratio_and_distance(s, r[0])
            string_ratio.append(ratio)
            ratio = self.levenshtein_ratio_and_distance(s, r[1])
            string_ratio.append(ratio)
            x = re.split(r"\s", r[2])
            for r1 in x:
                ratio = self.levenshtein_ratio_and_distance(s, r1)
                string_ratio.append(ratio)
            ratio = max(string_ratio)
            ratio_list.append(ratio)
        if 1.0 in ratio_list:
            idxs = [i for i, j in enumerate(ratio_list) if j == 1.0]
        else:
            idxs = [i for i, j in enumerate(ratio_list) if j > 0.7]
            print(s, "not found. Here's a close match:")

        isEmpty = len(idxs) == 0
        if isEmpty:
            return "No such element present. Please look through the catalog table."
        R = [r for i, r in enumerate(self.catalog.data) if i in idxs]

        return CatalogTable(data=R, columns=self.catalog.columns)

    def search_mode(self, s: str) -> CatalogTable:
        """
        Search the mode of training.

        Parameters
        ----------
        s : str
            The mode of training to search for.

        Returns
        -------
        CatalogTable
            The catalog table containing the models with the specified mode of training.
        """
        idxs = set()
        s = s.lower()
        data = self.catalog.data
        for i, r in enumerate(data):
            r_ = r[3].lower()
            if s in r_:
                idxs.update([i])
        isEmpty = len(idxs) == 0
        if isEmpty:
            return "No such element present. Look through the catalog table."

        R = [r for i, r in enumerate(self.catalog.data) if i in idxs]

        return CatalogTable(data=R, columns=self.catalog.columns)
