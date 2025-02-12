import os
import sqlite3

from ... import ErsiliaBase

SLUGDB_FILE = ".slug.db"


class SlugDb(ErsiliaBase):
    """
    Manages slug database operations for models.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for initializing the slug database.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.file_path = os.path.join(self.eos_dir, SLUGDB_FILE)
        self._table = "slugs"
        self.create_table()

    def _connect(self):
        return sqlite3.connect(self.file_path)

    def create_table(self):
        """
        Creates the slugs table in the database if it does not exist.
        """
        if self._table is None:
            return
        sql = """
        CREATE TABLE IF NOT EXISTS {0} (
            model_id text,
            slug text,
            PRIMARY KEY (model_id, slug)
        );
        """.format(self._table)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def insert(self, model_id, slug):
        """
        Inserts a model ID and slug into the database.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        slug : str
            Slug associated with the model.
        """
        if self._table is None:
            return
        sql = """
        INSERT OR IGNORE INTO {0} (model_id, slug) VALUES ('{1}', '{2}')
        """.format(self._table, model_id, slug)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def delete_by_model_id(self, model_id):
        """
        Deletes entries from the database by model ID.

        Parameters
        ----------
        model_id : str
            Identifier of the model to delete.
        """
        if self._table is None:
            return
        sql = """
        DELETE FROM {0}
            WHERE model_id = '{1}'
        """.format(self._table, model_id)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def delete_by_slug(self, slug):
        """
        Deletes entries from the database by slug.

        Parameters
        ----------
        slug : str
            Slug to delete.
        """
        if self._table is None:
            return
        sql = """
        DELETE FROM {0}
            WHERE slug = '{1}'
        """.format(self._table, slug)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def delete(self, model_id, slug):
        """
        Deletes a specific entry from the database by model ID and slug.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        slug : str
            Slug associated with the model.
        """
        if self._table is None:
            return
        sql = """
        DELETE FROM {0}
            WHERE model_id = '{1}' AND slug = '{2}'
        """.format(self._table, model_id, slug)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def models_of_slug(self, slug):
        """
        Retrieves model IDs associated with a specific slug.

        Parameters
        ----------
        slug : str
            Slug to search for.

        Returns
        -------
        set
            Set of model IDs associated with the slug.
        """
        sql = """
        SELECT model_id FROM {0}
            WHERE slug = '{1}'
        """.format(self._table, slug)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = {x[0] for x in c.fetchall()}
        conn.close()
        return res

    def slugs_of_model(self, model_id):
        """
        Retrieves slugs associated with a specific model ID.

        Parameters
        ----------
        model_id : str
            Identifier of the model to search for.

        Returns
        -------
        set
            Set of slugs associated with the model ID.
        """
        sql = """
        SELECT slug FROM {0}
            WHERE model_id = '{1}'
        """.format(self._table, model_id)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = {x[0] for x in c.fetchall()}
        conn.close()
        return res

    def clean(self):
        """
        Cleans the database by deleting all entries.
        """
        if self._table is None:
            return
        sql = """
        DELETE FROM {0}
        """.format(self._table)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()
