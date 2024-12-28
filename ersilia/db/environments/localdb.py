import os
import sqlite3

from ... import ErsiliaBase

ENVIRONMENTDB_FILE = ".environment.db"


class EnvironmentDb(ErsiliaBase):
    """
    Manages the storage and retrieval of model environments.

    Specifically it is essential for keeping track of the environments associated with each model.
    It provides a database where models and their corresponding environments (such as Docker images or conda environments) are stored and managed.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for initializing the environment database.


    Examples
    --------
    >>> env_db = EnvironmentDb(config_json)
    >>> env_db.table = "conda"
    >>> env_db.insert("model_id", "venv_name")
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.file_path = os.path.join(self.eos_dir, ENVIRONMENTDB_FILE)
        self._table = None

    @property
    def table(self):
        """
        Gets the current table name.

        Returns
        -------
        str
            The current table name.
        """
        return self._table

    @table.setter
    def table(self, table):
        """
        Sets the table name and creates the table if it does not exist.

        Parameters
        ----------
        table : str
            The name of the table to set.
        """
        self._table = table
        self.create_table()

    @table.deleter
    def table(self):
        """
        Deletes the table name.
        """
        del self._table

    def _connect(self):
        return sqlite3.connect(self.file_path)

    def create_table(self):
        """
        Create table if it does not exist.
        """
        if self._table is None:
            return
        sql = """
        CREATE TABLE IF NOT EXISTS {0} (
            model_id text,
            env text,
            PRIMARY KEY (model_id, env)
        );
        """.format(self._table)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def _fetch_tables(self):
        if self._table is None:
            return
        conn = self._connect()
        c = conn.cursor()
        c.execute('SELECT name FROM sqlite_master WHERE type = "table"')
        res = {x[0] for x in list(c.fetchall())}
        conn.close()
        return res

    def insert(self, model_id, env):
        """
        Inserts a model ID and environment into the database.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        env : str
            Environment associated with the model.
        """
        if self._table is None:
            return
        sql = """
        INSERT OR IGNORE INTO {0} (model_id, env) VALUES ('{1}', '{2}')
        """.format(self._table, model_id, env)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def delete(self, model_id, env):
        """
        Deletes a specific entry from the database by model ID and environment.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        env : str
            Environment associated with the model.
        """
        if self._table is None:
            return
        sql = """
        DELETE FROM {0}
            WHERE model_id = '{1}' AND env = '{2}'
        """.format(self._table, model_id, env)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def envs_of_model(self, model_id):
        """
        Retrieves environments associated with a specific model ID.

        Parameters
        ----------
        model_id : str
            Identifier of the model to search for.

        Returns
        -------
        set
            Set of environments associated with the model ID.
        """
        if self._table is None:
            return
        sql = """
        SELECT env FROM {0}
            WHERE model_id = '{1}'
        """.format(self._table, model_id)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = {x[0] for x in c.fetchall()}
        conn.close()
        return res

    def models_of_env(self, env):
        """
        Retrieves model IDs associated with a specific environment.

        Parameters
        ----------
        env : str
            Environment to search for.

        Returns
        -------
        set
            Set of model IDs associated with the environment.
        """
        if self._table is None:
            return
        sql = """
        SELECT model_id FROM {0}
            WHERE env = '{1}'
        """.format(self._table, env)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = {x[0] for x in c.fetchall()}
        conn.close()
        return res

    def models_with_same_env(self, model_id):
        """
        Retrieves model IDs that share the same environment as the given model ID.

        Parameters
        ----------
        model_id : str
            Identifier of the model to search for.

        Returns
        -------
        set
            Set of model IDs that share the same environment.
        """
        if self._table is None:
            return
        sql = """
        SELECT model_id FROM {0}
            WHERE env IN (SELECT env FROM {0} WHERE model_id = '{1}')
        """.format(self._table, model_id)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = {x[0] for x in c.fetchall()}
        conn.close()
        return res

    def envs_with_same_model(self, env):
        """
        Retrieves environments that share the same model as the given environment.

        Parameters
        ----------
        env : str
            Environment to search for.

        Returns
        -------
        set
            Set of environments that share the same model.
        """
        if self._table is None:
            return
        sql = """
        SELECT env FROM {0}
            WHERE model_id IN (SELECT model_id FROM {0} WHERE env = '{1}')
        """.format(self._table, env)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = {x[0] for x in c.fetchall()}
        conn.close()
        return res

    def fetchall(self):
        """
        Retrieves all entries from the current table.

        Returns
        -------
        list
            List of all entries in the table.
        """
        if self._table is None:
            return
        sql = """
        SELECT * FROM {0}
        """.format(self._table)
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = list(c.fetchall())
        conn.close()
        return res

    def clean(self):
        """
        Cleans the database by deleting all entries from the current table.
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
