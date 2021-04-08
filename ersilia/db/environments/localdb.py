import sqlite3
import os
from ... import ErsiliaBase

ENVIRONMENTDB_FILE = ".environment.db"


class EnvironmentDb(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.file_path = os.path.join(self.eos_dir, ENVIRONMENTDB_FILE)
        self._table = None

    @property
    def table(self):
        return self._table

    @table.setter
    def table(self, table):
        self._table = table
        self.create_table()

    @table.deleter
    def table(self):
        del self._table

    def _connect(self):
        return sqlite3.connect(self.file_path)

    def create_table(self):
        if self._table is None:
            return
        sql = """
        CREATE TABLE IF NOT EXISTS {0} (
            model_id text,
            env text,
            PRIMARY KEY (model_id, env)
        );
        """.format(
            self._table
        )
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
        if self._table is None:
            return
        sql = """
        INSERT OR IGNORE INTO {0} (model_id, env) VALUES ('{1}', '{2}')
        """.format(
            self._table, model_id, env
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def delete(self, model_id, env):
        if self._table is None:
            return
        sql = """
        DELETE FROM {0}
            WHERE model_id = '{1}' AND env = '{2}'
        """.format(
            self._table, model_id, env
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def envs_of_model(self, model_id):
        if self._table is None:
            return
        sql = """
        SELECT env FROM {0}
            WHERE model_id = '{1}'
        """.format(
            self._table, model_id
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = {x[0] for x in c.fetchall()}
        conn.close()
        return res

    def models_of_env(self, env):
        if self._table is None:
            return
        sql = """
        SELECT model_id FROM {0}
            WHERE env = '{1}'
        """.format(
            self._table, env
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = {x[0] for x in c.fetchall()}
        conn.close()
        return res

    def models_with_same_env(self, model_id):
        if self._table is None:
            return
        sql = """
        SELECT model_id FROM {0}
            WHERE env IN (SELECT env FROM {0} WHERE model_id = '{1}')
        """.format(
            self._table, model_id
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = {x[0] for x in c.fetchall()}
        conn.close()
        return res

    def envs_with_same_model(self, env):
        if self._table is None:
            return
        sql = """
        SELECT env FROM {0}
            WHERE model_id IN (SELECT model_id FROM {0} WHERE env = '{1}')
        """.format(
            self._table, env
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = {x[0] for x in c.fetchall()}
        conn.close()
        return res

    def fetchall(self):
        if self._table is None:
            return
        sql = """
        SELECT * FROM {0}
        """.format(
            self._table
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = list(c.fetchall())
        conn.close()
        return res

    def clean(self):
        if self._table is None:
            return
        sql = """
        DELETE FROM {0}
        """.format(
            self._table
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()
