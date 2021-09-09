import sqlite3
import os
from ... import ErsiliaBase

SLUGDB_FILE = ".slug.db"


class SlugDb(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.file_path = os.path.join(self.eos_dir, SLUGDB_FILE)
        self._table = "slugs"
        self.create_table()

    def _connect(self):
        return sqlite3.connect(self.file_path)

    def create_table(self):
        if self._table is None:
            return
        sql = """
        CREATE TABLE IF NOT EXISTS {0} (
            model_id text,
            slug text,
            PRIMARY KEY (model_id, slug)
        );
        """.format(
            self._table
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def insert(self, model_id, slug):
        if self._table is None:
            return
        sql = """
        INSERT OR IGNORE INTO {0} (model_id, slug) VALUES ('{1}', '{2}')
        """.format(
            self._table, model_id, slug
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def delete_by_model_id(self, model_id):
        if self._table is None:
            return
        sql = """
        DELETE FROM {0}
            WHERE model_id = '{1}'
        """.format(
            self._table, model_id
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def delete_by_slug(self, slug):
        if self._table is None:
            return
        sql = """
        DELETE FROM {0}
            WHERE slug = '{1}'
        """.format(
            self._table, slug
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def delete(self, model_id, slug):
        if self._table is None:
            return
        sql = """
        DELETE FROM {0}
            WHERE model_id = '{1}' AND slug = '{2}'
        """.format(
            self._table, model_id, slug
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        conn.commit()
        conn.close()

    def models_of_slug(self, slug):
        sql = """
        SELECT model_id FROM {0}
            WHERE slug = '{1}'
        """.format(
            self._table, slug
        )
        conn = self._connect()
        c = conn.cursor()
        c.execute(sql)
        res = {x[0] for x in c.fetchall()}
        conn.close()
        return res

    def slugs_of_model(self, model_id):
        sql = """
        SELECT slug FROM {0}
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
