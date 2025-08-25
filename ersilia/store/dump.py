import hashlib
import json

from redis import ConnectionError, Redis

from ..default import REDIS_EXPIRATION, REDIS_HOST, REDIS_PORT

# ruff: noqa: D101, D102


class DumpLocalCache:
    def conn_redis(self):
        return Redis(host=REDIS_HOST, port=REDIS_PORT, decode_responses=True)

    def init_redis(self):
        global redis_client
        try:
            redis_client = self.conn_redis()
            return True
        except ConnectionError:
            redis_client = None
            return False

    def generate_redis_key(self, raw_string):
        return hashlib.md5(raw_string.encode()).hexdigest()

    def fetch_cached_results(self, model_id, data, dim):
        hash_key = f"cache:{model_id}"
        fields = [
            item["input"] if isinstance(item, dict) and "input" in item else item
            for item in data
        ]
        raw = redis_client.hmget(hash_key, fields)
        results, missing = [], []
        for item, val in zip(data, raw):
            if val:
                results.append(json.loads(val))
            else:
                results.append([""] * dim)
                missing.append(item)
        return results, missing

    def cache_missing_results(self, model_id, missing_inputs, computed_results):
        hash_key = f"cache:{model_id}"
        pipe = redis_client.pipeline()
        for item, result in zip(missing_inputs, computed_results):
            field = (
                item["input"] if isinstance(item, dict) and "input" in item else item
            )
            pipe.hset(hash_key, field, json.dumps(result))
        pipe.expire(hash_key, REDIS_EXPIRATION)
        pipe.execute()

    def fetch_all_cached(self, model_id, dtype, cols):
        def dict_to_lists(d):
            keys = list(d.keys())
            values = list(d.values())
            return keys, values

        hash_key = f"cache:{model_id}"
        header = self.fetch_or_cache_header(model_id)
        header = header or cols
        assert header is not None, (
            "Headers can not be empty! This might happened either the header is not cached or resolved from model schema."
        )
        raw = redis_client.hgetall(hash_key)
        results = {field: json.loads(val) for field, val in raw.items()}
        inputs, results = dict_to_lists(results)
        if inputs and not isinstance(inputs[0], dict):
            inputs = [
                {"input": input, "key": hashlib.md5(input.encode()).hexdigest()}
                for input in inputs
            ]
        results = self.orient_to_json(results, header, inputs, "records", dtype)
        return results, inputs

    def fetch_or_cache_header(self, model_id, computed_headers=None):
        header_key = f"{model_id}:header"
        cached = redis_client.get(header_key)
        if cached:
            return json.loads(cached) if isinstance(cached, str) else cached
        if computed_headers is not None:
            redis_client.setex(
                header_key, REDIS_EXPIRATION, json.dumps(computed_headers)
            )
            return computed_headers
        return None

    def get_cached(self, model_id, data, dtype, cols, computed_headers=None):
        if not self.init_redis():
            return [], None

        results, missing = self.fetch_cached_results(model_id, data, len(cols))
        header = self.fetch_or_cache_header(model_id, computed_headers)
        header = header or cols
        assert header is not None, (
            "Headers can not be empty! This might happened either the header is not cached or resolved from model schema."
        )
        results = self.orient_to_json(results, header, data, "records", dtype)
        return results, missing

    def list_cached_inputs(self, model_id):
        hash_key = f"cache:{model_id}"
        return redis_client.hkeys(hash_key)

    def make_hashable(self, obj):
        if isinstance(obj, list):
            return tuple(self.make_hashable(x) for x in obj)
        elif isinstance(obj, dict):
            return tuple(sorted((k, self.make_hashable(v)) for k, v in obj.items()))
        return obj

    def orient_to_json(self, values, columns, index, orient, output_type):
        if len(output_type) > 1:
            output_type = "string"
        else:
            output_type = output_type[0].lower()

        def convert_value(x):
            if output_type == "string":
                return str(x)
            elif output_type == "float":
                if isinstance(x, str):
                    return float(x) if x != "" else None
                elif isinstance(x, (int, float)):
                    return float(x)
                elif isinstance(x, list) and x:
                    return float(x[0])
                else:
                    return None
            elif output_type == "integer":
                if isinstance(x, str):
                    return int(x) if x != "" else None
                elif isinstance(x, (int, float)):
                    return int(x)
                elif isinstance(x, list) and x:
                    return int(x[0])
                else:
                    return None
            return x

        if values and isinstance(values[0], list):
            serialized = [[convert_value(x) for x in row] for row in values]
        else:
            serialized = [convert_value(x) for x in values]

        if orient == "split":
            return {"columns": columns, "index": index, "data": serialized}
        elif orient == "records":
            return [dict(zip(columns, row)) for row in serialized]
        elif orient == "index":
            return {idx: dict(zip(columns, row)) for idx, row in zip(index, serialized)}
        elif orient == "columns":
            data = {}
            for j, col in enumerate(columns):
                col_data = {}
                for i, idx in enumerate(index):
                    col_data[self.make_hashable(idx)] = serialized[i][j]
                data[col] = col_data
            return data
        elif orient == "values":
            return serialized
        return None

    def _standardize_output(self, input_data, results, output, meta, n_samples=None):
        def _parse(v):
            if isinstance(v, str):
                try:
                    return int(v) if v.isdigit() else float(v)
                except ValueError:
                    return v
            return v

        _results = []
        try:
            keys = meta.get("outcome") if meta and meta.get("outcome") else ["outcome"]
            if isinstance(results, tuple):
                results, _ = results
            if n_samples and results:
                results = results[:n_samples]
            for inp, out in zip(input_data, results):
                vals = [_parse(v) for v in out.values()]
                if output.endswith(".json"):
                    if len(keys) == len(vals):
                        out_dict = dict(zip(keys, vals))
                    else:
                        out_dict = {
                            keys[0]: vals
                            if len(vals) > 1
                            else (vals[0] if vals else None)
                        }
                else:
                    out_dict = {"outcome": vals[0] if len(vals) == 1 else vals}
                _results.append({"input": inp, "output": out_dict})
        except Exception as e:
            raise Exception(f"Something went wring when standardizing the output: {e}")
        return _results
