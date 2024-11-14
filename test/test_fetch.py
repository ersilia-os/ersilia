import unittest
from ersilia.hub.fetch.fetch import ModelFetcher

class TestModelFetcher(unittest.TestCase):
    
    def test_from_dir_priority_when_repo_path_is_set(self):
        fetcher = ModelFetcher(from_dir="test_repo_path")
        self.assertEqual(fetcher.repo_path, "test_repo_path", "repo_path should match the value provided to from_dir")

    def test_from_dir_fallback_when_repo_path_is_none(self):
        fetcher = ModelFetcher(from_dir="fallback_dir")
        self.assertEqual(fetcher.repo_path, "fallback_dir", "from_dir should be used as the primary path")

    def test_no_repo_path_or_from_dir(self):
        fetcher = ModelFetcher(from_dir=None)
        self.assertIsNone(fetcher.repo_path, "repo_path should be None if from_dir is None")

if __name__ == '__main__':
    unittest.main()
