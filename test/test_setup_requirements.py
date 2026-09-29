from unittest.mock import patch

import pytest

from ersilia.setup.requirements.git import GitLfsRequirement


def test_git_lfs_not_installed_raises_error(capsys):
    with patch(
        "ersilia.setup.requirements.git.run_command_check_output",
        return_value="command not found",
    ):
        req = GitLfsRequirement()
        with pytest.raises(SystemExit):
            req.is_installed(install_if_necessary=False)
        captured = capsys.readouterr()
        assert "Git LFS is not installed" in captured.err


def test_git_lfs_installed_returns_true():
    with patch(
        "ersilia.setup.requirements.git.run_command_check_output",
        return_value="git-lfs/3.4.0 (GitHub; linux amd64; go 1.21.0)",
    ):
        req = GitLfsRequirement()
        result = req.is_installed(install_if_necessary=False)
        assert result is True
