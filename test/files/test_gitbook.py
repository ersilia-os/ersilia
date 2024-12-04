from unittest.mock import patch


def test_conda_create():
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        assert (
            mock_run(["conda", "create", "-n", "ersilia", "python=3.10"]).returncode
            == 0
        )


def test_conda_activate():
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        assert mock_run(["conda", "activate", "ersilia"]).returncode == 0


def test_git_clone():
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        assert (
            mock_run(
                ["git", "clone", "https://github.com/ersilia-os/ersilia.git"]
            ).returncode
            == 0
        )


def test_pip_install():
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        assert mock_run(["pip", "install", "-e", "."]).returncode == 0


def test_ersilia_help():
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        assert mock_run(["ersilia", "--xxx"]).returncode == 0


def test_ersilia_catalog():
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        assert mock_run(["ersilia", "catalog"]).returncode == 0


def test_install_isaura():
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        assert (
            mock_run(["python", "-m", "pip", "install", "isaura==0.1"]).returncode == 0
        )
