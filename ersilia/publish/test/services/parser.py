import os
import yaml
import re
import warnings
import shutil
import shlex
from pathlib import Path
import subprocess

from ....utils.terminal import run_command
from ....utils.logging import logger
from ....default import _CONDA_BOOTSTRAP

BASE = "base"
TEST_BASE = "test"
PYTHON_EXE = "python"


def eval_conda_prefix() -> str:
    logger.info(f"[eval_conda_prefix] start")
    conda_exe = os.environ.get("CONDA_EXE") or shutil.which("conda")
    logger.info(f"[eval_conda_prefix] CONDA_EXE={os.environ.get('CONDA_EXE')}")
    logger.info(f"[eval_conda_prefix] which conda={shutil.which('conda')}")
    logger.info(f"[eval_conda_prefix] conda_exe={conda_exe}")

    if not conda_exe:
        env_prefix = os.environ.get("CONDA_PREFIX")
        logger.info(f"[eval_conda_prefix] CONDA_PREFIX={env_prefix}")
        if env_prefix:
            base_candidate = Path(env_prefix).parents[1]
            candidate = base_candidate / "bin" / "conda"
            logger.info(f"[eval_conda_prefix] candidate={candidate}")
            if candidate.exists():
                conda_exe = str(candidate)
                logger.info(f"[eval_conda_prefix] using candidate conda_exe={conda_exe}")

    if not conda_exe:
        logger.info(f"[eval_conda_prefix] no conda executable found -> ''")
        return ""

    logger.info(f"[eval_conda_prefix] running: {conda_exe} info --base")
    p = subprocess.run([conda_exe, "info", "--base"], capture_output=True, text=True)
    logger.info(f"[eval_conda_prefix] returncode={p.returncode}")
    logger.info(f"[eval_conda_prefix] stdout={(p.stdout or '').strip()}")
    logger.info(f"[eval_conda_prefix] stderr={(p.stderr or '').strip()}")

    if p.returncode != 0:
        logger.info(f"[eval_conda_prefix] conda info failed -> ''")
        return ""

    out = (p.stdout or "").strip()
    logger.info(f"[eval_conda_prefix] base={out}")
    return out


def get_conda_source(env):
    logger.info(f"[get_conda_source] env={env} bootstrap_len={len(_CONDA_BOOTSTRAP or '')}")
    lines = [
        "set -euo pipefail",
        _CONDA_BOOTSTRAP.strip(),
        'command -v conda >/dev/null 2>&1 || { echo "ERROR: conda not available after bootstrap" >&2; exit 127; }',
        "unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER || true",
        "unset CONDA_SHLVL || true",
        r"""while IFS='=' read -r k _; do
  if [[ "$k" =~ ^CONDA_PREFIX_[0-9]+$ ]]; then
    unset "$k"
  fi
done < <(env)""",
        "export CONDA_NO_PLUGINS=true",
        f'conda activate "{env}"',
    ]
    logger.info(f"[get_conda_source] lines_count={len(lines)} last_line={lines[-1]}")
    return lines


def run_bash(script):
    logger.info(f"[run_bash] script_len={len(script or '')}")
    proc = subprocess.run(["bash", "-lc", script], text=True, capture_output=True)
    logger.info(f"[run_bash] returncode={proc.returncode}")
    logger.info(f"[run_bash] stdout_tail={((proc.stdout or '')[-800:]).strip()}")
    logger.info(f"[run_bash] stderr_tail={((proc.stderr or '')[-800:]).strip()}")
    return proc


def build_python_executable_script(env):
    logger.info(f"[build_python_executable_script] env={env}")
    lines = get_conda_source(env) + [r'python -c "import sys; print(sys.executable)"']
    script = "\n".join(lines)
    logger.info(f"[build_python_executable_script] script_len={len(script)} head={script[:400]}")
    return script


def get_python_executable_from_env(env):
    logger.info(f"[get_python_executable_from_env] env={env}")
    proc = run_bash(build_python_executable_script(env))
    if proc.returncode != 0:
        logger.info(f"[get_python_executable_from_env] FAILED env={env}")
        return "", proc
    out = (proc.stdout or "").strip()
    logger.info(f"[get_python_executable_from_env] python_exe={out}")
    return out, proc


def conda_python_executable(env=None):
    env = env or BASE
    logger.info(f"[conda_python_executable] start env={env}")

    python_path, proc = get_python_executable_from_env(env)
    logger.info(f"[conda_python_executable] primary python_path={python_path}")
    if python_path:
        return python_path

    python_path, proc = get_python_executable_from_env(TEST_BASE)
    logger.info(f"[conda_python_executable] fallback env={TEST_BASE} python_path={python_path}")
    if python_path:
        return python_path

    logger.info(
        f"[conda_python_executable] failed env='{env}' stdout={(proc.stdout or '').strip()} "
        f"stderr={(proc.stderr or '').strip()}"
    )
    return PYTHON_EXE


NATIVE = frozenset(
    {
        "apt",
        "apt-get",
        "apt-cache",
        "apt-key",
        "curl",
        "wget",
        "ls",
        "cd",
        "cat",
        "echo",
        "touch",
        "mkdir",
        "rm",
        "cp",
        "export",
        "mv",
        "bash",
        "sh",
        "sudo",
        "python",
        "python3",
        "pip",
        "conda",
        "from",
        "workdir",
        "copy",
        "maintainer",
    }
)


def get_native():
    return NATIVE


class InstallParser:
    def __init__(self, file_name, conda_env_name=None):
        self.conda_env_name = conda_env_name
        self.file_name = file_name
        logger.info(f"[InstallParser.__init__] file_name={file_name} conda_env_name={conda_env_name}")
        self.python_version = self._get_python_version()
        logger.info(f"[InstallParser.__init__] python_version={self.python_version}")

    def _get_python_version(self):
        raise NotImplementedError("Implement this in subclass")

    def _get_commands(self):
        raise NotImplementedError("Implement this in subclass")

    def get_python_exe(self):
        logger.info(f"[InstallParser.get_python_exe] conda_env_name={self.conda_env_name}")
        return conda_python_executable(self.conda_env_name)

    @staticmethod
    def _has_conda(commands):
        for cmd in commands:
            if isinstance(cmd, list):
                if cmd and str(cmd[0]).strip().lower() == "conda":
                    return True
                s = " ".join(map(str, cmd)).lstrip().lower()
                if s.startswith("conda "):
                    return True
            else:
                s = str(cmd).lstrip().lower()
                if s.startswith("conda "):
                    return True
        return False

    def _is_valid_url(self, url):
        pattern = re.compile(r"^(git\+https://|git\+ssh://|https://).*")
        ok = bool(pattern.match(url))
        logger.info(f"[InstallParser._is_valid_url] url={url} ok={ok}")
        return ok

    def _convert_pip_entry_to_bash(self, command):
        logger.info(f"[InstallParser._convert_pip_entry_to_bash] command={command}")

        if isinstance(command, str):
            s = command.strip()
            parts = shlex.split(s)
            if len(parts) >= 2 and parts[0] in ("pip", "pip3") and parts[1] == "install":
                out = s if parts[0] == "pip" else "pip " + " ".join(parts[1:])
                logger.info(f"[InstallParser._convert_pip_entry_to_bash] out={out}")
                return out
            logger.info(f"[InstallParser._convert_pip_entry_to_bash] out={s}")
            return s

        if len(command) == 2:
            pkg = command[1]
            if pkg.startswith("git+"):
                if not self._is_valid_url(pkg):
                    raise ValueError("Invalid Git URL provided")
                out = f"pip install {pkg}"
                logger.info(f"[InstallParser._convert_pip_entry_to_bash] out={out}")
                return out
            raise ValueError("pip install entry must include version")

        pkg, ver = command[1], command[2]
        spec = pkg if (ver == "" or "git+" in pkg) else f"{pkg}=={ver}"
        flags = command[2:] if "git+" in pkg else command[3:]
        out = f"pip install {spec}" + (" " + " ".join(flags) if flags else "")
        logger.info(f"[InstallParser._convert_pip_entry_to_bash] out={out}")
        return out

    def _convert_conda_entry_to_bash(self, command):
        logger.info(f"[InstallParser._convert_conda_entry_to_bash] command={command}")

        if isinstance(command, str):
            out = command.strip()
            logger.info(f"[InstallParser._convert_conda_entry_to_bash] out={out}")
            return out

        if len(command) >= 4 and command[1] != "install":
            _, pkg, ver, *rest = command
            channels = [x for x in rest if x not in ("-y",)]
            flags = [x for x in rest if x == "-y"]
            if not channels:
                warnings.warn(f"No channel specified for conda package '{pkg}'")
                channels = ["default"]
            channel_flags = []
            for ch in channels:
                channel_flags.extend(["-c", ch])
            cmd = ["conda", "install"] + flags + channel_flags + [f"{pkg}={ver}"]
            if "-y" not in flags:
                cmd.append("-y")
            out = " ".join(cmd)
            logger.info(f"[InstallParser._convert_conda_entry_to_bash] out={out}")
            return out

        parts = command[1:]
        cmd = ["conda", "install"]
        channels = []
        flags = []
        pkg_spec = None
        i = 0

        while i < len(parts):
            p = parts[i]
            if p in ("-c", "--channel") and i + 1 < len(parts):
                channels += [parts[i], parts[i + 1]]
                i += 2
            elif p == "-y":
                flags.append("-y")
                i += 1
            else:
                pkg_spec = p
                i += 1

        if not pkg_spec:
            raise ValueError("No conda package specified")

        cmd += flags + channels + [pkg_spec]
        if "-y" not in flags:
            cmd.append("-y")

        out = " ".join(cmd)
        logger.info(f"[InstallParser._convert_conda_entry_to_bash] out={out}")
        return out

    def _prefix_unknown(self, raw):
        s = str(raw).strip()
        if not s:
            return ""
        head = s.split()[0]
        native = get_native()
        out = s if head.lower() in native else s
        logger.info(f"[InstallParser._prefix_unknown] raw={raw} out={out}")
        return out

    @staticmethod
    def _head_of_string(cmd_str):
        try:
            parts = shlex.split(cmd_str.strip())
        except ValueError:
            parts = cmd_str.strip().split()
        return parts[0].lower() if parts else ""

    def _convert_commands_to_bash_script(self):
        commands = self._get_commands()
        has_conda = self._has_conda(commands)
        env = self.conda_env_name or "base"
        python_exe = self.get_python_exe()

        logger.info(f"[InstallParser._convert_commands_to_bash_script] has_conda={has_conda} env={env}")
        logger.info(f"[InstallParser._convert_commands_to_bash_script] python_exe={python_exe}")
        logger.info(f"[InstallParser._convert_commands_to_bash_script] commands_count={len(commands)}")

        lines = get_conda_source(env) if has_conda else []

        def add_pip(x):
            pip_cmd = f"{python_exe} -m {self._convert_pip_entry_to_bash(x)}"
            lines.append(pip_cmd)
            logger.info(f"[InstallParser._convert_commands_to_bash_script] add_pip={pip_cmd}")

        for cmd in commands:
            if isinstance(cmd, list):
                if not cmd:
                    continue
                head = str(cmd[0]).lower()
                if head == "pip":
                    add_pip(cmd)
                elif head == "conda":
                    out = self._convert_conda_entry_to_bash(cmd)
                    lines.append(out)
                    logger.info(f"[InstallParser._convert_commands_to_bash_script] add_conda_list={out}")
                else:
                    out = self._prefix_unknown(" ".join(map(str, cmd)))
                    lines.append(out)
                    logger.info(f"[InstallParser._convert_commands_to_bash_script] add_other_list={out}")
                continue

            s = str(cmd).strip()
            if not s:
                continue
            head = self._head_of_string(s)
            if head in ("pip", "pip3"):
                add_pip(s)
            elif head == "conda":
                out = self._convert_conda_entry_to_bash(s)
                lines.append(out)
                logger.info(f"[InstallParser._convert_commands_to_bash_script] add_conda_str={out}")
            else:
                out = self._prefix_unknown(s)
                lines.append(out)
                logger.info(f"[InstallParser._convert_commands_to_bash_script] add_other_str={out}")

        script = os.linesep.join(lines)
        logger.info(f"[InstallParser._convert_commands_to_bash_script] script_len={len(script)}")
        logger.info(f"[InstallParser._convert_commands_to_bash_script] script_head={script[:500]}")
        return script

    def write_bash_script(self, file_name=None):
        if file_name is None:
            file_name = os.path.splitext(self.file_name)[0] + ".sh"
        logger.info(f"[InstallParser.write_bash_script] file_name={file_name}")
        data = self._convert_commands_to_bash_script()
        with open(file_name, "w") as f:
            f.write(data)
        logger.info(f"[InstallParser.write_bash_script] wrote_bytes={len(data)}")

    def _install_packages(self, file_path):
        cmd = f"bash {file_path}"
        logger.info(f"[InstallParser._install_packages] cmd={cmd}")
        run_command(cmd)
        logger.info(f"[InstallParser._install_packages] done cmd={cmd}")

    def check_file_exists(self):
        ok = os.path.exists(self.file_name)
        logger.info(f"[InstallParser.check_file_exists] file={self.file_name} exists={ok}")
        return ok


INSTALL_YML = "install.yml"


class YAMLInstallParser(InstallParser):
    def __init__(self, file_dir, conda_env_name=None):
        self.file_type = INSTALL_YML
        self.file_name = os.path.join(file_dir, self.file_type)
        logger.info(f"[YAMLInstallParser.__init__] file_name={self.file_name} conda_env_name={conda_env_name}")
        self.data = self._load_yaml()
        super().__init__(self.file_name, conda_env_name)

    def _load_yaml(self):
        logger.info(f"[YAMLInstallParser._load_yaml] reading {self.file_name}")
        with open(self.file_name, "r") as file:
            data = yaml.safe_load(file)
        logger.info(f"[YAMLInstallParser._load_yaml] loaded_type={type(data)} keys={list(data.keys()) if isinstance(data, dict) else None}")
        return data

    def _get_python_version(self):
        if "python" not in self.data or not isinstance(self.data["python"], str):
            raise ValueError("Python version must be a string")
        v = self.data["python"]
        logger.info(f"[YAMLInstallParser._get_python_version] python={v}")
        return v

    def _get_commands(self):
        if "commands" not in self.data:
            raise KeyError("Missing 'commands' key in YAML file")
        cmds = self.data["commands"]
        logger.info(f"[YAMLInstallParser._get_commands] commands_count={len(cmds) if isinstance(cmds, list) else 'n/a'}")
        return cmds


DOCKERFILE = "Dockerfile"
FILE_TYPE = "Dockerfile"


class DockerfileInstallParser(InstallParser):
    def __init__(self, file_dir, conda_env_name=None):
        file_name = os.path.join(file_dir, FILE_TYPE)
        logger.info(f"[DockerfileInstallParser.__init__] file_name={file_name} conda_env_name={conda_env_name}")
        super().__init__(file_name, conda_env_name)

    def _get_python_version(self):
        with open(self.file_name) as f:
            for line in f:
                if line.startswith("FROM"):
                    match = re.search(r"py(\d+\.\d+|\d{2,3})", line)
                    if match:
                        v = match.group(1)
                        if "." not in v:
                            v = f"{v[0]}.{v[1:]}"
                        logger.info(f"[DockerfileInstallParser._get_python_version] python={v}")
                        return v
        raise ValueError("Python version not found")

    @staticmethod
    def _is_flag(tok):
        return tok.startswith("-")

    @staticmethod
    def _is_git_spec(tok):
        return tok.startswith("git+")

    @staticmethod
    def _validate_pip_pkg_spec(tok):
        if DockerfileInstallParser._is_git_spec(tok):
            return
        if "==" not in tok:
            raise ValueError("pip install must pin versions (use '=='): got '{0}'".format(tok))
        pkg, ver = tok.split("==", 1)
        if not pkg or not ver:
            raise ValueError("Invalid pip pin: '{0}'".format(tok))

    @staticmethod
    def _validate_conda_pkg_spec(tok):
        if "=" not in tok:
            raise ValueError("conda install must pin versions (use '='): got '{0}'".format(tok))
        pkg, ver = tok.split("=", 1)
        if pkg.endswith("="):
            pkg = pkg[:-1]
            if ver.startswith("="):
                ver = ver[1:]
        if not pkg or not ver:
            raise ValueError("Invalid conda pin: '{0}'".format(tok))
        if not re.match(r"^[A-Za-z0-9._-]+$", pkg):
            raise ValueError("Invalid conda package name: '{0}'".format(pkg))

    @staticmethod
    def _validate_pip_install(parts):
        if len(parts) < 3:
            raise ValueError("pip install must include at least one package spec")
        i = 2
        while i < len(parts):
            tok = parts[i]
            if DockerfileInstallParser._is_flag(tok):
                if tok in ("-i", "--index-url", "--extra-index-url", "-f", "--find-links") and i + 1 < len(parts):
                    i += 2
                else:
                    i += 1
                continue
            DockerfileInstallParser._validate_pip_pkg_spec(tok)
            i += 1

    @staticmethod
    def _validate_conda_install(parts):
        if len(parts) < 3:
            raise ValueError("conda install must include at least one package spec")
        i = 2
        saw_pkg = False
        while i < len(parts):
            tok = parts[i]
            if tok in ("-c", "--channel"):
                if i + 1 >= len(parts):
                    raise ValueError("conda install: -c/--channel must have a value")
                i += 2
                continue
            if tok in ("-y", "--yes", "--quiet", "-q", "--freeze-installed"):
                i += 1
                continue
            if DockerfileInstallParser._is_flag(tok):
                i += 1
                continue
            saw_pkg = True
            DockerfileInstallParser._validate_conda_pkg_spec(tok)
            i += 1
        if not saw_pkg:
            raise ValueError("conda install must include at least one pinned package spec")

    @staticmethod
    def _maybe_validate_install_command(cmd):
        parts = shlex.split(cmd)
        if len(parts) >= 2 and parts[0] == "pip" and parts[1] == "install":
            DockerfileInstallParser._validate_pip_install(parts)
            return
        if len(parts) >= 2 and parts[0] == "conda" and parts[1] == "install":
            DockerfileInstallParser._validate_conda_install(parts)
            return

    def _get_commands(self):
        cmds = []
        with open(self.file_name) as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                if s.startswith("RUN "):
                    cmd = s[4:].strip()
                    self._maybe_validate_install_command(cmd)
                    cmds.append(cmd)
        logger.info(f"[DockerfileInstallParser._get_commands] commands_count={len(cmds)}")
        return cmds
