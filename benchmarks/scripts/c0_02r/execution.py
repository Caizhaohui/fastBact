from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


class GateError(RuntimeError):
    pass


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


@dataclass(frozen=True, slots=True)
class Command:
    argv: tuple[str, ...]
    cwd: str
    exit_code: int
    wall_seconds: float
    peak_rss_kb: int
    stdout: str
    stderr: str
    stdout_sha256: str
    stderr_sha256: str


@dataclass(frozen=True, slots=True)
class Source:
    path: str
    commit: str
    git_tree: str
    source_tree_sha256: str
    dirty_status: str


def source_identity(path: Path) -> Source:
    def git(*args: str) -> str:
        return subprocess.check_output(
            ("git", "-C", str(path), *args), text=True
        ).strip()

    status = git("status", "--porcelain=v1", "--untracked-files=all")
    identity = hashlib.sha256()
    for filename in sorted(git("ls-files", "-z").split("\0")):
        if filename:
            identity.update(
                filename.encode() + b"\0" + digest(path / filename).encode() + b"\0"
            )
    return Source(
        str(path),
        git("rev-parse", "HEAD"),
        git("rev-parse", "HEAD^{tree}"),
        identity.hexdigest(),
        status,
    )


def require_compute() -> None:
    if not os.environ.get("SLURM_JOB_ID") or not os.environ.get("SLURMD_NODENAME"):
        raise GateError(
            "Run under a Slurm compute allocation; login-node execution is forbidden"
        )


@dataclass(frozen=True, slots=True)
class Executor:
    directory: Path
    history: list[Command]

    def run(
        self, argv: tuple[str, ...], label: str, cwd: Path | None = None
    ) -> Command:
        require_compute()
        prefix = self.directory / label
        prefix.parent.mkdir(parents=True, exist_ok=True)
        stdout, stderr, rss = (
            Path(str(prefix) + ".stdout"),
            Path(str(prefix) + ".stderr"),
            Path(str(prefix) + ".rss"),
        )
        if stdout.exists() or stderr.exists():
            raise GateError(f"Refusing to overwrite command evidence: {prefix}")
        started = time.monotonic()
        with stdout.open("w") as out, stderr.open("w") as err:
            result = subprocess.run(
                ("/usr/bin/time", "-f", "%M", "-o", str(rss), *argv),
                cwd=cwd,
                stdout=out,
                stderr=err,
                check=False,
            )
        receipt = Command(
            argv,
            str(cwd or Path.cwd()),
            result.returncode,
            time.monotonic() - started,
            int(rss.read_text().strip().splitlines()[-1]) if rss.exists() else 0,
            str(stdout),
            str(stderr),
            digest(stdout),
            digest(stderr),
        )
        self.history.append(receipt)
        Path(str(prefix) + ".json").write_text(
            json.dumps(asdict(receipt), indent=2) + "\n"
        )
        if receipt.exit_code != 0:
            raise GateError(
                f"Command failed ({receipt.exit_code}): {argv}; see {stderr}"
            )
        return receipt
