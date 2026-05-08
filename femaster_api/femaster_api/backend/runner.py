"""High-level FEMaster subprocess execution."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from femaster_api.backend.result_reader import Result, ResultReader
from femaster_api.export.femaster_writer import FEMasterWriter
from femaster_api.model.model import Model


@dataclass(frozen=True, slots=True)
class RunResult:
    """Low-level subprocess result for advanced callers."""

    input_path: Path
    result_path: Path
    returncode: int
    stdout: str
    stderr: str

    @property
    def ok(self) -> bool:
        """Return true when FEMaster exited successfully."""

        return self.returncode == 0


class FEMaster:
    """Run FEMaster for a model and return parsed Python results."""

    def __init__(
        self,
        model: Model | None = None,
        *,
        executable: str | Path | None = None,
        threads: int = 1,
        workdir: str | Path | None = None,
        input_path: str | Path | None = None,
        result_path: str | Path | None = None,
        temporary_files: bool = True,
        redirect_stdout: bool = True,
        extra_args: Sequence[str] = (),
    ) -> None:
        self.model           = model
        self.executable      = find_femaster_executable(executable)
        self.threads         = int(threads)
        self.workdir         = None if workdir is None else Path(workdir)
        self.input_path      = None if input_path is None else Path(input_path)
        self.result_path     = None if result_path is None else Path(result_path)
        self.temporary_files = bool(temporary_files)
        self.redirect_stdout = bool(redirect_stdout)
        self.extra_args      = tuple(extra_args)

        self.last_run: RunResult | None = None

    def run(
        self,
        model: Model | None = None,
        *,
        validate: bool = True,
        check: bool = True,
    ) -> Result:
        """Write, solve, read, and return a parsed result object."""

        active_model = model or self.model
        if active_model is None:
            raise ValueError("FEMaster.run() needs a model")
        if validate:
            active_model.validate(raise_on_error=True)

        if self._uses_temporary_directory:
            with tempfile.TemporaryDirectory(prefix="femaster_") as directory:
                return self._run_in_directory(active_model, Path(directory), check=check)

        directory = self.workdir or Path.cwd()
        directory.mkdir(parents=True, exist_ok=True)
        return self._run_in_directory(active_model, directory, check=check)

    def run_deck(
        self,
        input_path: str | Path,
        *,
        result_path: str | Path | None = None,
        cwd: str | Path | None = None,
        check: bool = True,
    ) -> RunResult:
        """Run an existing deck and return the low-level subprocess result."""

        inp = Path(input_path)
        res = Path(result_path) if result_path is not None else inp.with_suffix(".res")
        command = self._command(inp, res)

        completed = subprocess.run(
            command,
            cwd=None if cwd is None else str(cwd),
            text=True,
            capture_output=self.redirect_stdout,
            check=False,
        )

        stdout = completed.stdout if completed.stdout is not None else ""
        stderr = completed.stderr if completed.stderr is not None else ""
        run_result = RunResult(inp, res, completed.returncode, stdout, stderr)
        self.last_run = run_result

        if check and not run_result.ok:
            raise RuntimeError(_failure_message(run_result))
        return run_result

    def run_model(
        self,
        model: Model,
        input_path: str | Path,
        *,
        result_path: str | Path | None = None,
        cwd: str | Path | None = None,
        validate: bool = True,
        check: bool = True,
    ) -> Result:
        """Compatibility wrapper that writes explicit files and returns results."""

        if validate:
            model.validate(raise_on_error=True)
        FEMasterWriter(model).write(input_path)
        run_result = self.run_deck(input_path, result_path=result_path, cwd=cwd, check=check)
        return ResultReader().read(run_result.result_path)

    @property
    def _uses_temporary_directory(self) -> bool:
        return self.temporary_files and (self.input_path is None or self.result_path is None)

    def _run_in_directory(self, model: Model, directory: Path, *, check: bool) -> Result:
        inp = self.input_path or directory / f"{_file_stem(model.name)}.inp"
        res = self.result_path or directory / f"{_file_stem(model.name)}.res"

        inp.parent.mkdir(parents=True, exist_ok=True)
        res.parent.mkdir(parents=True, exist_ok=True)

        FEMasterWriter(model).write(inp)
        run_result = self.run_deck(inp, result_path=res, cwd=self.workdir, check=check)
        return ResultReader().read(run_result.result_path)

    def _command(self, input_path: Path, result_path: Path) -> list[str]:
        return [
            str(self.executable),
            str(input_path),
            "--output",
            str(result_path),
            "--ncpus",
            str(self.threads),
            *self.extra_args,
        ]


FEMasterRunner = FEMaster


def find_femaster_executable(executable: str | Path | None = None) -> Path:
    """Resolve a FEMaster executable from an explicit path, env var, local bin, or PATH."""

    candidates: list[str | Path] = []
    if executable is not None:
        resolved = _resolve_executable(executable)
        if resolved is not None:
            return resolved
        raise FileNotFoundError(f"Could not find the FEMaster executable at: {executable}")
    for variable in ("FEMASTER_EXECUTABLE", "FEMASTER_PATH"):
        value = os.environ.get(variable)
        if value:
            candidates.append(value)
    candidates.extend(_local_bin_candidates())
    candidates.extend(("FEMaster", "femaster", "FEMaster.exe", "femaster.exe"))

    for candidate in candidates:
        resolved = _resolve_executable(candidate)
        if resolved is not None:
            return resolved

    raise FileNotFoundError(
        "Could not find the FEMaster executable. Pass executable=..., set "
        "FEMASTER_EXECUTABLE, or put FEMaster on PATH."
    )


def _resolve_executable(candidate: str | Path) -> Path | None:
    path = Path(candidate)
    if path.is_dir():
        for name in ("FEMaster", "femaster", "FEMaster.exe", "femaster.exe"):
            executable = path / name
            if executable.exists():
                return executable.resolve()
    if path.exists():
        return path.resolve()
    found = shutil.which(str(candidate))
    return None if found is None else Path(found).resolve()


def _local_bin_candidates() -> tuple[Path, ...]:
    roots = [Path.cwd(), Path(__file__).resolve()]
    candidates: list[Path] = []
    for root in roots:
        for parent in (root, *root.parents):
            candidates.extend(
                (
                    parent / "bin" / "FEMaster",
                    parent / "bin" / "femaster",
                    parent / "bin" / "FEMaster.exe",
                    parent / "bin" / "femaster.exe",
                )
            )
    return tuple(candidates)


def _failure_message(result: RunResult) -> str:
    lines = [f"FEMaster failed with code {result.returncode}"]
    if result.stderr:
        lines.append(result.stderr)
    if result.stdout:
        lines.append(result.stdout)
    return "\n".join(lines)


def _file_stem(name: str) -> str:
    clean = "".join(char if char.isalnum() or char in ("-", "_") else "_" for char in name.strip())
    return clean or "model"
