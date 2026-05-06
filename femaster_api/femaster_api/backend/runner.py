"""Subprocess runner for FEMaster jobs."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from femaster_api.export.femaster_writer import FEMasterWriter
from femaster_api.model.model import Model


@dataclass(frozen=True, slots=True)
class RunResult:
    input_path: Path
    output_path: Path
    returncode: int
    stdout: str
    stderr: str

    @property
    def ok(self) -> bool:
        return self.returncode == 0


class FEMasterRunner:
    """Write a deck and run a FEMaster executable."""

    def __init__(self, executable: str | Path = "FEMaster", *, ncpus: int = 1) -> None:
        self.executable = str(executable)
        self.ncpus = int(ncpus)

    def run_deck(
        self,
        input_path: str | Path,
        *,
        output_path: str | Path | None = None,
        cwd: str | Path | None = None,
        extra_args: Sequence[str] = (),
        check: bool = True,
    ) -> RunResult:
        inp = Path(input_path)
        out = Path(output_path) if output_path is not None else inp.with_suffix(".res")
        command = [
            self.executable,
            str(inp),
            "--output",
            str(out),
            "--ncpus",
            str(self.ncpus),
            *extra_args,
        ]
        completed = subprocess.run(
            command,
            cwd=None if cwd is None else str(cwd),
            text=True,
            capture_output=True,
            check=False,
        )
        result = RunResult(inp, out, completed.returncode, completed.stdout, completed.stderr)
        if check and not result.ok:
            raise RuntimeError(f"FEMaster failed with code {result.returncode}\n{result.stderr}")
        return result

    def run_model(
        self,
        model: Model,
        input_path: str | Path,
        *,
        output_path: str | Path | None = None,
        cwd: str | Path | None = None,
        validate: bool = True,
        check: bool = True,
    ) -> RunResult:
        if validate:
            model.validate(raise_on_error=True)
        FEMasterWriter(model).write(input_path)
        return self.run_deck(input_path, output_path=output_path, cwd=cwd, check=check)
