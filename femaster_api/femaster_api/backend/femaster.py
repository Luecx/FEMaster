"""High-level FEMaster process facade."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from femaster_api.backend.result_reader import Result, ResultReader
from femaster_api.export.femaster_writer import FEMasterWriter
from femaster_api.model.model import Model


@dataclass(frozen=True, slots=True)
class ProcessResult:
    input_path: Path
    output_path: Path
    returncode: int
    stdout: str
    stderr: str

    @property
    def ok(self) -> bool:
        return self.returncode == 0


class FEMaster:
    """Write, run, and read a model through a FEMaster executable."""

    def __init__(self, executable: str | Path = "FEMaster", *, ncpus: int = 1, model: Model | None = None) -> None:
        self.executable = str(executable)
        self.ncpus = int(ncpus)
        self.model = model
        self.last_process: ProcessResult | None = None

    def set_model(self, model: Model) -> "FEMaster":
        if not isinstance(model, Model):
            raise TypeError("model must be a Model")
        self.model = model
        return self

    def write(self, path: str | Path, *, validate: bool = True) -> Path:
        model = self._require_model()
        if validate:
            model.validate(raise_on_error=True)
        output = Path(path)
        FEMasterWriter(model).write(output)
        return output

    def read_res(self, path: str | Path) -> Result:
        return ResultReader().read(path)

    def run(
        self,
        input_path: str | Path | None = None,
        *,
        output_path: str | Path | None = None,
        cwd: str | Path | None = None,
        extra_args: Sequence[str] = (),
        validate: bool = True,
        check: bool = True,
    ) -> Result:
        model = self._require_model()
        inp = Path(input_path) if input_path is not None else Path(f"{model.name}.inp")
        out = Path(output_path) if output_path is not None else inp.with_suffix(".res")
        self.write(inp, validate=validate)
        self.last_process = self._run_deck(inp, out, cwd=cwd, extra_args=extra_args, check=check)
        return self.read_res(out)

    def _run_deck(
        self,
        input_path: Path,
        output_path: Path,
        *,
        cwd: str | Path | None,
        extra_args: Sequence[str],
        check: bool,
    ) -> ProcessResult:
        command = [
            self.executable,
            str(input_path),
            "--output",
            str(output_path),
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
        result = ProcessResult(input_path, output_path, completed.returncode, completed.stdout, completed.stderr)
        if check and not result.ok:
            raise RuntimeError(f"FEMaster failed with code {result.returncode}\n{result.stderr}")
        return result

    def _require_model(self) -> Model:
        if self.model is None:
            raise ValueError("no model has been set")
        return self.model
