from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Optional, Union, TYPE_CHECKING
from os import PathLike
import subprocess
import tempfile
import uuid

from .solution import Solution

if TYPE_CHECKING:
    from .model import Model

PathInput = Union[str, Path, PathLike[str]]


class Runner:
    class Engine(Enum):
        FEMASTER = "FEMASTER"
        ASAMI = "ASAMI"

    class Option(Enum):
        NO_TEMP_FILES = "no_temp_files"

    """Coordinates execution of a model with a selected engine."""

    def __init__(self) -> None:
        self._engine: Optional[Runner.Engine] = None
        self._model: Optional["Model"] = None
        self._engine_path: Optional[Path] = None
        self._options: set[Runner.Option] = set()

    def set_engine(self, engine: Runner.Engine, *, path: Optional[PathInput] = None) -> "Runner":
        """Select the engine used for the next run. Provide a path for FEMASTER."""
        if engine is Runner.Engine.FEMASTER and path is None and self._engine_path is None:
            raise ValueError("FEMASTER engine requires an executable path.")

        if path is not None:
            self._engine_path = self._normalize_path(path)

        self._engine = engine
        return self

    def set_engine_path(self, path: PathInput) -> "Runner":
        """Configure or update the engine's executable path."""
        self._engine_path = self._normalize_path(path)
        return self

    def set_option(self, option: Option, enabled: bool = True) -> "Runner":
        """Enable or disable runtime options such as keeping temp files."""
        if enabled:
            self._options.add(option)
        else:
            self._options.discard(option)
        return self

    def set_model(self, model: "Model") -> "Runner":
        """Attach the model that should be run."""
        self._model = model
        return self

    def run(
        self,
        engine: Optional[Engine] = None,
        model: Optional["Model"] = None,
    ) -> Solution:
        """Execute the selected engine with the provided model and return a solution."""
        if engine is not None:
            self.set_engine(engine)
        if model is not None:
            self.set_model(model)

        if self._engine is None:
            raise ValueError("No engine selected. Call set_engine() or pass an engine to run().")
        if self._model is None:
            raise ValueError("No model provided. Call set_model() or pass a model to run().")
        if self._engine is Runner.Engine.FEMASTER and self._engine_path is None:
            raise ValueError(
                "FEMASTER engine requires an executable path. "
                "Use set_engine(path=...) or set_engine_path() before run()."
            )

        if self._engine is Runner.Engine.FEMASTER:
            return self._run_femaster()

        raise NotImplementedError(f"Engine {self._engine.value} is not implemented yet.")

    def _run_femaster(self) -> Solution:
        assert self._model is not None
        assert self._engine_path is not None

        model_input = self._model.to_femaster()

        keep_files = Runner.Option.NO_TEMP_FILES in self._options
        tmp_context = None
        if keep_files:
            tmp_path = Path.cwd()
        else:
            tmp_context = tempfile.TemporaryDirectory(prefix="fembridge_femaster_")
            tmp_path = Path(tmp_context.name)

        try:
            base_name = f"model_{uuid.uuid4().hex}"
            inp_path = tmp_path / f"{base_name}.inp"
            inp_path.write_text(model_input, encoding="utf-8")

            cmd = [str(self._engine_path), str(inp_path)]
            completed = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)
            if completed.returncode != 0:
                stdout = completed.stdout.strip()
                stderr = completed.stderr.strip()
                message_parts = [part for part in (stdout, stderr) if part]
                message = "\n".join(message_parts)
                raise RuntimeError(
                    "FEMASTER execution failed with exit code "
                    f"{completed.returncode}.\n{message}"
                )

            res_path = inp_path.with_suffix(".res")
            if not res_path.exists():
                raise FileNotFoundError(
                    "FEMASTER execution did not create the expected .res file "
                    f"at {res_path}."
                )

            solution = Solution()
            solution.read_res(res_path, self._model.steps)
            if keep_files:
                solution.files = {"inp": inp_path.resolve(), "res": res_path.resolve()}
            return solution
        finally:
            if not keep_files and tmp_context is not None:
                tmp_context.cleanup()

    @staticmethod
    def _normalize_path(path: PathInput) -> Path:
        normalized = Path(path).expanduser()
        return normalized.resolve(strict=False)
