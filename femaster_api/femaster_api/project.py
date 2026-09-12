"""Top-level FEMaster Project container and execution entry point.

Project mirrors the semantic scope of FEMaster input: reusable Parts own local
topology, while shared definitions, assembly regions/surfaces, constraints,
collectors, features and analysis steps live at project level. The class also
provides the single write/run boundary requested by the public Python API.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Iterable

from ._format import block, blocks, keyword
from .fields import FieldRepository
from .model.assembly import InstanceRepository, PartRepository
from .model.boundary import AmplitudeRepository, LoadCollectorRepository, SupportCollectorRepository
from .model.constraints import ConstraintRepository
from .model.coordinates import CoordinateSystemRepository
from .model.features import FeatureRepository
from .model.materials import MaterialRepository, ProfileRepository
from .model.mesh import Surface
from .model.regions import RegionRepository
from .model.steps import StepRepository
from .repository import NamedRepository


class Project:
    """Complete editable FEMaster project.

    The default Part lives exclusively inside ``parts`` at position zero. Project
    keeps no second default-part pointer or duplicate mesh repository. All
    part-local topology is therefore reached through ``project.parts``.
    """

    def __init__(self, name: str = "model") -> None:
        self.name = str(name)

        # ------------------------------------------------------------------
        # Semantic topology and assembly
        # ------------------------------------------------------------------

        self.parts     = PartRepository()
        self.instances = InstanceRepository()
        self.regions   = RegionRepository()
        self.surfaces  = NamedRepository[Surface]()

        # ------------------------------------------------------------------
        # Shared global definitions
        # ------------------------------------------------------------------

        self.materials          = MaterialRepository()
        self.profiles           = ProfileRepository()
        self.coordinate_systems = CoordinateSystemRepository()
        self.amplitudes         = AmplitudeRepository()
        self.fields             = FieldRepository()
        self.features           = FeatureRepository()

        # ------------------------------------------------------------------
        # Constraints and boundary-condition collectors
        # ------------------------------------------------------------------

        self.constraints        = ConstraintRepository()
        self.load_collectors    = LoadCollectorRepository()
        self.support_collectors = SupportCollectorRepository()

        # ------------------------------------------------------------------
        # Analysis definition
        # ------------------------------------------------------------------

        self.steps = StepRepository()

        # Reader fallback. Unknown keyword blocks are retained for diagnostics
        # but are not exported automatically because their scope/dependencies
        # cannot be reconstructed safely without a semantic implementation.
        self.unparsed_blocks: list[object] = []

    def to_femaster(self) -> str:
        """Return the complete FEMaster input deck in dependency order."""

        # ------------------------------------------------------------------
        # Build the assembly scope
        # ------------------------------------------------------------------

        assembly_body = blocks((
            self.instances.to_femaster(),
            self.regions.to_femaster(),
            "\n\n".join(surface.to_femaster() for surface in self.surfaces),
        ))

        assembly = ""
        if assembly_body:
            assembly = block([
                keyword("ASSEMBLY"),
                assembly_body,
                keyword("ENDASSEMBLY"),
            ])

        # ------------------------------------------------------------------
        # Emit definitions in semantic dependency order
        # ------------------------------------------------------------------

        return blocks((
            keyword("MODEL", NAME=self.name),

            # Shared definitions are independent of part instantiation.
            self.coordinate_systems.to_femaster(),
            self.materials.to_femaster(),
            self.profiles.to_femaster(),
            self.amplitudes.to_femaster(),

            # The default Part is written in root scope; explicit Parts receive
            # PART/ENDPART wrappers. Assembly instances materialize them later.
            self.parts.to_femaster(),
            assembly,

            # Remaining definitions operate on compiled assembly topology.
            self.fields.to_femaster(),
            self.features.to_femaster(),
            self.constraints.to_femaster(),

            # Collectors own individual BC/load entries. Steps only reference
            # collector names and therefore follow the definitions here.
            self.support_collectors.to_femaster(),
            self.load_collectors.to_femaster(),
            self.steps.to_femaster(),

            keyword("END"),
        )) + "\n"

    def write(self, path: str | Path) -> Path:
        """Write the complete project as a UTF-8 FEMaster input file."""

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_femaster(), encoding="utf-8")
        return path

    def run(
        self,
        *,
        executable: str | Path = "FEMaster",
        directory: str | Path = ".",
        arguments: Iterable[str] = (),
        input_name: str | None = None,
        check: bool = True,
    ) -> subprocess.CompletedProcess[str]:
        """Write the project and execute FEMaster synchronously.

        Additional command-line arguments are passed through unchanged instead
        of mirroring FEMaster's complete CLI in Python. The generated input path
        is always the final positional argument.
        """

        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)

        filename   = input_name or f"{self.name}.inp"
        input_path = self.write(directory / filename)

        command = [
            str(executable),
            "--format",
            "femaster",
            *[str(argument) for argument in arguments],
            str(input_path.resolve()),
        ]

        return subprocess.run(
            command,
            cwd=directory,
            text=True,
            check=check,
            capture_output=False,
        )

    @classmethod
    def read(cls, path: str | Path) -> "Project":
        """Read a FEMaster/Abaqus-like INP file into the Python object model."""

        from .io.inp import InpReader

        return InpReader().read(path)
