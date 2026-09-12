"""Top-level editable FEMaster project."""

from __future__ import annotations

import subprocess
from collections.abc import Iterable
from pathlib import Path

from ..amplitude.amplitude_repository import AmplitudeRepository
from ..load.load_collector_repository import LoadCollectorRepository
from ..support.support_collector_repository import SupportCollectorRepository
from ..common.format import block, blocks, keyword
from ..common.named_repository import NamedRepository
from ..constraint.constraint_repository import ConstraintRepository
from ..coordinate_system.coordinate_system_repository import CoordinateSystemRepository
from ..feature.feature_repository import FeatureRepository
from ..field.field_repository import FieldRepository
from ..instance.instance_repository import InstanceRepository
from ..material.material_repository import MaterialRepository
from ..mesh.surface import Surface
from ..part.part_repository import PartRepository
from ..profile.profile_repository import ProfileRepository
from ..region.region_repository import RegionRepository
from ..section.assembly_section_repository import AssemblySectionRepository
from ..step.step_repository import StepRepository


class Project:
    """Complete editable FEMaster project.

    The implicit default Part exists only as ``parts[0]``. Part-local mesh data
    never lives directly on Project. Shared definitions and assembly-level
    definitions are owned here.
    """

    def __init__(self, name: str = "model") -> None:
        self.name = str(name)

        # Semantic topology and assembly.
        self.parts = PartRepository()
        self.instances = InstanceRepository()
        self.regions = RegionRepository()
        self.surfaces = NamedRepository[Surface]()
        self.sections = AssemblySectionRepository()

        # Shared definitions.
        self.materials = MaterialRepository()
        self.profiles = ProfileRepository()
        self.coordinate_systems = CoordinateSystemRepository()
        self.amplitudes = AmplitudeRepository()
        self.fields = FieldRepository()
        self.features = FeatureRepository()

        # Constraints and collectors.
        self.constraints = ConstraintRepository()
        self.load_collectors = LoadCollectorRepository()
        self.support_collectors = SupportCollectorRepository()

        # Analyses.
        self.steps = StepRepository()

        # Unsupported imported blocks are retained explicitly.
        self.unparsed_blocks: list[object] = []

    def export(self) -> str:
        """Export the complete native FEMaster input deck."""

        assembly_body = blocks((
            self.instances.export(),
            self.regions.export(),
            "\n\n".join(surface.export() for surface in self.surfaces),
            self.sections.export(),
        ))

        assembly = ""
        if assembly_body:
            assembly = block([
                keyword("ASSEMBLY"),
                assembly_body,
                keyword("ENDASSEMBLY"),
            ])

        return blocks((
            keyword("MODEL", NAME=self.name),
            self.coordinate_systems.export(),
            self.materials.export(),
            self.profiles.export(),
            self.amplitudes.export(),
            self.parts.export(),
            assembly,
            self.fields.export(),
            self.features.export(),
            self.constraints.export(),
            self.support_collectors.export(),
            self.load_collectors.export(),
            self.steps.export(),
            keyword("END"),
        )) + "\n"

    def write(self, path: str | Path) -> Path:
        """Export this Project to one UTF-8 input file."""

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.export(), encoding="utf-8")
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
        """Export the project and execute FEMaster synchronously."""

        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)

        filename = input_name or f"{self.name}.inp"
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
    def import_file(cls, path: str | Path) -> "Project":
        """Import a FEMaster/Abaqus-like input file."""

        from ...io.inp_importer import InpImporter

        return InpImporter().import_file(path)
