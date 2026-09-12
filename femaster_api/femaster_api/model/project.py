"""Editable FEMaster project root and native input-deck reader.

``Project`` owns every editable model concept that is not local to a ``Part``.
More importantly, the Python object graph uses real object relationships: an
``Instance`` stores a ``Part``, a load stores its target / amplitude /
orientation objects, sections store their region/material/profile objects, and
analysis steps store their collector objects.  Native IDs and semantic names are
therefore serialization tokens, not the in-memory representation of references.

``read_inp`` / ``read_inp_text`` are the boundary where textual tokens are
resolved.  A parsed name such as ``AMPLITUDE=RAMP`` is looked up immediately and
the resulting ``Amplitude`` object is passed to the model constructor.  Unknown
or unresolved relationships fail explicitly rather than leaving a dangling
string inside an otherwise valid project.  Unsupported keyword blocks remain in
``unparsed_blocks`` so syntax is never silently discarded.

Native export works in the opposite direction: concrete model classes own their
``export`` representation, repositories own grouping/order, and ``Project`` only
composes scopes and dependencies.  No generic entity-reference wrapper or
serializer registry is involved.
"""

from __future__ import annotations

import subprocess
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from .amplitude.amplitude import Amplitude
from .amplitude.amplitude_interpolation import AmplitudeInterpolation
from .amplitude.amplitude_repository import AmplitudeRepository
from .common.format import block, blocks, keyword
from .constraint.constraint_connector import Connector
from .constraint.constraint_coupling import Coupling
from .constraint.constraint_coupling_type import CouplingType
from .constraint.constraint_equation import Equation
from .constraint.constraint_repository import ConstraintRepository
from .constraint.constraint_rigid_body import RigidBodyConstraint
from .constraint.constraint_tie import Tie
from .coordinate_system.coordinate_system import CoordinateSystem
from .coordinate_system.coordinate_system_cylindrical import CylindricalCoordinateSystem
from .coordinate_system.coordinate_system_rectangular import RectangularCoordinateSystem
from .coordinate_system.coordinate_system_repository import CoordinateSystemRepository
from .element.element import Element
from .element.element_types import ELEMENT_TYPES
from .feature.feature_point_mass import PointMass
from .feature.feature_repository import FeatureRepository
from .field.field import Field
from .field.field_domain import FieldDomain
from .field.field_repository import FieldRepository
from .instance.instance import Instance
from .instance.instance_repository import InstanceRepository
from .load.load_collector import LoadCollector
from .load.load_collector_repository import LoadCollectorRepository
from .load.load_inertial import InertialLoad
from .load.load_nodal_force import NodalForce
from .load.load_pressure import PressureLoad
from .load.load_surface_traction import SurfaceTraction
from .load.load_thermal import ThermalLoad
from .load.load_volume import VolumeLoad
from .material.material import Material
from .material.material_elasticity_abd import ABDElasticity
from .material.material_elasticity_generalized_isotropic import GeneralizedIsotropicElasticity
from .material.material_elasticity_isotropic import IsotropicElasticity
from .material.material_elasticity_orthotropic import OrthotropicElasticity
from .material.material_repository import MaterialRepository
from .node.node import Node
from .part.part import Part
from .part.part_repository import PartRepository
from .profile.profile import Profile
from .profile.profile_repository import ProfileRepository
from .region.region_element import ElementRegion
from .region.region_node import NodeRegion
from .region.region_repository import RegionRepository
from .region.region_surface import SurfaceRegion
from .section.section_beam import BeamSection
from .section.section_mass import MassSection
from .section.section_repository import SectionRepository
from .section.section_repository_assembly import AssemblySectionRepository
from .section.section_rotary_inertia import RotaryInertiaSection
from .section.section_shell import ShellSection
from .section.section_shell_abd import ABDShellSection
from .section.section_solid import SolidSection
from .section.section_spring import SpringSection
from .section.section_truss import TrussSection
from .step.step_buckling import BucklingStep
from .step.step_modal import ModalStep
from .step.step_nonlinear_static import NonlinearStaticStep
from .step.step_repository import StepRepository
from .step.step_static import StaticStep
from .step.step_transient import TransientStep
from .step.util.constraint_method import ConstraintMethod
from .step.util.newmark_control import NewmarkControl
from .step.util.rayleigh_damping import RayleighDamping
from .step.util.solver_control import SolverControl
from .step.util.solver_device import SolverDevice
from .step.util.solver_method import SolverMethod
from .step.util.time_control import TimeControl
from .support.support import Support
from .support.support_collector import SupportCollector
from .support.support_collector_repository import SupportCollectorRepository
from .surface.surface import Surface
from .surface.surface_element import ElementSurface
from .surface.surface_node import NodeSurface
from .surface.surface_repository import SurfaceRepository


_Block = dict[str, Any]


class Project:
    """Complete editable FEMaster project and input-deck import boundary."""

    def __init__(self, name: str = "model") -> None:
        self.name = str(name)

        # Reusable topology and assembly.  The default Part is created and
        # retained solely by PartRepository at position zero.
        self.parts = PartRepository()
        self.instances = InstanceRepository()
        self.regions = RegionRepository()
        self.surfaces = SurfaceRepository()
        self.sections = AssemblySectionRepository()

        # Shared global definitions.  Other model objects retain references to
        # these objects directly rather than copying their semantic names.
        self.materials = MaterialRepository()
        self.profiles = ProfileRepository()
        self.coordinate_systems = CoordinateSystemRepository()
        self.amplitudes = AmplitudeRepository()
        self.fields = FieldRepository()
        self.features = FeatureRepository()

        # Assembly constraints and collector-owned boundary definitions.
        self.constraints = ConstraintRepository()
        self.load_collectors = LoadCollectorRepository()
        self.support_collectors = SupportCollectorRepository()

        # Ordered analysis procedures.
        self.steps = StepRepository()

        # Unsupported imported blocks remain visible for diagnostics / future
        # support.  They are not re-exported because their semantics are unknown.
        self.unparsed_blocks: list[_Block] = []

    # ------------------------------------------------------------------
    # Native export and execution
    # ------------------------------------------------------------------

    def export(self) -> str:
        """Export the complete native FEMaster deck in dependency order."""

        assembly_body = blocks((
            self.instances.export(),
            self.regions.export(),
            self.surfaces.export(),
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
        """Write this project as one UTF-8 native input deck."""

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
        """Export this model and execute FEMaster synchronously."""

        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)

        input_path = self.write(
            directory / (input_name or f"{self.name}.inp")
        )
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

    # ------------------------------------------------------------------
    # Public INP readers
    # ------------------------------------------------------------------

    @classmethod
    def read_inp(cls, path: str | Path) -> "Project":
        """Read a FEMaster/Abaqus-like keyword deck into public model objects."""

        return cls.read_inp_text(Path(path).read_text(encoding="utf-8"))

    @classmethod
    def read_inp_text(cls, text: str) -> "Project":
        """Parse input text and resolve every supported relation to an object."""

        parsed = cls._blocks(text)
        project = cls(cls._model_name(parsed))

        current_part = project.parts.default()
        current_material: Material | None = None
        current_step: Any = None
        in_assembly = False

        for item in parsed:
            name = item["name"]

            # --------------------------------------------------------------
            # Scope control
            # --------------------------------------------------------------
            if name in {"MODEL", "END"}:
                continue

            if name == "PART":
                current_part = project.parts.add(
                    Part(cls._required(item, "NAME"))
                )
                current_material = None
                current_step = None
                in_assembly = False
                continue

            if name == "ENDPART":
                current_part = project.parts.default()
                current_material = None
                continue

            if name == "ASSEMBLY":
                in_assembly = True
                current_part = project.parts.default()
                current_material = None
                continue

            if name == "ENDASSEMBLY":
                in_assembly = False
                continue

            if name == "INSTANCE":
                project.instances.add(cls._read_instance(project, item))
                continue

            # --------------------------------------------------------------
            # Analysis procedures and child commands
            # --------------------------------------------------------------
            if name == "LOADCASE":
                current_step = cls._create_step(item)
                project.steps.add(current_step)
                continue

            if current_step is not None and cls._read_step_subcommand(
                project,
                current_step,
                item,
            ):
                continue

            # --------------------------------------------------------------
            # Part/assembly topology and named regions
            # --------------------------------------------------------------
            regions = project.regions if in_assembly else current_part.regions
            surfaces = project.surfaces if in_assembly else current_part.surfaces

            if name == "NODE":
                nodes: list[Node] = []
                for row in item["data"]:
                    values = [value for value in row if value != ""]
                    if len(values) < 3:
                        raise ValueError(
                            f"NODE at line {item['line']} requires id, x, y"
                        )
                    node = current_part.nodes.add(
                        Node(
                            int(values[0]),
                            float(values[1]),
                            float(values[2]),
                            float(values[3]) if len(values) > 3 else 0.0,
                        )
                    )
                    nodes.append(node)

                nset = cls._key(item, "NSET")
                if nset:
                    cls._node_region(regions, nset).add(*nodes)
                continue

            if name == "ELEMENT":
                type_name = cls._required(item, "TYPE").upper()
                element_class = ELEMENT_TYPES.get(type_name)
                elements: list[Element] = []

                for row in item["data"]:
                    values = [value for value in row if value != ""]
                    if len(values) < 2:
                        raise ValueError(
                            f"ELEMENT at line {item['line']} requires "
                            "id and connectivity"
                        )

                    element_id = int(values[0])
                    nodes = tuple(
                        current_part.nodes[int(value)]
                        for value in values[1:]
                    )

                    if element_class is None:
                        element = Element(element_id, nodes)
                        element.type_name = type_name
                    else:
                        element = element_class(element_id, nodes)

                    current_part.elements.add(element)
                    elements.append(element)

                elset = cls._key(item, "ELSET")
                if elset:
                    cls._element_region(regions, elset).add(*elements)
                continue

            if name == "NSET":
                region_name = cls._key(item, "NAME") or cls._key(item, "NSET")
                if not region_name:
                    raise ValueError(
                        f"NSET at line {item['line']} requires NAME"
                    )
                cls._node_region(regions, region_name).add(
                    *cls._resolve_node_members(
                        current_part,
                        regions,
                        cls._flat_tokens(item),
                    )
                )
                continue

            if name == "ELSET":
                region_name = cls._key(item, "NAME") or cls._key(item, "ELSET")
                if not region_name:
                    raise ValueError(
                        f"ELSET at line {item['line']} requires NAME"
                    )
                cls._element_region(regions, region_name).add(
                    *cls._resolve_element_members(
                        current_part,
                        regions,
                        cls._flat_tokens(item),
                    )
                )
                continue

            if name == "SFSET":
                region_name = cls._key(item, "SFSET") or cls._key(item, "NAME")
                if not region_name:
                    raise ValueError(
                        f"SFSET at line {item['line']} requires SFSET or NAME"
                    )
                region = (
                    regions.surfaces[region_name]
                    if region_name in regions.surfaces
                    else regions.surfaces.add(SurfaceRegion(region_name))
                )
                region.add(
                    *cls._resolve_surface_members(
                        surfaces,
                        regions,
                        cls._flat_tokens(item),
                    )
                )
                continue

            if name == "SURFACE":
                surface_name = (
                    cls._key(item, "NAME")
                    or cls._key(item, "SFSET")
                    or cls._key(item, "SURFACE")
                )
                if not surface_name:
                    raise ValueError(
                        f"SURFACE at line {item['line']} requires NAME or SFSET"
                    )

                surface_type = (
                    cls._key(item, "TYPE", "ELEMENT") or "ELEMENT"
                ).upper()

                if surface_type == "NODE":
                    surface = NodeSurface(surface_name)
                    for row in item["data"]:
                        if row and row[0]:
                            surface.add(
                                cls._resolve_node_target(
                                    project,
                                    current_part,
                                    cls._token(row[0]),
                                )
                            )
                    surfaces.add(surface)
                    continue

                surface = ElementSurface(surface_name)
                for row in item["data"]:
                    if len(row) >= 2 and row[0] and row[1]:
                        surface.add(
                            cls._resolve_element_target(
                                project,
                                current_part,
                                cls._token(row[0]),
                            ),
                            cls._surface_side(row[1]),
                        )
                surfaces.add(surface)
                continue

            # --------------------------------------------------------------
            # Shared global definitions
            # --------------------------------------------------------------
            if name == "MATERIAL":
                current_material = Material(cls._required(item, "NAME"))
                project.materials.add(current_material)
                continue

            if name == "ELASTIC" and current_material is not None:
                values = cls._flat_floats(item)
                type_name = (
                    cls._key(item, "TYPE", "ISOTROPIC") or "ISOTROPIC"
                ).upper()

                if type_name == "ISOTROPIC":
                    current_material.elasticity = IsotropicElasticity(
                        values[0], values[1]
                    )
                elif type_name == "GENISO":
                    current_material.elasticity = GeneralizedIsotropicElasticity(
                        *values[:3]
                    )
                elif type_name == "ENGINEERINGCONSTANTS":
                    current_material.elasticity = OrthotropicElasticity(
                        *values[:9]
                    )
                elif type_name == "ABD":
                    current_material.elasticity = ABDElasticity(values)
                else:
                    project.unparsed_blocks.append(item)
                continue

            if name == "DENSITY" and current_material is not None:
                current_material.density = cls._flat_floats(item)[0]
                continue

            if name == "THERMALEXPANSION" and current_material is not None:
                current_material.thermal_expansion = cls._flat_floats(item)[0]
                continue

            if name == "PROFILE":
                values = cls._flat_floats(item)
                values.extend([0.0] * max(0, 9 - len(values)))
                project.profiles.add(
                    Profile(cls._required(item, "NAME"), *values[:9])
                )
                continue

            if name == "ORIENTATION":
                values = cls._flat_floats(item)
                orientation_name = cls._required(item, "NAME")
                type_name = (
                    cls._key(item, "TYPE", "RECTANGULAR") or "RECTANGULAR"
                ).upper()

                if type_name == "CYLINDRICAL":
                    project.coordinate_systems.add(
                        CylindricalCoordinateSystem(
                            orientation_name,
                            values[0:3],
                            values[3:6],
                            values[6:9],
                        )
                    )
                else:
                    project.coordinate_systems.add(
                        RectangularCoordinateSystem(
                            orientation_name,
                            values[0:3],
                            values[3:6] if len(values) >= 6 else None,
                            values[6:9] if len(values) >= 9 else None,
                        )
                    )
                continue

            # --------------------------------------------------------------
            # Sections and concentrated point-element properties
            # --------------------------------------------------------------
            if name in {
                "SOLIDSECTION",
                "SHELLSECTION",
                "BEAMSECTION",
                "TRUSSSECTION",
            }:
                cls._read_part_section(project, current_part, item)
                continue

            if name in {"MASS", "ROTARYINERTIA", "SPRING"}:
                target_sections = (
                    project.sections if in_assembly else current_part.sections
                )
                target_regions = (
                    project.regions if in_assembly else current_part.regions
                )
                cls._read_point_section(target_sections, target_regions, item)
                continue

            # --------------------------------------------------------------
            # Fields, amplitudes and other global features
            # --------------------------------------------------------------
            if name == "AMPLITUDE":
                interpolation_name = (
                    cls._key(item, "TYPE", "LINEAR") or "LINEAR"
                ).upper()
                interpolation = AmplitudeInterpolation[interpolation_name]
                amplitude = Amplitude(
                    cls._required(item, "NAME"),
                    interpolation=interpolation,
                )
                values = cls._flat_floats(item)
                for index in range(0, len(values) - 1, 2):
                    amplitude.add(values[index], values[index + 1])
                project.amplitudes.add(amplitude)
                continue

            if name == "FIELD":
                project.fields.add(
                    cls._read_model_field(project, current_part, item)
                )
                continue

            if name == "POINTMASS":
                values = cls._flat_floats(item)
                values.extend([0.0] * max(0, 10 - len(values)))
                project.features.add(
                    PointMass(
                        cls._resolve_node_region(
                            project,
                            current_part,
                            cls._required(item, "NSET"),
                        ),
                        values[0],
                        inertia=values[1:4],
                        spring=values[4:7],
                        rotational_spring=values[7:10],
                    )
                )
                continue

            # --------------------------------------------------------------
            # Supports and loads
            # --------------------------------------------------------------
            if name == "SUPPORT":
                collector = cls._support_collector(
                    project,
                    cls._required(item, "SUPPORT_COLLECTOR"),
                )
                orientation = cls._resolve_coordinate_system(
                    project,
                    cls._key(item, "ORIENTATION"),
                )

                for row in item["data"]:
                    if not row:
                        continue
                    values = [
                        None if value == "" else float(value)
                        for value in row[1:]
                    ]
                    collector.add(
                        Support(
                            cls._resolve_node_target(
                                project,
                                current_part,
                                cls._token(row[0]),
                            ),
                            values,
                            orientation=orientation,
                        )
                    )
                continue

            if name in {
                "CLOAD",
                "DLOAD",
                "PLOAD",
                "VLOAD",
                "TLOAD",
                "INERTIALOAD",
            }:
                cls._read_load(project, current_part, item)
                continue

            # --------------------------------------------------------------
            # Assembly-level constraints
            # --------------------------------------------------------------
            if name == "RBM":
                project.constraints.add(
                    RigidBodyConstraint(
                        cls._resolve_element_region(
                            project,
                            current_part,
                            cls._required(item, "ELSET"),
                        )
                    )
                )
                continue

            if name == "COUPLING":
                surface_name = (
                    cls._key(item, "SURFACE")
                    or cls._key(item, "SFSET")
                )
                slave_name = surface_name or cls._key(item, "SLAVE")
                if not slave_name:
                    raise ValueError(
                        f"COUPLING at line {item['line']} requires a slave target"
                    )

                slave = (
                    cls._resolve_surface_target(
                        project,
                        current_part,
                        surface_name,
                    )
                    if surface_name is not None
                    else cls._resolve_node_target(
                        project,
                        current_part,
                        cls._token(slave_name),
                    )
                )
                dofs = (
                    [int(value) for value in item["data"][0] if value != ""]
                    if item["data"]
                    else (1, 1, 1, 1, 1, 1)
                )
                project.constraints.add(
                    Coupling(
                        cls._resolve_node_target(
                            project,
                            current_part,
                            cls._token(cls._required(item, "MASTER")),
                        ),
                        slave,
                        type=CouplingType[
                            (
                                cls._key(item, "TYPE", "KINEMATIC")
                                or "KINEMATIC"
                            ).upper()
                        ],
                        dofs=dofs,
                    )
                )
                continue

            if name == "CONNECTOR":
                project.constraints.add(
                    Connector(
                        Connector.Type[cls._required(item, "TYPE").upper()],
                        cls._resolve_node_region(
                            project,
                            current_part,
                            cls._required(item, "NSET1"),
                        ),
                        cls._resolve_node_region(
                            project,
                            current_part,
                            cls._required(item, "NSET2"),
                        ),
                        cls._require_coordinate_system(
                            project,
                            cls._required(item, "COORDINATESYSTEM"),
                        ),
                    )
                )
                continue

            if name == "TIE":
                project.constraints.add(
                    Tie(
                        cls._resolve_surface_target(
                            project,
                            current_part,
                            cls._required(item, "MASTER"),
                        ),
                        cls._resolve_surface_target(
                            project,
                            current_part,
                            cls._required(item, "SLAVE"),
                        ),
                        adjust=(
                            cls._key(item, "ADJUST", "YES") or "YES"
                        ).upper() == "YES",
                        distance=(
                            float(cls._key(item, "DISTANCE"))
                            if cls._key(item, "DISTANCE") is not None
                            else None
                        ),
                    )
                )
                continue

            if name == "EQUATION":
                tokens = [
                    value
                    for row in item["data"]
                    for value in row
                    if value != ""
                ]
                if tokens:
                    term_count = int(tokens[0])
                    equation = Equation()
                    values = tokens[1:]
                    for index in range(term_count):
                        start = 3 * index
                        equation.add(
                            cls._resolve_node_target(
                                project,
                                current_part,
                                cls._token(values[start]),
                            ),
                            int(values[start + 1]),
                            float(values[start + 2]),
                        )
                    project.constraints.add(equation)
                continue

            project.unparsed_blocks.append(item)

        return project

    # ------------------------------------------------------------------
    # Input syntax helpers
    # ------------------------------------------------------------------

    @classmethod
    def _blocks(cls, text: str) -> list[_Block]:
        """Split input text into simple keyword dictionaries."""

        result: list[_Block] = []
        current: _Block | None = None

        for line_number, raw in enumerate(text.splitlines(), start=1):
            stripped = raw.strip()
            if not stripped or stripped.startswith("**"):
                continue

            if stripped.startswith("*"):
                current = cls._parse_keyword(stripped, line_number)
                result.append(current)
                continue

            if current is not None:
                current["data"].append(
                    [token.strip() for token in raw.split(",")]
                )

        return result

    @staticmethod
    def _parse_keyword(line: str, line_number: int) -> _Block:
        tokens = [token.strip() for token in line[1:].split(",")]
        name = tokens[0].upper().replace(" ", "")
        keys: dict[str, str] = {}

        for token in tokens[1:]:
            if "=" in token:
                key, value = token.split("=", 1)
                keys[key.strip().upper().replace(" ", "")] = value.strip()

        return {
            "name": name,
            "keys": keys,
            "data": [],
            "line": line_number,
        }

    @staticmethod
    def _key(
        item: _Block,
        name: str,
        default: str | None = None,
    ) -> str | None:
        return item["keys"].get(name.upper().replace(" ", ""), default)

    @classmethod
    def _required(cls, item: _Block, key: str) -> str:
        value = cls._key(item, key)
        if value is None or value == "":
            raise ValueError(
                f"{item['name']} at line {item['line']} requires {key}"
            )
        return value

    @classmethod
    def _model_name(cls, items: list[_Block]) -> str:
        for item in items:
            if item["name"] == "MODEL":
                return cls._key(item, "NAME", "model") or "model"
        return "model"

    @classmethod
    def _flat_tokens(cls, item: _Block) -> list[int | str]:
        return [
            cls._token(value)
            for row in item["data"]
            for value in row
            if value != ""
        ]

    @staticmethod
    def _flat_floats(item: _Block) -> list[float]:
        return [
            float(value)
            for row in item["data"]
            for value in row
            if value != ""
        ]

    @staticmethod
    def _token(value: str) -> int | str:
        value = value.strip()
        try:
            return int(value)
        except ValueError:
            return value

    @staticmethod
    def _surface_side(value: str) -> int:
        value = value.strip().upper()
        if value == "SPOS":
            return 1
        if value == "SNEG":
            return 2
        if value.startswith("S"):
            value = value[1:]
        return int(value)

    # ------------------------------------------------------------------
    # Native-token -> model-object resolution
    # ------------------------------------------------------------------

    @staticmethod
    def _resolve_node_members(
        part: Part,
        regions: RegionRepository,
        tokens: Iterable[int | str],
    ) -> list[Node]:
        result: list[Node] = []
        for token in tokens:
            if isinstance(token, int):
                result.append(part.nodes[token])
            elif token in regions.nodes:
                result.extend(regions.nodes[token].members)
            else:
                raise ValueError(f"cannot resolve node-region member: {token}")
        return result

    @staticmethod
    def _resolve_element_members(
        part: Part,
        regions: RegionRepository,
        tokens: Iterable[int | str],
    ) -> list[Element]:
        result: list[Element] = []
        for token in tokens:
            if isinstance(token, int):
                result.append(part.elements[token])
            elif token in regions.elements:
                result.extend(regions.elements[token].members)
            else:
                raise ValueError(f"cannot resolve element-region member: {token}")
        return result

    @staticmethod
    def _resolve_surface_members(
        surfaces: SurfaceRepository,
        regions: RegionRepository,
        tokens: Iterable[int | str],
    ) -> list[Surface]:
        result: list[Surface] = []
        for token in tokens:
            if not isinstance(token, str):
                raise ValueError("surface-region members must be named surfaces")
            if token in surfaces:
                result.append(surfaces[token])
            elif token in regions.surfaces:
                result.extend(regions.surfaces[token].members)
            else:
                raise ValueError(f"cannot resolve surface-region member: {token}")
        return result

    @classmethod
    def _resolve_node_target(
        cls,
        project: "Project",
        part: Part,
        token: int | str,
    ) -> Node | NodeRegion:
        if isinstance(token, int):
            return part.nodes[token]
        if token in part.regions.nodes:
            return part.regions.nodes[token]
        if token in project.regions.nodes:
            return project.regions.nodes[token]
        raise ValueError(f"cannot resolve node target: {token}")

    @classmethod
    def _resolve_element_target(
        cls,
        project: "Project",
        part: Part,
        token: int | str,
    ) -> Element | ElementRegion:
        if isinstance(token, int):
            return part.elements[token]
        if token in part.regions.elements:
            return part.regions.elements[token]
        if token in project.regions.elements:
            return project.regions.elements[token]
        raise ValueError(f"cannot resolve element target: {token}")

    @classmethod
    def _resolve_surface_target(
        cls,
        project: "Project",
        part: Part,
        name: str,
    ) -> Surface | SurfaceRegion:
        if name in part.surfaces:
            return part.surfaces[name]
        if name in part.regions.surfaces:
            return part.regions.surfaces[name]
        if name in project.surfaces:
            return project.surfaces[name]
        if name in project.regions.surfaces:
            return project.regions.surfaces[name]
        raise ValueError(f"cannot resolve surface target: {name}")

    @classmethod
    def _resolve_node_region(
        cls,
        project: "Project",
        part: Part,
        name: str,
    ) -> NodeRegion:
        if name in part.regions.nodes:
            return part.regions.nodes[name]
        if name in project.regions.nodes:
            return project.regions.nodes[name]
        raise ValueError(f"cannot resolve node region: {name}")

    @classmethod
    def _resolve_element_region(
        cls,
        project: "Project",
        part: Part,
        name: str,
    ) -> ElementRegion:
        if name in part.regions.elements:
            return part.regions.elements[name]
        if name in project.regions.elements:
            return project.regions.elements[name]
        raise ValueError(f"cannot resolve element region: {name}")

    @staticmethod
    def _resolve_coordinate_system(
        project: "Project",
        name: str | None,
    ) -> CoordinateSystem | None:
        return None if name is None else project.coordinate_systems[name]

    @staticmethod
    def _require_coordinate_system(
        project: "Project",
        name: str,
    ) -> CoordinateSystem:
        return project.coordinate_systems[name]

    @staticmethod
    def _resolve_amplitude(
        project: "Project",
        name: str | None,
    ) -> Amplitude | None:
        return None if name is None else project.amplitudes[name]

    # ------------------------------------------------------------------
    # Object construction helpers used by the input reader
    # ------------------------------------------------------------------

    @classmethod
    def _read_instance(cls, project: "Project", item: _Block) -> Instance:
        translation = None
        rotation = None

        if item["data"] and len(
            [value for value in item["data"][0] if value != ""]
        ) == 3:
            translation = tuple(float(value) for value in item["data"][0])

        rotation_row = None
        if len(item["data"]) > 1:
            rotation_row = item["data"][1]
        elif item["data"] and len(
            [value for value in item["data"][0] if value != ""]
        ) == 7:
            rotation_row = item["data"][0]

        if rotation_row is not None:
            values = tuple(
                float(value)
                for value in rotation_row
                if value != ""
            )
            rotation = (values[0:3], values[3:6], values[6])

        return Instance(
            cls._required(item, "NAME"),
            project.parts[cls._required(item, "PART")],
            translation=translation,
            rotation=rotation,
        )

    @classmethod
    def _read_part_section(
        cls,
        project: "Project",
        part: Part,
        item: _Block,
    ) -> None:
        name = item["name"]
        section_name = cls._key(item, "NAME") or f"{name}_{len(part.sections)}"
        element_region = part.regions.elements[cls._required(item, "ELSET")]

        material_name = cls._key(item, "MATERIAL")
        material = project.materials[material_name] if material_name else None
        orientation_name = cls._key(item, "ORIENTATION")
        orientation = cls._resolve_coordinate_system(project, orientation_name)

        if name == "SOLIDSECTION":
            if material is None:
                raise ValueError("SOLIDSECTION requires MATERIAL")
            part.sections.add(
                SolidSection(
                    section_name,
                    element_region,
                    material,
                    orientation,
                )
            )
            return

        if name == "TRUSSSECTION":
            if material is None:
                raise ValueError("TRUSSSECTION requires MATERIAL")
            part.sections.add(
                TrussSection(
                    section_name,
                    element_region,
                    material,
                    cls._flat_floats(item)[0],
                )
            )
            return

        if name == "BEAMSECTION":
            if material is None:
                raise ValueError("BEAMSECTION requires MATERIAL")
            beam_orientation = cls._flat_floats(item)
            part.sections.add(
                BeamSection(
                    section_name,
                    element_region,
                    material,
                    project.profiles[cls._required(item, "PROFILE")],
                    beam_orientation[:3],
                )
            )
            return

        section_type = (
            cls._key(item, "TYPE", "INTEGRATED") or "INTEGRATED"
        ).upper()
        csys_axis = int(cls._key(item, "CSYSAXIS", "1") or 1)

        if section_type == "ABD":
            values = cls._flat_floats(item)
            if len(values) != 40:
                raise ValueError(
                    f"SHELLSECTION TYPE=ABD at line {item['line']} "
                    "requires 40 values"
                )
            part.sections.add(
                ABDShellSection(
                    section_name,
                    element_region,
                    values[:36],
                    values[36:40],
                    thickness=float(
                        cls._key(item, "THICKNESS", "1.0") or 1.0
                    ),
                    material=material,
                    orientation=orientation,
                    csys_axis=csys_axis,
                )
            )
            return

        if material is None:
            raise ValueError("integrated SHELLSECTION requires MATERIAL")
        thickness_values = cls._flat_floats(item)
        thickness = (
            thickness_values[0]
            if thickness_values
            else float(cls._key(item, "THICKNESS", "1.0") or 1.0)
        )
        part.sections.add(
            ShellSection(
                section_name,
                element_region,
                material,
                thickness,
                orientation,
                csys_axis,
            )
        )

    @classmethod
    def _read_point_section(
        cls,
        sections: SectionRepository,
        regions: RegionRepository,
        item: _Block,
    ) -> None:
        section_name = (
            cls._key(item, "NAME")
            or f"{item['name']}_{len(sections)}"
        )
        element_region = regions.elements[cls._required(item, "ELSET")]
        values = cls._flat_floats(item)

        if item["name"] == "MASS":
            sections.add(MassSection(section_name, element_region, values[0]))
            return

        if item["name"] == "ROTARYINERTIA":
            values.extend([0.0] * max(0, 6 - len(values)))
            if any(value != 0.0 for value in values[3:6]):
                raise ValueError(
                    "ROTARY INERTIA products I12, I13 and I23 must be zero"
                )
            sections.add(
                RotaryInertiaSection(section_name, element_region, values[:3])
            )
            return

        if len(values) < 2:
            raise ValueError(
                f"SPRING at line {item['line']} requires DOF and stiffness"
            )
        sections.add(
            SpringSection(section_name, element_region, int(values[0]), values[1])
        )

    @classmethod
    def _read_model_field(
        cls,
        project: "Project",
        part: Part,
        item: _Block,
    ) -> Field:
        """Read one model field and resolve its entity address to model objects."""

        domain_name = (
            cls._key(item, "TYPE")
            or cls._key(item, "DOMAIN")
            or "UNKNOWN"
        ).upper().replace("_", "")

        aliases = {
            "NODE": FieldDomain.NODE,
            "ELEMENT": FieldDomain.ELEMENT,
            "ELEMENTNODAL": FieldDomain.ELEMENT_NODAL,
            "ELEMENTIP": FieldDomain.ELEMENT_IP,
            "IP": FieldDomain.ELEMENT_IP,
            "ELEMENTMP": FieldDomain.ELEMENT_MP,
            "MP": FieldDomain.ELEMENT_MP,
        }
        domain = aliases.get(domain_name, FieldDomain.UNKNOWN)
        if domain is FieldDomain.UNKNOWN:
            raise ValueError(
                f"FIELD at line {item['line']} has unsupported domain: "
                f"{domain_name}"
            )

        cols = int(cls._key(item, "COLS", "0") or 0)
        result = Field(
            cls._required(item, "NAME"),
            domain,
            [f"C{index + 1}" for index in range(cols)],
        )

        # The first address column always denotes a model entity and is resolved
        # immediately.  Remaining columns are zero-based positions local to that
        # element and therefore remain intrinsic integers.
        for row in item["data"]:
            if not row or row[0] == "":
                continue

            if domain is FieldDomain.NODE:
                target = cls._resolve_node_target(
                    project,
                    part,
                    cls._token(row[0]),
                )
                if not isinstance(target, Node):
                    raise ValueError(
                        "NODE field rows must address a concrete Node"
                    )
                key = target
                value_start = 1

            elif domain is FieldDomain.ELEMENT:
                target = cls._resolve_element_target(
                    project,
                    part,
                    cls._token(row[0]),
                )
                if not isinstance(target, Element):
                    raise ValueError(
                        "ELEMENT field rows must address a concrete Element"
                    )
                key = target
                value_start = 1

            elif domain in {
                FieldDomain.ELEMENT_NODAL,
                FieldDomain.ELEMENT_IP,
            }:
                if len(row) < 2 or row[1] == "":
                    raise ValueError(
                        f"{domain.value} field row requires element and local index"
                    )
                target = cls._resolve_element_target(
                    project,
                    part,
                    cls._token(row[0]),
                )
                if not isinstance(target, Element):
                    raise ValueError(
                        f"{domain.value} field rows must address a concrete Element"
                    )
                key = (target, int(row[1]))
                value_start = 2

            else:
                if len(row) < 3 or row[1] == "" or row[2] == "":
                    raise ValueError(
                        "ELEMENT_MP field row requires element, local IP and local MP"
                    )
                target = cls._resolve_element_target(
                    project,
                    part,
                    cls._token(row[0]),
                )
                if not isinstance(target, Element):
                    raise ValueError(
                        "ELEMENT_MP field rows must address a concrete Element"
                    )
                key = (target, int(row[1]), int(row[2]))
                value_start = 3

            result.set(
                key,
                (
                    float(value)
                    for value in row[value_start:]
                    if value != ""
                ),
            )

        return result

    @classmethod
    def _read_load(
        cls,
        project: "Project",
        part: Part,
        item: _Block,
    ) -> None:
        collector = cls._load_collector(
            project,
            cls._required(item, "LOAD_COLLECTOR"),
        )
        amplitude = cls._resolve_amplitude(project, cls._key(item, "AMPLITUDE"))
        orientation = cls._resolve_coordinate_system(
            project,
            cls._key(item, "ORIENTATION"),
        )

        if item["name"] == "CLOAD":
            for row in item["data"]:
                collector.add(
                    NodalForce(
                        cls._resolve_node_target(
                            project,
                            part,
                            cls._token(row[0]),
                        ),
                        (float(value) for value in row[1:7]),
                        orientation=orientation,
                        amplitude=amplitude,
                    )
                )
            return

        if item["name"] == "DLOAD":
            for row in item["data"]:
                collector.add(
                    SurfaceTraction(
                        cls._resolve_surface_target(project, part, row[0]),
                        (float(value) for value in row[1:4]),
                        orientation=orientation,
                        amplitude=amplitude,
                    )
                )
            return

        if item["name"] == "PLOAD":
            for row in item["data"]:
                collector.add(
                    PressureLoad(
                        cls._resolve_surface_target(project, part, row[0]),
                        float(row[1]),
                        amplitude=amplitude,
                    )
                )
            return

        if item["name"] == "VLOAD":
            for row in item["data"]:
                collector.add(
                    VolumeLoad(
                        cls._resolve_element_target(
                            project,
                            part,
                            cls._token(row[0]),
                        ),
                        (float(value) for value in row[1:4]),
                        orientation=orientation,
                        amplitude=amplitude,
                    )
                )
            return

        if item["name"] == "TLOAD":
            collector.add(
                ThermalLoad(
                    project.fields[cls._required(item, "TEMPERATUREFIELD")],
                    float(
                        cls._key(item, "REFERENCETEMPERATURE", "0") or 0
                    ),
                )
            )
            return

        row = item["data"][0]
        values = [float(value) for value in row[1:] if value != ""]
        values.extend([0.0] * max(0, 12 - len(values)))
        collector.add(
            InertialLoad(
                cls._resolve_element_target(
                    project,
                    part,
                    cls._token(row[0]),
                ),
                center=values[0:3],
                center_acceleration=values[3:6],
                omega=values[6:9],
                alpha=values[9:12],
                consider_point_masses=(
                    cls._key(item, "CONSIDER_POINT_MASSES", "0") == "1"
                ),
            )
        )

    @classmethod
    def _create_step(cls, item: _Block) -> Any:
        type_name = cls._required(item, "TYPE").upper()
        name = cls._key(item, "NAME") or f"STEP_{item['line']}"

        if type_name == "LINEARSTATIC":
            return StaticStep(name)
        if type_name == "EIGENFREQ":
            return ModalStep(name, 0)
        if type_name == "LINEARBUCKLING":
            return BucklingStep(name, 0)
        if type_name == "NONLINEARSTATIC":
            return NonlinearStaticStep(name)
        if type_name == "LINEARTRANSIENT":
            return TransientStep(name, TimeControl(0.0, 0.0, 0.0))

        raise ValueError(
            f"unsupported LOADCASE TYPE at line {item['line']}: {type_name}"
        )

    @classmethod
    def _read_step_subcommand(
        cls,
        project: "Project",
        step: Any,
        item: _Block,
    ) -> bool:
        name = item["name"]

        if name == "SUPPORTS":
            step.supports = tuple(
                project.support_collectors[value]
                for row in item["data"]
                for value in row
                if value
            )
            return True

        if name == "LOADS":
            step.loads = tuple(
                project.load_collectors[value]
                for row in item["data"]
                for value in row
                if value
            )
            return True

        if name == "SOLVER":
            step.solver = SolverControl(
                SolverDevice[
                    (cls._key(item, "DEVICE", "CPU") or "CPU").upper()
                ],
                SolverMethod[
                    (cls._key(item, "METHOD", "DIRECT") or "DIRECT").upper()
                ],
            )
            return True

        if name == "CONSTRAINTMETHOD":
            step.constraint_method = ConstraintMethod[
                cls._required(item, "TYPE").upper()
            ]
            return True

        if name == "NUMEIGENVALUES" and hasattr(step, "number_of_modes"):
            step.number_of_modes = int(item["data"][0][0])
            return True

        if name == "SIGMA" and hasattr(step, "sigma"):
            step.sigma = float(item["data"][0][0])
            return True

        if name == "NONLINEAR" and isinstance(step, NonlinearStaticStep):
            step.control = (
                cls._key(item, "CONTROL", step.control) or step.control
            ).upper()
            step.increments = cls._int_key(item, "INCREMENTS")
            step.max_increments = cls._int_key(item, "MAX_INCREMENTS")
            step.initial_increment = cls._float_key(
                item,
                "INITIAL_INCREMENT",
            )
            step.minimum_increment = cls._float_key(
                item,
                "MINIMUM_INCREMENT",
            )
            step.maximum_increment = cls._float_key(
                item,
                "MAXIMUM_INCREMENT",
            )
            step.max_iterations = cls._int_key(item, "MAXITER")
            step.tolerance = cls._float_key(item, "TOL")
            return True

        if name == "TIME" and isinstance(step, TransientStep):
            step.time = TimeControl(*tuple(map(float, item["data"][0][:3])))
            return True

        if name == "NEWMARK" and isinstance(step, TransientStep):
            step.newmark = NewmarkControl(
                *tuple(map(float, item["data"][0][:2]))
            )
            return True

        if name == "DAMPING" and isinstance(step, TransientStep):
            step.damping = RayleighDamping(
                *tuple(map(float, item["data"][0][:2]))
            )
            return True

        if name == "WRITEEVERY" and isinstance(step, TransientStep):
            step.write_every = int(item["data"][0][0])
            return True

        if name == "INERTIARELIEF" and isinstance(step, StaticStep):
            step.inertia_relief = True
            return True

        if name == "REBALANCELOADS" and isinstance(step, StaticStep):
            step.rebalance_loads = True
            return True

        return False

    @staticmethod
    def _node_region(regions: RegionRepository, name: str) -> NodeRegion:
        if name in regions.nodes:
            return regions.nodes[name]
        return regions.nodes.add(NodeRegion(name))

    @staticmethod
    def _element_region(
        regions: RegionRepository,
        name: str,
    ) -> ElementRegion:
        if name in regions.elements:
            return regions.elements[name]
        return regions.elements.add(ElementRegion(name))

    @staticmethod
    def _load_collector(project: "Project", name: str) -> LoadCollector:
        if name in project.load_collectors:
            return project.load_collectors[name]
        return project.load_collectors.add(LoadCollector(name))

    @staticmethod
    def _support_collector(
        project: "Project",
        name: str,
    ) -> SupportCollector:
        if name in project.support_collectors:
            return project.support_collectors[name]
        return project.support_collectors.add(SupportCollector(name))

    @classmethod
    def _int_key(cls, item: _Block, key: str) -> int | None:
        value = cls._key(item, key)
        return None if value is None else int(value)

    @classmethod
    def _float_key(cls, item: _Block, key: str) -> float | None:
        value = cls._key(item, key)
        return None if value is None else float(value)
