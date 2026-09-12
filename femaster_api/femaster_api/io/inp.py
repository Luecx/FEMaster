"""Reader for FEMaster/Abaqus-like keyword input decks.

The reader first separates source text into syntax-level KeywordBlock objects and
then populates the same public Project model that users construct manually.
Unsupported blocks are retained on ``Project.unparsed_blocks`` rather than being
silently discarded.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..fields import Field, FieldDomain
from ..model.assembly import Instance, Part
from ..model.boundary import (
    Amplitude,
    AmplitudeInterpolation,
    InertialLoad,
    LoadCollector,
    NodalForce,
    PressureLoad,
    Support,
    SupportCollector,
    SurfaceTraction,
    ThermalLoad,
    VolumeLoad,
)
from ..model.constraints import Connector, Coupling, CouplingType, Equation, RigidBodyConstraint, Tie
from ..model.coordinates import CylindricalCoordinateSystem, RectangularCoordinateSystem
from ..model.features import PointMass
from ..model.materials import (
    ABDElasticity,
    GeneralizedIsotropicElasticity,
    IsotropicElasticity,
    Material,
    OrthotropicElasticity,
    Profile,
)
from ..model.mesh import ELEMENT_TYPES, Element, Node, Surface
from ..model.regions import ElementRegion, NodeRegion, RegionRepository
from ..model.sections import (
    ABDShellSection,
    BeamSection,
    MassSection,
    RotaryInertiaSection,
    SectionRepository,
    ShellSection,
    SolidSection,
    SpringSection,
    TrussSection,
)
from ..model.steps import (
    BucklingStep,
    ConstraintMethod,
    ModalStep,
    NewmarkControl,
    NonlinearStaticStep,
    RayleighDamping,
    SolverControl,
    SolverDevice,
    SolverMethod,
    StaticStep,
    TimeControl,
    TransientStep,
)
from ..project import Project


@dataclass(slots=True)
class KeywordBlock:
    """One parsed keyword line and all following data rows."""

    name: str
    keys: dict[str, str] = field(default_factory=dict)
    data: list[list[str]] = field(default_factory=list)
    line: int = 0

    def key(self, name: str, default: str | None = None) -> str | None:
        return self.keys.get(name.upper().replace(" ", ""), default)


class InpReader:
    """Read FEMaster keyword decks into the public Project object model."""

    def read(self, path: str | Path) -> Project:
        """Read one UTF-8 input file."""

        return self.parse(Path(path).read_text(encoding="utf-8"))

    def parse(self, text: str) -> Project:
        """Parse input text into a new Project."""

        parsed = self.blocks(text)
        project = Project(self._model_name(parsed))

        current_part = project.parts.default()
        current_material: Material | None = None
        current_step = None
        in_assembly = False

        for item in parsed:
            name = item.name

            # -----------------------------------------------------------------
            # Scope control
            # -----------------------------------------------------------------

            if name in {"MODEL", "END"}:
                continue

            if name == "PART":
                current_part = project.parts.add(Part(self._required(item, "NAME")))
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
                project.instances.add(self._read_instance(item))
                continue

            # -----------------------------------------------------------------
            # Analysis steps and their child commands
            # -----------------------------------------------------------------

            if name == "LOADCASE":
                current_step = self._create_step(item)
                project.steps.add(current_step)
                continue

            if current_step is not None and self._read_step_subcommand(current_step, item):
                continue

            # -----------------------------------------------------------------
            # Part-local topology and named regions
            # -----------------------------------------------------------------

            regions = project.regions if in_assembly else current_part.regions

            if name == "NODE":
                node_ids: list[int] = []
                for row in item.data:
                    values = [value for value in row if value != ""]
                    if len(values) < 3:
                        raise ValueError(f"NODE at line {item.line} requires id, x, y")

                    node = Node(
                        int(values[0]),
                        float(values[1]),
                        float(values[2]),
                        float(values[3]) if len(values) > 3 else 0.0,
                    )
                    current_part.nodes.add(node)
                    node_ids.append(node.id)

                if item.key("NSET"):
                    self._node_region(regions, item.key("NSET") or "").add(*node_ids)
                continue

            if name == "ELEMENT":
                type_name = self._required(item, "TYPE").upper()
                element_class = ELEMENT_TYPES.get(type_name)
                element_ids: list[int] = []

                for row in item.data:
                    values = [value for value in row if value != ""]
                    if len(values) < 2:
                        raise ValueError(f"ELEMENT at line {item.line} requires id and connectivity")

                    element_id = int(values[0])
                    connectivity = tuple(int(value) for value in values[1:])

                    if element_class is None:
                        element = Element(element_id, connectivity)
                        element.type_name = type_name
                    else:
                        element = element_class(element_id, connectivity)

                    current_part.elements.add(element)
                    element_ids.append(element_id)

                if item.key("ELSET"):
                    self._element_region(regions, item.key("ELSET") or "").add(*element_ids)
                continue

            if name == "NSET":
                region_name = item.key("NAME") or item.key("NSET")
                if not region_name:
                    raise ValueError(f"NSET at line {item.line} requires NAME")
                self._node_region(regions, region_name).add(*self._flat_tokens(item))
                continue

            if name == "ELSET":
                region_name = item.key("NAME") or item.key("ELSET")
                if not region_name:
                    raise ValueError(f"ELSET at line {item.line} requires NAME")
                self._element_region(regions, region_name).add(*self._flat_tokens(item))
                continue

            if name == "SURFACE":
                surface_name = item.key("NAME") or item.key("SFSET") or item.key("SURFACE")
                if not surface_name:
                    raise ValueError(f"SURFACE at line {item.line} requires NAME or SFSET")

                surface_type = (item.key("TYPE", "ELEMENT") or "ELEMENT").upper()
                if surface_type == "NODE":
                    region = self._node_region(regions, surface_name)
                    for row in item.data:
                        if row and row[0]:
                            region.add(self._token(row[0]))
                    continue

                surface = Surface(surface_name)
                for row in item.data:
                    if len(row) >= 2 and row[0] and row[1]:
                        surface.add(self._token(row[0]), self._surface_side(row[1]))

                target = project.surfaces if in_assembly else current_part.surfaces
                target.add(surface)
                continue

            # -----------------------------------------------------------------
            # Shared material, profile and coordinate-system definitions
            # -----------------------------------------------------------------

            if name == "MATERIAL":
                current_material = Material(self._required(item, "NAME"))
                project.materials.add(current_material)
                continue

            if name == "ELASTIC" and current_material is not None:
                values = self._flat_floats(item)
                type_name = (item.key("TYPE", "ISOTROPIC") or "ISOTROPIC").upper()

                if type_name == "ISOTROPIC":
                    current_material.elasticity = IsotropicElasticity(values[0], values[1])
                elif type_name == "GENISO":
                    current_material.elasticity = GeneralizedIsotropicElasticity(*values[:3])
                elif type_name == "ENGINEERINGCONSTANTS":
                    current_material.elasticity = OrthotropicElasticity(*values[:9])
                elif type_name == "ABD":
                    current_material.elasticity = ABDElasticity(tuple(values))
                else:
                    project.unparsed_blocks.append(item)
                continue

            if name == "DENSITY" and current_material is not None:
                current_material.density = self._flat_floats(item)[0]
                continue

            if name == "THERMALEXPANSION" and current_material is not None:
                current_material.thermal_expansion = self._flat_floats(item)[0]
                continue

            if name == "PROFILE":
                values = self._flat_floats(item)
                values.extend([0.0] * max(0, 9 - len(values)))
                project.profiles.add(Profile(self._required(item, "NAME"), *values[:9]))
                continue

            if name == "ORIENTATION":
                values = self._flat_floats(item)
                orientation_name = self._required(item, "NAME")
                type_name = (item.key("TYPE", "RECTANGULAR") or "RECTANGULAR").upper()

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

            # -----------------------------------------------------------------
            # Sections and point-element properties
            # -----------------------------------------------------------------

            if name in {"SOLIDSECTION", "SHELLSECTION", "BEAMSECTION", "TRUSSSECTION"}:
                self._read_part_section(project, current_part.sections, item)
                continue

            if name in {"MASS", "ROTARYINERTIA", "SPRING"}:
                target_sections = project.sections if in_assembly else current_part.sections
                self._read_point_section(target_sections, item)
                continue

            # -----------------------------------------------------------------
            # Global fields, amplitudes and non-topological features
            # -----------------------------------------------------------------

            if name == "AMPLITUDE":
                interpolation_name = (item.key("TYPE", "LINEAR") or "LINEAR").upper()
                interpolation = AmplitudeInterpolation[interpolation_name]
                amplitude = Amplitude(self._required(item, "NAME"), interpolation=interpolation)
                values = self._flat_floats(item)
                for index in range(0, len(values) - 1, 2):
                    amplitude.add(values[index], values[index + 1])
                project.amplitudes.add(amplitude)
                continue

            if name == "FIELD":
                project.fields.add(self._read_model_field(item))
                continue

            if name == "POINTMASS":
                values = self._flat_floats(item)
                values.extend([0.0] * max(0, 10 - len(values)))
                project.features.add(
                    PointMass(
                        self._required(item, "NSET"),
                        values[0],
                        inertia=values[1:4],
                        spring=values[4:7],
                        rotational_spring=values[7:10],
                    )
                )
                continue

            # -----------------------------------------------------------------
            # Supports and loads
            # -----------------------------------------------------------------

            if name == "SUPPORT":
                collector = self._support_collector(project, self._required(item, "SUPPORT_COLLECTOR"))
                orientation = item.key("ORIENTATION")

                for row in item.data:
                    if not row:
                        continue
                    values = [None if value == "" else float(value) for value in row[1:]]
                    collector.add(Support(self._token(row[0]), values, orientation=orientation))
                continue

            if name in {"CLOAD", "DLOAD", "PLOAD", "VLOAD", "TLOAD", "INERTIALOAD"}:
                self._read_load(project, item)
                continue

            # -----------------------------------------------------------------
            # Assembly-level constraints
            # -----------------------------------------------------------------

            if name == "RBM":
                project.constraints.add(RigidBodyConstraint(self._required(item, "ELSET")))
                continue

            if name == "COUPLING":
                slave = item.key("SFSET") or item.key("SLAVE")
                if not slave:
                    raise ValueError(f"COUPLING at line {item.line} requires SLAVE or SFSET")

                dofs = [int(value) for value in item.data[0]] if item.data else (1, 1, 1, 1, 1, 1)
                project.constraints.add(
                    Coupling(
                        self._required(item, "MASTER"),
                        slave,
                        type=CouplingType[(item.key("TYPE", "KINEMATIC") or "KINEMATIC").upper()],
                        dofs=dofs,
                        slave_is_surface=item.key("SFSET") is not None,
                    )
                )
                continue

            if name == "CONNECTOR":
                project.constraints.add(
                    Connector(
                        self._required(item, "TYPE"),
                        self._required(item, "NSET1"),
                        self._required(item, "NSET2"),
                        self._required(item, "COORDINATESYSTEM"),
                    )
                )
                continue

            if name == "TIE":
                project.constraints.add(
                    Tie(
                        self._required(item, "MASTER"),
                        self._required(item, "SLAVE"),
                        adjust=(item.key("ADJUST", "YES") or "YES").upper() == "YES",
                        distance=float(item.key("DISTANCE")) if item.key("DISTANCE") is not None else None,
                    )
                )
                continue

            if name == "EQUATION":
                tokens = [value for row in item.data for value in row if value != ""]
                if tokens:
                    term_count = int(tokens[0])
                    equation = Equation()
                    values = tokens[1:]
                    for index in range(term_count):
                        start = 3 * index
                        equation.add(
                            self._token(values[start]),
                            int(values[start + 1]),
                            float(values[start + 2]),
                        )
                    project.constraints.add(equation)
                continue

            project.unparsed_blocks.append(item)

        return project

    def blocks(self, text: str) -> list[KeywordBlock]:
        """Split source text into keyword blocks without semantic interpretation."""

        result: list[KeywordBlock] = []
        current: KeywordBlock | None = None

        for line_number, raw in enumerate(text.splitlines(), start=1):
            stripped = raw.strip()
            if not stripped or stripped.startswith("**"):
                continue

            if stripped.startswith("*"):
                current = self._parse_keyword(stripped, line_number)
                result.append(current)
                continue

            if current is not None:
                current.data.append([token.strip() for token in raw.split(",")])

        return result

    # -------------------------------------------------------------------------
    # Model-definition helpers
    # -------------------------------------------------------------------------

    @staticmethod
    def _model_name(items: list[KeywordBlock]) -> str:
        for item in items:
            if item.name == "MODEL":
                return item.key("NAME", "model") or "model"
        return "model"

    def _read_instance(self, item: KeywordBlock) -> Instance:
        translation = None
        rotation = None

        if item.data and len([value for value in item.data[0] if value != ""]) == 3:
            translation = tuple(float(value) for value in item.data[0])

        rotation_row = None
        if len(item.data) > 1:
            rotation_row = item.data[1]
        elif item.data and len([value for value in item.data[0] if value != ""]) == 7:
            rotation_row = item.data[0]

        if rotation_row is not None:
            values = tuple(float(value) for value in rotation_row if value != "")
            rotation = (values[0:3], values[3:6], values[6])

        return Instance(
            self._required(item, "NAME"),
            self._required(item, "PART"),
            translation=translation,
            rotation=rotation,
        )

    def _read_part_section(self, project: Project, sections: SectionRepository, item: KeywordBlock) -> None:
        name = item.name
        section_name = item.key("NAME") or f"{name}_{len(sections)}"
        elset = self._required(item, "ELSET")

        if name == "SOLIDSECTION":
            sections.add(
                SolidSection(
                    section_name,
                    elset,
                    self._required(item, "MATERIAL"),
                    item.key("ORIENTATION"),
                )
            )
            return

        if name == "TRUSSSECTION":
            sections.add(
                TrussSection(
                    section_name,
                    elset,
                    self._required(item, "MATERIAL"),
                    self._flat_floats(item)[0],
                )
            )
            return

        if name == "BEAMSECTION":
            orientation = self._flat_floats(item)
            sections.add(
                BeamSection(
                    section_name,
                    elset,
                    self._required(item, "MATERIAL"),
                    self._required(item, "PROFILE"),
                    orientation[:3],
                )
            )
            return

        section_type = (item.key("TYPE", "INTEGRATED") or "INTEGRATED").upper()
        csys_axis = int(item.key("CSYSAXIS", "1") or 1)

        if section_type == "ABD":
            values = self._flat_floats(item)
            if len(values) != 40:
                raise ValueError(f"SHELLSECTION TYPE=ABD at line {item.line} requires 40 values")

            sections.add(
                ABDShellSection(
                    section_name,
                    elset,
                    values[:36],
                    values[36:40],
                    thickness=float(item.key("THICKNESS", "1.0") or 1.0),
                    material=item.key("MATERIAL"),
                    orientation=item.key("ORIENTATION"),
                    csys_axis=csys_axis,
                )
            )
            return

        thickness_values = self._flat_floats(item)
        thickness = thickness_values[0] if thickness_values else float(item.key("THICKNESS", "1.0") or 1.0)
        sections.add(
            ShellSection(
                section_name,
                elset,
                self._required(item, "MATERIAL"),
                thickness,
                item.key("ORIENTATION"),
                csys_axis,
            )
        )

    def _read_point_section(self, sections: SectionRepository, item: KeywordBlock) -> None:
        section_name = item.key("NAME") or f"{item.name}_{len(sections)}"
        elset = self._required(item, "ELSET")
        values = self._flat_floats(item)

        if item.name == "MASS":
            sections.add(MassSection(section_name, elset, values[0]))
            return

        if item.name == "ROTARYINERTIA":
            values.extend([0.0] * max(0, 6 - len(values)))
            if any(value != 0.0 for value in values[3:6]):
                raise ValueError("ROTARY INERTIA products I12, I13 and I23 must be zero")
            sections.add(RotaryInertiaSection(section_name, elset, values[:3]))
            return

        if len(values) < 2:
            raise ValueError(f"SPRING at line {item.line} requires DOF and stiffness")
        sections.add(SpringSection(section_name, elset, int(values[0]), values[1]))

    def _read_model_field(self, item: KeywordBlock) -> Field:
        domain_name = (item.key("TYPE") or item.key("DOMAIN") or "UNKNOWN").upper()
        try:
            domain = FieldDomain[domain_name]
        except KeyError:
            domain = FieldDomain.UNKNOWN

        cols = int(item.key("COLS", "0") or 0)
        result = Field(
            self._required(item, "NAME"),
            domain,
            [f"C{index + 1}" for index in range(cols)],
        )

        index_cols = 2 if domain in {
            FieldDomain.ELEMENT_NODAL,
            FieldDomain.ELEMENT_IP,
            FieldDomain.ELEMENT_MP,
        } else 1

        for row in item.data:
            if len(row) < index_cols:
                continue
            indices = tuple(self._token(value) for value in row[:index_cols])
            key = indices[0] if len(indices) == 1 else indices
            result.set(key, (float(value) for value in row[index_cols:] if value != ""))

        return result

    def _read_load(self, project: Project, item: KeywordBlock) -> None:
        collector = self._load_collector(project, self._required(item, "LOAD_COLLECTOR"))
        amplitude = item.key("AMPLITUDE")
        orientation = item.key("ORIENTATION")

        if item.name == "CLOAD":
            for row in item.data:
                collector.add(
                    NodalForce(
                        self._token(row[0]),
                        (float(value) for value in row[1:7]),
                        orientation=orientation,
                        amplitude=amplitude,
                    )
                )
            return

        if item.name == "DLOAD":
            for row in item.data:
                collector.add(
                    SurfaceTraction(
                        self._token(row[0]),
                        (float(value) for value in row[1:4]),
                        orientation=orientation,
                        amplitude=amplitude,
                    )
                )
            return

        if item.name == "PLOAD":
            for row in item.data:
                collector.add(PressureLoad(self._token(row[0]), float(row[1]), amplitude=amplitude))
            return

        if item.name == "VLOAD":
            for row in item.data:
                collector.add(
                    VolumeLoad(
                        self._token(row[0]),
                        (float(value) for value in row[1:4]),
                        orientation=orientation,
                        amplitude=amplitude,
                    )
                )
            return

        if item.name == "TLOAD":
            collector.add(
                ThermalLoad(
                    self._required(item, "TEMPERATUREFIELD"),
                    float(item.key("REFERENCETEMPERATURE", "0") or 0),
                )
            )
            return

        row = item.data[0]
        values = [float(value) for value in row[1:] if value != ""]
        values.extend([0.0] * max(0, 12 - len(values)))
        collector.add(
            InertialLoad(
                row[0],
                center=values[0:3],
                center_acceleration=values[3:6],
                omega=values[6:9],
                alpha=values[9:12],
                consider_point_masses=(item.key("CONSIDER_POINT_MASSES", "0") == "1"),
            )
        )

    # -------------------------------------------------------------------------
    # Analysis helpers
    # -------------------------------------------------------------------------

    def _create_step(self, item: KeywordBlock):
        type_name = self._required(item, "TYPE").upper()
        name = item.key("NAME") or f"STEP_{item.line}"

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

        raise ValueError(f"unsupported LOADCASE TYPE at line {item.line}: {type_name}")

    def _read_step_subcommand(self, step, item: KeywordBlock) -> bool:
        name = item.name

        if name == "SUPPORTS":
            step.supports = tuple(value for row in item.data for value in row if value)
            return True
        if name == "LOADS":
            step.loads = tuple(value for row in item.data for value in row if value)
            return True
        if name == "SOLVER":
            step.solver = SolverControl(
                SolverDevice[(item.key("DEVICE", "CPU") or "CPU").upper()],
                SolverMethod[(item.key("METHOD", "DIRECT") or "DIRECT").upper()],
            )
            return True
        if name == "CONSTRAINTMETHOD":
            step.constraint_method = ConstraintMethod[self._required(item, "TYPE").upper()]
            return True
        if name == "NUMEIGENVALUES" and hasattr(step, "number_of_modes"):
            step.number_of_modes = int(item.data[0][0])
            return True
        if name == "SIGMA" and hasattr(step, "sigma"):
            step.sigma = float(item.data[0][0])
            return True
        if name == "NONLINEAR" and isinstance(step, NonlinearStaticStep):
            step.control           = (item.key("CONTROL", step.control) or step.control).upper()
            step.increments        = self._int_key(item, "INCREMENTS")
            step.max_increments    = self._int_key(item, "MAX_INCREMENTS")
            step.initial_increment = self._float_key(item, "INITIAL_INCREMENT")
            step.minimum_increment = self._float_key(item, "MINIMUM_INCREMENT")
            step.maximum_increment = self._float_key(item, "MAXIMUM_INCREMENT")
            step.max_iterations    = self._int_key(item, "MAXITER")
            step.tolerance         = self._float_key(item, "TOL")
            return True
        if name == "TIME" and isinstance(step, TransientStep):
            step.time = TimeControl(*tuple(map(float, item.data[0][:3])))
            return True
        if name == "NEWMARK" and isinstance(step, TransientStep):
            step.newmark = NewmarkControl(*tuple(map(float, item.data[0][:2])))
            return True
        if name == "DAMPING" and isinstance(step, TransientStep):
            step.damping = RayleighDamping(*tuple(map(float, item.data[0][:2])))
            return True
        if name == "WRITEEVERY" and isinstance(step, TransientStep):
            step.write_every = int(item.data[0][0])
            return True
        if name == "INERTIARELIEF" and isinstance(step, StaticStep):
            step.inertia_relief = True
            return True
        if name == "REBALANCELOADS" and isinstance(step, StaticStep):
            step.rebalance_loads = True
            return True

        return False

    # -------------------------------------------------------------------------
    # Repository and syntax helpers
    # -------------------------------------------------------------------------

    @staticmethod
    def _node_region(regions: RegionRepository, name: str) -> NodeRegion:
        if name in regions.nodes:
            return regions.nodes[name]
        return regions.nodes.add(NodeRegion(name))

    @staticmethod
    def _element_region(regions: RegionRepository, name: str) -> ElementRegion:
        if name in regions.elements:
            return regions.elements[name]
        return regions.elements.add(ElementRegion(name))

    @staticmethod
    def _load_collector(project: Project, name: str) -> LoadCollector:
        if name in project.load_collectors:
            return project.load_collectors[name]
        return project.load_collectors.add(LoadCollector(name))

    @staticmethod
    def _support_collector(project: Project, name: str) -> SupportCollector:
        if name in project.support_collectors:
            return project.support_collectors[name]
        return project.support_collectors.add(SupportCollector(name))

    @staticmethod
    def _parse_keyword(line: str, line_number: int) -> KeywordBlock:
        tokens = [token.strip() for token in line[1:].split(",")]
        name = tokens[0].upper().replace(" ", "")
        keys: dict[str, str] = {}

        for token in tokens[1:]:
            if "=" in token:
                key, value = token.split("=", 1)
                keys[key.strip().upper().replace(" ", "")] = value.strip()

        return KeywordBlock(name, keys, [], line_number)

    @staticmethod
    def _required(item: KeywordBlock, key: str) -> str:
        value = item.key(key)
        if value is None or value == "":
            raise ValueError(f"{item.name} at line {item.line} requires {key}")
        return value

    @classmethod
    def _flat_tokens(cls, item: KeywordBlock) -> list[int | str]:
        return [cls._token(value) for row in item.data for value in row if value != ""]

    @staticmethod
    def _flat_floats(item: KeywordBlock) -> list[float]:
        return [float(value) for row in item.data for value in row if value != ""]

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

    @staticmethod
    def _int_key(item: KeywordBlock, key: str) -> int | None:
        value = item.key(key)
        return None if value is None else int(value)

    @staticmethod
    def _float_key(item: KeywordBlock, key: str) -> float | None:
        value = item.key(key)
        return None if value is None else float(value)
