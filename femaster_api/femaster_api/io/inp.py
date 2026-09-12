"""Reader for FEMaster/Abaqus-like keyword input decks.

The reader first separates the file into syntax-level KeywordBlock objects and
then maps supported blocks into the public Python model classes. Unknown blocks
are retained on Project.unparsed_blocks so unsupported syntax is visible instead
of being silently discarded.
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
from ..model.constraints import (
    Connector,
    Coupling,
    CouplingType,
    Equation,
    RigidBodyConstraint,
    Tie,
)
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
from ..model.sections import BeamSection, ShellSection, SolidSection, TrussSection
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
from ..repository import NamedRepository


@dataclass(slots=True)
class KeywordBlock:
    """One parsed keyword line and all following data rows."""

    name: str
    keys: dict[str, str] = field(default_factory=dict)
    data: list[list[str]] = field(default_factory=list)
    line: int = 0

    def key(self, name: str, default: str | None = None) -> str | None:
        return self.keys.get(name.upper(), default)


class InpReader:
    """Read FEMaster keyword decks into the public Project object model."""

    def read(self, path: str | Path) -> Project:
        """Read and parse one UTF-8 input file."""

        return self.parse(Path(path).read_text(encoding="utf-8"))

    def parse(self, text: str) -> Project:
        """Parse input text into a new Project."""

        parsed = self.blocks(text)
        model_name = "model"
        for item in parsed:
            if item.name == "MODEL":
                model_name = item.key("NAME", "model") or "model"
                break

        project = Project(model_name)
        current_part = project.parts.default()
        current_material: Material | None = None
        current_step = None
        in_assembly = False

        for item in parsed:
            name = item.name

            # --------------------------------------------------------------
            # Scope control
            # --------------------------------------------------------------

            if name in {"MODEL", "END"}:
                continue

            if name == "ASSEMBLY":
                in_assembly = True
                current_part = project.parts.default()
                continue

            if name == "ENDASSEMBLY":
                in_assembly = False
                continue

            if name == "PART":
                part_name = self._required(item, "NAME")
                current_part = project.parts.add(Part(part_name))
                current_material = None
                current_step = None
                in_assembly = False
                continue

            if name == "ENDPART":
                current_part = project.parts.default()
                current_material = None
                continue

            if name == "INSTANCE":
                part_name = self._required(item, "PART")
                instance_name = self._required(item, "NAME")
                translation = None
                rotation = None
                if item.data and len(item.data[0]) == 3:
                    translation = tuple(map(float, item.data[0]))
                rotation_row = item.data[1] if len(item.data) > 1 else (
                    item.data[0] if item.data and len(item.data[0]) == 7 else None
                )
                if rotation_row is not None:
                    values = tuple(map(float, rotation_row))
                    rotation = (values[0:3], values[3:6], values[6])
                project.instances.add(
                    Instance(instance_name, part_name, translation=translation, rotation=rotation)
                )
                continue

            # --------------------------------------------------------------
            # Analysis-step subcommands are interpreted before root commands
            # with the same spelling, such as LOADS and SUPPORTS.
            # --------------------------------------------------------------

            if name == "LOADCASE":
                current_step = self._create_step(item)
                project.steps.add(current_step)
                continue

            if current_step is not None and self._read_step_subcommand(current_step, item):
                continue

            # --------------------------------------------------------------
            # Part-local topology
            # --------------------------------------------------------------

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

                nset = item.key("NSET")
                if nset:
                    regions = project.regions if in_assembly else current_part.regions
                    self._node_region(regions, nset).add(*node_ids)
                continue

            if name == "ELEMENT":
                type_name = self._required(item, "TYPE").upper()
                element_class = ELEMENT_TYPES.get(type_name)
                element_ids: list[int] = []

                for row in item.data:
                    values = [value for value in row if value != ""]
                    if len(values) < 2:
                        raise ValueError(f"ELEMENT at line {item.line} requires id and connectivity")
                    id = int(values[0])
                    nodes = tuple(int(value) for value in values[1:])
                    if element_class is None:
                        element = Element(id, nodes)
                        element.type_name = type_name
                    else:
                        element = element_class(id, nodes)
                    current_part.elements.add(element)
                    element_ids.append(id)

                elset = item.key("ELSET")
                if elset:
                    regions = project.regions if in_assembly else current_part.regions
                    self._element_region(regions, elset).add(*element_ids)
                continue

            if name == "NSET":
                region_name = item.key("NAME") or item.key("NSET")
                if region_name is None:
                    raise ValueError(f"NSET at line {item.line} requires NAME")
                regions = project.regions if in_assembly else current_part.regions
                region = self._node_region(regions, region_name)
                region.add(*(self._token(value) for row in item.data for value in row if value != ""))
                continue

            if name == "ELSET":
                region_name = item.key("NAME") or item.key("ELSET")
                if region_name is None:
                    raise ValueError(f"ELSET at line {item.line} requires NAME")
                regions = project.regions if in_assembly else current_part.regions
                region = self._element_region(regions, region_name)
                region.add(*(self._token(value) for row in item.data for value in row if value != ""))
                continue

            if name == "SURFACE":
                surface_name = item.key("NAME") or item.key("SFSET") or item.key("SURFACE")
                if surface_name is None:
                    raise ValueError(f"SURFACE at line {item.line} requires NAME or SFSET")
                surface = Surface(surface_name)
                for row in item.data:
                    if len(row) < 2:
                        continue
                    surface.add(self._token(row[0]), self._surface_side(row[1]))
                surfaces: NamedRepository[Surface] = project.surfaces if in_assembly else current_part.surfaces
                surfaces.add(surface)
                continue

            # --------------------------------------------------------------
            # Shared material/profile/coordinate definitions
            # --------------------------------------------------------------

            if name == "MATERIAL":
                current_material = Material(self._required(item, "NAME"))
                project.materials.add(current_material)
                continue

            if name == "ELASTIC" and current_material is not None:
                values = [float(value) for row in item.data for value in row if value != ""]
                type_name = (item.key("TYPE", "ISOTROPIC") or "ISOTROPIC").upper()
                if type_name == "ISOTROPIC":
                    current_material.elasticity = IsotropicElasticity(values[0], values[1])
                elif type_name == "GENISO":
                    current_material.elasticity = GeneralizedIsotropicElasticity(*values[:3])
                elif type_name == "ENGINEERINGCONSTANTS":
                    current_material.elasticity = OrthotropicElasticity(*values[:9])
                elif type_name == "ABD":
                    current_material.elasticity = ABDElasticity(values)
                else:
                    project.unparsed_blocks.append(item)
                continue

            if name == "DENSITY" and current_material is not None:
                current_material.density = float(item.data[0][0])
                continue

            if name == "THERMALEXPANSION" and current_material is not None:
                current_material.thermal_expansion = float(item.data[0][0])
                continue

            if name == "PROFILE":
                values = [float(value) for value in item.data[0]]
                values.extend([0.0] * (9 - len(values)))
                project.profiles.add(Profile(self._required(item, "NAME"), *values[:9]))
                continue

            if name == "ORIENTATION":
                values = [float(value) for value in item.data[0]]
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
                    y_axis = values[3:6] if len(values) >= 6 else None
                    z_axis = values[6:9] if len(values) >= 9 else None
                    project.coordinate_systems.add(
                        RectangularCoordinateSystem(orientation_name, values[0:3], y_axis, z_axis)
                    )
                continue

            # --------------------------------------------------------------
            # Part-local sections
            # --------------------------------------------------------------

            if name in {"SOLIDSECTION", "SHELLSECTION", "BEAMSECTION", "TRUSSSECTION"}:
                section_name = item.key("NAME") or f"{name}_{len(current_part.sections)}"
                elset = self._required(item, "ELSET")
                material = self._required(item, "MATERIAL")

                if name == "SOLIDSECTION":
                    section = SolidSection(section_name, elset, material, item.key("ORIENTATION"))
                elif name == "SHELLSECTION":
                    section = ShellSection(
                        section_name,
                        elset,
                        material,
                        float(item.data[0][0]),
                        item.key("ORIENTATION"),
                    )
                elif name == "BEAMSECTION":
                    orientation = tuple(map(float, item.data[0])) if item.data else None
                    section = BeamSection(
                        section_name,
                        elset,
                        material,
                        self._required(item, "PROFILE"),
                        orientation,
                    )
                else:
                    section = TrussSection(section_name, elset, material, float(item.data[0][0]))

                current_part.sections.add(section)
                continue

            # --------------------------------------------------------------
            # Global fields, amplitudes and features
            # --------------------------------------------------------------

            if name == "AMPLITUDE":
                interpolation = AmplitudeInterpolation[(item.key("TYPE", "LINEAR") or "LINEAR").upper()]
                amplitude = Amplitude(self._required(item, "NAME"), interpolation=interpolation)
                for row in item.data:
                    values = [float(value) for value in row if value != ""]
                    for start in range(0, len(values), 2):
                        if start + 1 < len(values):
                            amplitude.add(values[start], values[start + 1])
                project.amplitudes.add(amplitude)
                continue

            if name == "FIELD":
                domain_name = (item.key("TYPE") or item.key("DOMAIN") or "UNKNOWN").upper()
                try:
                    domain = FieldDomain[domain_name]
                except KeyError:
                    domain = FieldDomain.UNKNOWN
                cols = int(item.key("COLS", "0") or 0)
                field_item = Field(self._required(item, "NAME"), domain, [f"C{i + 1}" for i in range(cols)])
                index_cols = 2 if domain in {FieldDomain.ELEMENT_NODAL, FieldDomain.ELEMENT_IP, FieldDomain.ELEMENT_MP} else 1
                for row in item.data:
                    if len(row) < index_cols:
                        continue
                    keys = tuple(self._token(value) for value in row[:index_cols])
                    key = keys[0] if len(keys) == 1 else keys
                    field_item.set(key, (float(value) for value in row[index_cols:]))
                project.fields.add(field_item)
                continue

            if name == "POINTMASS":
                values = [float(value) for value in item.data[0]]
                values.extend([0.0] * (10 - len(values)))
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

            # --------------------------------------------------------------
            # Boundary-condition and load collectors
            # --------------------------------------------------------------

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
                collector = self._load_collector(project, self._required(item, "LOAD_COLLECTOR"))
                amplitude = item.key("AMPLITUDE")
                orientation = item.key("ORIENTATION")

                if name == "CLOAD":
                    for row in item.data:
                        collector.add(NodalForce(self._token(row[0]), map(float, row[1:7]), orientation=orientation, amplitude=amplitude))
                elif name == "DLOAD":
                    for row in item.data:
                        collector.add(SurfaceTraction(self._token(row[0]), map(float, row[1:4]), orientation=orientation, amplitude=amplitude))
                elif name == "PLOAD":
                    for row in item.data:
                        collector.add(PressureLoad(self._token(row[0]), float(row[1]), amplitude=amplitude))
                elif name == "VLOAD":
                    for row in item.data:
                        collector.add(VolumeLoad(self._token(row[0]), map(float, row[1:4]), orientation=orientation, amplitude=amplitude))
                elif name == "TLOAD":
                    collector.add(ThermalLoad(
                        self._required(item, "TEMPERATUREFIELD"),
                        float(item.key("REFERENCETEMPERATURE", "0") or 0),
                    ))
                else:
                    values = item.data[0]
                    numbers = [float(value) for value in values[1:]]
                    numbers.extend([0.0] * (12 - len(numbers)))
                    collector.add(InertialLoad(
                        values[0],
                        center=numbers[0:3],
                        center_acceleration=numbers[3:6],
                        omega=numbers[6:9],
                        alpha=numbers[9:12],
                        consider_point_masses=(item.key("CONSIDER_POINT_MASSES", "0") == "1"),
                    ))
                continue

            # --------------------------------------------------------------
            # Assembly constraints
            # --------------------------------------------------------------

            if name == "RBM":
                project.constraints.add(RigidBodyConstraint(self._required(item, "ELSET")))
                continue

            if name == "COUPLING":
                slave = item.key("SFSET") or item.key("SLAVE")
                if slave is None:
                    raise ValueError(f"COUPLING at line {item.line} requires SLAVE or SFSET")
                dofs = [int(value) for value in item.data[0]] if item.data else (1, 1, 1, 1, 1, 1)
                project.constraints.add(Coupling(
                    self._required(item, "MASTER"),
                    slave,
                    type=CouplingType[(item.key("TYPE", "KINEMATIC") or "KINEMATIC").upper()],
                    dofs=dofs,
                    slave_is_surface=item.key("SFSET") is not None,
                ))
                continue

            if name == "CONNECTOR":
                project.constraints.add(Connector(
                    self._required(item, "TYPE"),
                    self._required(item, "NSET1"),
                    self._required(item, "NSET2"),
                    self._required(item, "COORDINATESYSTEM"),
                ))
                continue

            if name == "TIE":
                project.constraints.add(Tie(
                    self._required(item, "MASTER"),
                    self._required(item, "SLAVE"),
                    adjust=(item.key("ADJUST", "YES") or "YES").upper() == "YES",
                    distance=float(item.key("DISTANCE")) if item.key("DISTANCE") is not None else None,
                ))
                continue

            if name == "EQUATION":
                flat = [value for row in item.data for value in row if value != ""]
                if flat:
                    count = int(flat[0])
                    equation = Equation()
                    values = flat[1:]
                    for index in range(count):
                        start = 3 * index
                        equation.add(self._token(values[start]), int(values[start + 1]), float(values[start + 2]))
                    project.constraints.add(equation)
                continue

            project.unparsed_blocks.append(item)

        return project

    def blocks(self, text: str) -> list[KeywordBlock]:
        """Split raw deck text into keyword blocks without semantic interpretation."""

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

            if current is None:
                continue

            current.data.append([token.strip() for token in raw.split(",")])

        return result

    # ------------------------------------------------------------------
    # Analysis helpers
    # ------------------------------------------------------------------

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
            values = tuple(map(float, item.data[0][:3]))
            step.time = TimeControl(*values)
            return True
        if name == "NEWMARK" and isinstance(step, TransientStep):
            values = tuple(map(float, item.data[0][:2]))
            step.newmark = NewmarkControl(*values)
            return True
        if name == "DAMPING" and isinstance(step, TransientStep):
            values = tuple(map(float, item.data[0][:2]))
            step.damping = RayleighDamping(*values)
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

    # ------------------------------------------------------------------
    # Repository helpers
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # Syntax helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_keyword(line: str, line_number: int) -> KeywordBlock:
        tokens = [token.strip() for token in line[1:].split(",")]
        name = tokens[0].upper()
        keys: dict[str, str] = {}
        for token in tokens[1:]:
            if "=" in token:
                key, value = token.split("=", 1)
                keys[key.strip().upper()] = value.strip()
        return KeywordBlock(name, keys, [], line_number)

    @staticmethod
    def _required(item: KeywordBlock, key: str) -> str:
        value = item.key(key)
        if value is None or value == "":
            raise ValueError(f"{item.name} at line {item.line} requires {key}")
        return value

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
