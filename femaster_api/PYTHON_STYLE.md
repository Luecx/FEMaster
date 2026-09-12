# FEMaster Python API Style

These rules intentionally mirror the FEMaster C++ style: preserve finite-element
semantics first, make ownership explicit, and prefer readable local code over
framework abstractions.

## Source structure

Every model concept lives below `femaster_api/model/`. Do not create parallel
model hierarchies under generic packages such as `io`, `export`, `importers` or
`mesh`. `Project` lives directly in `model/project.py`. Nodes, elements,
surfaces, regions, materials, sections, loads, supports, constraints, steps,
fields and results each have their own semantic package.

A Python source file may define **at most one class**. This includes base classes,
concrete classes, repositories, enums and helper/control classes. Pure functions,
type aliases and lookup tables may share a module when they naturally belong to
that module and do not hide an additional class.

Use grouped filenames inside semantic packages. Examples are `step_modal.py`,
`step_static.py`, `region_node.py`, `surface_element.py`, `load_pressure.py` and
`section_shell.py`. Concrete FEM element classes use one file each, such as
`element_c3d8.py` and `element_mitc4.py`.

Step implementation details that are not analysis steps themselves belong below
`step/util/`, for example `solver_control.py`, `time_control.py` and
`rayleigh_damping.py`.

## Documentation

Every Python file starts with a substantial module docstring. It should explain
what the module represents, where the concept is owned, how identifiers are
interpreted, and any important native FEMaster export/import behavior. Tiny
one-line module docstrings are not acceptable for model modules.

Every public class and every non-trivial public method has a docstring. Comments
inside functions should mark logical phases or explain FEM semantics, parsing
state, identifier spaces, numerical meaning or non-obvious format behavior.
Do not add comments that merely restate an assignment.

## Ownership

`Project` owns shared/global definitions and assembly-level concepts. `Part`
owns part-local nodes, elements, regions, surfaces and sections. The implicit
default part exists only as `project.parts[0]`; `parts.default()` returns that
same object and it cannot be removed.

Loads and supports belong directly to their collectors. Do not keep a second
global copy. Persistent references between named objects use semantic names, not
repository positions.

Node and element IDs are native FEMaster IDs and must never be silently
renumbered. Integer subscription on an ID repository addresses the FEMaster ID;
explicit positional access uses `.at(index)`.

## Export

Every independently exportable class owns an `export()` method. Repositories own
`export()` when grouping/order is part of their responsibility. `Project.export()`
only composes these pieces in native scope/dependency order. Do not introduce a
serializer registry, visitor hierarchy, decorator-based dispatch or metaclass.

## Input reading

Input reading belongs to `Project` through `Project.read_inp()` and
`Project.read_inp_text()`. Parsing populates the same public model users construct
manually. Temporary syntax state should use ordinary private data structures,
not a second public DTO object hierarchy.

Successfully parsed but unsupported keyword blocks remain visible as unparsed
data. Never silently discard syntax that the public model cannot represent.

## Results

The canonical post-processing hierarchy is:

`Result -> Solution -> Frame -> Field`

`Solution` represents one solver/loadcase result. `Frame` owns a numeric ID, an
optional physical `value`, and fields. `FieldDomain` defines storage location and
`FieldType` defines semantic meaning. RES and FRD parsing must normalize into
these shared classes rather than define format-specific result models.

Result readers are class methods on `Result`: `read()`, `read_res()`,
`read_res_text()`, `read_frd()` and `read_frd_text()`.

## Validation

Validate invariants where invalid data first becomes meaningful: duplicate names
or IDs, element connectivity size, vector dimension, immutable semantic names,
field row width/address width, protected default-part removal and required native
keyword data. Do not add speculative validation that changes FEMaster semantics.

## Dependencies and abstractions

Prefer ordinary Python and the standard library. Add an abstraction only when it
represents a real FEMaster concept or substantial repeated behavior. Do not add
Pydantic, attrs, serialization frameworks, dependency injection or generic
repository infrastructure beyond the small repository primitives already used.

Tests must enforce the structural rules above. Do not claim tests pass unless
they were actually executed.
