# FEMaster Python API Style

These rules mirror the FEMaster C++ model: preserve finite-element semantics,
make ownership and relationships explicit, and prefer readable local code over
framework abstractions.

## Source structure

Every model concept lives below `femaster_api/model/`. Do not create parallel
model hierarchies under generic packages such as `io`, `export`, `importers` or
`mesh`. `Project` lives directly in `model/project.py`.

A Python source file defines at most one **top-level** model class. Tightly scoped
implementation types may be nested when they have no independent model lifecycle;
`Equation.Term` and `Connector.Type` are the intended examples. Do not split such
helpers into public modules merely to satisfy a mechanical one-class rule.

Use grouped filenames inside semantic packages, e.g. `step_modal.py`,
`region_node.py`, `surface_element.py`, `load_pressure.py` and
`section_shell.py`. Concrete FEM elements use one file each such as
`element_c3d8.py`.

Step controls that are not analyses themselves belong below `step/util/`.

## Documentation

Every Python file starts with a substantial module docstring explaining ownership,
relationship semantics and relevant native FEMaster behavior. Public classes and
non-trivial methods have docstrings. Inline comments should explain phases,
scopes, parsing resolution or non-obvious FEM semantics rather than restating
assignments.

## Ownership and object relationships

`Project` owns shared/global definitions and assembly concepts. `Part` owns
part-local nodes, elements, regions, surfaces and sections. The implicit default
part exists only as `project.parts[0]`.

**A reference to another model concept is always the actual Python object.**
Never accept or retain a name/ID string as a surrogate object reference. Write
the legal concrete types directly at the use site instead of defining generic
aliases. Examples:

```python
Node | NodeRegion
Element | ElementRegion
Surface | SurfaceRegion
Amplitude | None
CoordinateSystem | None
```

Consequently:

- elements contain `Node` objects, not node IDs;
- regions contain their entity objects;
- surfaces contain node/element objects or the corresponding region objects;
- loads/supports store concrete targets and shared amplitude/orientation objects;
- sections store `ElementRegion`, `Material`, `Profile` and orientation objects;
- instances store `Part` objects;
- constraints store concrete region/surface/node objects;
- steps store collector objects.

Intrinsic identities (`Node.id`, `Material.name`, `Region.name`, etc.) remain
native FEMaster identities and repository keys. Import/export may read/write those
values, but editable model constructors must not use them in place of objects.
There must be no `EntityReference`, `NodeReference` or `ElementReference` alias.

Loads and supports belong directly to their collectors; do not keep a second
global copy.

## Export

Every independently exportable class owns `export()`. Serialization is the point
where object relationships are deliberately reduced to native IDs/names.
Repositories own export only when grouping/order is part of their responsibility.
`Project.export()` composes the pieces in native dependency order. Do not add a
serializer registry, visitor hierarchy, decorator dispatch or metaclass.

## Input reading

Input reading belongs to `Project.read_inp()` / `Project.read_inp_text()`. The
reader may temporarily see textual IDs/names because that is what the file stores,
but each supported token must be resolved immediately to the corresponding model
object before constructing the public object graph. An unresolved supported
reference is an error; do not preserve it as a fallback string.

Unsupported keyword blocks remain visible in `project.unparsed_blocks`.

## Results

The canonical post-processing hierarchy is:

`Result -> Solution -> Frame -> Field`

Result-file IDs are output addresses rather than editable model relationships and
may remain semantic identifiers from the file. RES and FRD normalize into the
shared result classes.

## Validation

Validate invalid relationships at construction time: wrong object types,
duplicate names/IDs, element connectivity size, vector dimensions, field widths,
protected default-part removal and required native data. Runtime checks should
make accidental string/ID surrogates fail immediately.

## Dependencies and abstractions

Prefer ordinary Python and the standard library. Add abstractions only for real
FEMaster concepts or substantial repeated behavior. Do not add Pydantic, attrs,
serialization frameworks, dependency injection or generic reference wrappers.
Tests must enforce these structural rules.
