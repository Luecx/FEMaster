# FEMaster Python API

The Python API mirrors the semantic FEMaster model as a connected object graph.
All editable FEM concepts live below `femaster_api.model`; the package root only
re-exports the public API for concise user code.

## Ownership and references

`Project` owns shared definitions, assembly definitions, collectors and analysis
steps. `Part` owns part-local nodes, elements, regions, surfaces and sections.
The implicit default part is permanently stored as `project.parts[0]` and is
returned by `project.parts.default()`.

Cross-object relationships always store **objects, never names or IDs**. For
example:

- element connectivity contains `Node` objects,
- `NodeRegion` contains `Node` objects,
- a `NodalForce` target is `Node | NodeRegion`,
- a load amplitude is `Amplitude | None`,
- a load/support orientation is `CoordinateSystem | None`,
- a section stores its `ElementRegion`, `Material` and (where applicable)
  `Profile` / `CoordinateSystem` objects,
- an `Instance` stores its `Part`,
- a `Step` stores `LoadCollector` / `SupportCollector` objects.

IDs and names remain stable identities and native file tokens. They are used for
repository lookup and during import/export, but they are not substitutes for
objects inside model constructors. There are deliberately no generic
`EntityReference`, `NodeReference` or `ElementReference` aliases.

## Example

```python
from femaster_api import (
    Amplitude,
    ElementRegion,
    IsotropicElasticity,
    LoadCollector,
    Material,
    Node,
    NodeRegion,
    NodalForce,
    Project,
    StaticStep,
    Support,
    SupportCollector,
    T3,
    TrussSection,
)

project = Project("bar")
part = project.parts.default()

n1 = part.nodes.add(Node(1, 0.0, 0.0, 0.0))
n2 = part.nodes.add(Node(2, 1.0, 0.0, 0.0))
bar = part.elements.add(T3(1, (n1, n2)))

root = part.regions.add(NodeRegion("ROOT", (n1,)))
tip = part.regions.add(NodeRegion("TIP", (n2,)))
bar_region = part.regions.add(ElementRegion("BAR", (bar,)))

steel = project.materials.add(
    Material("STEEL", elasticity=IsotropicElasticity(210000.0, 0.3))
)
part.sections.add(TrussSection("BAR_SECTION", bar_region, steel, 100.0))

ramp = project.amplitudes.add(
    Amplitude("RAMP").add(0.0, 0.0).add(1.0, 1.0)
)

bc = project.support_collectors.add(SupportCollector("BC"))
bc.add(Support(root, (0.0, 0.0, 0.0)))

loads = project.load_collectors.add(LoadCollector("LOAD"))
loads.add(NodalForce(tip, (1000.0, 0.0, 0.0, 0.0, 0.0, 0.0), amplitude=ramp))

project.steps.add(StaticStep("STATIC", loads=(loads,), supports=(bc,)))
project.write("bar.inp")
```

Input reading resolves native IDs/names immediately back to those same public
objects:

```python
project = Project.read_inp("bar.inp")
```

Unsupported keyword blocks are retained on `project.unparsed_blocks` rather than
silently discarded.

## Package layout

The source is split by FEM concept and permits at most one **top-level** class per
Python file. Tightly scoped implementation types may be nested when they have no
independent model lifecycle, e.g. `Equation.Term` and `Connector.Type`.

```text
model/
  project.py
  node/
  element/
  surface/
  region/
  material/
  section/
  load/
  support/
  constraint/
  step/
    util/
  field/
  result/
```

There is deliberately no generic `mesh/`, `io/` or `model/project/` package.

## Results

Post-processing uses one format-independent hierarchy:

```text
Result -> Solution -> Frame -> Field
```

`Result.read_res()` and `Result.read_frd()` normalize native result formats into
that common representation. Result entity identifiers are output addresses, not
editable model-object relationships, so they remain the semantic IDs encoded by
the result file.
