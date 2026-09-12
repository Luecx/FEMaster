# FEMaster Python API

The Python API mirrors the semantic FEMaster model instead of exposing a flat
collection of loosely related objects. All editable FEM concepts live below
`femaster_api.model`; the package root only re-exports the public API for concise
user code.

## Ownership

`Project` owns shared definitions, assembly definitions, collectors and analysis
steps. `Part` owns part-local nodes, elements, regions, surfaces and sections.
The implicit default part is permanently stored as `project.parts[0]` and is
returned by `project.parts.default()`; it is never duplicated on `Project`.

Nodes and elements preserve their native sparse FEMaster IDs. Named repositories
preserve insertion order and support both positional and name-based lookup, while
repository positions are never used as persistent model identifiers.

## Package layout

The source is intentionally split by FEM concept and follows a strict rule of at
most one class per Python file. Concrete variants are grouped by filename prefix:

```text
model/
  project.py
  node/
    node.py
    node_repository.py
  element/
    element.py
    element_repository.py
    element_c3d8.py
    element_s4.py
    ...
  surface/
    surface.py
    surface_element.py
    surface_node.py
    surface_repository.py
  region/
    region.py
    region_node.py
    region_element.py
    region_surface.py
    region_line.py
    region_repository.py
  step/
    step.py
    step_static.py
    step_modal.py
    step_buckling.py
    step_nonlinear_static.py
    step_transient.py
    step_repository.py
    util/
      solver_control.py
      solver_device.py
      solver_method.py
      constraint_method.py
      time_control.py
      newmark_control.py
      rayleigh_damping.py
  field/
    field.py
    field_domain.py
    field_type.py
    field_repository.py
  result/
    result.py
    result_solution.py
    result_frame.py
```

The same convention is used for materials, sections, loads, supports,
constraints, amplitudes, profiles, coordinate systems, features, parts and
instances. There is deliberately no generic `mesh/`, `io/` or `model/project/`
package.

## Export and input reading

Every object that has a native FEMaster representation implements `export()`.
Repositories implement `export()` only where they own ordering or grouping.
`Project.export()` therefore only orchestrates global dependency order and
scope.

```python
from femaster_api import Node, Project, T3D2

project = Project("model")
part = project.parts.default()

part.nodes.add(Node(1, 0.0, 0.0, 0.0))
part.nodes.add(Node(2, 1.0, 0.0, 0.0))
part.elements.add(T3D2(1, (1, 2)))

project.write("model.inp")
```

Input reading belongs to the model root rather than a parallel importer object:

```python
project = Project.read_inp("model.inp")
project = Project.read_inp_text(text)
```

Unsupported keyword blocks are retained on `project.unparsed_blocks` instead of
being silently discarded.

## Results

Post-processing uses one format-independent hierarchy:

```text
Result -> Solution -> Frame -> Field
```

A `Solution` represents one solver/loadcase result. A `Frame` has a discrete
`id`, an optional physical `value` and a collection of fields. The meaning of
`value` depends on the analysis and file format: for example time, frequency,
buckling factor or another frame coordinate.

```python
from femaster_api import Result

result = Result.read("job.res")
result = Result.read_res("job.res")
result = Result.read_frd("job.frd")

solution = result.solution(1)
frame = solution.frame(1)
displacement = frame.field("DISPLACEMENT")
```

FEMaster FRD output stores the physical frame value in its `100CL` record, so it
is retained as `Frame.value`. The current native RES writer does not persist the
`frame_value` argument; in that case `Frame.value` is intentionally `None`
unless an explicit `FRAME` record is present. Semantic entity identifiers such
as `17` and `bolt.17` are preserved by result readers.
