# FEMaster API

This package provides a small, object-based Python modeling interface for
FEMaster input decks and result files.

The in-memory model is independent from FEMaster deck syntax. Model objects do
not store export IDs; `FEMasterWriter` assigns local IDs only while writing.

Implemented scope:

- object repositories for nodes, elements, sets, surfaces, orientations,
  materials, sections, fields, features, loads, supports, constraints, and
  loadcases
- FEMaster deck writer for the currently registered DSL keywords
- validation diagnostics for object ownership and cross-repository references
- subprocess runner for FEMaster
- text `.res` result reader with `Result -> LoadcaseResult -> Frame -> fields`
- optional VTK XML and VTKHDF mesh export helpers under
  `femaster_api.tools.tovtk`

Example:

```python
from femaster_api import FEMasterWriter, Model
from femaster_api.model import (
    Element,
    ElementSet,
    IsotropicElasticity,
    LinearStaticLoadcase,
    LoadCollector,
    Material,
    NodalForce,
    Node,
    NodeSet,
    ShellSection,
    Support,
    SupportCollector,
)
from femaster_api.model.elements import S4

model = Model("cantilever")

n1 = model.nodes.add(Node(0.0, 0.0, 0.0))
n2 = model.nodes.add(Node(1.0, 0.0, 0.0))
n3 = model.nodes.add(Node(1.0, 1.0, 0.0))
n4 = model.nodes.add(Node(0.0, 1.0, 0.0))

plate = model.elements.add(Element(S4, [n1, n2, n3, n4]))

fixed        = model.node_sets.add(NodeSet("FIXED", [n1, n4]))
loaded       = model.node_sets.add(NodeSet("LOADED", [n2, n3]))
plate_set    = model.element_sets.add(ElementSet("PLATE_ELEMENTS", [plate]))

steel = model.materials.add(
    Material(
        "STEEL",
        elasticity=IsotropicElasticity(
            youngs_modulus=210000.0,
            poisson_ratio=0.3,
        ),
    )
)

model.sections.add(
    ShellSection(
        "PLATE",
        material    = steel,
        element_set = plate_set,
        thickness   = 1.0,
    )
)

bcs = model.support_collectors.add(
    SupportCollector(
        "BCS",
        [
            Support(fixed, ux=0.0, uy=0.0, uz=0.0),
        ],
    )
)

loads = model.load_collectors.add(
    LoadCollector(
        "LOADS",
        [
            NodalForce(loaded, fz=-1000.0),
        ],
    )
)

model.loadcases.add(
    LinearStaticLoadcase(
        "LC1",
        supports = [bcs],
        loads    = [loads],
    )
)

model.validate().raise_for_errors()
FEMasterWriter(model).write("cantilever.inp")
```

Run FEMaster and read results directly:

```python
from femaster_api import FEMaster

solver = FEMaster(
    model,
    threads         = 8,
    executable      = None,   # auto-detect from FEMASTER_EXECUTABLE or PATH
    temporary_files = True,   # keep .inp/.res in a temporary directory
    redirect_stdout = True,   # capture solver stdout/stderr
)

result = solver.run()

frame = result.loadcases["LC1"].frames[0]
disp  = frame.fields["DISPLACEMENT"]
```

Use explicit files when you want to keep the generated input and result:

```python
solver = FEMaster(
    model,
    input_path      = "cantilever.inp",
    result_path     = "cantilever.res",
    temporary_files = False,
)

result = solver.run()
```

VTK mesh/result export:

```python
from femaster_api.tools.tovtk import write_vtk, write_vtkhdf

write_vtk(model, "cantilever.vtk", result=result)
write_vtkhdf(model, "cantilever.vtkhdf")
```
