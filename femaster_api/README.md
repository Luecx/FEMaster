# FEMaster API

This package is a new, independent Python interface for FEMaster.

The design keeps the in-memory model separate from the FEMaster input deck.
Model objects are plain data objects and repositories; the FEMaster-specific
text format is generated only by `FEMasterWriter`.

Implemented scope:

- neutral model repositories for nodes, elements, sets, surfaces, orientations,
  materials, sections, fields, features, loads, supports, constraints, and
  loadcases
- FEMaster deck writer for the currently registered DSL keywords
- model validation diagnostics
- subprocess runner for FEMaster
- text `.res` result reader for both documented `LOADCASE` blocks and the
  current writer format using `LC` plus `FIELD, NAME=...`

Example:

```python
from femaster_api import S4, FEMasterWriter, Model

model = Model("cantilever")
n1 = model.add_node(1, 0, 0)
n2 = model.add_node(2, 0, 0)
n3 = model.add_node(2, 1, 0)
n4 = model.add_node(1, 1, 0)

plate = model.add_element(S4, n1, n2, n3, n4)
fixed = model.sets.nodes("FIXED", [n1, n4])
loaded = model.sets.nodes("LOADED", [n2, n3])
model.sets.elements("PLATE_ELEMENTS", [plate])

model.materials.isotropic("STEEL", youngs_modulus=210000.0, poisson_ratio=0.3)
model.sections.shell("PLATE", material="STEEL", element_set="PLATE_ELEMENTS", thickness=1.0)
model.boundaries.fix(fixed, ux=0, uy=0, uz=0, name="FIXED_SUPP")
model.loads.force(loaded, fz=-1000, name="FORCE_Z")
model.steps.static("LOADCASE_1", loads=["FORCE_Z"], boundaries=["FIXED_SUPP"])

model.validate(raise_on_error=True)
FEMasterWriter(model).write("cantilever.inp")
```
