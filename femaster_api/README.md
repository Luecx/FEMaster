# FEMaster Python API

The Python API mirrors the semantic FEMaster model and exports native FEMaster
keyword input directly from the model objects.

```python
from femaster_api import Node, Project, T3D2

project = Project("model")
part = project.parts.default()

part.nodes.add(Node(1, 0.0, 0.0, 0.0))
part.nodes.add(Node(2, 1.0, 0.0, 0.0))
part.elements.add(T3D2(1, (1, 2)))

text = project.export()
project.write("model.inp")
```

## Repository access

Named repositories support both positional and semantic-name access:

```python
project.parts[0]
project.parts["WING"]
project.parts.default()
```

The default Part is permanently stored at position zero. It is not duplicated on
`Project` and cannot be deleted.

Nodes and elements are different: their integer keys are FEMaster IDs rather
than repository positions.

```python
part.nodes[100]
part.elements[250]
```

Use `repository.at(index)` only when insertion-order positional access is really
required for ID-based repositories.

## Importing input and results

```python
from femaster_api import InpImporter, import_result

project = InpImporter().import_file("model.inp")
result = import_result("model.res")
field = result.field("DISP", loadcase=1, frame=0)
```

For in-memory text, importers also provide `import_text()`.

RES and FEMaster-generated nodal FRD output use the common
`Result -> LoadCase -> Frame -> Field` hierarchy and the central `FieldDomain`
and `FieldType` enums.

See `PYTHON_STYLE.md` for the architectural and documentation conventions used
by the package.
