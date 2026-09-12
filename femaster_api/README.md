# FEMaster Python API

The Python API mirrors FEMaster's semantic model. All model objects live below
`femaster_api.model`; importers are the only public object layer outside the
model package.

## Structure

The package is intentionally split by semantic concept:

```text
femaster_api/
    model/
        project/
        part/
        instance/
        mesh/
            elements/
        region/
        material/
        profile/
        coordinate_system/
        section/
        amplitude/
        load/
        support/
        constraint/
        feature/
        step/
        field/
        result/
        common/
    io/
```

Every Python source file contains **at most one class definition**. Repositories,
enums, controls and concrete element types therefore each have their own file.

`model/__init__.py` and the package root re-export the public classes, so normal
use remains compact:

```python
from femaster_api import Node, Project, T3D2

project = Project("model")
part = project.parts.default()

part.nodes.add(Node(1, 0.0, 0.0, 0.0))
part.nodes.add(Node(2, 1.0, 0.0, 0.0))
part.elements.add(T3D2(1, (1, 2)))

project.write("model.inp")
```

## Repository access

Named repositories support both positional and semantic-name access:

```python
project.parts[0]
project.parts["WING"]
project.parts.default()
```

The default Part is permanently stored at position zero. It is not duplicated
on `Project` and cannot be deleted.

Node and element repositories preserve FEMaster IDs. Their integer subscription
is semantic ID lookup rather than positional lookup:

```python
part.nodes[100]
part.elements[250]
```

Use `repository.at(index)` only for explicit insertion-order access.

## Export and import

Every model object with an independent native representation implements
`export()`. Repositories implement `export()` only when they own grouping or
ordering.

```python
text = project.export()
```

Input and result parsing lives under `femaster_api.io`:

```python
from femaster_api import InpImporter, ResImporter, FrdImporter

project = InpImporter().import_file("model.inp")
result = ResImporter().import_file("model.res")
```

`Field`, `FieldDomain`, `FieldType`, `Frame`, `LoadCase` and `Result` are model
objects and therefore live under `femaster_api.model.field` and
`femaster_api.model.result`, not beside the model package.
