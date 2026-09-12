# FEMaster Python Code Style

These rules mirror the C++ code philosophy: correctness and explicit FEM
semantics take precedence over compactness or framework abstractions.

## 1. Package structure

- Every model object lives below `femaster_api/model/`.
- Model concepts are grouped into semantic subpackages (`part`, `mesh`, `field`,
  `result`, `section`, `load`, ...).
- `io` contains only parsing/import infrastructure and syntax-level records.
- Do not create model classes at package root.
- Do not reintroduce flat collection modules such as `mesh.py`, `fields.py`,
  `materials.py` or `steps.py`.

## 2. One class per file

**Every Python source file contains at most one class definition.**

This applies equally to:

- public model classes,
- repositories,
- enums,
- control/value classes,
- abstract/base classes,
- parser/importer classes,
- private helper classes.

A module may contain helper functions and type aliases when no class is defined
or when those helpers are local to the single class in that file.

`__init__.py` files only re-export names and must not define classes.

## 3. Model architecture

- `Project` is the top-level model owner.
- `Part` owns part-local nodes, elements, regions, surfaces and sections.
- The default Part exists only as `project.parts[0]`.
- `project.parts.default()` returns exactly `project.parts[0]`.
- The default Part cannot be deleted or replaced.
- Named repositories support `repository[index]` and `repository[name]`.
- Repository positions are convenience indices, never persistent model IDs.
- Persistent references between named objects use names.
- Node and element repositories preserve FEMaster IDs.
- Loads and supports belong directly to their collectors.
- Field and result objects are model objects and live under `model/field` and
  `model/result`.

## 4. Export

Every model class with an independent native representation implements
`export()` itself.

Repositories implement `export()` only where they own ordering or grouping
required by the input syntax. `Project.export()` orchestrates dependency order;
it must not become a type-switching serializer registry.

Do not introduce visitors, serializer registries, decorators or metaclasses for
ordinary native export.

## 5. Import

Importers populate the same public model classes users construct manually.
Do not maintain a second DTO hierarchy.

Use:

- `InpImporter.import_file()` / `import_text()`
- `ResImporter.import_file()` / `import_text()`
- `FrdImporter.import_file()` / `import_text()`

Successfully parsed but unsupported input blocks are retained explicitly rather
than silently discarded.

## 6. Fields and results

`FieldDomain` and `FieldType` are the central format-independent definitions.
RES and FRD importers map into:

```text
Result -> LoadCase -> Frame -> Field
```

Format-specific readers must not create competing field classes or enums.

## 7. Code organization

Prefer direct, locally understandable Python. Keep mathematical and semantic
ownership explicit. Use small abstractions only when they represent a genuine
concept or meaningful repeated behavior.

Validate invariants at the object boundary where invalid data first becomes
meaningful. Never silently renumber node or element IDs.

Do not claim tests passed unless they were actually executed.
