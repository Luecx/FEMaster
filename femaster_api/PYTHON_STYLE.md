# FEMaster Python Code Style

These rules apply to the FEMaster Python API. They intentionally mirror the C++
style: correctness and explicit FEM semantics take precedence over compactness,
framework abstractions, or formatter-driven code shape.

## 1. Priorities

When several implementations are possible, prefer in this order:

1. Correct FEMaster semantics.
2. Clear ownership and scope.
3. Direct, locally understandable Python.
4. Readable top-to-bottom control flow.
5. Stable semantic names and FEMaster identifiers.
6. Reuse of existing public model classes.
7. Small abstractions only where they remove real repeated behavior.

Do not introduce a framework merely to reduce line count.

## 2. Model architecture

- `Project` is the top-level owner of global definitions and repositories.
- `Part` owns part-local nodes, elements, regions, surfaces and sections.
- The default Part exists only as `project.parts[0]`.
- `project.parts.default()` must return exactly `project.parts[0]`.
- The default Part cannot be deleted or replaced.
- Named repositories support `repository[index]` and `repository[name]`.
- Repository indices are positional conveniences, never persistent model IDs.
- Persistent references between named objects use names.
- Node and element repositories preserve FEMaster IDs; `nodes[id]` and
  `elements[id]` address those semantic IDs.
- Loads and supports belong directly to their collectors. Do not maintain a
  second global copy of collector entries.

## 3. Export

Every class that has an independent FEMaster representation implements
`to_femaster()` itself.

Repositories implement `to_femaster()` when they own ordering or grouping
required by the input format. `Project.to_femaster()` only orchestrates global
output order; it must not become a type-switching serializer registry.

Use the small `_format.py` helpers for keyword lines, CSV rows and block joining.
Do not introduce visitors, decorators, serializer registries or metaclasses for
ordinary FEMaster output.

## 4. Readers

Readers populate the same public object model that users construct manually.
Do not create a second DTO hierarchy for parsed files.

Syntax that is parsed successfully but has no semantic implementation must be
retained explicitly as unsupported/unparsed data rather than silently ignored.

Result readers map format-specific data into the common hierarchy:

`Result -> LoadCase -> Frame -> Field`

`FieldDomain` and `FieldType` are the central definitions for field semantics.
Format readers must not maintain independent competing field enums.

## 5. Files and modules

Each module starts with a descriptive module docstring explaining:

- what the module defines,
- where the objects live in FEMaster scope,
- ownership responsibilities,
- important identifier or export conventions.

Keep related small value classes together when they form one subsystem. Split a
module when it becomes difficult to understand as one semantic unit, not merely
because it exceeds an arbitrary line count.

## 6. Classes

Every public class receives a docstring describing responsibility and scope.
Class contents should normally follow this order:

1. Class constants.
2. Construction and persistent definition data.
3. Public modification/access operations.
4. FEMaster export.
5. Container protocol methods.
6. Private helpers.

Avoid properties that only hide a trivial public attribute. Properties are
appropriate for invariants such as immutable semantic names.

## 7. Functions

Non-trivial functions should read as a sequence of documented phases. Use
section comments for substantial phases:

```python
# ------------------------------------------------------------------
# Part-local topology
# ------------------------------------------------------------------
```

Within algorithms, comments should explain FEM semantics, state transitions,
identifier spaces, parsing decisions or numerical meaning. Do not comment every
obvious assignment individually.

Long functions are acceptable when they remain a clear ordered semantic pass.
Do not extract one-use helpers solely to shorten a function.

## 8. Data classes and dependencies

Prefer ordinary Python classes and small standard-library types. `dataclass` is
acceptable for passive records, but it is not required for every model object.

Do not add Pydantic, attrs, validation frameworks, serializer frameworks or
other third-party dependencies for functionality that is straightforward with
the standard library.

## 9. Formatting

- Use four spaces.
- Use type annotations on public APIs.
- Keep short expressions on one line.
- Align closely related assignments when it improves visual structure.
- Prefer explicit names over abbreviations except established FEM terminology.
- Avoid formatter rules that destroy useful mathematical or structural layout.

## 10. Validation

Validate invariants at the object boundary where invalid data first becomes
meaningful, for example:

- duplicate repository names or IDs,
- element connectivity size,
- vector dimension,
- immutable names,
- required keyword keys,
- field row width,
- protected default Part removal.

Do not silently renumber user node or element IDs.

## 11. Scope of changes

Do not retain compatibility abstractions for an obsolete API unless explicitly
required. When an API is intentionally redesigned, keep one clear implementation
rather than a new model plus adapters for every previous representation.

Do not claim tests have passed unless they were actually executed.
