"""Structural regression tests for the deliberately explicit package layout."""

from __future__ import annotations

import ast
from pathlib import Path

import femaster_api


PACKAGE_ROOT = Path(femaster_api.__file__).parent
MODEL_ROOT = PACKAGE_ROOT / "model"


def test_every_python_file_defines_at_most_one_class():
    """Prevent future drift back to large multi-class aggregation modules."""

    offenders: list[str] = []

    for path in PACKAGE_ROOT.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        classes = [
            node
            for node in tree.body
            if isinstance(node, ast.ClassDef)
        ]
        if len(classes) > 1:
            offenders.append(
                f"{path.relative_to(PACKAGE_ROOT)}: "
                + ", ".join(node.name for node in classes)
            )

    assert not offenders, "\n".join(offenders)


def test_model_files_have_descriptive_module_docstrings():
    """Require real file-level documentation, not empty aggregation modules."""

    offenders: list[str] = []

    for path in PACKAGE_ROOT.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        docstring = ast.get_docstring(tree) or ""
        if len(docstring.strip()) < 120:
            offenders.append(str(path.relative_to(PACKAGE_ROOT)))

    assert not offenders, "\n".join(offenders)


def test_domain_layout_has_no_io_mesh_or_project_subpackage():
    """Keep I/O behavior on Project/Result and topology in domain folders."""

    assert not (PACKAGE_ROOT / "io").exists()
    assert not (MODEL_ROOT / "mesh").exists()
    assert not (MODEL_ROOT / "project").exists()

    assert (MODEL_ROOT / "project.py").is_file()
    assert (MODEL_ROOT / "node").is_dir()
    assert (MODEL_ROOT / "element").is_dir()
    assert (MODEL_ROOT / "surface").is_dir()
    assert (MODEL_ROOT / "region").is_dir()
    assert (MODEL_ROOT / "result").is_dir()


def test_concrete_element_set_matches_canonical_femaster_types():
    """Expose only the canonical internal element names requested by the API."""

    expected = {
        "B33",
        "C3D4",
        "C3D5",
        "C3D6",
        "C3D8",
        "C3D8R",
        "C3D10",
        "C3D15",
        "C3D20",
        "C3D20R",
        "S3",
        "S4",
        "S6",
        "S8",
        "T3",
    }

    assert set(femaster_api.ELEMENT_TYPES) == expected

    element_root = MODEL_ROOT / "element"
    forbidden_modules = {
        "element_t3d2.py",
        "element_mitc4.py",
        "element_mitc8.py",
        "element_qspt.py",
        "element_mitc3frt.py",
        "element_mitc4frt.py",
        "element_mitc6frt.py",
        "element_mitc8frt.py",
        "element_mass.py",
        "element_rotary_inertia.py",
        "element_spring.py",
    }

    assert not any((element_root / name).exists() for name in forbidden_modules)
    assert not hasattr(femaster_api, "T3D2")
    assert not hasattr(femaster_api, "MITC4")
    assert not hasattr(femaster_api, "QSPT")
