import ast
from pathlib import Path

import femaster_api


def test_every_source_file_contains_at_most_one_class():
    package_root = Path(femaster_api.__file__).parent

    for path in package_root.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        classes = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.ClassDef)
        ]
        assert len(classes) <= 1, (
            f"{path.relative_to(package_root)} defines "
            f"{len(classes)} classes"
        )


def test_model_classes_are_not_defined_at_package_root():
    package_root = Path(femaster_api.__file__).parent
    root_modules = [
        path
        for path in package_root.glob("*.py")
        if path.name != "__init__.py"
    ]
    assert root_modules == []
