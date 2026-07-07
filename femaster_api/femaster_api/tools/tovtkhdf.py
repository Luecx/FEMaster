"""Compatibility entry point for VTK-HDF export."""

from __future__ import annotations

from femaster_api.tools.vtk_export.cli import main as _main


def main(argv: list[str] | None = None) -> None:
    args = list(argv or [])
    if "--writer" not in args:
        args = ["--writer", "vtkhdf", *args]
    _main(args)


if __name__ == "__main__":
    main()
