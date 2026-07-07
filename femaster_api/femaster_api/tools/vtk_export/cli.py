"""Command-line converter for FEMaster decks and results."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable

from femaster_api.backend import ResultReader
from femaster_api.importers import FEMasterReader

from .writer_vtk import write_vtk
from .writer_vtkhdf import write_vtkhdf


def gather_inp_files(inputs: Iterable[str], *, recursive: bool = False) -> list[Path]:
    files: list[Path] = []
    for item in inputs:
        path = Path(item)
        if path.is_dir():
            files.extend(sorted(path.rglob("*.inp") if recursive else path.glob("*.inp")))
        elif path.is_file() and path.suffix.lower() == ".inp":
            files.append(path)
    return list(dict.fromkeys(files))


def convert_one(
    inp_path: Path,
    *,
    writer: str,
    solution_ext: str,
    export_sets: bool,
    export_couplings: bool,
    export_shell_thickness: bool,
    pack_sets: bool,
    log_mode: str,
    log_format: str,
) -> tuple[bool, str]:
    try:
        model = FEMasterReader().read(inp_path)
        result_path = inp_path.with_suffix(solution_ext)
        result = ResultReader().read(result_path) if result_path.exists() else None
        outputs: list[str] = []
        for writer_name in (("vtk", "vtkhdf") if writer == "both" else (writer,)):
            if writer_name == "vtk":
                out_path = inp_path.with_suffix(".vtk")
                write_vtk(
                    model,
                    result,
                    out_path,
                    export_sets=export_sets,
                    pack_sets=pack_sets,
                    export_shell_thickness=export_shell_thickness,
                    export_couplings=export_couplings,
                    log_mode=log_mode,
                    log_format=log_format,
                )
            elif writer_name == "vtkhdf":
                out_path = inp_path.with_suffix(".vtkhdf")
                write_vtkhdf(
                    model,
                    result,
                    out_path,
                    export_sets=export_sets,
                    pack_sets=pack_sets,
                    export_shell_thickness=export_shell_thickness,
                    export_couplings=export_couplings,
                    log_mode=log_mode,
                    log_format=log_format,
                )
            else:
                raise ValueError(f"unsupported writer {writer_name!r}")
            outputs.append(str(out_path))
        return True, ", ".join(outputs)
    except Exception as exc:
        return False, f"{inp_path}: {exc}"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert FEMaster .inp decks and optional .res files to VTK outputs.")
    parser.add_argument("inputs", nargs="+", help="One or more .inp files or directories.")
    parser.add_argument("--writer", choices=["vtk", "vtkhdf", "both"], default="vtk")
    parser.add_argument("--solution-ext", "--solution_ext", dest="solution_ext", default=".res")
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--export-sets", action="store_true")
    parser.add_argument("--pack-sets", action="store_true")
    parser.add_argument("--export-couplings", action="store_true")
    parser.add_argument("--export-shell-thickness", action="store_true")
    parser.add_argument("--log", dest="log_mode", choices=["quiet", "summary", "verbose"], default="summary")
    parser.add_argument("--format", dest="log_format", choices=["table", "json"], default="table")
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    files = gather_inp_files(args.inputs, recursive=args.recursive)
    if not files:
        print("[WARNING] No .inp files found.")
        return

    kwargs = dict(
        writer=args.writer,
        solution_ext=args.solution_ext,
        export_sets=args.export_sets,
        export_couplings=args.export_couplings,
        export_shell_thickness=args.export_shell_thickness,
        pack_sets=args.pack_sets,
        log_mode=args.log_mode,
        log_format=args.log_format,
    )
    ok = 0
    failed = 0
    print(f"[INFO] Found {len(files)} file(s).")
    if args.workers <= 1:
        for path in files:
            success, message = convert_one(path, **kwargs)
            ok += int(success)
            failed += int(not success)
            print("[INFO] " + ("OK  -> " if success else "ERR -> ") + message)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(convert_one, path, **kwargs) for path in files]
            for future in as_completed(futures):
                success, message = future.result()
                ok += int(success)
                failed += int(not success)
                print("[INFO] " + ("OK  -> " if success else "ERR -> ") + message)
    print(f"\nDone. ok={ok}, failed={failed}.")


if __name__ == "__main__":
    main()
