from femaster_api.tools.vtk_export.cli import convert_one, gather_inp_files


def test_cli_gathers_and_converts_inp(tmp_path) -> None:
    case = tmp_path / "case.inp"
    case.write_text(
        """
        *NODE
        1, 0, 0, 0
        2, 1, 0, 0
        3, 0, 1, 0
        4, 0, 0, 1
        *ELEMENT, TYPE=C3D4
        1, 1, 2, 3, 4
        """,
        encoding="utf-8",
    )

    assert gather_inp_files([str(tmp_path)]) == [case]
    ok, message = convert_one(
        case,
        writer="vtk",
        solution_ext=".res",
        export_sets=False,
        export_couplings=False,
        export_shell_thickness=False,
        pack_sets=False,
        log_mode="quiet",
        log_format="table",
    )

    assert ok, message
    assert case.with_suffix(".vtk").exists()
