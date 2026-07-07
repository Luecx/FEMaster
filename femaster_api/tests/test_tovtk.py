import pytest
import vtk

from femaster_api import C3D4, CouplingConstraint, Element, ElementSet, Field, FieldDomain, Frame, LoadCase, Material, Model, Node, NodeSet, Result, ShellSection, tovtk


def test_tovtk_writes_model_fields_without_result(tmp_path) -> None:
    model = _tet_model()
    temperature = Field("TEMPERATURE", FieldDomain.NODE, 1).set(0, (20.0,)).set(1, (21.0,))
    model.fields.add(temperature)

    path = tmp_path / "model.vtu"
    written = tovtk(model, path)

    grid = _read_vtu(written)
    assert grid.GetNumberOfPoints() == 4
    assert grid.GetNumberOfCells() == 1
    assert grid.GetPointData().GetArray("TEMPERATURE") is not None


def test_tovtk_writes_result_fields_and_derived_fields(tmp_path) -> None:
    model = _tet_model()
    result = Result(
        [
            LoadCase(
                1,
                frames=[
                    Frame(
                        0,
                        fields={
                            "DISPLACEMENT": Field("DISPLACEMENT", FieldDomain.NODE, 3)
                            .set(0, (0.1, 0.0, 0.0))
                            .set(1, (0.0, 0.2, 0.0)),
                            "STRESS": Field("STRESS", FieldDomain.ELEMENT, 6).set(0, (100.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
                            "BUCKLING_MODE": Field("BUCKLING_MODE", FieldDomain.ELEMENT_NODAL, 3)
                            .set((0, 0), (1.0, 0.0, 0.0))
                            .set((0, 1), (0.0, 1.0, 0.0))
                            .set((0, 2), (0.0, 0.0, 1.0)),
                        },
                    )
                ],
            )
        ]
    )

    path = tmp_path / "result.vtu"
    tovtk(model, result, path)

    grid = _read_vtu(path)
    point_data = grid.GetPointData()
    cell_data = grid.GetCellData()
    assert point_data.GetArray("DISPLACEMENT") is not None
    assert point_data.GetArray("DISPLACEMENT_XYZ") is not None
    assert point_data.GetArray("BUCKLING_MODE") is not None
    assert point_data.GetArray("BUCKLING_MODE_XYZ") is not None
    assert cell_data.GetArray("BUCKLING_MODE_ELEMENT_NODAL") is not None
    assert cell_data.GetArray("STRESS") is not None
    assert cell_data.GetArray("MISES").GetTuple1(0) == pytest.approx(100.0)
    assert cell_data.GetArray("STRESS_MISES").GetTuple1(0) == pytest.approx(100.0)


def test_tovtk_exports_sets_and_shell_thickness(tmp_path) -> None:
    model = _tet_model()
    node_set = model.sets.add(NodeSet("SUPPORT", (model.nodes[0], model.nodes[1])))
    element_set = model.sets.add(ElementSet("SHELLS", (model.elements[0],)))
    material = model.materials.add(Material("MAT").set_isotropic_elasticity(E=1.0, nu=0.3))
    model.sections.add(ShellSection("SEC", material, element_set, 0.25))

    path = tmp_path / "sets.vtu"
    tovtk(model, path, export_sets=True, export_shell_thickness=True)

    grid = _read_vtu(path)
    assert grid.GetPointData().GetArray(f"SET_N_{node_set.name}") is not None
    assert grid.GetCellData().GetArray(f"SET_E_{element_set.name}") is not None
    assert grid.GetCellData().GetArray("SHELL_THICKNESS").GetTuple1(0) == pytest.approx(0.25)


def test_tovtk_writes_vtkhdf_when_vtk_supports_it(tmp_path) -> None:
    if not hasattr(vtk, "vtkHDFWriter"):
        pytest.skip("VTK build has no vtkHDFWriter")

    path = tmp_path / "model.vtkhdf"
    tovtk(_tet_model(), path)

    assert path.exists()
    assert path.stat().st_size > 0


def test_tovtk_writes_vtkhdf_time_steps(tmp_path) -> None:
    if not hasattr(vtk, "vtkHDFWriter") or not hasattr(vtk, "vtkHDFReader"):
        pytest.skip("VTK build has no VTK-HDF reader/writer")

    result = Result(
        [
            LoadCase(
                1,
                frames=[
                    Frame(0, fields={"TIME": Field("TIME", FieldDomain.UNKNOWN, 1).set(0, (0.0,))}),
                    Frame(1, fields={"TIME": Field("TIME", FieldDomain.UNKNOWN, 1).set(0, (0.5,))}),
                ],
            )
        ]
    )
    path = tmp_path / "series.vtkhdf"
    tovtk(_tet_model(), result, path)

    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(path))
    reader.UpdateInformation()
    info = reader.GetOutputInformation(0)
    key = vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS()
    assert info.Length(key) == 2


def test_tovtk_exports_kinematic_couplings(tmp_path) -> None:
    model = _tet_model()
    master = model.sets.add(NodeSet("MASTER", (model.nodes[0],)))
    slave = model.sets.add(NodeSet("SLAVE", (model.nodes[1], model.nodes[2])))
    model.constraints.add(CouplingConstraint(master, slave, (True, True, False, False, False, False)))

    path = tmp_path / "couplings.vtu"
    tovtk(model, path, export_couplings=True)

    grid = _read_vtu(path)
    assert grid.GetNumberOfCells() == 3
    assert grid.GetCellData().GetArray("COUPLING_ID") is not None
    assert grid.GetCellData().GetArray("COUPLING_DOF_MASK") is not None


def _tet_model() -> Model:
    model = Model("tet")
    n0 = model.nodes.add(Node(0.0, 0.0, 0.0))
    n1 = model.nodes.add(Node(1.0, 0.0, 0.0))
    n2 = model.nodes.add(Node(0.0, 1.0, 0.0))
    n3 = model.nodes.add(Node(0.0, 0.0, 1.0))
    model.elements.add(Element(C3D4, (n0, n1, n2, n3)))
    return model


def _read_vtu(path):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()
