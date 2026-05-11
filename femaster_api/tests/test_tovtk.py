import pytest
import vtk

from femaster_api import C3D4, Element, Field, FieldDomain, Frame, LoadCase, Model, Node, Result, tovtk


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
