import pytest

from femaster_api import (
    C3D4,
    EntityType,
    FEMasterReader,
    FieldDomain,
    Model,
    NonlinearStaticStep,
    StaticStep,
    load_model_from_inp,
)
from femaster_api.export import FEMasterWriter


def test_importer_maps_deck_ids_to_repository_ids() -> None:
    model = load_model_from_inp(
        """
        *NODE
        10, 0, 0, 0
        20, 1, 0, 0
        30, 0, 1, 0
        40, 0, 0, 1
        *ELEMENT, TYPE=C3D4, ELSET=TETS
        5, 10, 20, 30, 40
        *FIELD, NAME=TEMP, TYPE=NODE, COLS=1
        20, 42.0
        *MATERIAL, NAME=MAT
        *ELASTIC, TYPE=ISOTROPIC
        1.0, 0.3
        *SOLIDSECTION, ELSET=TETS, MATERIAL=MAT
        """.splitlines()
    )

    assert len(model.nodes) == 4
    assert len(model.elements) == 1
    assert model.elements[0].topology is C3D4
    assert model.elements[0].node_ids == (0, 1, 2, 3)
    assert model.sets.get(EntityType.ELEMENT, "TETS").members[0].id == 0
    assert model.fields.get("TEMP").row(1) == (42.0,)


def test_importer_reads_orthotropic_engineering_constants() -> None:
    model = load_model_from_inp(
        """
        *MATERIAL, NAME=ORTHO
        *ELASTIC, TYPE=ENGINEERINGCONSTANTS
        100.0, 200.0, 300.0, 0.12, 0.13, 0.23, 12.0, 13.0, 23.0
        """.splitlines()
    )

    elasticity = model.materials.get("ORTHO").elasticity
    assert elasticity.e1 == pytest.approx(100.0)
    assert elasticity.e2 == pytest.approx(200.0)
    assert elasticity.e3 == pytest.approx(300.0)
    assert elasticity.nu12 == pytest.approx(0.12)
    assert elasticity.nu13 == pytest.approx(0.13)
    assert elasticity.nu23 == pytest.approx(0.23)
    assert elasticity.g12 == pytest.approx(12.0)
    assert elasticity.g13 == pytest.approx(13.0)
    assert elasticity.g23 == pytest.approx(23.0)


def test_importer_converts_orthotropic_stiffness_coefficients() -> None:
    model = load_model_from_inp(
        """
        *MATERIAL, NAME=ORTHO
        *ELASTIC, TYPE=ORTHOTROPIC
        120.0, 40.0, 120.0, 40.0, 40.0, 120.0, 40.0, 40.0, 40.0
        """.splitlines()
    )

    elasticity = model.materials.get("ORTHO").elasticity
    assert elasticity.e1 == pytest.approx(100.0)
    assert elasticity.e2 == pytest.approx(100.0)
    assert elasticity.e3 == pytest.approx(100.0)
    assert elasticity.nu12 == pytest.approx(0.25)
    assert elasticity.nu13 == pytest.approx(0.25)
    assert elasticity.nu23 == pytest.approx(0.25)
    assert elasticity.g12 == pytest.approx(40.0)
    assert elasticity.g13 == pytest.approx(40.0)
    assert elasticity.g23 == pytest.approx(40.0)


def test_importer_reads_static_and_nonlinear_loadcases() -> None:
    model = load_model_from_inp(
        """
        *NODE
        1, 0, 0, 0
        *NSET, NSET=FIXED
        1
        *SUPPORT, SUPPORT_COLLECTOR=BCS
        FIXED, 0, 0, 0, 0, 0, 0
        *LOADCASE, TYPE=LINEARSTATIC, NAME=static
        *SUPPORTS
        BCS
        *CONSTRAINTMETHOD, TYPE=NULLSPACE
        *LOADCASE, TYPE=NONLINEARSTATIC, NAME=nl
        *SUPPORTS
        BCS
        *NONLINEAR, CONTROL=ARC_LENGTH, MAX_INCREMENTS=12, ADAPTIVE=OFF, TOL=1e-7
        """.splitlines()
    )

    assert isinstance(model.steps[0], StaticStep)
    assert isinstance(model.steps[1], NonlinearStaticStep)
    assert model.steps[1].control == "ARC_LENGTH"
    assert model.steps[1].max_increments == 12
    assert model.steps[1].adaptive is False
    assert model.steps[1].tolerance == pytest.approx(1e-7)


def test_writer_output_roundtrips_through_importer(tmp_path) -> None:
    original = Model("roundtrip")
    n0 = original.nodes.add(__import__("femaster_api").Node(0.0, 0.0, 0.0))
    n1 = original.nodes.add(__import__("femaster_api").Node(1.0, 0.0, 0.0))
    n2 = original.nodes.add(__import__("femaster_api").Node(0.0, 1.0, 0.0))
    n3 = original.nodes.add(__import__("femaster_api").Node(0.0, 0.0, 1.0))
    original.elements.add(__import__("femaster_api").Element(C3D4, (n0, n1, n2, n3)))

    path = tmp_path / "roundtrip.inp"
    FEMasterWriter(original, include_header=False).write(path)

    imported = FEMasterReader().read(path)
    assert len(imported.nodes) == 4
    assert len(imported.elements) == 1
    assert imported.elements[0].node_ids == (0, 1, 2, 3)
