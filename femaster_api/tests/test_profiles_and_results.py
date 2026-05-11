from femaster_api import (
    B33,
    BeamSection,
    Element,
    ElementSet,
    FieldDomain,
    Material,
    Model,
    Node,
    Profile,
    ResultReader,
)


def test_profiles_are_model_owned_outside_sections() -> None:
    model = Model("profiles")
    n1 = model.nodes.add(Node(0.0, 0.0, 0.0))
    n2 = model.nodes.add(Node(1.0, 0.0, 0.0))
    element = model.elements.add(Element(B33, (n1, n2)))
    element_set = model.sets.add(ElementSet("BEAM", (element,)))
    material = model.materials.add(Material("STEEL").set_isotropic_elasticity(E=1.0, nu=0.3))
    profile = model.profiles.add(Profile("BOX", area=1.0, iy=1.0, iz=1.0, j=1.0))

    model.sections.add(BeamSection("SEC", material, element_set, profile))

    assert model.validate().errors == ()


def test_result_reader_maps_loadcases_frames_and_element_nodal_fields() -> None:
    result = ResultReader().parse(
        """
        LOADCASE 1 STATIC
        FRAME 2 converged
        FIELD DISPLACEMENT ROWS=2 COLS=5 INDEX_COLS=2 VALUE_COLS=3
        1, 10, 0.1, 0.2, 0.3
        1, 11, 0.4, 0.5, 0.6
        END FIELD
        FIELD ENERGY ROWS=1 COLS=1
        42.0
        END FIELD
        """
    )

    frame = result.loadcase(1).frame(2)
    displacement = frame.field("DISPLACEMENT")
    energy = frame.field("ENERGY")

    assert displacement.domain is FieldDomain.ELEMENT_NODAL
    assert displacement.row((1, 10)) == (0.1, 0.2, 0.3)
    assert energy.domain is FieldDomain.UNKNOWN
    assert energy.row(0) == (42.0,)
