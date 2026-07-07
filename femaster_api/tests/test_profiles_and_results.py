from femaster_api import (
    B33,
    BeamSection,
    Element,
    ElementSet,
    FEMasterWriter,
    FieldDomain,
    Material,
    ModalStep,
    Model,
    Node,
    BucklingStep,
    NonlinearStaticStep,
    Profile,
    ResultReader,
    TimeControl,
    TransientStep,
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
        FIELD DISPLACEMENT TYPE=ELEMENT_NODAL ROWS=2 INDEX_COLS=2 VALUE_COLS=3
        1, 10, 0.1, 0.2, 0.3
        1, 11, 0.4, 0.5, 0.6
        END FIELD
        FIELD ENERGY TYPE=ELEMENT ROWS=1 COLS=1
        42.0
        END FIELD
        FIELD IP_STRESS TYPE=ELEMENT_IP ROWS=1 COLS=6
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0
        END FIELD
        """
    )

    frame = result.loadcase(1).frame(2)
    displacement = frame.field("DISPLACEMENT")
    energy = frame.field("ENERGY")
    ip_stress = frame.field("IP_STRESS")

    assert displacement.domain is FieldDomain.ELEMENT_NODAL
    assert displacement.row((1, 10)) == (0.1, 0.2, 0.3)
    assert energy.domain is FieldDomain.ELEMENT
    assert energy.row(0) == (42.0,)
    assert ip_stress.domain is FieldDomain.ELEMENT_IP
    assert ip_stress.row(0) == (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)


def test_nonlinear_static_step_writes_current_femaster_keywords() -> None:
    model = Model("nonlinear")
    model.steps.add(
        NonlinearStaticStep(
            "arc",
            control="ARC_LENGTH",
            initial_increment=0.05,
            minimum_increment=1e-6,
            maximum_increment=0.1,
            max_increments=300,
            adaptive=False,
            max_iterations=35,
            tolerance=1e-7,
            regularize_zero_rows=False,
            constraint_summary=True,
        )
    )

    deck = FEMasterWriter(model, include_header=False).render()

    assert "*LOADCASE, TYPE=NONLINEARSTATIC, NAME=arc" in deck
    assert "*NONLINEAR, CONTROL=ARC_LENGTH" in deck
    assert "INITIAL_INCREMENT=0.05" in deck
    assert "MINIMUM_INCREMENT=1e-06" in deck
    assert "MAX_INCREMENTS=300" in deck
    assert "ADAPTIVE=OFF" in deck
    assert "MAXITER=35" in deck
    assert "TOL=1e-07" in deck
    assert "REGULARIZE_ZERO_ROWS=OFF" in deck
    assert "*CONSTRAINTSUMMARY" in deck


def test_writer_does_not_emit_unsupported_constraintmethod_for_modal_buckling_or_transient() -> None:
    model = Model("constraintmethod")
    model.steps.add(ModalStep("modes"))
    model.steps.add(BucklingStep("buckling", loads=()))
    model.steps.add(TransientStep("transient", loads=(), time=TimeControl(0.0, 1.0, 0.1)))

    deck = FEMasterWriter(model, include_header=False).render()

    assert "*CONSTRAINTMETHOD" not in deck
    assert "*LOADCASE, TYPE=EIGENFREQ, NAME=modes" in deck
    assert "*LOADCASE, TYPE=LINEARBUCKLING, NAME=buckling" in deck
    assert "*LOADCASE, TYPE=LINEARTRANSIENT, NAME=transient" in deck
