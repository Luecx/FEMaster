from femaster_api import InpReader


DECK = """\
*MODEL, NAME=READ_TEST
*NODE, NSET=NALL
10, 0.0, 0.0, 0.0
20, 1.0, 0.0, 0.0
*ELEMENT, TYPE=T3D2, ELSET=BAR
50, 10, 20
*NSET, NAME=ROOT
10
*MATERIAL, NAME=STEEL
*ELASTIC, TYPE=ISOTROPIC
210000.0, 0.3
*TRUSSSECTION, ELSET=BAR, MATERIAL=STEEL
100.0
*SUPPORT, SUPPORT_COLLECTOR=BC
ROOT, 0.0, 0.0, 0.0
*LOADCASE, TYPE=EIGENFREQ, NAME=MODES
*SUPPORTS
BC
*NUMEIGENVALUES
3
*END
"""


def test_inp_reader_builds_public_model():
    project = InpReader().parse(DECK)
    part = project.parts.default()

    assert project.name == "READ_TEST"
    assert part.nodes[10].x == 0.0
    assert part.nodes[20].x == 1.0
    assert part.elements[50].nodes == (10, 20)
    assert part.regions.nodes["ROOT"].members == [10]
    assert project.materials["STEEL"].elasticity.youngs_modulus == 210000.0
    assert project.support_collectors["BC"].supports[0].target == "ROOT"
    assert project.steps["MODES"].number_of_modes == 3
