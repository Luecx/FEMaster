"""Checks for the public Result -> Solution -> Frame -> Field hierarchy."""

from femaster_api import FieldDomain, FieldType, Result


RES = """\
LC 1 STATIC
FRAME, ID=3, VALUE=0.25
FIELD, NAME=DISPLACEMENT, TYPE=NODE, INDEX_COLS=1, VALUE_COLS=3, ROWS=2
17 1.0 2.0 3.0
bolt.17 4.0 5.0 6.0
END FIELD
FIELD, NAME=STRESS, TYPE=ELEMENT_IP, INDEX_COLS=2, VALUE_COLS=2, ROWS=1
wing.4 2 100.0 200.0
END FIELD
"""


FRD = """\
    1PSTEP                         1           2           7
  100CL  101  5.00000E+00           1                     2    1MODAL      1
 -4  DISP        3    1
 -5  D1          1    1    1    0
 -5  D2          1    1    2    0
 -5  D3          1    1    3    0
 -1         17 1.0 2.0 3.0
 -3
"""


def test_read_res_builds_solution_frame_and_semantic_fields():
    result = Result.read_res_text(RES)
    solution = result.solution(1)
    frame = solution.frame(3)

    assert solution.name == "STATIC"
    assert frame.value == 0.25

    displacement = frame.field("DISPLACEMENT")
    assert displacement.domain is FieldDomain.NODE
    assert displacement.type is FieldType.DISPLACEMENT
    assert displacement[17] == (1.0, 2.0, 3.0)
    assert displacement["bolt.17"] == (4.0, 5.0, 6.0)

    stress = frame.field(FieldType.STRESS)
    assert stress[("wing.4", 2)] == (100.0, 200.0)


def test_read_frd_preserves_step_frame_and_physical_frame_value():
    result = Result.read_frd_text(FRD)
    solution = result.solution(7)
    frame = solution.frame(2)

    assert solution.type == "MODAL"
    assert frame.value == 5.0
    assert frame.field(FieldType.DISPLACEMENT)[17] == (1.0, 2.0, 3.0)
