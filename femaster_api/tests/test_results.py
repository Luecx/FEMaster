from femaster_api import FieldDomain, FieldType, ResReader


RES = """\
LC 1 STATIC
FRAME 0
FIELD, NAME=DISPLACEMENT, TYPE=NODE, INDEX_COLS=1, VALUE_COLS=3, ROWS=2
               17    1.000000e+00    2.000000e+00    3.000000e+00
          bolt.17    4.000000e+00    5.000000e+00    6.000000e+00
END FIELD
FIELD, NAME=STRESS, TYPE=ELEMENT_IP, INDEX_COLS=2, VALUE_COLS=2, ROWS=1
           wing.4               2    1.000000e+02    2.000000e+02
END FIELD
"""


def test_res_reader_preserves_semantic_entity_ids():
    result = ResReader().parse(RES)
    displacement = result.field("DISPLACEMENT", loadcase=1, frame=0)

    assert displacement.domain is FieldDomain.NODE
    assert displacement.type is FieldType.DISPLACEMENT
    assert displacement[17] == (1.0, 2.0, 3.0)
    assert displacement["bolt.17"] == (4.0, 5.0, 6.0)

    stress = result.field("STRESS", loadcase=1, frame=0)
    assert stress.domain is FieldDomain.ELEMENT_IP
    assert stress[("wing.4", 2)] == (100.0, 200.0)
