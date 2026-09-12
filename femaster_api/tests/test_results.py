from femaster_api import FieldDomain, FieldType, ResReader


RES = """\
LC 1 STATIC
FRAME 0
FIELD DISP ROWS=2 COLS=4 INDEX_COLS=1 VALUE_COLS=3 TYPE=NODE
1 1.0 2.0 3.0
2 4.0 5.0 6.0
"""


def test_res_reader_uses_common_field_types():
    result = ResReader().parse(RES)
    field = result.field("DISP", loadcase=1, frame=0)

    assert field.domain is FieldDomain.NODE
    assert field.type is FieldType.DISPLACEMENT
    assert field[1] == (1.0, 2.0, 3.0)
    assert field[2] == (4.0, 5.0, 6.0)
