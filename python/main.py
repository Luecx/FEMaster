from geometry import *
from topopt import *

# outer = SegmentGroup([StraightSegment([0, 0], [5, 0], subdivisions=40),
#                       StraightSegment([5, 0], [5, 5], subdivisions=40),
#                       StraightSegment([5, 5], [0, 5], subdivisions=40),
#                       StraightSegment([0, 5], [0, 0], subdivisions=40)])
#
# inner = SegmentGroup([StraightSegment([2, 2], [3, 2], subdivisions=10),
#                       StraightSegment([3, 2], [3, 3], subdivisions=10),
#                       CurvedSegment  ([3, 3], [2, 3], mid=[2.5, 3.5], subdivisions=10),
#                       StraightSegment([2, 3], [2, 2], subdivisions=10, name="MY_CUSTOM_NAME_1")])
#
# # from 4 to 4.5
# inner2 = SegmentGroup([StraightSegment([4, 4], [4.5, 4], subdivisions=10, name="HERE_I_WANT_LOADS"),
#                        StraightSegment([4.5, 4], [4.5, 4.5], subdivisions=10),
#                        StraightSegment([4.5, 4.5], [4, 4.5], subdivisions=10),
#                        StraightSegment([4, 4.5], [4, 4], subdivisions=10)])
#
# geom = Geometry.mesh_interior([outer, inner, inner2])
# geom = geom.extruded(5, spacing=0.1)
# geom = geom.as_second_order()
# geom.write_input_deck("my_deck.inp")
