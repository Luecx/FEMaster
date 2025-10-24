from .c3d8 import C3D8

class C3D5(C3D8):
    num_nodes = 5
    def __init__(self, element_id, node_ids):
        # Expand to 8 with repeated 5th as in original code
        super().__init__(element_id, [node_ids[0], node_ids[1], node_ids[2], node_ids[3],
                                      node_ids[4], node_ids[4], node_ids[4], node_ids[4]])
