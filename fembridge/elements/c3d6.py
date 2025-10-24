from .element import Element

class C3D6(Element):
    num_nodes = 6
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'C3D6')

    def to_second_order(self, new_node_ids):
        return self  # Already second order

    def connectivity(self):
        return [(self.node_ids[0], self.node_ids[1]),
                (self.node_ids[1], self.node_ids[2]),
                (self.node_ids[2], self.node_ids[0]),
                (self.node_ids[0], self.node_ids[3]),
                (self.node_ids[1], self.node_ids[4]),
                (self.node_ids[2], self.node_ids[5]),
                (self.node_ids[3], self.node_ids[4]),
                (self.node_ids[4], self.node_ids[5]),
                (self.node_ids[5], self.node_ids[3])]

    def subdivide(self, edge_nodes, geometry, only_quads=False):
        raise NotImplementedError("SOME ERROR IN THIS CODE")
        # (original experimental code kept intentionally, per your note)

    def mirror_ids(self):
        self.node_ids = [self.node_ids[2], self.node_ids[1], self.node_ids[0],
                         self.node_ids[5], self.node_ids[4], self.node_ids[3]]
