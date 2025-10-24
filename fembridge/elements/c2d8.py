from .element import Element

class C2D8(Element):
    num_nodes = 8
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'C2D8')

    def to_second_order(self, new_node_ids):
        return self  # Already second order

    def connectivity(self):
        return [(self.node_ids[0], self.node_ids[1]),
                (self.node_ids[1], self.node_ids[2]),
                (self.node_ids[2], self.node_ids[3]),
                (self.node_ids[3], self.node_ids[0]),
                (self.node_ids[0], self.node_ids[4]),
                (self.node_ids[1], self.node_ids[5]),
                (self.node_ids[2], self.node_ids[6]),
                (self.node_ids[3], self.node_ids[7])]

    def mirror_ids(self):
        self.node_ids = [self.node_ids[3], self.node_ids[2], self.node_ids[1], self.node_ids[0],
                         self.node_ids[6], self.node_ids[5], self.node_ids[4], self.node_ids[7]]
