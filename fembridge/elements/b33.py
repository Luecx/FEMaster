from .element import Element

class B33(Element):
    num_nodes = 2

    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'B33')

    def to_second_order(self, new_node_ids):
        return self

    def connectivity(self):
        return [(self.node_ids[0], self.node_ids[1])]

    def subdivide(self, edge_nodes, geometry, only_quads=False):
        n1, n2 = self.node_ids
        nmid = edge_nodes[(n1, n2)]
        geometry.add_element(B33(self.element_id, [n1, nmid]))
        geometry.add_element(B33(None, [nmid, n2]))

    def mirror_ids(self):
        self.node_ids = self.node_ids[::-1]
