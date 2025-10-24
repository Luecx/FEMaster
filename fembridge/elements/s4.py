from .element import Element

class S4(Element):
    num_nodes = 4
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'S4')

    def to_second_order(self, new_node_ids):
        from .s8 import S8
        n1, n2, n3, n4 = self.node_ids
        n5 = new_node_ids[(n1, n2)]
        n6 = new_node_ids[(n2, n3)]
        n7 = new_node_ids[(n3, n4)]
        n8 = new_node_ids[(n4, n1)]
        return S8(self.element_id, [n1, n2, n3, n4, n5, n6, n7, n8])

    def connectivity(self):
        return [(self.node_ids[0], self.node_ids[1]),
                (self.node_ids[1], self.node_ids[2]),
                (self.node_ids[2], self.node_ids[3]),
                (self.node_ids[3], self.node_ids[0])]

    def subdivide(self, edge_nodes, geometry, only_quads=False):
        n1, n2, n3, n4 = self.node_ids
        n5 = edge_nodes[(n1, n2)]
        n6 = edge_nodes[(n2, n3)]
        n7 = edge_nodes[(n3, n4)]
        n8 = edge_nodes[(n4, n1)]
        center = self.compute_center(geometry)
        nm = geometry.add_node(x=center[0], y=center[1], z=center[2])
        geometry.add_element(element_id=self.element_id,
                             element_type='S4', node_ids=[n1, n5, nm, n8])
        geometry.add_element(element_type='S4', node_ids=[n5, n2, n6, nm])
        geometry.add_element(element_type='S4', node_ids=[nm, n6, n3, n7])
        geometry.add_element(element_type='S4', node_ids=[n8, nm, n7, n4])

    def mirror_ids(self):
        # Explicit mapping; do not use default reverse
        self.node_ids = self.node_ids[::-1]
