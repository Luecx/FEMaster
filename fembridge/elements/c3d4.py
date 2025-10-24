from .element import Element

class C3D4(Element):
    num_nodes = 4
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'C3D4')

    def to_second_order(self, new_node_ids):
        from .c3d10 import C3D10
        n1, n2, n3, n4 = self.node_ids
        n5 = new_node_ids[(n1, n2)]
        n6 = new_node_ids[(n2, n3)]
        n7 = new_node_ids[(n3, n1)]
        n8 = new_node_ids[(n1, n4)]
        n9 = new_node_ids[(n2, n4)]
        n10 = new_node_ids[(n3, n4)]
        return C3D10(self.element_id, [n1, n2, n3, n4, n5, n6, n7, n8, n9, n10])

    def connectivity(self):
        return [(self.node_ids[0], self.node_ids[1]),
                (self.node_ids[1], self.node_ids[2]),
                (self.node_ids[2], self.node_ids[0]),
                (self.node_ids[0], self.node_ids[3]),
                (self.node_ids[1], self.node_ids[3]),
                (self.node_ids[2], self.node_ids[3])]

    def subdivide(self, edge_nodes, geometry, only_quads=False):
        n1, n2, n3, n4 = self.node_ids
        n5 = edge_nodes[(n1, n2)]
        n6 = edge_nodes[(n2, n3)]
        n7 = edge_nodes[(n3, n1)]
        n8 = edge_nodes[(n1, n4)]
        n9 = edge_nodes[(n2, n4)]
        n10 = edge_nodes[(n3, n4)]
        geometry.add_element(self.element_id,
                             element_type='C3D4', node_ids=[n1, n5, n7, n8])
        geometry.add_element(element_type='C3D4', node_ids=[n2, n6, n5, n9])
        geometry.add_element(element_type='C3D4', node_ids=[n3, n7, n6, n10])
        geometry.add_element(element_type='C3D4', node_ids=[n4, n8, n9, n10])
        geometry.add_element(element_type='C3D4', node_ids=[n5, n6, n7, n10])
        geometry.add_element(element_type='C3D4', node_ids=[n5, n6, n9, n10])
        geometry.add_element(element_type='C3D4', node_ids=[n5, n7, n8, n10])
        geometry.add_element(element_type='C3D4', node_ids=[n7, n6, n10, n5])

    def mirror_ids(self):
        self.node_ids = [self.node_ids[2], self.node_ids[3], self.node_ids[1], self.node_ids[0]]
