from .element import Element

class C3D8(Element):
    num_nodes = 8
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'C3D8')

    def to_second_order(self, new_node_ids):
        from .c3d20 import C3D20
        n1, n2, n3, n4, n5, n6, n7, n8 = self.node_ids
        n9 = new_node_ids[(n1, n2)]
        n10 = new_node_ids[(n2, n3)]
        n11 = new_node_ids[(n3, n4)]
        n12 = new_node_ids[(n4, n1)]
        n13 = new_node_ids[(n5, n6)]
        n14 = new_node_ids[(n6, n7)]
        n15 = new_node_ids[(n7, n8)]
        n16 = new_node_ids[(n8, n5)]
        n17 = new_node_ids[(n1, n5)]
        n18 = new_node_ids[(n2, n6)]
        n19 = new_node_ids[(n3, n7)]
        n20 = new_node_ids[(n4, n8)]
        return C3D20(self.element_id, [n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20])

    def connectivity(self):
        return [(self.node_ids[0], self.node_ids[1]),
                (self.node_ids[1], self.node_ids[2]),
                (self.node_ids[2], self.node_ids[3]),
                (self.node_ids[3], self.node_ids[0]),
                (self.node_ids[4], self.node_ids[5]),
                (self.node_ids[5], self.node_ids[6]),
                (self.node_ids[6], self.node_ids[7]),
                (self.node_ids[7], self.node_ids[4]),
                (self.node_ids[0], self.node_ids[4]),
                (self.node_ids[1], self.node_ids[5]),
                (self.node_ids[2], self.node_ids[6]),
                (self.node_ids[3], self.node_ids[7])]

    def subdivide(self, edge_nodes, geometry, only_quads=False):
        n1, n2, n3, n4, n5, n6, n7, n8 = self.node_ids
        n9 = edge_nodes[(n1, n2)]
        n10 = edge_nodes[(n2, n3)]
        n11 = edge_nodes[(n3, n4)]
        n12 = edge_nodes[(n4, n1)]
        n13 = edge_nodes[(n5, n6)]
        n14 = edge_nodes[(n6, n7)]
        n15 = edge_nodes[(n7, n8)]
        n16 = edge_nodes[(n8, n5)]
        n17 = edge_nodes[(n1, n5)]
        n18 = edge_nodes[(n2, n6)]
        n19 = edge_nodes[(n3, n7)]
        n20 = edge_nodes[(n4, n8)]

        center = self.compute_center(geometry)
        n_center = geometry.add_node(x=center[0], y=center[1], z=center[2])

        geometry.add_element(element_type='C3D8', node_ids=[n1, n9, n_center, n12, n17, n_center, n_center, n16])
        geometry.add_element(element_type='C3D8', node_ids=[n9, n2, n10, n_center, n_center, n18, n_center, n_center])
        geometry.add_element(element_type='C3D8', node_ids=[n_center, n10, n3, n11, n_center, n_center, n19, n_center])
        geometry.add_element(element_type='C3D8', node_ids=[n12, n_center, n11, n4, n16, n_center, n_center, n20])
        geometry.add_element(element_type='C3D8', node_ids=[n17, n_center, n_center, n16, n5, n13, n_center, n_center])
        geometry.add_element(element_type='C3D8', node_ids=[n_center, n18, n_center, n_center, n13, n6, n14, n_center])
        geometry.add_element(element_type='C3D8', node_ids=[n_center, n_center, n19, n_center, n_center, n14, n7, n15])
        geometry.add_element(element_type='C3D8', node_ids=[n16, n_center, n_center, n20, n_center, n_center, n15, n8])

    def mirror_ids(self):
        self.node_ids = [self.node_ids[3], self.node_ids[2], self.node_ids[1], self.node_ids[0],
                         self.node_ids[7], self.node_ids[6], self.node_ids[5], self.node_ids[4]]
