from .element import Element

class S3(Element):
    num_nodes = 3

    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'S3')

    def to_second_order(self, new_node_ids):
        from .s6 import S6
        n1, n2, n3 = self.node_ids
        n4 = new_node_ids[(n1, n2)]
        n5 = new_node_ids[(n2, n3)]
        n6 = new_node_ids[(n3, n1)]
        return S6(self.element_id, [n1, n2, n3, n4, n5, n6])

    def connectivity(self):
        return [(self.node_ids[0], self.node_ids[1]),
                (self.node_ids[1], self.node_ids[2]),
                (self.node_ids[2], self.node_ids[0])]

    def subdivide(self, edge_nodes, geometry, only_quads=False):
        n1, n2, n3 = self.node_ids
        n4 = edge_nodes[(n1, n2)]
        n5 = edge_nodes[(n2, n3)]
        n6 = edge_nodes[(n3, n1)]

        if only_quads:
            nm = geometry.add_node(x=(geometry.nodes[n1][0] + geometry.nodes[n2][0] + geometry.nodes[n3][0]) / 3,
                                   y=(geometry.nodes[n1][1] + geometry.nodes[n2][1] + geometry.nodes[n3][1]) / 3,
                                   z=(geometry.nodes[n1][2] + geometry.nodes[n2][2] + geometry.nodes[n3][2]) / 3)

            geometry.add_element(element_id=self.element_id,
                                 element_type='S4', node_ids=[n1, n4, nm, n6])
            geometry.add_element(element_type='S4', node_ids=[n2, n5, nm, n4])
            geometry.add_element(element_type='S4', node_ids=[n3, n6, nm, n5])
        else:
            geometry.add_element(element_id=self.element_id,
                                 element_type='S3', node_ids=[n1, n4, n6])
            geometry.add_element(element_type='S3', node_ids=[n4, n2, n5])
            geometry.add_element(element_type='S3', node_ids=[n4, n5, n6])
            geometry.add_element(element_type='S3', node_ids=[n3, n6, n5])

    def mirror_ids(self):
        self.node_ids = self.node_ids[::-1]
