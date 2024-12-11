class Element:
    def __init__(self, element_id, node_ids, elem_type):
        self.element_id = element_id
        self.node_ids = node_ids
        self.elem_type = elem_type

    def to_second_order(self, new_node_ids):
        raise NotImplementedError("This method should be implemented by subclasses")

    def compute_center(self, geometry):
        nodes = [geometry.nodes[node_id] for node_id in self.node_ids]
        x_center = sum(x for x, y, z in nodes) / len(nodes)
        y_center = sum(y for x, y, z in nodes) / len(nodes)
        z_center = sum(z for x, y, z in nodes) / len(nodes)
        return (x_center, y_center, z_center)

    def connectivity(self):
        raise NotImplementedError("This method should be implemented by subclasses")

    def subdivide(self, edge_nodes, geometry, only_quads=False):
        raise NotImplementedError("This method should be implemented by subclasses")

    def mirror_ids(self):
        raise NotImplementedError("This method should be implemented by subclasses")

class B33(Element):
    num_nodes = 2
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'B33')

    def to_second_order(self, new_node_ids):
        return self

    def connectivity(self):
        return [(self.node_ids[0], self.node_ids[1])]

    def subdivide(self, edge_nodes, geometry, only_quads=False):
        raise NotImplementedError("This method should be implemented by subclasses")

    def mirror_ids(self):
        self.node_ids = self.node_ids[::-1]

class Point(Element):
    num_nodes = 1
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'Point')

    def to_second_order(self, new_node_ids):
        return self

    def connectivity(self):
        return []

    def subdivide(self, edge_nodes, geometry, only_quads=False):
        raise NotImplementedError("This method should be implemented by subclasses")

    def mirror_ids(self):
        pass

class C2D3(Element):
    num_nodes = 3

    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'C2D3')

    def to_second_order(self, new_node_ids):
        n1, n2, n3 = self.node_ids
        n4 = new_node_ids[(n1, n2)]
        n5 = new_node_ids[(n2, n3)]
        n6 = new_node_ids[(n3, n1)]
        return C2D6(self.element_id, [n1, n2, n3, n4, n5, n6])

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
                                 element_type='C2D4', node_ids=[n1, n4, nm, n6])
            geometry.add_element(element_type='C2D4', node_ids=[n2, n5, nm, n4])
            geometry.add_element(element_type='C2D4', node_ids=[n3, n6, nm, n5])
        else:
            geometry.add_element(element_id=self.element_id,
                                 element_type='C2D3', node_ids=[n1, n4, n6])
            geometry.add_element(element_type='C2D3', node_ids=[n4, n2, n5])
            geometry.add_element(element_type='C2D3', node_ids=[n4, n5, n6])
            geometry.add_element(element_type='C2D3', node_ids=[n3, n6, n5])

    def mirror_ids(self):
        self.node_ids = self.node_ids[::-1]

class C2D4(Element):
    num_nodes = 4
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'C2D4')

    def to_second_order(self, new_node_ids):
        n1, n2, n3, n4 = self.node_ids
        n5 = new_node_ids[(n1, n2)]
        n6 = new_node_ids[(n2, n3)]
        n7 = new_node_ids[(n3, n4)]
        n8 = new_node_ids[(n4, n1)]
        return C2D8(self.element_id, [n1, n2, n3, n4, n5, n6, n7, n8])

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
                             element_type='C2D4', node_ids=[n1, n5, nm, n8])
        geometry.add_element(element_type='C2D4', node_ids=[n5, n2, n6, nm])
        geometry.add_element(element_type='C2D4', node_ids=[nm, n6, n3, n7])
        geometry.add_element(element_type='C2D4', node_ids=[n8, nm, n7, n4])

    def mirror_ids(self):
        self.node_ids = self.node_ids[::-1]

class C3D4(Element):
    num_nodes = 4
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'C3D4')

    def to_second_order(self, new_node_ids):
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

class C2D6(Element):
    num_nodes = 6
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'C2D6')

    def to_second_order(self, new_node_ids):
        return self  # Already second order

    def connectivity(self):
        return [(self.node_ids[0], self.node_ids[3]),
                (self.node_ids[3], self.node_ids[1]),
                (self.node_ids[1], self.node_ids[4]),
                (self.node_ids[4], self.node_ids[2]),
                (self.node_ids[2], self.node_ids[5]),
                (self.node_ids[5], self.node_ids[0])]


    def mirror_ids(self):
        self.node_ids = [self.node_ids[2], self.node_ids[1], self.node_ids[0],
                         self.node_ids[4], self.node_ids[3], self.node_ids[5]]

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
        # TODO: Implement this method
        n1, n2, n3, n4, n5, n6 = self.node_ids
        n7 = edge_nodes[(n1, n2)]
        n8 = edge_nodes[(n2, n3)]
        n9 = edge_nodes[(n3, n1)]
        n10 = edge_nodes[(n4, n5)]
        n11 = edge_nodes[(n5, n6)]
        n12 = edge_nodes[(n6, n4)]
        n13 = edge_nodes[(n1, n4)]
        n14 = edge_nodes[(n2, n5)]
        n15 = edge_nodes[(n3, n6)]

        center1 = ((geometry.nodes[n1][0] + geometry.nodes[n2][0] + geometry.nodes[n3][0]) / 3,
                   (geometry.nodes[n1][1] + geometry.nodes[n2][1] + geometry.nodes[n3][1]) / 3,
                   (geometry.nodes[n1][2] + geometry.nodes[n2][2] + geometry.nodes[n3][2]) / 3)
        center2 = ((geometry.nodes[n4][0] + geometry.nodes[n5][0] + geometry.nodes[n6][0]) / 3,
                   (geometry.nodes[n4][1] + geometry.nodes[n5][1] + geometry.nodes[n6][1]) / 3,
                   (geometry.nodes[n4][2] + geometry.nodes[n5][2] + geometry.nodes[n6][2]) / 3)

        n_center1 = geometry.add_node(x=center1[0], y=center1[1], z=center1[2])
        n_center2 = geometry.add_node(x=center2[0], y=center2[1], z=center2[2])

        geometry.add_element(element_type='C3D6', node_ids=[n1, n7, n_center1, n4, n10, n_center2])
        geometry.add_element(element_type='C3D6', node_ids=[n7, n2, n8, n10, n5, n11])
        geometry.add_element(element_type='C3D6', node_ids=[n_center1, n8, n3, n_center2, n11, n6])
        geometry.add_element(element_type='C3D6', node_ids=[n9, n_center1, n3, n12, n_center2, n6])
        geometry.add_element(element_type='C3D6', node_ids=[n1, n9, n_center1, n4, n12, n_center2])
        geometry.add_element(element_type='C3D6', node_ids=[n7, n8, n_center1, n10, n11, n_center2])
        geometry.add_element(element_type='C3D6', node_ids=[n8, n9, n_center1, n11, n12, n_center2])
        geometry.add_element(element_type='C3D6', node_ids=[n9, n7, n_center1, n12, n10, n_center2])

    def mirror_ids(self):
        self.node_ids = [self.node_ids[2], self.node_ids[1], self.node_ids[0],
                         self.node_ids[5], self.node_ids[4], self.node_ids[3]]

class C3D8(Element):
    num_nodes = 8
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'C3D8')

    def to_second_order(self, new_node_ids):
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

class C3D5(C3D8):
    num_nodes = 5
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, [node_ids[0], node_ids[1], node_ids[2], node_ids[3]
                                            , node_ids[4], node_ids[4], node_ids[4], node_ids[4]])


class C3D10(Element):
    num_nodes = 10
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'C3D10')

    def to_second_order(self, new_node_ids):
        return self  # Already second order

    def connectivity(self):
        return [(self.node_ids[0], self.node_ids[1]),
                (self.node_ids[1], self.node_ids[2]),
                (self.node_ids[2], self.node_ids[0]),
                (self.node_ids[0], self.node_ids[3]),
                (self.node_ids[1], self.node_ids[3]),
                (self.node_ids[2], self.node_ids[3]),
                (self.node_ids[0], self.node_ids[4]),
                (self.node_ids[1], self.node_ids[5]),
                (self.node_ids[2], self.node_ids[6]),
                (self.node_ids[3], self.node_ids[7]),
                (self.node_ids[4], self.node_ids[5]),
                (self.node_ids[5], self.node_ids[6]),
                (self.node_ids[6], self.node_ids[7]),
                (self.node_ids[7], self.node_ids[4])]

    def mirror_ids(self):
        # mirroring using: 3 2 1 4 6 5 7 10 9 8
        self.node_ids = [self.node_ids[2], self.node_ids[1], self.node_ids[0], self.node_ids[3],
                         self.node_ids[5], self.node_ids[4], self.node_ids[6], self.node_ids[9],
                         self.node_ids[8], self.node_ids[7]]

class C3D15(Element):
    num_nodes = 15
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'C3D15')

    def to_second_order(self, new_node_ids):
        return self  # Already second order

    def connectivity(self):
        return [(self.node_ids[0], self.node_ids[1]),
                (self.node_ids[1], self.node_ids[2]),
                (self.node_ids[2], self.node_ids[0]),
                (self.node_ids[3], self.node_ids[4]),
                (self.node_ids[4], self.node_ids[5]),
                (self.node_ids[5], self.node_ids[3]),
                (self.node_ids[0], self.node_ids[3]),
                (self.node_ids[1], self.node_ids[4]),
                (self.node_ids[2], self.node_ids[5]),
                (self.node_ids[0], self.node_ids[6]),
                (self.node_ids[1], self.node_ids[7]),
                (self.node_ids[2], self.node_ids[8]),
                (self.node_ids[3], self.node_ids[9]),
                (self.node_ids[4], self.node_ids[10]),
                (self.node_ids[5], self.node_ids[11]),
                (self.node_ids[6], self.node_ids[7]),
                (self.node_ids[7], self.node_ids[8]),
                (self.node_ids[8], self.node_ids[9]),
                (self.node_ids[9], self.node_ids[10]),
                (self.node_ids[10], self.node_ids[11]),
                (self.node_ids[11], self.node_ids[6])]

    def mirror_ids(self):
        # mirroring using: 3 2 1 6 5 4, 8 7 9 11 10 12, 15 14 13
        self.node_ids = [self.node_ids[2], self.node_ids[1], self.node_ids[0], self.node_ids[5],
                            self.node_ids[4], self.node_ids[3], self.node_ids[7], self.node_ids[6],
                            self.node_ids[8], self.node_ids[10], self.node_ids[9], self.node_ids[11],
                            self.node_ids[14], self.node_ids[13], self.node_ids[12]]

class C3D20(Element):
    num_nodes = 20
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'C3D20')

    def to_second_order(self, new_node_ids):
        return self  # Already second order

    def connectivity(self):
        return [(self.node_ids[0], self.node_ids[1]),
                (self.node_ids[1], self.node_ids[2]),
                (self.node_ids[2], self.node_ids[3]),
                (self.node_ids[3], self.node_ids[0]),
                (self.node_ids[4], self.node_ids[5]),
                (self.node_ids[5], self.node_ids[6]),
                (self.node_ids[6], self.node_ids[7]),
                (self.node_ids[7], self.node_ids[4]),
                (self.node_ids[0], self.node_ids[8]),
                (self.node_ids[1], self.node_ids[9]),
                (self.node_ids[2], self.node_ids[10]),
                (self.node_ids[3], self.node_ids[11]),
                (self.node_ids[4], self.node_ids[12]),
                (self.node_ids[5], self.node_ids[13]),
                (self.node_ids[6], self.node_ids[14]),
                (self.node_ids[7], self.node_ids[15]),
                (self.node_ids[8], self.node_ids[9]),
                (self.node_ids[9], self.node_ids[10]),
                (self.node_ids[10], self.node_ids[11]),
                (self.node_ids[11], self.node_ids[12]),
                (self.node_ids[12], self.node_ids[13]),
                (self.node_ids[13], self.node_ids[14]),
                (self.node_ids[14], self.node_ids[15])]

    def mirror_ids(self):
        # mirroring using: 4, 3, 2, 1, 8, 7, 6, 5 | 11, 10, 9, 12, 15, 14, 13, 16 | 20, 19, 18, 17
        self.node_ids = [self.node_ids[3], self.node_ids[2], self.node_ids[1], self.node_ids[0],
                            self.node_ids[7], self.node_ids[6], self.node_ids[5], self.node_ids[4],
                            self.node_ids[10], self.node_ids[9], self.node_ids[8], self.node_ids[11],
                            self.node_ids[14], self.node_ids[13], self.node_ids[12], self.node_ids[15],
                            self.node_ids[19], self.node_ids[18], self.node_ids[17], self.node_ids[16]]