from .element import Element

class Point(Element):
    num_nodes = 1
    def __init__(self, element_id, node_ids):
        super().__init__(element_id, node_ids, 'P')

    def to_second_order(self, new_node_ids):
        return self

    def connectivity(self):
        return []

    def subdivide(self, edge_nodes, geometry, only_quads=False):
        raise NotImplementedError("This method should be implemented by subclasses")

    def mirror_ids(self):
        pass
