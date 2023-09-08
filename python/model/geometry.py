class Geometry:
    def __init__(self):
        self.nodes = []
        self.elements = []

        self.node_sets = {}
        self.elem_sets = {}

    def add_node(self, node_id, x, y, z):
        # Resize nodes list if necessary
        if node_id >= len(self.nodes):
            self.nodes.extend([None] * (node_id - len(self.nodes) + 1))
        self.nodes[node_id] = (x, y, z)

    def add_element(self, element_id, element_type, node_ids):
        # Resize elements list if necessary
        if element_id >= len(self.elements):
            self.elements.extend([None] * (element_id - len(self.elements) + 1))
        self.elements[element_id] = {'type': element_type, 'nodes': node_ids}

    def add_node_set(self, name):
        if name not in self.node_sets:
            self.node_sets[name] = []

    def add_element_set(self, name):
        if name not in self.elem_sets:
            self.elem_sets[name] = []

    @staticmethod
    def read_input_deck(filename):
        geometry = Geometry()
        with open(filename, 'r') as file:
            lines = file.readlines()

        def parse_line(line, type=float):
            if line:
                s = line.split(',')
            return int(s[0]), map(type, s[1:])

        key_word = None
        elem_type = None

        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith('**'):
                continue
            if line.startswith('*'):
                params = line.split(',')
                key_word = params[0].strip()

                # parse element type
                if key_word == "*ELEMENT":
                    for v in params[1:]:
                        if 'TYPE' in v:
                            elem_type = v.split("=")[1]
            else:
                if key_word == "*NODE":
                    id, coords = parse_line(line, float)
                    geometry.add_node(id, *coords)

                if key_word == "*ELEMENT":
                    id, ids = parse_line(line, int)
                    geometry.add_element(id, elem_type, list(ids))
        return geometry


