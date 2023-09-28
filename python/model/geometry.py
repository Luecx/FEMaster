class Geometry:
    def __init__(self):
        self.nodes = []
        self.elements = []

        self.node_sets = {"NALL": []}
        self.elem_sets = {"EALL": []}

    def add_node(self, node_id, x, y, z=0):
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

    def compute_element_midpoints(self):
        element_midpoints = []
        for element in self.elements:
            if element is not None:
                node_ids = element['nodes']
                nodes = [self.nodes[node_id] for node_id in node_ids]
                x_mid = sum(x for x, y, z in nodes) / len(nodes)
                y_mid = sum(y for x, y, z in nodes) / len(nodes)
                z_mid = sum(z for x, y, z in nodes) / len(nodes)
                element_midpoints.append((x_mid, y_mid, z_mid))
            else:
                element_midpoints.append(None)
        return element_midpoints

    def __str__(self):
        ret_str = ""
        ret_str += "Nodes    : {}\n".format(len(self.nodes))
        ret_str += "Elements : {}\n".format(len(self.elements))
        ret_str += "Node Sets:\n"
        for name, ids in self.node_sets.items():
            ret_str += "  {}: {} nodes\n".format(name, len(ids))
        ret_str += "Element Sets:\n"
        for name, ids in self.elem_sets.items():
            ret_str += "  {}: {} elements\n".format(name, len(ids))
        return ret_str

    def extrude(self, k, spacing=1):
        new_geom = Geometry()

        max_node_id = len(self.nodes)
        max_elem_id = len(self.elements)

        # create / copy nodes (k+1) times
        for i in range(k+1):
            # create nodes
            for id, node in enumerate(self.nodes):
                if node is not None:
                    x, y, z = node
                    new_geom.add_node(i * max_node_id + id, x, y, z + i * spacing)

        for i in range(k):
            # create elements
            for id, ids in enumerate(self.elements):
                if ids is not None:
                    ids1 = [m + i     * max_node_id for m in ids['nodes']]
                    ids2 = [m + (i+1) * max_node_id for m in ids['nodes']]

                    ids  = ids1 + ids2
                    new_geom.add_element(i * max_elem_id + id, element_type=f"C3D{len(ids)}", node_ids=ids)

        # copy all the old sets and adjust the nodes
        for nset in self.node_sets:
            new_geom.add_node_set(nset)
            for n in self.node_sets[nset]:
                for i in range(k+1):
                    new_geom.node_sets[nset].append(n + i * max_node_id)

        for elset in self.elem_sets:
            new_geom.add_element_set(elset)
            for n in self.elem_sets[elset]:
                for i in range(k):
                    new_geom.elem_sets[elset].append(n + i * max_elem_id)
        return new_geom

    @staticmethod
    def read_input_deck(filename):
        geometry = Geometry()
        with open(filename, 'r') as file:
            lines = file.readlines()

        def parse_line(line, type=float):
            if line:
                s = line.split(',')
                s = [x for x in s if x]  # remove empty entries
            return int(s[0]), list(map(type, s[1:]))

        def parse_command(lines):
            result = {}
            for item in lines:
                item = item.replace(' ', '')
                if '=' in item:
                    k, v = item.split('=')
                    result[k] = v
            return result

        key_word = None
        elem_type = None
        nset_name = None
        elset_name = None

        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith('**'):
                continue
            if line.startswith('*'):
                params = line.split(',')
                key_word = params[0].strip()
                keys = parse_command(params[1:])

                if key_word == "*NODE":
                    if 'NSET' in keys:
                        nset_name = keys['NSET']
                        if nset_name != "NALL":
                            geometry.add_node_set(nset_name)
                        else:
                            nset_name = None
                    else:
                        nset_name = None

                if key_word == "*ELEMENT":
                    if 'ELSET' in keys:
                        elset_name = keys['ELSET']
                        if elset_name != "EALL":
                            geometry.add_element_set(elset_name)
                        else:
                            elset_name = None
                    else:
                        elset_name = None
                    elem_type = keys['TYPE']

                if key_word == "*NSET":
                    nset_name = keys['NAME']
                    geometry.add_node_set(nset_name)

                if key_word == "*ELSET":
                    elset_name = keys['NAME']
                    geometry.add_element_set(elset_name)
            else:
                if key_word == "*NODE":
                    id, coords = parse_line(line, float)
                    geometry.add_node(id, *coords)
                    if nset_name:
                        geometry.node_sets[nset_name].append(id)
                    geometry.node_sets["NALL"].append(id)

                if key_word == "*ELEMENT":
                    id, ids = parse_line(line, int)
                    geometry.add_element(id, elem_type, list(ids))
                    if elset_name:
                        geometry.elem_sets[elset_name].append(id)
                    geometry.elem_sets["EALL"].append(id)

                if key_word == "*NSET":
                    id, ids = parse_line(line, int)
                    ids = list(ids)
                    if ids:
                        geometry.node_sets[nset_name].extend([id, *ids])
                    else:
                        geometry.node_sets[nset_name].append(id)

                if key_word == "*ELSET":
                    id, ids = parse_line(line, int)
                    ids = list(ids)
                    if ids:
                        geometry.elem_sets[elset_name].extend([id, *ids])
                    else:
                        geometry.elem_sets[elset_name].append(id)

        return geometry

    def write_input_deck(self, filename):
        with open(filename, 'w') as file:
            # write nodes
            file.write("*NODE, NSET=NALL\n")
            for id, node in enumerate(self.nodes):
                if node is not None:
                    file.write(f"{id + 1},{node[0]},{node[1]},{node[2]}\n")

            # write elements
            file.write("*ELEMENT, TYPE=C3D8, ELSET=EALL\n")
            for id, element in enumerate(self.elements):
                if element is not None:
                    nodes = ','.join(map(str, [n + 1 for n in element['nodes']]))
                    file.write(f"{id + 1},{nodes}\n")

            # write node sets
            for name, ids in self.node_sets.items():
                if name != "NALL":
                    file.write(f"*NSET, NAME={name}\n")
                    file.write(',\n'.join(map(str, [id + 1 for id in ids])) + '\n')

            # write element sets
            for name, ids in self.elem_sets.items():
                if name != "EALL":
                    file.write(f"*ELSET, NAME={name}\n")
                    file.write(',\n'.join(map(str, [id + 1 for id in ids])) + '\n')

geom = Geometry.read_input_deck("../../naca_flat.txt")
geom2 =geom.extrude(1, spacing=10)
geom2.write_input_deck("../../naca_extruded_mini.inp")