
def write_input_deck(self, filename):
    from .fem_geometry import Geometry

    with open(filename, 'w') as file:
        # Write nodes
        file.write("*NODE, NSET=NALL\n")
        for id, node in enumerate(self.nodes):
            if node is not None:
                file.write(f"{id},{node[0]},{node[1]},{node[2]}\n")

        # Write elements by type
        element_classes = {k : [] for k in Geometry.element_classes.keys()}

        for element in self.elements:
            if element is not None:
                element_classes[element.elem_type].append(element)

        for elem_type, elements in element_classes.items():
            if elements:
                file.write(f"*ELEMENT, TYPE={elem_type}, ELSET=EALL\n")
                for element in elements:
                    nodes = ','.join(map(str, [n for n in element.node_ids]))
                    file.write(f"{element.element_id},{nodes}\n")

        # Write node sets
        for name, ids in self.node_sets.items():
            if name != "NALL":
                file.write(f"*NSET, NAME={name}\n")
                file.write(',\n'.join(map(str, [id for id in ids])) + '\n')

        # Write element sets
        for name, ids in self.elem_sets.items():
            if name != "EALL":
                file.write(f"*ELSET, NAME={name}\n")
                file.write(',\n'.join(map(str, [id for id in ids])) + '\n')