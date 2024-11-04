def extruded_geometry(self, n, spacing=1):

    # Lazy import to avoid circular dependency
    from .fem_geometry import Geometry

    # Determine element order
    element_order = self.determine_element_order()
    if element_order == 0:
        raise ValueError("Element order could not be determined or no elements are present.")

    geometry = Geometry(3)  # New 3D geometry
    max_node_id = len(self.nodes)

    # Copy nodes and node sets
    for i in range(n * element_order + 1):
        for node_id, node in enumerate(self.nodes):
            if node is not None:
                x, y, z = node
                new_node_id = i * max_node_id + node_id
                geometry.add_node(new_node_id, x, y, z + i * spacing / (2 if element_order == 2 else 1))

    new_element_id = 1
    # Copy elements and element sets
    for i in range(n):
        for element in self.elements:
            if element is not None:
                if element_order == 1:
                    node_ids_layer1 = [i * max_node_id + node_id for node_id in element.node_ids]
                    node_ids_layer2 = [(i + 1) * max_node_id + node_id for node_id in element.node_ids]

                    new_ids = node_ids_layer1 + node_ids_layer2

                else:
                    node_ids_layer1 = [i * max_node_id + node_id for node_id in element.node_ids]
                    node_ids_layer2 = [(i + 2) * max_node_id + node_id for node_id in element.node_ids]
                    node_ids_layerm = [(i + 1) * max_node_id + node_id for node_id in element.node_ids]

                    # remove the second half of the mid layer
                    node_ids_layerm = node_ids_layerm[:len(node_ids_layerm) // 2]

                    # split the layer1 into two parts of equal size
                    node_ids_layer1a = node_ids_layer1[:len(node_ids_layer1) // 2]
                    node_ids_layer1b = node_ids_layer1[len(node_ids_layer1) // 2:]

                    # split the layer2 into two parts
                    node_ids_layer2a = node_ids_layer2[:len(node_ids_layer2) // 2]
                    node_ids_layer2b = node_ids_layer2[len(node_ids_layer2) // 2:]

                    new_ids = node_ids_layer1a + node_ids_layer2a + node_ids_layer1b + node_ids_layer2b + node_ids_layerm


                geometry.add_element(new_element_id, "C3D"+str(len(new_ids)), new_ids)
                new_element_id += 1
    # Copy sets
    for name, ids in self.node_sets.items():
        if name == "NALL":
            continue
        geometry.add_node_set(name)
        for node_id in ids:
            for i in range(n + 1):
                geometry.node_sets[name].append(i * max_node_id + node_id)

    for name, ids in self.elem_sets.items():
        if name == "EALL":
            continue
        geometry.add_element_set(name)
        for elem_id in ids:
            for i in range(n):
                geometry.elem_sets[name].append(i * len(self.elements) + elem_id)

    return geometry