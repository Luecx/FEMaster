
def read_input_deck(filename):
    from .fem_geometry import Geometry
    import re
    import tqdm

    geometry = Geometry()
    with open(filename, 'r') as file:
        lines = file.readlines()

    def preprocess_line(line):
        """Convert the line to upper case, remove spaces, and replace them with commas."""
        return re.sub(r'\s+', ',', line.strip().upper())

    def is_command(line):
        """Check if the line starts with '*' indicating a command."""
        return line.startswith('*')

    def require(key_names, parts, default=None, options=None):
        """Extract a key value from the parts based on the provided key names."""
        for key in key_names:
            if key in parts:
                key_index = parts.index(key)
                if key_index + 1 < len(parts):
                    if options is not None and parts[key_index + 1] not in options:
                        raise ValueError(f"Invalid value for {key}: {parts[key_index + 1]} \n\tExpected: {options}")
                    return parts[key_index + 1]
        return default

    def parse_line(line, dtype=float):
        """Parse a line of data into a list of values, converting to a specified type."""
        line = re.sub(r',\s*$', '', line)  # Remove trailing commas and spaces
        return list(map(dtype, re.split(r'[,\s]+', line.strip())))

    key_word        = None
    current_nset    = None
    current_elset   = None
    index = 0

    # init progress bar
    pbar = tqdm.tqdm(total=len(lines), desc="Reading input deck")

    while index < len(lines):
        # since we skip a few indices, set the pbar value to the index
        pbar.update(index - pbar.n)

        line = preprocess_line(lines[index])

        if not line or line.startswith('**'):  # Skip comments or empty lines
            index += 1
            continue

        if is_command(line):
            parts = re.split(r'[,=\s]+', line)
            key_word = parts[0].upper()

            # Process different commands
            if key_word == '*NODE':
                current_nset = require(["NSET", "NAME"], parts, default='NALL')
                geometry.add_node_set(current_nset)

            elif key_word == '*ELEMENT':
                elem_type = require(["TYPE"], parts, options=Geometry.element_classes.keys())
                current_elset = require(["ELSET", "NAME"], parts, default='EALL')
                geometry.add_element_set(current_elset)

            elif key_word == '*NSET':
                current_nset = require(["NAME", "NSET"], parts)
                geometry.add_node_set(current_nset)

            elif key_word == '*ELSET':
                current_elset = require(["NAME", "ELSET"], parts)
                geometry.add_element_set(current_elset)

            else:
                key_word = None  # Unknown command
            index += 1

        else:
            if key_word is None:
                index += 1
                continue

            data = parse_line(line)

            # Process node data
            if key_word == '*NODE':
                node_id, *coords = data
                geometry.add_node(int(node_id), *coords)
                if current_nset:
                    geometry.add_node_to_set(current_nset, int(node_id))

            # Process element data
            elif key_word == '*ELEMENT' and elem_type in Geometry.element_classes:
                required_nodes = Geometry.element_classes[elem_type].num_nodes
                element_data = data

                # Read more lines if the number of nodes is not enough
                while len(element_data) - 1 < required_nodes:
                    index += 1
                    extra_data = parse_line(preprocess_line(lines[index]))
                    element_data.extend(extra_data)

                element_id, *node_ids = map(int, element_data)
                geometry.add_element(int(element_id), elem_type, [nid for nid in node_ids])
                if current_elset:
                    geometry.add_element_to_set(current_elset, int(element_id))

            # Process NSET data
            elif key_word == '*NSET' and current_nset:
                for nid in data:
                    geometry.add_node_to_set(current_nset, int(nid))

            # Process ELSET data
            elif key_word == '*ELSET' and current_elset:
                for eid in data:
                    geometry.add_element_to_set(current_elset, int(eid))

            index += 1

    pbar.close()
    return geometry