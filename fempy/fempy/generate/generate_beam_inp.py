import yaml
import argparse

from ..geometry import Geometry

def parse_config(config_file):
    """Parse the YAML configuration file."""
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)
    raise ValueError("Invalid configuration file.")

def generate_geom(sections):
    x = 0.0
    geom = Geometry()
    geom.add_node(x=x)


    for section in sections:
        length  = section['length']
        n_elems = section['elements']
        name    = section['name']

        dx = length / n_elems

        for i in range(n_elems):
            node_id = geom.add_node(x=x+dx)
            elem_id = geom.add_element(element_type='B33', node_ids=[node_id-1, node_id])
            geom.add_element_to_set(name, elem_id)
            x += dx

    return geom


def write_input_deck(config, output_file="output.inp"):
    geom = generate_geom(config['sections'])
    geom.write_input_deck(output_file)

    def get(dict, val, default=0):
        return dict[val] if val in dict else default

    with open(output_file, 'a') as out:
        # material
        material = config['material']
        out.write(f"*MATERIAL, NAME=MAT\n")
        out.write(f"*ELASTIC, TYPE=ISOTROPIC\n{material['E']}, {get(material, 'nu')}\n")
        out.write(f"*DENSITY\n{get(material,'rho')}\n\n")
        # sections
        for section in config['sections']:
            out.write(f"*PROFILE, NAME={section['name']}\n")
            out.write(f"{section['area']}, {section['Ix']}, {section['Iy']}, {section['It']}\n")
            out.write(f"*BEAM SECTION, ELSET={section['name']}, MATERIAL=MAT, PROFILE={section['name']}\n")
            out.write(f"0, 1, 0\n")

        if 'bc' in config:
            out.write("*SUPPORT, SUPPORT_COLLECTOR=SUPPS\n")
            if 'left' in config['bc']:
                out.write(f"0, {config['bc']['left']}\n")
            if 'right' in config['bc']:
                out.write(f"{len(geom.nodes)-1}, {config['bc']['right']}\n")

        if 'loads' in config:
            out.write("*CLOAD, LOAD_COLLECTOR=LOADS\n")
            if 'left' in config['loads']:
                out.write(f"1, {config['loads']['left']}\n")
            if 'right' in config['loads']:
                out.write(f"{len(geom.nodes)-1}, {config['loads']['right']}\n")

        out.write(f'*LOADCASE, TYPE={config["scenario"]}\n')
        if 'bc' in config:
            out.write("*SUPPORT\nSUPPS\n")
        if 'loads' in config:
            out.write("*LOAD\nLOADS\n")
        out.write(f'*END')

def main():

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate input deck for a beam model.")
    parser.add_argument("config_file", help="Path to the configuration file.")
    parser.add_argument("--output_file", help="Path to the output file.", default=None)
    args = parser.parse_args()

    if args.output_file is None:
        args.output_file = args.config_file.replace(".yaml", ".inp")

    config = parse_config(args.config_file)
    write_input_deck(config, args.output_file)

if __name__ == "__main__":
    main()
