import yaml
import argparse
from ..geometry import Geometry


def parse_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)


def generate_geom(sections):
    x = 0.0
    geom = Geometry()

    left_id = geom.add_node(x=x)
    geom.add_node_to_set("LEFT", left_id)

    for section in sections:
        length = section['length']
        n_elems = section['elements']
        name = section['name']

        dx = length / n_elems

        for i in range(n_elems):
            x += dx
            node_id = geom.add_node(x=x)
            if i == n_elems - 1 and section == sections[-1]:
                geom.add_node_to_set("RIGHT", node_id)

            elem_id = geom.add_element(element_type='B33', node_ids=[node_id - 1, node_id])
            geom.add_element_to_set(name, elem_id)

    return geom


def add_optional_cube(geom, size=1e4):
    base_x = -size / 2
    base_y = -size / 2
    base_z = -size / 2
    dx, dy, dz = size / 2, size, size

    n = lambda dx=0, dy=0, dz=0: geom.add_node(x=base_x + dx, y=base_y + dy, z=base_z + dz)
    cube_nodes = [n(0, 0, 0), n(dx, 0, 0), n(dx, dy, 0), n(0, dy, 0),
                  n(0, 0, dz), n(dx, 0, dz), n(dx, dy, dz), n(0, dy, dz)]

    elem_id = geom.add_element(element_type="C3D8", node_ids=cube_nodes)
    geom.add_element_to_set("CUBE", elem_id)

    return cube_nodes

def write_input_deck(config, output_file="output.inp"):
    geom = generate_geom(config['sections'])

    if config.get("show_cube", False):
        cube_nodes = add_optional_cube(geom)
        geom.add_material("CUBEMAT", 1, 0.3, 1)
        geom.add_section_solid("CUBEMAT", "CUBE")
        geom.add_node_set("CUBE_FIX")
        for nid in cube_nodes:
            geom.add_node_to_set("CUBE_FIX", nid)
        geom.add_supp("CUBE_FIX", set="CUBE_FIX", fx=0, fy=0, fz=0)

    # Add material
    mat = config['material']
    geom.add_material("MAT", mat['E'], mat.get('nu', 0), mat.get('rho', 0))

    # Add all sections
    for section in config['sections']:
        name = section['name']
        geom.add_profile(
            name=name,
            area=section['area'],
            Iy=section['Iy'],
            Iz=section['Iz'],
            It=section['It']
        )
        geom.add_section_beam(
            elset=name,
            material="MAT",
            profile=name,
            orientation=section.get("orientation", (0, 1, 0))
        )

    # Add point masses
    if 'point_masses' in config:
        for entry in config['point_masses']:
            elset = entry['set']
            name = entry.get('name', elset)  # use set name if no explicit name

            if elset not in geom.node_sets:
                raise ValueError(f"Point mass set '{elset}' not found in node sets.")

            # Create point element(s) for each node in the set
            node_ids = geom.node_sets[elset]
            new_set_name = f"PM_{name}"

            geom.add_element_set(new_set_name)
            for nid in node_ids:
                elem_id = geom.add_element(element_type='P', node_ids=[nid])
                geom.add_element_to_set(new_set_name, elem_id)

            # Add point mass section on this element set
            mass = entry.get('mass')
            if mass is not None and not isinstance(mass, (int, float)):
                raise ValueError("mass must be a scalar (float or int).")

            inertia = entry.get('inertia')
            spring = entry.get('spring')
            rotary_spring = entry.get('rotary_spring')

            geom.add_section_pointmass(
                elset=new_set_name,
                mass=mass,
                inertia=inertia,
                spring=spring,
                rotary_spring=rotary_spring
            )



    # Add supports
    if 'bc' in config:
        for setname, dofs in config['bc'].items():
            if dofs is None or dofs == []:
                continue
            if setname not in geom.node_sets:
                raise ValueError(f"Support set '{setname}' not found in node sets.")
            geom.add_supp(setname, set=setname, **parse_dofs(dofs))

    # Add loads
    if 'loads' in config:
        for setname, dofs in config['loads'].items():
            if dofs is None or dofs == []:
                continue
            if setname not in geom.node_sets:
                raise ValueError(f"Load set '{setname}' not found in node sets.")
            geom.add_load(setname, set=setname, **parse_dofs(dofs))

    # Add step (based on scenario)
    scenario = config.get("scenario", "LINEAR STATIC").upper()
    if scenario == "LINEAR STATIC":
        geom.add_step(
            type="LINEAR STATIC",
            name="Step-1",
            loads=list(geom.loads.keys()),
            supps=list(geom.supps.keys())
        )
    elif scenario == "EIGENFREQ":
        geom.add_step(
            type="EIGENFREQ",
            name="Step-1",
            supps=list(geom.supps.keys()),
            numeigenvalues=config.get("numeigenvalues", 10)
        )
    else:
        raise ValueError(f"Unsupported scenario type: {scenario}")

    geom.write_input_deck(output_file)
    print(f"Written input deck to {output_file}")


def parse_dofs(lst):
    keys = ['fx', 'fy', 'fz', 'mx', 'my', 'mz']
    return {k: v for k, v in zip(keys, lst) if v is not None}


def generate_example_yaml(path="example_beam.yaml"):
    example = {
        'scenario': "LINEAR STATIC",
        'show_cube': True,
        'material': {
            'E': 210000,
            'nu': 0.3,
            'rho': 8.5e-9
        },
        'sections': [{
            'name': f'SECTION{i + 1}',
            'length': 13800.0,
            'elements': 10,
            'area': 400000 + i * 10000,
            'Ix': 1.0e12 + i * 1e11,
            'Iy': 1.0e12 + i * 1e11,
            'It': 2.0e12 + i * 2e11
        } for i in range(10)],
        'bc': {
            'LEFT': [0, 0, 0, null, 0, 0],
            'MID': []
        },
        'loads': {
            'RIGHT': [0, -1000, 0, 0, 0, 0]
        }
    }
    with open(path, 'w') as f:
        yaml.dump(example, f, sort_keys=False)
    print(f"Sample input deck written to: {path}")


def main():
    parser = argparse.ArgumentParser(description="Generate input deck for a beam model.")
    parser.add_argument("config_file", nargs="?", help="Path to the YAML configuration file.")
    parser.add_argument("--output_file", help="Path to the output file.", default=None)
    parser.add_argument("--example", action="store_true", help="Write a sample input YAML.")
    args = parser.parse_args()

    if args.example:
        generate_example_yaml()
        return

    if not args.config_file:
        print("Error: No config file provided.")
        return

    if args.output_file is None:
        args.output_file = args.config_file.replace(".yaml", ".inp")

    config = parse_config(args.config_file)
    write_input_deck(config, args.output_file)


if __name__ == "__main__":
    main()
