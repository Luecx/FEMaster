def write_input_deck(self, filename, format='femaster'):
    from .fem_geometry import Geometry


    with open(filename, 'w') as file:


        def spacing(n=3):
            file.write('**\n' * n)

        # Write nodes
        file.write("*NODE, NSET=NALL\n")
        for id, node in enumerate(self.nodes):
            if node is not None:
                file.write(f"{id},{node[0]},{node[1]},{node[2]}\n")

        # Write elements by type
        element_classes = {k: [] for k in Geometry.element_classes.keys()}
        for element in self.elements:
            if element is not None:
                element_classes[element.elem_type].append(element)

        for elem_type, elements in element_classes.items():
            if elements:
                # adjust S8 -> S8R and S6 -> STRI65
                if format == 'abaqus':
                    if elem_type == 'S8':
                        elem_type = 'S8R'
                    elif elem_type == 'S6':
                        elem_type = 'STRI65'

                file.write(f"*ELEMENT, TYPE={elem_type}, ELSET=EALL\n")
                for element in elements:
                    nodes = ','.join(map(str, element.node_ids))
                    file.write(f"{element.element_id},{nodes}\n")


        # Write node sets
        for name, ids in self.node_sets.items():
            if name != "NALL":
                file.write(f"*NSET, NSET={name}\n")
                file.write(',\n'.join(map(str, ids)) + '\n')

        # Write element sets
        for name, ids in self.elem_sets.items():
            if name != "EALL":
                file.write(f"*ELSET, ELSET={name}\n")
                file.write(',\n'.join(map(str, ids)) + '\n')

        # Write material properties
        spacing()
        for name, material in self.materials.items():
            file.write(f"*MATERIAL, NAME={name}\n")
            file.write(f"*ELASTIC, TYPE=ISOTROPIC\n")
            file.write(f"{material['young']}, {material['poisson']}\n")

        # Write sections
        for section in self.sections:
            spacing()
            type = section['type']
            if type == 'SOLID':
                file.write(f"*SOLID SECTION, ELSET={section['elset']}, MATERIAL={section['material']}\n")
            if type == 'SHELL':
                file.write(f"*SHELL SECTION, ELSET={section['elset']}, MATERIAL={section['material']}\n")
                file.write(f"{section['thickness']}\n")

        spacing()
        if format == 'femaster':
            # Write boundary conditions and loads before steps
            spacing()
            for name in self.supps:
                bc = self.supps[name]
                set_name = bc['set']
                bc_data  = bc['data']
                file.write(f"*SUPPORT, SUPPORT_COLLECTOR={name}\n")
                file.write(f"{set_name}, {', '.join([str(k) if k is not None else ' ' for k in bc_data])}\n")

            spacing()
            for name in self.loads:
                load = self.loads[name]
                set_name  = load['set']
                load_data = load['data']
                file.write(f"*CLOAD, LOAD_COLLECTOR={name}\n")
                file.write(f"{set_name}, {', '.join([str(k) if k is not None else ' ' for k in bc_data])}\n")

            # Write steps
            for step in self.steps:
                spacing()
                file.write(f"*LOADCASE, TYPE=LINEAR STATIC\n")
                file.write(f"*SUPPORT\n")
                for support_name in step['supports']:
                    file.write(f"{support_name}\n")
                file.write(f"*LOAD\n")
                for load_name in step['loads']:
                    file.write(f"{load_name}\n")
                file.write("*END")


        elif format == 'abaqus':
            # write assembly
            spacing()
            # file.write("*ASSEMBLY\n")
            # file.write("*INSTANCE, NAME=PART-1, PART=PART-1\n")
            # file.write("*END INSTANCE\n")

            # Write steps with boundary conditions and loads
            for step in self.steps:
                spacing()
                file.write(f"*STEP, NAME={step['name']}\n")
                file.write(f"*STATIC\n")

                # Write loads for this step
                spacing(1)
                file.write(f"*CLOAD\n")
                for load_name in step['loads']:
                    load = self.loads[load_name]
                    set_name  = load['set']
                    load_data = load['data']
                    # file.write(f"{set_name}, {', '.join(map(str, load_data))}\n")
                    for idx, dat in enumerate(load_data):
                        if dat is not None and dat != 0:
                            file.write(f"{set_name}, {idx+1}, {dat}\n")

                # Write supports for this step
                spacing(1)
                file.write(f"*BOUNDARY\n")
                for support_name in step['supports']:
                    support = self.supps[support_name]
                    set_name = support['set']
                    bc_data = support['data']
                    for idx, dat in enumerate(bc_data):
                        if dat is not None:
                            file.write(f"{set_name}, {idx+1}, {idx+1}, {float(dat)}\n")

                file.write(f"*END STEP\n")
        else:
            raise ValueError(f"Unsupported format {format}. Supported formats are 'femaster' and 'abaqus'.")
