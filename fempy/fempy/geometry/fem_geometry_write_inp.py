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
            file.write(f"*DENSITY\n")
            file.write(f"{material['density']}\n")

        # Write profiles
        if hasattr(self, "profiles"):
            for profile in self.profiles:
                spacing()
                file.write(f"*PROFILE, NAME={profile['name']}\n")
                file.write(f"{profile['area']}, {profile['Iy']}, {profile['Iz']}, {profile['It']}\n")

        # Write sections
        for section in self.sections:
            spacing()
            stype = section['type']
            if stype == 'SOLID':
                file.write(f"*SOLID SECTION, ELSET={section['elset']}, MATERIAL={section['material']}\n")
            elif stype == 'SHELL':
                file.write(f"*SHELL SECTION, ELSET={section['elset']}, MATERIAL={section['material']}\n")
                file.write(f"{section['thickness']}\n")
            elif stype == 'POINTMASS':
                file.write(f"*POINT MASS SECTION, ELSET={section['elset']}\n")
                file.write(f"{section['mass'] if section['mass'] is not None else 0.0}\n")
                file.write(f"{', '.join(map(str, section['inertia'] or [0.0]*3))}\n")
                file.write(f"{', '.join(map(str, section['spring'] or [0.0]*3))}\n")
                file.write(f"{', '.join(map(str, section['rotary_spring'] or [0.0]*3))}\n")
            elif stype == 'BEAM':
                file.write(f"*BEAM SECTION, ELSET={section['elset']}, MATERIAL={section['material']}, PROFILE={section['profile']}\n")
                orientation = section.get('orientation', (0, 1, 0))
                file.write(f"{orientation[0]}, {orientation[1]}, {orientation[2]}\n")

        spacing()
        if format == 'femaster':
            # Write boundary conditions and loads before steps
            spacing()
            for name in self.supps:
                bc = self.supps[name]
                set_name = bc['set']
                bc_data = bc['data']
                file.write(f"*SUPPORT, SUPPORT_COLLECTOR={name}\n")
                file.write(f"{set_name}, {', '.join([str(k) if k is not None else ' ' for k in bc_data])}\n")

            spacing()
            for name in self.loads:
                load = self.loads[name]
                set_name = load['set']
                load_data = load['data']
                file.write(f"*CLOAD, LOAD_COLLECTOR={name}\n")
                file.write(f"{set_name}, {', '.join([str(k) if k is not None else '0.0' for k in load_data])}\n")

            # Write steps
            for step in self.steps:
                spacing()
                step_type = step.get('type', 'LINEAR STATIC').upper()
                file.write(f"*LOADCASE, TYPE={step_type}\n")

                if step_type == "LINEAR STATIC":
                    file.write(f"*SUPPORT\n")
                    for support_name in step['supps']:
                        file.write(f"{support_name}\n")
                    file.write(f"*LOAD\n")
                    for load_name in step['loads']:
                        file.write(f"{load_name}\n")

                elif step_type == "EIGENFREQ":
                    file.write(f"*SUPPORT\n")
                    for support_name in step['supps']:
                        file.write(f"{support_name}\n")
                    file.write(f"*NUMEIGENVALUES\n{step['numeigenvalues']}\n")

                file.write("*END")

        elif format == 'abaqus':
            # Write steps with boundary conditions and loads
            for step in self.steps:
                spacing()
                step_type = step.get("type", "LINEAR STATIC").upper()
                file.write(f"*STEP, NAME={step['name']}\n")

                if step_type == "LINEAR STATIC":
                    file.write(f"*STATIC\n")
                elif step_type == "EIGENFREQ":
                    file.write(f"*FREQUENCY\n{step['numeigenvalues']}\n")
                else:
                    raise ValueError(f"Unsupported step type for Abaqus export: {step_type}")

                spacing(1)
                if "loads" in step and step["loads"]:
                    file.write(f"*CLOAD\n")
                    for load_name in step['loads']:
                        load = self.loads[load_name]
                        set_name = load['set']
                        load_data = load['data']
                        for idx, dat in enumerate(load_data):
                            if dat is not None and dat != 0:
                                file.write(f"{set_name}, {idx+1}, {dat}\n")

                spacing(1)
                if "supps" in step and step["supps"]:
                    file.write(f"*BOUNDARY\n")
                    for support_name in step['supps']:
                        support = self.supps[support_name]
                        set_name = support['set']
                        bc_data = support['data']
                        for idx, dat in enumerate(bc_data):
                            if dat is not None:
                                file.write(f"{set_name}, {idx+1}, {idx+1}, {float(dat)}\n")

                file.write(f"*END STEP\n")

        else:
            raise ValueError(f"Unsupported format {format}. Supported formats are 'femaster' and 'abaqus'.")
