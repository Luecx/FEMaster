import model

def generate_model(model, densities):
    output = []
    output.append("*NODE")
    for node_id, (x, y, z) in enumerate(model.nodes):
        output.append(f"{node_id}, {x}, {y}, {z}")  # Adding 1 because node_id starts from 0

    output.append("\n*ELEMENT, TYPE=C3D8")
    for elem_id, node_ids in enumerate(model.elements):
        output.append(f"{elem_id}, " + ", ".join(map(str, [nid for nid in node_ids])))  # Adding 1 because indexing starts from 0

    output.append(f"\n*MATERIAL, NAME=MAT1")
    output.append(f"*ELASTIC, TYPE=ISO")
    output.append(f"{model.youngs},{model.nu}")
    output.append(f"*SOLID SECTION, ELSET=EALL, MAT=MAT1")

    output.append("\n*SUPPORT, SUPPORT_COLLECTOR=END_SUPPORTS")
    for node_id, values in model.supports.items():
        output.append(f"{node_id}"
                      f", {' ' if model.supports[node_id][0] is None else model.supports[node_id][0]}"
                      f", {' ' if model.supports[node_id][1] is None else model.supports[node_id][1]}"
                      f", {' ' if model.supports[node_id][2] is None else model.supports[node_id][2]}")

    output.append("\n*CLOAD, LOAD_COLLECTOR=END_LOADS")
    for node_id, values in model.loads.items():
        output.append(f"{node_id}"
                      f", {' ' if model.loads[node_id][0] is None else model.loads[node_id][0]}"
                      f", {' ' if model.loads[node_id][1] is None else model.loads[node_id][1]}"
                      f", {' ' if model.loads[node_id][2] is None else model.loads[node_id][2]}")

    output.append(f"\n*LOADCASE, TYPE= LINEAR STATIC TOPO")
    output.append(f"*SUPPORT \nEND_SUPPORTS")
    output.append(f"*LOAD \nEND_LOADS")

    output.append("\n*DENSITY")
    for i, density in enumerate(densities):
        output.append(f"{i},{density}")

    output.append(f"\n*EXPONENT\n{model.exponent}")
    output.append(f"\n*SOLVER, TYPE={'DIRECT' if model.direct else 'INDIRECT'}, DEVICE={'CPU' if model.cpu else 'GPU'}\n*END")

    return '\n'.join(output)
