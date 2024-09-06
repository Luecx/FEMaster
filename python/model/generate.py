import copy

import numpy as np

from spring_sim import SpringSimulation

from geometry import Geometry
from solution import Solution
from scipy.optimize import curve_fit


def generate_beam(length, height, elem_length, elem_height):
    geom = Geometry()
    n_elem_x = int(length // elem_length)
    n_elem_y = int(height // elem_height)

    geom.add_node_set('LEFT')
    geom.add_node_set('RIGHT')
    geom.add_node_set('BOTTOM')
    geom.add_node_set('TOP')
    geom.add_node_set('CENTER')

    # Generate nodes
    node_id = 1
    for i in range(n_elem_x + 1):
        for j in range(n_elem_y + 1):
            x = i * elem_length
            y = j * elem_height - height / 2
            geom.add_node(node_id, x, y)

            # Add to relevant sets
            if i == 0:
                geom.add_node_to_set('LEFT', node_id)
            elif i == n_elem_x:
                geom.add_node_to_set('RIGHT', node_id)
            if j == 0:
                geom.add_node_to_set('BOTTOM', node_id)
            elif j == n_elem_y:
                geom.add_node_to_set('TOP', node_id)
            if j == n_elem_y // 2:
                geom.add_node_to_set('CENTER', node_id)

            node_id += 1

    # Generate C2D4 elements
    elem_id = 1
    for i in range(n_elem_x):
        for j in range(n_elem_y):
            bottom_left = i * (n_elem_y + 1) + j
            bottom_right = (i + 1) * (n_elem_y + 1) + j
            top_left = i * (n_elem_y + 1) + j + 1
            top_right = (i + 1) * (n_elem_y + 1) + j + 1

            # Add as C2D4
            node_ids = [bottom_left + 1, bottom_right + 1, top_right + 1, top_left + 1]
            geom.add_element(elem_id, 'C2D4', node_ids)
            elem_id += 1

    return geom


def generate_arc(width=60, offset=100, R_inner=30, R_outer=5, angle=140, num_elements_width=10,
                           num_elements_inner_arc=30, num_elements_outer_arc=10, num_elements_inner_straight=10,
                           num_elements_outer_straight=10):
    # make sure angle is < 180 degrees
    if angle >= 180:
        raise ValueError("Angle must be less than 180 degrees")


    # 8 arrays for the 8 segments
    end_a = []
    end_b = []
    inner_straight_a = []
    inner_straight_b = []
    inner_arc = []
    outer_straight_a = []
    outer_straight_b = []
    outer_arc = []

    # computer the center of the outer arc
    dir_angle       = np.array([np.cos(np.radians(angle    )), np.sin(np.radians(angle    ))])
    dir_angle_90    = np.array([-dir_angle[1], dir_angle[0]])
    dir_angle_half  = np.array([np.cos(np.radians(angle / 2)), np.sin(np.radians(angle / 2))])
    dx = (width + R_inner - R_outer)
    dy = dx * dir_angle_half[1] / dir_angle_half[0]

    outer_center = [dx, dy]
    outer_arc_p1 = [dx + R_outer * dir_angle[0], dy + R_outer * dir_angle[1]]
    outer_arc_p2 = [dx + R_outer               , dy]

    # Generate boundary points for the inner arc (counter clock wise) at radius R_inner
    for i in range(num_elements_inner_arc + 1):
        theta = i * angle / num_elements_inner_arc
        x = R_inner * np.cos(np.radians(theta))
        y = R_inner * np.sin(np.radians(theta))
        inner_arc.append([x, y])
    inner_arc = np.array(inner_arc)

    # Generate boundary points for the inner straight_b
    inner_straight_b = np.linspace(inner_arc[-1],
                                   inner_arc[-1] + (dir_angle_90 * offset), num_elements_inner_straight + 1)

    # generate points for the end_b perpendicular by turning clockwise
    end_b = np.linspace(inner_straight_b[-1], inner_straight_b[-1] + (dir_angle * width), num_elements_width + 1)

    # turn right again and fill to outer_arc_p1
    outer_straight_b = np.linspace(end_b[-1],
                                   [outer_arc_p1[0], outer_arc_p1[1]], num_elements_outer_straight + 1)

    # boudnary points for the outer arc (clock wise) at radius R_outer with base point end_b[-1]
    for i in range(num_elements_outer_arc + 1):
        theta = angle - i * angle / num_elements_outer_arc
        x = R_outer * np.cos(np.radians(theta)) + outer_center[0]
        y = R_outer * np.sin(np.radians(theta)) + outer_center[1]
        outer_arc.append([x, y])
    outer_arc = np.array(outer_arc)

    # outer straight a just going vertically down to the x-axis (width, -offset)
    outer_straight_a = np.linspace([outer_arc[-1][0], outer_arc[-1][1]], [width + R_inner, -offset], num_elements_outer_straight + 1)

    # base point for the inner straight a is the end of the outer straight a
    end_a = np.linspace([outer_straight_a[-1][0], outer_straight_a[-1][1]], [R_inner, -offset], num_elements_width + 1)

    # back to origin
    inner_straight_a = np.linspace([end_a[-1][0], end_a[-1][1]], [R_inner, 0], num_elements_inner_straight + 1)

    # Combine all boundary points to form a closed loop
    boundary_points = np.vstack([inner_arc, inner_straight_b, end_b, outer_straight_b, outer_arc, outer_straight_a, end_a, inner_straight_a])

    # remove points that are too close to each other
    for i in range(len(boundary_points) - 1, 0, -1):
        if np.linalg.norm(boundary_points[i] - boundary_points[i - 1]) < 1e-6:
            boundary_points = np.delete(boundary_points, i, axis=0)
    # create a triangle mesh
    geom = Geometry.mesh_interior(boundary_points, force_quads=True)
    return geom

if __name__ == '__main__':

    def find_bottom_nodes(geom):
        nodes = np.array(geom.nodes)
        # find all node ids that are near the bottom -> find min y value and take all within 1e-6
        min_y = np.min(nodes[:,1])
        return np.where(np.abs(nodes[:,1] - min_y) < 1e-6)

    def find_left_nodes(geom):
        nodes = np.array(geom.nodes)
        # find all node ids that are near the left -> find min x value and take all within 1e-6
        min_x = np.min(nodes[:,0])
        return np.where(np.abs(nodes[:,0] - min_x) < 1e-6)

    def write_file(file_name, geom):
        geom = geom.to_second_order()
        geom = geom.extrude(1,spacing=-10)
        geom.write_input_deck(file_name)

        with open(file_name, "a") as f:
            f.write("""
*MATERIAL, NAME=MAT1
*ELASTIC, TYPE=ISO
210000,0.3
*SOLID SECTION, ELSET=EALL, MAT=MAT1

*SUPPORT, SUPPORT_COLLECTOR=SUPPS
LEFT, 0, 0, 0
*CLOAD, LOAD_COLLECTOR=LOADS
RIGHT, 0, -1, 0

*LOADCASE, TYPE= LINEAR STATIC
*SUPPORT
SUPPS
*LOAD
LOADS
*SOLVER, METHOD=DIRECT, DEVICE=CPU
*END""")


    geom = generate_arc(angle=130,
                                  num_elements_inner_arc=10,
                                  num_elements_width=4,
                                  num_elements_outer_arc=5,
                                  offset=100,
                                  num_elements_inner_straight=4,
                                  num_elements_outer_straight=7,
                                  R_inner=20,
                                  R_outer=30)
    # geom = generate_beam(10, 1, 0.1, 0.1)

    write_file("./python/model/arc.inp", copy.deepcopy(geom))

    exit(0)


    for angle in [120,50,75,90,120,150]:
        for offset in [100, 10,20,30,40,50]:
            for R_inner in [10,20,30,40,50]:

                # generate geometry
                error = True
                while error:
                    try:
                        geom = generate_irregular_arc(angle=angle,
                                      num_elements_inner_arc=20,
                                      num_elements_width=15,
                                      num_elements_outer_arc=6,
                                      offset=offset,
                                      num_elements_inner_straight=3,
                                      num_elements_outer_straight=10,
                                      R_inner=R_inner,
                                      R_outer=30)
                        error = False
                    except ValueError:
                        error = True


                file_name_0 = "test_{}_{}_{}_0.inp".format(angle,offset,R_inner)
                file_name_a = "test_{}_{}_{}_a.inp".format(angle,offset,R_inner)
                file_name_b = "test_{}_{}_{}_b.inp".format(angle,offset,R_inner)

                write_file(file_name_0, copy.deepcopy(geom))
                # subdivide once
                geom = geom.subdivide(1)
                geom.equalize_node_spacing(100, 0.01)

                write_file(file_name_a, copy.deepcopy(geom))
                # subdivide once
                geom = geom.subdivide(1, only_quads=True)
                geom.write_input_deck(file_name_b)
                write_file(file_name_b, copy.deepcopy(geom))

                exit(0)


    # geom = generate_irregular_arc(angle=75,
    #                               num_elements_inner_arc=20,
    #                               num_elements_width=15,
    #                               num_elements_outer_arc=6,
    #                               offset=10,
    #                               num_elements_inner_straight=3,
    #                               num_elements_outer_straight=10,
    #                               R_inner=30,
    #                               R_outer=30)

    # subdivide once
    geom = geom.subdivide(1)
    geom.plot_2d()
    geom.equalize_node_spacing(100, 0.01)
    geom.plot_2d()
    geom.write_input_deck("test1.inp")
    geom = geom.subdivide(1)
    geom.plot_2d()
    geom.write_input_deck("test2.inp")
