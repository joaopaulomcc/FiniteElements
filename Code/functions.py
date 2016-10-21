from classes import Material
from classes import Property
from classes import Load
from classes import AppliedLoad
from classes import Constrain
from classes import Point
from classes import Line
from classes import Node
from classes import Element


# This file contains the functions used

def create_mesh(lines_array, points_array):
    # This function creates the nodes, elements and the applied loads
    # used to create the K, M and C matrices

    node_count = 0
    elem_count = 0

    counter = 0  # Counts how many applied loads were created from a load

    for i, point in enumerate(points_array):
        # Creates a node for every point defined
        Node(point.x, point.y, 0)
        node_count += 1

        if point.load.kind == "CONCENTRATED":
            # Creates an applied load if there is a load applied to the point
            AppliedLoad(point.load, counter, [point.number])
            counter += 1

    for line in lines_array:
        counter = 0
        local_nodes = line.mesh()  # array with coord of the points in the line

        if len(local_nodes) > 2:
            # Create first node of the line
            Node(local_nodes[1][0], local_nodes[1][1], 0)
            node_count += 1

            # Create first element of the line
            Element(line.prop,
                    Node.instances[line.point_0.number],
                    Node.instances[node_count - 1])
            elem_count += 1

            if line.load.kind == "DISTRIBUTED":
                # Creates an applied load for every element created
                node_0 = Element.instances[elem_count - 1].node_0
                node_1 = Element.instances[elem_count - 1].node_1
                AppliedLoad(line.load, counter, [node_0, node_1])
                counter += 1

            for i in range(1, line.n - 1):
                # Creates the internal nodes of the line
                Node(local_nodes[i + 1][0], local_nodes[i + 1][1], 0)
                node_count += 1

                # Creates the internal elements of the line in a sequence
                node_0 = Element.instances[elem_count - 1].node_1
                node_1 = Node.instances[node_count - 1]
                Element(line.prop,
                        node_0,
                        node_1)
                elem_count += 1

                if line.load.kind == "DISTRIBUTED":
                    node_0 = Element.instances[elem_count - 1].node_0
                    node_1 = Element.instances[elem_count - 1].node_1
                    AppliedLoad(line.load, counter, [node_0, node_1])
                    counter += 1

            # Create last element of the line
            node_0 = Element.instances[elem_count - 1].node_1
            node_1 = Node.instances[line.point_1.number]
            Element(line.prop,
                    node_0,
                    node_1)
            elem_count += 1

            if line.load.kind == "DISTRIBUTED":
                node_0 = Element.instances[elem_count - 1].node_0
                node_1 = Element.instances[elem_count - 1].node_1
                AppliedLoad(line.load, counter, [node_0, node_1])
                counter += 1

        else:
            # If a line has only one element, this creates it
            node_0 = Node.instances[line.point_0.number]
            node_1 = Node.instances[line.point_1.number]
            Element(line.prop,
                    node_0,
                    node_1)
            elem_count += 1

            if line.load.kind == "DISTRIBUTED":
                node_0 = Element.instances[elem_count - 1].node_0
                node_1 = Element.instances[elem_count - 1].node_1
                AppliedLoad(line.load, counter, [node_0, node_1])
                counter += 1

    nodes_array = Node.instances
    elements_array = Element.instances
    app_loads_array = AppliedLoad.instances

    # print("DEBUG")
    # print("\nnodes_array")
    # print(nodes_array)
    # print("\nelements_array")
    # print(elements_array)
    # print("\napp_loads_array")
    # print(app_loads_array)

    return nodes_array, elements_array, app_loads_array


def local2global(matrix_local, matrix_global, map_vector):
    # Uses the elements from the local matrix to build the
    # global matrix

    for i in range(len(map_vector)):
        for j in range(len(map_vector)):
            i_g = int(map_vector[i])
            j_g = int(map_vector[j])

            # Relate nodes with Degrees of Freedom
            matrix_global[i_g][j_g] += matrix_local[i][j]
