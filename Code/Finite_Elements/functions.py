from classes import Material
from classes import Property
from classes import Load
from classes import AppliedLoad
from classes import Constrain
from classes import Point
from classes import Line
from classes import Node
from classes import Element

from scipy import linalg
import numpy as np
import time


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


def get_normalized_eigs(stiff_matrix, mass_matrix):
    # This function finds the eigenvalues and eigenvectors of the problem,
    # sorts the eigenvectors in order of increasing eigenvalue and
    # normalizes the eigenvectors by the mass matrix

    # Solves the eigenvalue problem
    [eig_values, eig_vectors] = linalg.eig(stiff_matrix, mass_matrix)

    # Sorts the eingenvectors
    idx = eig_values.argsort()[::1]
    eig_values = eig_values[idx]
    eig_vectors = eig_vectors[:, idx]

    # Normalizes the eigenvectors
    factor_matrix = np.dot(np.dot(eig_vectors.transpose(), mass_matrix),
                           eig_vectors)
    factor_vector = factor_matrix.diagonal()
    factor_matrix = np.diag(factor_vector)  # eliminates numerical errors
    norm_eig_vectors = np.dot(eig_vectors, linalg.inv(np.sqrt(factor_matrix)))

    return eig_values, norm_eig_vectors


def damp_rayleigh(mass_matrix, stiff_matrix, freq_vector):
    # This function calculates the damping matrix using Rayleigh's method
    # user input is required

    mode_0 = int(input("Enter first mode number: "))
    damp_0 = float(input("Enter mode damping: "))

    mode_1 = int(input("Enter second mode number: "))
    damp_1 = float(input("Enter mode damping: "))

    freq_0 = freq_vector[mode_0]
    freq_1 = freq_vector[mode_1]

    alpha_beta = ((2 * freq_0 * freq_1) / (freq_1**2 - freq_0**2)) * np.array([[freq_1, -freq_0], [-1 / freq_1, 1 / freq_0]]) * np.array([[damp_0], [damp_1]])

    alpha = alpha_beta[0][0]
    beta = alpha_beta[1][0]

    damp_matrix = alpha * mass_matrix + beta * stiff_matrix

    damp_vector = np.zeros(len(freq_vector))

    for i in range(len(freq_vector)):
        damp_vector[i] = 0.5 * (alpha / freq_vector[i] + beta * freq_vector[i])

    return damp_matrix, damp_vector


def damp_bismarck(mass_matrix, eig_vectors, freq_vector):
    # This function calculates the damping matrix using Bismark's experimental
    # method user input is required

    damp_vector = np.zeros(len(freq_vector))
    damp_matrix = np.zeros((len(freq_vector), len(freq_vector)))

    for i in range(len(freq_vector)):
        damp_vector[i] = float(input("Enter damping for mode " + str(i) + " :"))
        eig_vector = eig_vectors[:, i]
        eig_vector_tp = eig_vector.transpose()
        m_ii = np.dot(np.dot(eig_vector_tp, mass_matrix), eig_vector)
        beta = 2 * damp_vector[i] * freq_vector[i] * m_ii
        damp_matrix[i][i] = 1 / beta

    return damp_matrix, damp_vector


def get_damp_matrix(mass_matrix, stiff_matrix, eig_vectors, freq_vector):

    # Calculates Damping Matrix
    damp_yn = input("\nCalculate Damping Matrix (Y/N): ")
    damp_yn = damp_yn.upper()

    print("\nCalculating Damp Matrix ...")
    star_time = time.clock()

    if damp_yn == 'Y':
        damp_method = input("Method (R for Rayleigh, B for Bismark)?: ")
        damp_method = damp_method.upper()

        if damp_method == 'R':
            damp_matrix, damp_vector = damp_rayleigh(mass_matrix,
                                                     stiff_matrix,
                                                     freq_vector)
        elif damp_method == 'B':
            damp_matrix, damp_vector = damp_bismarck(mass_matrix,
                                                     eig_vectors,
                                                     freq_vector)
        else:
            print("Method was not recognized, continuing with zero damping")
            damp_matrix = np.zeros((len(freq_vector), len(freq_vector)))
    else:
        damp_matrix = np.zeros((len(freq_vector), len(freq_vector)))

    print("Time used: " + str(round(time.clock() - star_time, 4)) + "s")

    return damp_matrix
