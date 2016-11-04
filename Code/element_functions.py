# This file contains the functions to calculate the element matrices

import numpy as np
import math
import sys


def local_stiff_matrix(element, nodes_array):

    element_kind = element.prop.kind

    x_0 = element.node_0.x
    y_0 = element.node_0.y
    x_1 = element.node_1.x
    y_1 = element.node_1.y

    node_0_index = nodes_array.index(element.node_0)
    node_1_index = nodes_array.index(element.node_1)

    l = math.sqrt((x_1 - x_0) ** 2 + (y_1 - y_0) ** 2)  # Element's length

    c = (x_1 - x_0) / l  # Element's cos
    s = (y_1 - y_0) / l  # Element's sin

    if element_kind == "BAR":
        # Calculates the stiffness matrix for a "BAR" element

        map_vector = np.zeros(4)
        map_vector[0] = node_0_index * 3
        map_vector[1] = node_0_index * 3 + 1
        map_vector[2] = node_1_index * 3
        map_vector[3] = node_1_index * 3 + 1

        e = element.prop.material.young
        a = element.prop.area

        # Element's stiffness matrix before rotation
        k_mat = np.array([[1, -1], [-1, 1]]) * (e * a / l)

        # Element's rotation matrix
        t_mat = np.array([[c, s, 0, 0], [0, 0, c, s]])

        stiff_local = np.dot(np.dot(np.transpose(t_mat), k_mat), t_mat)

        return stiff_local, map_vector

    elif element_kind == "BEAM4":
        # This function receives the information of a BEAM4
        # element and returns, it's stiffness matrix

        map_vector = np.zeros(4)
        map_vector[0] = node_0_index * 3 + 1
        map_vector[1] = node_0_index * 3 + 2
        map_vector[2] = node_1_index * 3 + 1
        map_vector[3] = node_1_index * 3 + 2

        e = element.prop.material.young  # Element's young modulus
        i = element.prop.inertia  # Element's moment of inertia

        # Element's stiffness matrix before rotation
        k = np.zeros((4, 4))
        k[0][0] = 12
        k[0][1] = 6 * l
        k[0][2] = -12
        k[0][3] = 6 * l
        k[1][0] = 6 * l
        k[1][1] = 4 * l ** 2
        k[1][2] = -6 * l
        k[1][3] = 2 * l ** 2
        k[2][0] = -12
        k[2][1] = -6 * l
        k[2][2] = 12
        k[2][3] = -6 * l
        k[3][0] = 6 * l
        k[3][1] = 2 * l ** 2
        k[3][2] = -6 * l
        k[3][3] = 4 * l ** 2

        k *= ((e * i) / (l ** 3))

        # Element's rotation matrix
        t = np.zeros((4, 4))
        t[0][0] = c
        t[0][1] = s
        t[1][0] = -s
        t[1][1] = c
        t[2][2] = c
        t[2][3] = s
        t[3][2] = -s
        t[3][3] = c

        stiff_local = np.dot(np.dot(np.transpose(t), k), t)

        return stiff_local, map_vector

    elif element_kind == "BEAM6":
        # This function receives the information of a BEAM6
        # element and returns, it's stiffness matrix

        map_vector = np.zeros(6)
        map_vector[0] = node_0_index * 3
        map_vector[1] = node_0_index * 3 + 1
        map_vector[2] = node_0_index * 3 + 2
        map_vector[3] = node_1_index * 3
        map_vector[4] = node_1_index * 3 + 1
        map_vector[5] = node_1_index * 3 + 2

        e = element.prop.material.young  # Element's young modulus
        a = element.prop.area
        i = element.prop.inertia  # Element's moment of inertia

        # Element's stiffness matrix before rotation
        k = np.zeros((6, 6))
        k[0][0] = e * a / l
        k[1][1] = 12 * e * i / (l ** 3)
        k[2][2] = 4 * e * i / l
        k[3][3] = e * a / l
        k[4][4] = 12 * e * i / (l ** 3)
        k[5][5] = 4 * e * i / l
        k[0][3] = -e * a / l
        k[3][0] = k[0][3]
        k[1][2] = 6 * e * i / (l ** 2)
        k[2][1] = k[1][2]
        k[1][4] = -12 * e * i / (l ** 3)
        k[4][1] = k[1][4]
        k[1][5] = 6 * e * i / (l ** 2)
        k[5][1] = k[1][5]
        k[2][4] = -6 * e * i / (l ** 2)
        k[4][2] = k[2][4]
        k[2][5] = 2 * e * i / l
        k[5][2] = k[2][5]
        k[5][4] = -6 * e * i / (l ** 2)
        k[4][5] = k[5][4]

        # Element's rotation matrix
        t = np.zeros((6, 6))
        t[0][0] = c
        t[0][1] = s
        t[1][0] = -s
        t[1][1] = c
        t[2][2] = 1
        t[3][3] = c
        t[3][4] = s
        t[4][3] = -s
        t[4][4] = c
        t[5][5] = 1

        stiff_local = np.dot(np.dot(np.transpose(t), k), t)

        return stiff_local, map_vector


def local_mass_matrix(element):
    element_kind = element.prop.kind

    x_0 = element.node_0.x
    y_0 = element.node_0.y
    x_1 = element.node_1.x
    y_1 = element.node_1.y

    l = math.sqrt((x_1 - x_0) ** 2 + (y_1 - y_0) ** 2)  # Element's length

    c = (x_1 - x_0) / l  # Element's cos
    s = (y_1 - y_0) / l  # Element's sin

    rho = element.prop.material.density  # Element's young modulus
    a = element.prop.area

    if element_kind == "BEAM4":
        # Element's stiffness matrix before rotation
        m = np.zeros((4, 4))
        m[0][0] = 156
        m[0][1] = 22 * l
        m[0][2] = 54
        m[0][3] = -13 * l
        m[1][0] = 22 * l
        m[1][1] = 4 * l ** 2
        m[1][2] = 13 * l
        m[1][3] = -3 * l ** 2
        m[2][0] = 54
        m[2][1] = 13 * l
        m[2][2] = 156
        m[2][3] = -22 * l
        m[3][0] = -13 * l
        m[3][1] = -3 * l ** 2
        m[3][2] = -22 * l
        m[3][3] = 4 * l ** 2

        m *= (rho * a * l / 420)

        # Element's rotation matrix
        t = np.zeros((4, 4))
        t[0][0] = c
        t[0][1] = s
        t[1][0] = -s
        t[1][1] = c
        t[2][2] = c
        t[2][3] = s
        t[3][2] = -s
        t[3][3] = c

        m_matrix = np.dot(np.dot(np.transpose(t), m), t)

        return m_matrix

    elif element_kind == "BEAM6":

        m = np.zeros((6, 6))
        mm = rho * a * l / 420
        ma = rho * a * l / 6

        m[0][0] = 2 * ma
        m[1][1] = 156 * mm
        m[2][2] = 4 * l ** 2 * mm
        m[3][3] = 2 * ma
        m[4][4] = 156 * mm
        m[5][5] = 4 * l ** 2 * mm
        m[0][3] = ma
        m[3][0] = m[0][3]
        m[1][2] = 22 * l * mm
        m[2][1] = m[1][2]
        m[1][4] = 54 * mm
        m[4][1] = m[1][4]
        m[1][5] = -13 * l * mm
        m[5][1] = m[1][5]
        m[2][4] = 13 * l * mm
        m[4][2] = m[2][4]
        m[2][5] = -3 * l ** 2 * mm
        m[5][2] = m[2][5]
        m[5][4] = -22 * l * mm
        m[4][5] = m[5][4]

        # Element's rotation matrix
        t = np.zeros((6, 6))
        t[0][0] = c
        t[0][1] = s
        t[1][0] = -s
        t[1][1] = c
        t[2][2] = 1
        t[3][3] = c
        t[3][4] = s
        t[4][3] = -s
        t[4][4] = c
        t[5][5] = 1

        mass_local = np.dot(np.dot(np.transpose(t), m), t)

        return mass_local


def global_force_vector(app_load, force_global, nodes_array):
    # This function creates the global force vector, it is actually 3 vectors
    # containing the load magnitude, frequency and phase, they are used
    # in transient analysis

    if app_load.kind == "CONCENTRATED":

        node = int(app_load.nodes[0])
        fx = app_load.fx
        fy = app_load.fy
        mz = app_load.mz
        freq = app_load.freq
        phase = app_load.phase
        force_global[node * 3][0] += fx
        force_global[node * 3][1] += freq
        force_global[node * 3][2] += phase
        force_global[node * 3 + 1][0] += fy
        force_global[node * 3 + 1][1] += freq
        force_global[node * 3 + 1][2] += phase
        force_global[node * 3 + 2][0] += mz
        force_global[node * 3 + 2][1] += freq
        force_global[node * 3 + 2][2] += phase

    elif app_load.kind == "DISTRIBUTED":

        node_0 = app_load.nodes[0]
        node_1 = app_load.nodes[1]
        fx = app_load.fx
        fy = app_load.fy
        freq = app_load.freq
        phase = app_load.phase

        node_0_index = nodes_array.index(node_0)
        node_1_index = nodes_array.index(node_1)

        l = math.sqrt((node_1.x - node_0.x)**2 + (node_1.y - node_0.y)**2)

        # Node 0
        force_global[node_0_index * 3][0] += (fx * l / 2)
        force_global[node_0_index * 3][1] = freq
        force_global[node_0_index * 3][2] = phase

        force_global[node_0_index * 3 + 1][0] += (fy * l / 2)
        force_global[node_0_index * 3 + 1][1] = freq
        force_global[node_0_index * 3 + 1][2] = phase

        force_global[node_0_index * 3 + 2][0] += (fy * l**2 / 12)
        force_global[node_0_index * 3 + 2][1] = freq
        force_global[node_0_index * 3 + 2][2] = phase

        # Node 1
        force_global[node_1_index * 3][0] += (fx * l / 2)
        force_global[node_1_index * 3][1] = freq
        force_global[node_1_index * 3][2] = phase

        force_global[node_1_index * 3 + 1][0] += (fy * l / 2)
        force_global[node_1_index * 3 + 1][1] = freq
        force_global[node_1_index * 3 + 1][2] = phase

        force_global[node_1_index * 3 + 2][0] += (fy * l ** 2 / 12)
        force_global[node_1_index * 3 + 2][1] = freq
        force_global[node_1_index * 3 + 2][2] = phase

    return

def local_geo_stiff_matrix(element):

    element_kind = element.prop.kind

    x_0 = element.node_0.x
    y_0 = element.node_0.y
    x_1 = element.node_1.x
    y_1 = element.node_1.y

    l = math.sqrt((x_1 - x_0) ** 2 + (y_1 - y_0) ** 2)  # Element's length

    c = (x_1 - x_0) / l  # Element's cos
    s = (y_1 - y_0) / l  # Element's sin

    if element_kind == "BEAM4":
        # This function receives the information of a BEAM4
        # element and returns, it's geometric stiffness matrix

        e = element.prop.material.young  # Element's young modulus
        i = element.prop.inertia  # Element's moment of inertia

        # Element's stiffness matrix before rotation
        k = np.zeros((4, 4))
        k[0][0] = 36
        k[0][1] = 3 * l
        k[0][2] = -36
        k[0][3] = 3 * l
        k[1][0] = 3 * l
        k[1][1] = 4 * l ** 2
        k[1][2] = -3 * l
        k[1][3] = -l ** 2
        k[2][0] = -36
        k[2][1] = -3 * l
        k[2][2] = 36
        k[2][3] = -3 * l
        k[3][0] = 3 * l
        k[3][1] = -l ** 2
        k[3][2] = -3 * l
        k[3][3] = 4 * l ** 2

        k *= 1 / (30 * l)

        # Element's rotation matrix
        t = np.zeros((4, 4))
        t[0][0] = c
        t[0][1] = s
        t[1][0] = -s
        t[1][1] = c
        t[2][2] = c
        t[2][3] = s
        t[3][2] = -s
        t[3][3] = c

        geo_stiff_local = np.dot(np.dot(np.transpose(t), k), t)

        return geo_stiff_local

    else:
        print("Buckling analysis is implemented for BAR4 elements only\n"
              "Quiting program")
        sys.exit()
