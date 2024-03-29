import numpy as np
from scipy import linalg
import math
import matplotlib.pyplot as plt
import time
from Transient import transient

###############################################################################
###############################################################################
###############################################################################

# Classes
# Classes created in order to store the data retrieved from the
# input file.


class Material:
    # stores the material name and physical properties
    # density, young's modulus and poisson's ratio

    def __init__(self, name, density, young, poisson):
        self.name = name
        self.density = density
        self.young = young
        self.poisson = poisson


class Property:
    # A property hAe a kind (BAR or BEAM), a material
    # and some parameters, for the BAR the parameters
    # are only the area, for the BEAM the parameters
    # are the area and the moment of inertia

    def __init__(self, name, kind, material, parameters):
        self.name = name
        self.kind = kind
        self.material = material
        self.parameters = parameters


class Load:
    # A load have a name, a kind, CONCENTRATED or DISTRIBUTED,
    # for the CONCENTRATED load the parameters are the node
    # where the force is applied and FX, FY and MZ, for the
    # DISTRIBUTED load the parameters are the two nodes
    # between witch the load is applied and QX and QY.

    def __init__(self, name, kind, parameters):
        self.name = name
        self.kind = kind
        self.parameters = parameters


class Constrain:
    # A constrain hAe a name, the node where it is applied
    # and it's values in the X axis, Y axis and rotation around
    # the Z axis, if a constrain does not exist, it's value is
    # the stting "FREE".

    def __init__(self, name, node, X, Y, RZ):
        self.name = name
        self.node = node
        self.X = X
        self.Y = Y
        self.RZ = RZ


class JobData:
    # job_data is te object that contains all the information
    # necessary to make the analysis

    def __init__(self, title, analysis, materials, properties,
                 loads, constrains, elements, nodes):
        self.title = title
        self.analysis = analysis
        self.materials = materials
        self.properties = properties
        self.loads = loads
        self.constrains = constrains
        self.elements = elements
        self.nodes = nodes


###############################################################################
###############################################################################
###############################################################################

# Functions

def read_input_file(input_file):
    # This function reads an .txt file and returns a job_data object, this
    # object contains all the information necessary in order to perform
    # the analysis

    # materials, properties, loads, and constrains are stored in dictonaries,
    # the key is the name of the material, propertie, etc. Elements and nodes
    # are stored in arrays
    materials = {}
    properties = {}
    loads = {}
    constrains = {}
    elements = []
    nodes = []

    flag = None

    for line in input_file:
        # This loop will through the whole input file, line by line. When
        # a line starts with a "#" it means that there's no information
        # to be stored in that line. #TITLE, #ANALYSIS, and so on signals
        # what information is stored in the next lines
        # With the exception of material, all information is stored as strings.

        line_s = line.split()

        if line[0] == "#":
            if line_s[0] == "#":
                flag = None
            elif line_s[0] == "#TITLE":
                flag = "title"
            elif line_s[0] == "#ANALYSIS":
                flag = "analysis"
            elif line_s[0] == "#MATERIALS":
                flag = "materials"
            elif line_s[0] == "#PROPERTIES":
                flag = "properties"
            elif line_s[0] == "#LOADS":
                flag = "loads"
            elif line_s[0] == "#CONSTRAINS":
                flag = "constrains"
            elif line_s[0] == "#ELEMENTS":
                flag = "elements"
            elif line_s[0] == "#NODES":
                flag = "nodes"

        else:

            if flag == "title":
                title = line.strip()

            elif flag == "analysis":
                analysis = line.strip()

            elif flag == "materials":
                materials[line_s[0]] = Material(line_s[0],
                                                float(line_s[1]),
                                                float(line_s[2]),
                                                float(line_s[3]))

            elif flag == "properties":
                properties[line_s[0]] = Property(line_s[0],
                                                 line_s[1],
                                                 line_s[2],
                                                 line_s[3:])

            elif flag == "loads":
                loads[line_s[0]] = Load(line_s[0],
                                        line_s[1],
                                        line_s[2:])

            elif flag == "constrains":
                constrains[line_s[0]] = Constrain(line_s[0],
                                                  int(line_s[1]),
                                                  line_s[2],
                                                  line_s[3],
                                                  line_s[4])

            elif flag == "elements":
                elements.append(line_s)

            elif flag == "nodes":
                nodes.append(line_s)

    return JobData(title, analysis, materials,
                   properties, loads, constrains, elements, nodes)


def K_bar(element, job):
    # This function receIes the information of a BAR element and returns,
    # it's stiffness matrix

    element_propertie_name = element[1]
    element_propertie = job.properties[element_propertie_name]
    element_material_name = element_propertie.material
    element_material = job.materials[element_material_name]

    E = element_material.young  # Element's young modulus
    A = float(element_propertie.parameters[0])  # Element's area

    node_0 = job.nodes[int(element[2])]
    node_1 = job.nodes[int(element[3])]

    x_0 = float(node_0[1])
    y_0 = float(node_0[2])

    x_1 = float(node_1[1])
    y_1 = float(node_1[2])

    L = math.sqrt((x_1 - x_0)**2 + (y_1 - y_0)**2)  # Element's length

    c = (x_1 - x_0) / L  # Element's cos
    s = (y_1 - y_0) / L  # Element's sin

    # Element's stiffness matrix before rotation
    K = np.array([[1, -1], [-1, 1]]) * (E * A / L)

    # Element's rotation matrix
    T = np.array([[c, s, 0, 0], [0, 0, c, s]])

    K_Local = np.dot(np.dot(np.transpose(T), K), T)

    return(K_Local)


def K_beam6(element, job):
    # This function receives the information of a BEAM6 element and returns,
    # it's stiffness matrix
    element_propertie_name = element[1]
    element_propertie = job.properties[element_propertie_name]
    element_material_name = element_propertie.material
    element_material = job.materials[element_material_name]

    E = element_material.young  # Element's young modulus
    A = float(element_propertie.parameters[0])  # Element's area
    I = float(element_propertie.parameters[1])  # Element's moment of inertia

    node_0 = job.nodes[int(element[2])]
    node_1 = job.nodes[int(element[3])]

    x_0 = float(node_0[1])
    y_0 = float(node_0[2])

    x_1 = float(node_1[1])
    y_1 = float(node_1[2])

    L = math.sqrt((x_1 - x_0)**2 + (y_1 - y_0)**2)  # Element's length

    c = (x_1 - x_0) / L  # Element's cos
    s = (y_1 - y_0) / L  # Element's sin

    # Element's stiffness matrix before rotation
    K = np.zeros((6, 6))
    K[0][0] = E * A / L
    K[1][1] = 12 * E * I / (L**3)
    K[2][2] = 4 * E * I / L
    K[3][3] = E * A / L
    K[4][4] = 12 * E * I / (L**3)
    K[5][5] = 4 * E * I / L
    K[0][3] = -E * A / L
    K[3][0] = K[0][3]
    K[1][2] = 6 * E * I / (L**2)
    K[2][1] = K[1][2]
    K[1][4] = -12 * E * I / (L**3)
    K[4][1] = K[1][4]
    K[1][5] = 6 * E * I / (L**2)
    K[5][1] = K[1][5]
    K[2][4] = -6 * E * I / (L**2)
    K[4][2] = K[2][4]
    K[2][5] = 2 * E * I / L
    K[5][2] = K[2][5]
    K[5][4] = -6 * E * I / (L**2)
    K[4][5] = K[5][4]

    # Element's rotation matrix
    T = np.zeros((6, 6))
    T[0][0] = c
    T[0][1] = -s
    T[1][0] = s
    T[1][1] = c
    T[2][2] = 1
    T[3][3] = c
    T[3][4] = -s
    T[4][3] = s
    T[4][4] = c
    T[5][5] = 1

    K_Local = np.dot(np.dot(np.transpose(T), K), T)

    return(K_Local)

def K_beam4(element, job):
    # This function receives the information of a BEAM4 element and returns,
    # it's stiffness matrix
    element_propertie_name = element[1]
    element_propertie = job.properties[element_propertie_name]
    element_material_name = element_propertie.material
    element_material = job.materials[element_material_name]

    E = element_material.young  # Element's young modulus
    A = float(element_propertie.parameters[0])  # Element's area
    I = float(element_propertie.parameters[1])  # Element's moment of inertia

    node_0 = job.nodes[int(element[2])]
    node_1 = job.nodes[int(element[3])]

    x_0 = float(node_0[1])
    y_0 = float(node_0[2])

    x_1 = float(node_1[1])
    y_1 = float(node_1[2])

    L = math.sqrt((x_1 - x_0)**2 + (y_1 - y_0)**2)  # Element's length

    c = (x_1 - x_0) / L  # Element's cos
    s = (y_1 - y_0) / L  # Element's sin

    # Element's stiffness matrix before rotation
    K = np.zeros((4, 4))
    K[0][0] = 12
    K[0][1] = 6 * L
    K[0][2] = -12
    K[0][3] = 6 * L
    K[1][0] = 6 * L
    K[1][1] = 4 * L**2
    K[1][2] = -6 * L
    K[1][3] = 2 * L**2
    K[2][0] = -12
    K[2][1] = -6 * L
    K[2][2] = 12
    K[2][3] = -6 * L
    K[3][0] = 6 * L
    K[3][1] = 2 * L**2
    K[3][2] = -6 * L
    K[3][3] = 4 * L**2

    K = ((E * I)/(L**3)) * K

    # Element's rotation matrix
    T = np.zeros((4, 4))
    T[0][0] = c
    T[0][1] = -s
    T[1][0] = s
    T[1][1] = c
    T[2][2] = c
    T[2][3] = -s
    T[3][2] = s
    T[3][3] = c

    K_Local = np.dot(np.dot(np.transpose(T), K), T)

    return(K_Local)

def M_beam4(element, job):
    # This function receives the information of a BEAM4 element and returns,
    # it's mass matrix
    element_propertie_name = element[1]
    element_propertie = job.properties[element_propertie_name]
    element_material_name = element_propertie.material
    element_material = job.materials[element_material_name]

    # E = element_material.young  # Element's young modulus
    rho = element_material.density  # Element's density
    A = float(element_propertie.parameters[0])  # Element's area
    # I = float(element_propertie.parameters[1])  # Element's moment of inertia

    node_0 = job.nodes[int(element[2])]
    node_1 = job.nodes[int(element[3])]

    x_0 = float(node_0[1])
    y_0 = float(node_0[2])

    x_1 = float(node_1[1])
    y_1 = float(node_1[2])

    L = math.sqrt((x_1 - x_0) ** 2 + (y_1 - y_0) ** 2)  # Element's length

    c = (x_1 - x_0) / L  # Element's cos
    s = (y_1 - y_0) / L  # Element's sin

    # Element's stiffness matrix before rotation
    M = np.zeros((4, 4))
    M[0][0] = 156
    M[0][1] = 22 * L
    M[0][2] = 54
    M[0][3] = -13 * L
    M[1][0] = 22 * L
    M[1][1] = 4 * L ** 2
    M[1][2] = 13 * L
    M[1][3] = -3 * L**2
    M[2][0] = 54
    M[2][1] = 13 * L
    M[2][2] = 156
    M[2][3] = -22 * L
    M[3][0] = -13 * L
    M[3][1] = -3 * L ** 2
    M[3][2] = -22 * L
    M[3][3] = 4 * L ** 2

    M = (rho * A * L / 420) * M

    # Element's rotation matrix
    T = np.zeros((4, 4))
    T[0][0] = c
    T[0][1] = -s
    T[1][0] = s
    T[1][1] = c
    T[2][2] = c
    T[2][3] = -s
    T[3][2] = s
    T[3][3] = c

    m_matrix = np.dot(np.dot(np.transpose(T), M), T)

    return m_matrix

def M_beam6(element, job):
    # This function receives the information of a BEAM6 element and returns,
    # it's mass matrix
    element_propertie_name = element[1]
    element_propertie = job.properties[element_propertie_name]
    element_material_name = element_propertie.material
    element_material = job.materials[element_material_name]

    # E = element_material.young  # Element's young modulus
    rho = element_material.density  # Element's density
    A = float(element_propertie.parameters[0])  # Element's area
    # I = float(element_propertie.parameters[1])  # Element's moment of inertia

    node_0 = job.nodes[int(element[2])]
    node_1 = job.nodes[int(element[3])]

    x_0 = float(node_0[1])
    y_0 = float(node_0[2])

    x_1 = float(node_1[1])
    y_1 = float(node_1[2])

    L = math.sqrt((x_1 - x_0) ** 2 + (y_1 - y_0) ** 2)  # Element's length

    c = (x_1 - x_0) / L  # Element's cos
    s = (y_1 - y_0) / L  # Element's sin


    # Element's stiffness matrix before rotation
    M = np.zeros((6, 6))
    mm = rho * A * L / 420
    ma = rho * A * L / 6
    leng = L

    M[0][0] = 2 * ma
    M[1][1] = 156 * mm
    M[2][2] = 4 * leng**2 * mm
    M[3][3] = 2 * ma
    M[4][4] = 156 * mm
    M[5][5] = 4 * leng**2 * mm
    M[0][3] = ma
    M[3][0] = M[0][3]
    M[1][2] = 22 * leng * mm
    M[2][1] = M[1][2]
    M[1][4] = 54 * mm
    M[4][1] = M[1][4]
    M[1][5] = -13 * leng * mm
    M[5][1] = M[1][5]
    M[2][4] = 13 * leng * mm
    M[4][2] = M[2][4]
    M[2][5] = -3 * leng**2 * mm
    M[5][2] = M[2][5]
    M[5][4] = -22 * leng * mm
    M[4][5] = M[5][4]

    # Element's rotation matrix
    T = np.zeros((6, 6))
    T[0][0] = c
    T[0][1] = -s
    T[1][0] = s
    T[1][1] = c
    T[2][2] = 1
    T[3][3] = c
    T[3][4] = -s
    T[4][3] = s
    T[4][4] = c
    T[5][5] = 1

    M_Local = np.dot(np.dot(np.transpose(T), M), T)

    return(M_Local)

def plot(title, node_numx, node_numy, results_global):
    scale_factor = float(input("Enter the Scale Factor for the displacements: "))

    if scale_factor == 0:
        print("Scale Factor cannot be zero, set to 1")
        scale_factor = 1

    node_res_x = np.zeros(len(node_numx))
    node_res_y = np.zeros(len(node_numx))
    node_res_rz = np.zeros(len(node_numx))

    for i in range(len(node_numx)):
        node_res_x[i] = node_numx[i] + scale_factor * results_global[i * 3]
        node_res_y[i] = node_numy[i] + scale_factor * results_global[i * 3 + 1]
        node_res_rz[i] = results_global[i * 3 + 2]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    for element in job.elements:
        x = [node_numx[int(element[2])], node_numx[int(element[3])]]
        y = [node_numy[int(element[2])], node_numy[int(element[3])]]
        x_def = [node_res_x[int(element[2])], node_res_x[int(element[3])]]
        y_def = [node_res_y[int(element[2])], node_res_y[int(element[3])]]
        ax1.plot(x, y, 'ko-')
        ax1.plot(x_def, y_def, 'bs-')

    plt.axis('equal')
    plt.title(title)
    plt.show()

def damp_rayleigh(mass_matrix, stiff_matrix, freq_vector):

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

def damp_bismark(mass_matrix, eig_vectors, freq_vector):

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
###############################################################################
###############################################################################
###############################################################################

# Execution code

# Read Input File

print("#####################################")
print("##         SOLVER - 2D             ##")
print("#####################################")

while True:
    try:
        input_file_name = input("Enter the input file name: ")
        input_file = open(input_file_name, 'r')
        break
    except OSError:
        print("ERROR: The file was not found\n")

print("\nReading input file ...")
star_time = time.clock()
input_file = open(input_file_name, 'r')
job = read_input_file(input_file)
input_file.close()
print("Time used: " + str(round(time.clock() - star_time, 4)) + "s\n")

N_DOF = 3 * len(job.nodes)
K_Global = np.zeros((N_DOF, N_DOF))
DOF_Global = np.zeros((N_DOF, 1))
F_Global = np.zeros((N_DOF, 1))

if job.analysis == "MODAL":
    M_Global = np.zeros((N_DOF, N_DOF))

active_DOFs = [False for x in range(N_DOF)]

star_time = time.clock()
print("Creating problem matrix ...")
total_elements = len(job.elements)

for element in job.elements:

    if job.properties[element[1]].kind == "BAR":

        K_Local = K_bar(element, job)
        K_MAP = np.zeros(4)
        K_MAP[0] = int(element[2]) * 3
        K_MAP[1] = int(element[2]) * 3 + 1
        K_MAP[2] = int(element[3]) * 3
        K_MAP[3] = int(element[3]) * 3 + 1

        for position in K_MAP:

            active_DOFs[int(position)] = True

        for i in range(len(K_Local)):
            for j in range(len(K_Local)):

                i_g = int(K_MAP[i])
                j_g = int(K_MAP[j])

                # Relate nodes with Degrees of Freedom
                K_Global[i_g][j_g] += K_Local[i][j]

    elif job.properties[element[1]].kind == "BEAM4":

        K_Local = K_beam4(element, job)
        K_MAP = np.zeros(4)
        K_MAP[0] = int(element[2]) * 3 + 1
        K_MAP[1] = int(element[2]) * 3 + 2
        K_MAP[2] = int(element[3]) * 3 + 1
        K_MAP[3] = int(element[3]) * 3 + 2

        for position in K_MAP:
            active_DOFs[int(position)] = True

        for i in range(len(K_Local)):
            for j in range(len(K_Local)):
                i_g = int(K_MAP[i])
                j_g = int(K_MAP[j])

                # Relate nodes with Degrees of Freedom
                K_Global[i_g][j_g] += K_Local[i][j]

    elif job.properties[element[1]].kind == "BEAM6":

        K_Local = K_beam6(element, job)
        K_MAP = np.zeros(6)
        K_MAP[0] = int(element[2]) * 3
        K_MAP[1] = int(element[2]) * 3 + 1
        K_MAP[2] = int(element[2]) * 3 + 2
        K_MAP[3] = int(element[3]) * 3
        K_MAP[4] = int(element[3]) * 3 + 1
        K_MAP[5] = int(element[3]) * 3 + 2

        for position in K_MAP:

            active_DOFs[int(position)] = True

        for i in range(len(K_Local)):
            for j in range(len(K_Local)):

                i_g = int(K_MAP[i])
                j_g = int(K_MAP[j])

                # Relate nodes with Degrees of Freedom
                K_Global[i_g][j_g] += K_Local[i][j]

    if job.analysis == "MODAL":

        if job.properties[element[1]].kind == "BEAM4":

            M_Local = M_beam4(element, job)
            M_MAP = np.zeros(4)
            M_MAP[0] = int(element[2]) * 3 + 1
            M_MAP[1] = int(element[2]) * 3 + 2
            M_MAP[2] = int(element[3]) * 3 + 1
            M_MAP[3] = int(element[3]) * 3 + 2

            for position in M_MAP:
                active_DOFs[int(position)] = True

            for i in range(len(M_Local)):
                for j in range(len(M_Local)):
                    i_g = int(M_MAP[i])
                    j_g = int(M_MAP[j])

                    # Relate nodes with Degrees of Freedom
                    M_Global[i_g][j_g] += M_Local[i][j]

        if job.properties[element[1]].kind == "BEAM6":

            M_Local = M_beam6(element, job)
            M_MAP = np.zeros(6)
            M_MAP[0] = int(element[2]) * 3
            M_MAP[1] = int(element[2]) * 3 + 1
            M_MAP[2] = int(element[2]) * 3 + 2
            M_MAP[3] = int(element[3]) * 3
            M_MAP[4] = int(element[3]) * 3 + 1
            M_MAP[5] = int(element[3]) * 3 + 2

            for position in M_MAP:
                active_DOFs[int(position)] = True

            for i in range(len(M_Local)):
                for j in range(len(M_Local)):
                    i_g = int(M_MAP[i])
                    j_g = int(M_MAP[j])

                    # Relate nodes with Degrees of Freedom
                    M_Global[i_g][j_g] += M_Local[i][j]

# Create F_Global Matrix
for key, value in job.loads.items():

    single_load = value

    if single_load.kind == "CONCENTRATED":

        node = int(single_load.parameters[0])
        FX = float(single_load.parameters[1])
        FY = float(single_load.parameters[2])
        MZ = float(single_load.parameters[3])
        F_Global[node * 3] += float(FX)
        F_Global[node * 3 + 1] += float(FY)
        F_Global[node * 3 + 2] += float(MZ)

    elif single_load.kind == "DISTRIBUTED":

        node0_number = int(single_load.parameters[0])
        node0_vector = job.nodes[node0_number]
        node0_x = float(node0_vector[1])
        node0_y = float(node0_vector[2])
        node1_number = int(single_load.parameters[1])
        node1_vector = job.nodes[node1_number]
        node1_x = float(node1_vector[1])
        node1_y = float(node1_vector[2])

        L = math.sqrt((node1_x - node0_x)**2 + (node1_y - node0_y)**2)
        QX = float(single_load.parameters[2])
        QY = float(single_load.parameters[3])

        F_Global[node0_number * 3] += float(QX * L / 2)
        F_Global[node0_number * 3 + 1] += float(QY * L / 2)
        F_Global[node0_number * 3 + 2] += float(QY * L**2 / 12)

        F_Global[node1_number * 3] += float(QX * L / 2)
        F_Global[node1_number * 3 + 1] += float(QY * L / 2)
        F_Global[node1_number * 3 + 2] += float(-QY * L**2 / 12)

inactive_DOFs = []

for i, value in enumerate(active_DOFs):

    if not value:

        inactive_DOFs.append(i)

K_Global = np.delete(K_Global, inactive_DOFs, 0)
K_Global = np.delete(K_Global, inactive_DOFs, 1)

if job.analysis == "MODAL":
    M_Global = np.delete(M_Global, inactive_DOFs, 0)
    M_Global = np.delete(M_Global, inactive_DOFs, 1)

F_Global = np.delete(F_Global, inactive_DOFs, 0)

DOF = ["FREE" for x in range(N_DOF)]

for key, value in job.constrains.items():

    single_constrain = value
    DOF[int(single_constrain.node) * 3] = single_constrain.X
    DOF[int(single_constrain.node) * 3 + 1] = single_constrain.Y
    DOF[int(single_constrain.node) * 3 + 2] = single_constrain.RZ


DOF = np.delete(DOF, inactive_DOFs, 0)

unknown = []
for i, D_DOF in enumerate(DOF):
    if D_DOF == 'FREE':
        unknown.append(i)

stiff_matrix = np.zeros((len(unknown), len(unknown)))
force_vector = np.zeros((len(unknown), 1))

if job.analysis == "MODAL":
    mass_matrix = np.zeros((len(unknown), len(unknown)))

for i, i_g in enumerate(unknown):
    force_vector[i][0] = F_Global[i_g][0]
    for j, j_g in enumerate(unknown):
        stiff_matrix[i][j] = K_Global[i_g][j_g]

        if job.analysis == "MODAL":
            mass_matrix[i][j] = M_Global[i_g][j_g]

print("Time used: " + str(round(time.clock() - star_time, 4)) + "s\n")

star_time = time.clock()
print("Solving global system ...")

if job.analysis == "MODAL":
    [eig_values, eig_vectors] = linalg.eig(stiff_matrix, mass_matrix)

else:
    displacements = linalg.solve(stiff_matrix, force_vector)

print("Time used: " + str(round(time.clock() - star_time, 4)) + "s\n")
# print(displacements)
###############################################################################
###############################################################################
###############################################################################

# Post Processing

node_numx = np.zeros(len(job.nodes))
node_numy = np.zeros(len(job.nodes))
results_global_active = np.zeros(len(DOF))

for i, value in enumerate(job.nodes):
    node_numx[i] = float(value[1])
    node_numy[i] = float(value[2])

if job.analysis == "MODAL":

    idx = eig_values.argsort()[::1]
    eig_values = eig_values[idx]
    eig_vectors = eig_vectors[:, idx]

    freq_vector = np.sqrt(abs(eig_values))

    print("Modes:")
    print("   #    Frequencies (rad/s)")

    for j, freq in enumerate(freq_vector):
        print('{0:4}    {1:4e}'.format(j, freq))

    # Normalizing eigenvectors by the mass matrix

    for k in range(len(eig_values)):
        eig_vector = eig_vectors[:, k]
        eig_vector_tp = eig_vector.transpose()
        m_ii = np.dot(np.dot(eig_vector_tp, mass_matrix), eig_vector)
        eig_vectors[:, k] = math.sqrt((1 / m_ii)) * eig_vectors[:, k]

    # Calculate damping matrix

    calc_yn = input("\nCalculate Damping Matrix (Y/N): ")
    calc_yn = calc_yn.upper()

    damp_matrix = np.zeros((len(eig_values), len(eig_values)))

    if calc_yn == 'Y':
        damp_method = input("Method (R for Rayleigh, B for Bismark)?: ")
        damp_method = damp_method.upper()

        if damp_method == 'R':
            damp_matrix, damp_vector = damp_rayleigh(mass_matrix,
                                                     stiff_matrix,
                                                     freq_vector)
        elif damp_method == 'B':
            damp_matrix, damp_vector = damp_bismark(mass_matrix,
                                                    eig_vectors,
                                                    freq_vector)
        else:
            print("Method was not recognized, continuing without calculating damping matrix.")

    calc_yn = input("\nCalculate Transient Response (Y/N): ")
    calc_yn = calc_yn.upper()

    if calc_yn == 'Y':
        t_step = math.pi / max(freq_vector)
        disp_0 = np.zeros(len(freq_vector))
        vel_0 = np.zeros(len(freq_vector))

        t_f = float(input("Enter simulation time: "))

        disp, vel, acc = transient(mass_matrix, stiff_matrix, damp_matrix,
                                   force_vector, disp_0, vel_0, 0.0001, 0,
                                   t_f)

        print(disp[:, 1])
        plt.plot(disp[:, 1])
        plt.show()


    while True:
        n_mode = int(input("\nChoose a mode to plot (-1 to close): "))

        if n_mode == -1:
            break

        displacements = eig_vectors[:, n_mode]

        for i, i_g in enumerate(unknown):
            results_global_active[i_g] = displacements[i]

        results_global = np.zeros(N_DOF)

        j = 0

        for i in range(N_DOF):
            if active_DOFs[i] == False:
                results_global[i] = 0
                j -= 1
            else:
                results_global[i] = results_global_active[j]
            j += 1

        plot(job.title, node_numx, node_numy, results_global)



else:

    for i, i_g in enumerate(unknown):
        results_global_active[i_g] = displacements[i]

    results_global = np.zeros(N_DOF)

    j = 0

    for i in range(N_DOF):
        if active_DOFs[i] == False:
            results_global[i] = 0
            j -= 1
        else:
            results_global[i] = results_global_active[j]
        j += 1

    print("Displacements and Rotations:")
    print("#       Delta_X         Delta_Y         Delta_RZ")
    for i in range(0, len(results_global), 3):
        print('{0:4d}    {1:4e}    {2:4e}    {3:4e}'.format(i // 3, results_global[i], results_global[i + 1], results_global[i + 2]))

    plot(job.title, node_numx, node_numy, results_global)

###############################################################################
###############################################################################
###############################################################################

# Plot







