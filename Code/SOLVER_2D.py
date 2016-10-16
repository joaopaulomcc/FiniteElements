import numpy as np
import math
import matplotlib.pyplot as plt
import time

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
        # what information is sotered in the next lines
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


def BAR(element, job):
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


def BEAM(element, job):
    # This function receIes the information of a BEAM element and returns,
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
    pass

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

active_DOFs = [False for x in range(N_DOF)]

star_time = time.clock()
print("Creating problem matrix ...")
total_elements = len(job.elements)

for element in job.elements:

    if job.properties[element[1]].kind == "BAR":

        K_Local = BAR(element, job)
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

    if job.properties[element[1]].kind == "BEAM":

        K_Local = BEAM(element, job)
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

Matrix_A = np.zeros((len(unknown), len(unknown)))
Vector_b = np.zeros((len(unknown), 1))

for i, i_g in enumerate(unknown):
    Vector_b[i][0] = F_Global[i_g][0]
    for j, j_g in enumerate(unknown):
        Matrix_A[i][j] = K_Global[i_g][j_g]

print("Time used: " + str(round(time.clock() - star_time, 4)) + "s\n")


star_time = time.clock()
print("Solving global system ...")
displacements = np.linalg.solve(Matrix_A, Vector_b)

print("Time used: " + str(round(time.clock() - star_time, 4)) + "s\n")

###############################################################################
###############################################################################
###############################################################################

# Post Processing

results_global_active = np.zeros(len(DOF))

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

print("Displacements and Rotations")
print("#       Delta_X         Delta_Y         Delta_RZ")
for i in range(0, len(results_global), 3):
    print('{0:4d}    {1:4e}    {2:4e}    {3:4e}'.format(i // 3, results_global[i], results_global[i + 1], results_global[i + 2]))

node_numx = np.zeros(len(job.nodes))
node_numy = np.zeros(len(job.nodes))

for i, value in enumerate(job.nodes):
    node_numx[i] = float(value[1])
    node_numy[i] = float(value[2])

###############################################################################
###############################################################################
###############################################################################

# Plot

while True:
    scale_factor = int(input("\nEnter the Scale Factor for the displacements (0 to close the program): "))

    if scale_factor == 0:
        break

    node_res_x = np.zeros(len(node_numx))
    node_res_y = np.zeros(len(node_numx))
    node_res_rz = np.zeros(len(node_numx))

    for i in range(len(node_numx)):
        node_res_x[i] = node_numx[i] + scale_factor * results_global[i * 3]
        node_res_y[i] = node_numy[i] + scale_factor * results_global[i * 3 + 1]
        node_res_rz[i] = results_global[i * 3 + 2]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(node_numx, node_numy, c='k', label="Original")
    ax1.scatter(node_res_x, node_res_y, c='r', label="Deformed")
    plt.axis('equal')
    plt.title(job.title)
    plt.legend()
    plt.show()




