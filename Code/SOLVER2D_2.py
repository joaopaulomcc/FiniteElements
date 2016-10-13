import numpy as np
import math

###############################################################################
###############################################################################
###############################################################################

# Classes
# Classes created in order to store the data retrieved from the
# input file.


class material:
    # stores the material name and physical properties
    # density, young's modulus and poisson's ratio

    def __init__(self, name, density, young, poisson, *args):
        self.name = name
        self.density = density
        self.young = young
        self.poisson = poisson


class propertie:
    # A propertie have a kind (BAR or BEAM), a material
    # and some parameters, for the BAR the parameters
    # are only the area, for the BEAM the parameters
    # are the area and the moment of inertia

    def __init__(self, name, kind, material, parameters, *args):
        self.name = name
        self.kind = kind
        self.material = material
        self.parameters = parameters


class load:
    # A load have a name, a kind, CONCENTRATED or DISTRIBUTED,
    # for the CONCENTRATED load the parameters are the node
    # where the force is applied and FX, FY and MZ, for the
    # DISTRIBUTED load the parameters are the two nodes
    # between wich the load is applied and QX and QY.

    def __init__(self, name, kind, parameters, *args):
        self.name = name
        self.kind = kind
        self.parameters = parameters


class constrain:
    # A constrain have a name, the node where it is applied
    # and it's values in the X axis, Y axis and rotation around
    # the Z axis, if a constrain does not exist, it's value is
    # the stting "FREE".

    def __init__(self, name, node, X, Y, RZ):
        self.name = name
        self.node = node
        self.X = X
        self.Y = Y
        self.RZ = RZ


class job_data:
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
    # object contains all the information necessary in order to performe
    # the analisys

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
                materials[line_s[0]] = material(line_s[0],
                                                float(line_s[1]),
                                                float(line_s[2]),
                                                float(line_s[3]))

            elif flag == "properties":
                properties[line_s[0]] = propertie(line_s[0],
                                                  line_s[1],
                                                  line_s[2],
                                                  line_s[3:])

            elif flag == "loads":
                loads[line_s[0]] = load(line_s[0],
                                        line_s[1],
                                        line_s[2:])

            elif flag == "constrains":
                constrains[line_s[0]] = constrain(line_s[0],
                                                  int(line_s[1]),
                                                  line_s[2],
                                                  line_s[3],
                                                  line_s[4])

            elif flag == "elements":
                elements.append(line_s)

            elif flag == "nodes":
                nodes.append(line_s)

    return job_data(title, analysis, materials,
                    properties, loads, constrains, elements, nodes)


def BAR(element, job):
    # This function receives the information of a BAR element and returns,
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
    # This function receives the information of a BEAM element and returns,
    # it's stiffness matrix
    pass

###############################################################################
###############################################################################
###############################################################################

# Execution code

# Read Input File


input_file = open("input2d_2.inp", 'r')
job = read_input_file(input_file)

N_DOF = 3 * len(job.nodes)
K_Global = np.zeros((N_DOF, N_DOF))
DOF_Global = np.zeros((N_DOF, 1))
F_Global = np.zeros((N_DOF, 1))

active_DOFs = [False for x in range(N_DOF)]

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

                # Relate nodes with Degrees of Freedon
                K_Global[i_g][j_g] += K_Local[i][j]

    if job.properties[element[1]].kind == "BEAM":

        K_Local = BAR(element, job)
        K_MAP = np.zeros(4)
        K_MAP[0] = int(element[2]) * 3
        K_MAP[1] = int(element[2]) * 3 + 1
        K_MAP[2] = int(element[2]) * 3 + 2
        K_MAP[3] = int(element[3]) * 3
        K_MAP[4] = int(element[3]) * 3 + 1
        K_MAP[4] = int(element[3]) * 3 + 2

        for position in K_MAP:

            active_DOFs[int(position)] = True

        for i in range(len(K_Local)):
            for j in range(len(K_Local)):

                i_g = int(K_MAP[i])
                j_g = int(K_MAP[j])

                # Relate nodes with Degrees of Freedon
                K_Global[i_g][j_g] += K_Local[i][j]


# Ceate F_Global Matriz
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
        node0_vector = job.nodes(node0_number)
        node0_x = node0_vector[1]
        node0_y = node0_vector[2]
        node1_number = int(single_load.parameters[1])
        node1_vector = job.nodes(node1_number)
        node1_x = node1_vector[1]
        node1_y = node1_vector[2]

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


displacements = np.linalg.solve(Matrix_A, Vector_b)

print(job.title)
print(job.analysis)
print("Displacements")
print(displacements)
