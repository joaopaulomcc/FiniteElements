import numpy as np 
import math
import scipy 

################################################################
#Classes

class material:
    def __init__(self, name, density, young, poisson, *args):
        self.name = name
        self.density = density
        self.young = young
        self.poisson = poisson

class propertie:
    def __init__(self, name, kind, material, parameters, *args):
        self.name = name
        self.kind = kind
        self.material = material
        self.parameters = parameters

class load:
    def __init__(self, name, node, FX, FY, MZ):
        self.name = name
        self.node = node
        self.FX = FX
        self.FY = FY
        self.MZ = MZ

class constrain:
    def __init__(self, name, node, X, Y, RZ):
        self.name = name
        self.node = node
        self.X = X
        self.Y = Y
        self.RZ = RZ

class job_data:
    def __init__(self, materials, properties, loads, constrains, elements, nodes):
        self.materials = materials
        self.properties = properties
        self.loads = loads
        self.constrains = constrains
        self.elements = elements
        self.nodes = nodes

##################################################################################
#Functions
def read_input_file(input_file):
    materials = {}
    properties = {}
    loads = {}
    constrains = {}
    elements = []
    nodes = []

    flag = None

    for line in input_file:

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
                title = line

            elif flag == "analysis":
                analysis = line

            elif flag == "materials":
                materials[line_s[0]] = material(line_s[0], float(line_s[1]), float(line_s[2]), float(line_s[3]))

            elif flag == "properties":
                properties[line_s[0]] = propertie(line_s[0], line_s[1], line_s[2], line_s[3:])

            elif flag == "loads":
                loads[line_s[0]] = load(line_s[0], int(line_s[1]), float(line_s[2]), float(line_s[3]), float(line_s[4]))

            elif flag == "constrains":
                constrains[line_s[0]] = constrain(line_s[0], int(line_s[1]), line_s[2], line_s[3], line_s[4])

            elif flag == "elements":
                elements.append(line_s)

            elif flag == "nodes":
                nodes.append(line_s)

    return job_data(materials, properties, loads, constrains, elements, nodes)

def BAR(element, job):

    # Função que calcula a matrix de rigidez para um elemento de barra

    element_propertie_name = element[1]
    element_propertie = job.properties[element_propertie_name]
    element_material_name = element_propertie.material
    element_material = job.materials[element_material_name]

    E = element_material.young
    A = float(element_propertie.parameters[0])

    # Problems for elements with different number of nodes!!!!
    Node1 = job.nodes[int(element[2]) - 1]
    Node2 = job.nodes[int(element[3]) - 1]

    X1 = float(Node1[1])
    Y1 = float(Node1[2])

    X2 = float(Node2[1])
    Y2 = float(Node2[2])

    L = math.sqrt((X2 - X1)**2 + (Y2 - Y1)**2)

    c = (X2 - X1)/L
    s = (Y2 - Y1)/L

    K = np.array([[1, -1], [-1, 1]])*(E*A/L)
    T = np.array([[c, s, 0, 0], [0, 0, c, s]])

    K_Local = np.dot(np.dot(np.transpose(T), K), T)

    # F_Local = np.array
    # TODO Think about the forces, how to apply then to the elements

    return K_Local
##################################################################################
# Read Input File

input_file = open("input2d.inp", 'r')
job = read_input_file(input_file)

N_DOF = 3*len(job.nodes)
K_Global = np.zeros((N_DOF, N_DOF))
DOF_Global = np.zeros((N_DOF, 1))
F_Global = np.zeros((N_DOF, 1))

active_DOFs = [False for x in range(N_DOF)]

for element in job.elements:	

    if job.properties[element[1]].kind == "BAR":

        K_Local = BAR(element, job)
        K_MAP = np.zeros(4)
        K_MAP[0] = int(element[2])*3 - 3
        K_MAP[1] = int(element[2])*3 - 2
        K_MAP[2] = int(element[3])*3 - 3
        K_MAP[3] = int(element[3])*3 - 2

        for position in K_MAP:

            active_DOFs[int(position)] = True


        for i in range(len(K_Local)):
            for j in range(len(K_Local)):

                i_g = int(K_MAP[i])
                j_g = int(K_MAP[j])

                K_Global[i_g][j_g] += K_Local[i][j]
                # Relate nodes with Degrees of Freedon


for key, value in job.loads.items():

    single_load = value
    F_Global[int(single_load.node)*3 - 3] = float(single_load.FX)
    F_Global[int(single_load.node)*3 - 2] = float(single_load.FY)
    F_Global[int(single_load.node)*3 - 1] = float(single_load.MZ)


inactive_DOFs = []
for i, value in enumerate(active_DOFs):
    if value:
        inactive_DOFs.append(i)

K_Global = np.delete(K_Global, inactive_DOFs, 0)
K_Global = np.delete(K_Global, inactive_DOFs, 1)
F_Global = np.delete(F_Global, inactive_DOFs, 0)

DOF = ["Free" for x in range(N_DOF)]

for key, value in job.constrains.items():

    single_constrain = value
    DOF[int(single_constrain.node)*3 - 3] = single_constrain.X
    DOF[int(single_constrain.node)*3 - 2] = single_constrain.Y
    DOF[int(single_constrain.node)*3 - 1] = single_constrain.RZ


DOF = np.delete(DOF, inactive_DOFs, 0)

unknown = []
for i, D_DOF in enumerate(DOF):
    if D_DOF == 'Free':
        unknown.append(i)

Matrix_A = np.zeros((len(unknown), len(unknown)))
Vector_b = np.zeros((len(unknown), 1))

for i, i_g in enumerate(unknown):
    Vector_b[i][0] = F_Global[i_g][0]
    for j, j_g in enumerate(unknown):
        Matrix_A[i][j] = K_Global[i_g][j_g]

print(Matrix_A)
print(Vector_b)

displacements = np.linalg.solve(Matrix_A, Vector_b)
print(displacements)

for i in range(10):
    print("Hello")