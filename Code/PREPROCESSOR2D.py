import numpy as np


class Point:
    def __init__(self, number, x, y, load):
        self.number = number
        self.x = x
        self.y = y
        self.load = load


class Line:
    def __init__(self, number, node_1, node_2, propertie, load, n):
        self.number = number
        self.node_1 = node_1
        self.node_2 = node_2
        self.propertie = propertie
        self.load = load
        self.n = n

    def mesh(self):
        x1 = self.node_1.x
        y1 = self.node_1.y
        x2 = self.node_2.x
        y2 = self.node_2.y
        n = self.n

        inclination = (y2 - y1) / (x2 - x1)
        delta_x = (x2 - x1) / n

        local_nodes = np.zeros((n + 1, 2))
        local_nodes[0][0] = x1
        local_nodes[0][1] = y1

        for i in range(1, n + 1):
            local_nodes[i][0] = local_nodes[i - 1][0] + delta_x
            local_nodes[i][1] = local_nodes[i - 1][1] + delta_x * inclination

        return local_nodes


class propertie:
    def __init__(self, name, kind, material, parameters, *args):
        self.name = name
        self.kind = kind
        self.material = material
        self.parameters = parameters


class material:
    def __init__(self, name, density, young, poisson, *args):
        self.name = name
        self.density = density
        self.young = young
        self.poisson = poisson


class load:
    def __init__(self, name, kind, parameters, *args):
        self.name = name
        self.kind = kind
        self.paramenters = parameters


class constrain:
    def __init__(self, name, node, X, Y, RZ):
        self.name = name
        self.node = node
        self.X = X
        self.Y = Y
        self.RZ = RZ


# Materials
AL7075_T6 = material("AL7075-T6", 2810, 71.7e9, 0.33)

# Properties
Bar_1 = propertie("Bar_1", "BAR", AL7075_T6, [1e4])
Bar_2 = propertie("Bar_2", "BAR", AL7075_T6, [2e4])


# Points
Points_array = [Point(1, 0, 0.5, 0),
                Point(2, 1, 0.5, 0),
                Point(3, 0, 0, 0)]

# Lines
Lines = [Line(1, Points_array[0], Points_array[1], Bar_1, 10, 10),
         Line(2, Points_array[2], Points_array[1], Bar_2, 10, 10)]




# Force_1 = load("Force_1", "CONCENTRATED", [Points_array[1], 0, -1000, 0])

# Distributed_1 = load("Distributed_1", "CONCENTRATED",
#                     [Points_array[0], Points_array[1], -100, 0])

# x = Lines[1]
# print(x.mesh())

n_nodes = 0
n_elem = 0
node_count = 0
elem_count = 0

for line in Lines:
    n_nodes += line.n + 1
    n_elem += line.n

node_matrix = np.zeros((n_nodes, 3))
elements_matrix = []

for i, point in enumerate(Points_array):
    node_matrix[i][0] = point.x
    node_matrix[i][1] = point.y
    node_matrix[i][2] = 0
    node_count += 1

for i, line in enumerate(Lines):
    elements_matrix.append([elem_count + 1,
                            line.propertie.name,
                            line.node_1.number,
                            node_count + 1])
    elem_count += 1

    local_nodes = line.mesh()

    for j in range(1, len((local_nodes)) - 1):
        node_matrix[node_count][0] = local_nodes[j][0]
        node_matrix[node_count][1] = local_nodes[j][1]
        node_matrix[node_count][2] = 0
        node_count += 1

    for k in range(1, line.n):
        node_1 = elements_matrix[k - 1][3]
        node_2 = 
        elements_matrix.append([elem_count + 1,
                                line.propertie.name,
                                node_1,
                                node_2])
        elem_count += 1

print(node_matrix)
print(elements_matrix)
