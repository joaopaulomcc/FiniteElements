import numpy as np
import weakref

###############################################################################
###############################################################################
###############################################################################

# Classes Definitions


class Point:

    instances = []

    def __init__(self, number, x, y, load, *args):
        self.__class__.instances.append(weakref.proxy(self))
        self.number = number
        self.x = x
        self.y = y
        self.load = load


class Line:

    instances = []

    def __init__(self, number, point_0, point_1, propertie, load, n, *args):
        self.__class__.instances.append(weakref.proxy(self))
        self.number = number
        self.point_0 = point_0
        self.point_1 = point_1
        self.propertie = propertie
        self.load = load
        self.n = n

    def mesh(self):
        x_0 = self.point_0.x
        y_0 = self.point_0.y
        x_1 = self.point_1.x
        y_1 = self.point_1.y
        n = self.n

        inclination = (y_1 - y_0) / (x_1 - x_0)
        delta_x = (x_1 - x_0) / n

        local_nodes = np.zeros((n + 1, 2))
        local_nodes[0][0] = x_0
        local_nodes[0][1] = y_0

        for i in range(1, n + 1):
            local_nodes[i][0] = local_nodes[i - 1][0] + delta_x
            local_nodes[i][1] = local_nodes[i - 1][1] + delta_x * inclination

        return local_nodes


class propertie:

    instances = []

    def __init__(self, name, kind, material, parameters, *args):
        self.__class__.instances.append(weakref.proxy(self))
        self.name = name
        self.kind = kind
        self.material = material
        self.parameters = parameters


class material:

    instances = []

    def __init__(self, name, density, young, poisson, *args):
        self.__class__.instances.append(weakref.proxy(self))
        self.name = name
        self.density = density
        self.young = young
        self.poisson = poisson


class load:

    instances = []

    def __init__(self, name, kind, parameters, *args):
        self.__class__.instances.append(weakref.proxy(self))
        self.name = name
        self.kind = kind
        self.parameters = parameters


class applied_load():

    instances = []

    def __init__(self, load, counter, nodes, *args):
        self.__class__.instances.append(weakref.proxy(self))
        self.name = load.name + "_" + str(counter)
        self.kind = load.kind
        self.nodes = nodes
        self.parameters = load.parameters


class constrain:

    instances = []

    def __init__(self, name, point, X, Y, RZ, *args):
        self.__class__.instances.append(weakref.proxy(self))
        self.name = name
        self.point = point
        self.X = X
        self.Y = Y
        self.RZ = RZ


###############################################################################
###############################################################################
###############################################################################

# Problem Inputs

# Title

title = "EX_002"

# Analysis Type

analysis_type = "STATIC"

# Materials
Aluminium = material("Aluminium", 2810, 80e9, 0.33)

# Properties
Beam_0 = propertie("Beam_0", "BEAM4", Aluminium, [0.0008, 1.06666667e-07])
Bar_0 = propertie("Bar_0", "BAR", Aluminium, [1.9634954084936207e-05])


# Loads
No_load = load("No_load", "NO_LOAD", [0])
Distributed_0 = load("Distributed_0", "DISTRIBUTED", [0, 10000])

# Points
Points_array = [Point(0, 0, 0.5, No_load),
                Point(1, 2, 0.5, No_load),
                Point(2, 0, 0, No_load)]

# Lines
Lines = [Line(0, Points_array[0], Points_array[1], Beam_0, Distributed_0, 1),
         Line(0, Points_array[2], Points_array[1], Bar_0, No_load, 1)]
# Lines = [Line(0, Points_array[0], Points_array[1], Bar_0, Distributed_0, 1),
#          Line(0, Points_array[2], Points_array[1], Bar_1, Distributed_0, 1)]

# Constrains
Fixed_0 = constrain("Fixed_0", Points_array[0], 0, 0, 0)
Fixed_1 = constrain("Fixed_1", Points_array[2], 0, 0, "FREE")

###############################################################################
###############################################################################
###############################################################################

# Mesh Creation

# Create Nodes
n_nodes = 0
n_elem = 0
node_count = 0
elem_count = 0

for line in Lines:
    n_nodes += line.n + 1
    n_elem += line.n

node_matrix = np.zeros((n_nodes, 3))
elements_matrix = []

counter = 0
applied_loads = []

for i, point in enumerate(Points_array):
    node_matrix[i][0] = point.x
    node_matrix[i][1] = point.y
    node_matrix[i][2] = 0
    node_count += 1

    if point.load.kind == "CONCENTRATED":
        applied_loads.append(applied_load(point.load, counter, [point.number]))
        counter += 1

for line in Lines:
    counter = 0

    local_nodes = line.mesh()

    if len(local_nodes) > 2:

        node_matrix[node_count][0] = local_nodes[1][0]
        node_matrix[node_count][1] = local_nodes[1][1]
        node_matrix[node_count][2] = 0
        node_count += 1

        elements_matrix.append([elem_count,
                                line.propertie.name,
                                line.point_0.number,
                                node_count - 1])
        elem_count += 1

        if line.load.kind == "DISTRIBUTED":
            first_node = elements_matrix[elem_count - 1][2]
            second_node = elements_matrix[elem_count - 1][3]
            applied_loads.append(applied_load(line.load, counter,
                                 [first_node, second_node]))
            counter += 1

        for i in range(1, line.n - 1):
            node_matrix[node_count][0] = local_nodes[i + 1][0]
            node_matrix[node_count][1] = local_nodes[i + 1][1]
            node_matrix[node_count][2] = 0
            node_count += 1

            node_0 = elements_matrix[elem_count - 1][3]
            node_1 = node_count - 1
            elements_matrix.append([elem_count,
                                    line.propertie.name,
                                    node_0,
                                    node_1])
            elem_count += 1

            if line.load.kind == "DISTRIBUTED":
                first_node = elements_matrix[elem_count - 1][2]
                second_node = elements_matrix[elem_count - 1][3]
                applied_loads.append(applied_load(line.load, counter,
                                     [first_node, second_node]))
                counter += 1

        elements_matrix.append([elem_count,
                                line.propertie.name,
                                elements_matrix[elem_count - 1][3],
                                line.point_1.number])

        elem_count += 1

        if line.load.kind == "DISTRIBUTED":
            first_node = elements_matrix[elem_count - 1][2]
            second_node = elements_matrix[elem_count - 1][3]
            applied_loads.append(applied_load(line.load, counter,
                                 [first_node, second_node]))
            counter += 1

    else:
        elements_matrix.append([elem_count,
                                line.propertie.name,
                                line.point_0.number,
                                line.point_1.number])

        if line.load.kind == "DISTRIBUTED":
            first_node = elements_matrix[elem_count - 1][2]
            second_node = elements_matrix[elem_count - 1][3]
            applied_loads.append(applied_load(line.load, counter,
                                 [first_node, second_node]))
            counter += 1

node_matrix = node_matrix[:node_count]

###############################################################################
###############################################################################
###############################################################################

# Write input file for the solver

file_name = title + ".input"
input_file = open(file_name, 'w')

with input_file as f:
    f.write("#TITLE\n")
    f.write(title + "\n")
    f.write("#\n")
    f.write("#ANALYSIS\n")
    f.write(analysis_type + "\n")
    f.write("#\n")
    f.write("#MATERIALS\n")

    for instance in material.instances:
        string_array = [instance.name,
                        str(instance.density),
                        str(instance.young),
                        str(instance.poisson)]
        f.write("    ".join(string_array) + "\n")

    f.write("#\n")
    f.write("#PROPERTIES\n")

    for instance in propertie.instances:
        string_array = [instance.name,
                        instance.kind,
                        instance.material.name]
        for parameter in instance.parameters:
            string_array.append(str(parameter))

        f.write("    ".join(string_array) + "\n")

    f.write("#\n")
    f.write("#LOADS\n")

    for instance in applied_loads:
        string_array = [instance.name,
                        instance.kind]
        if instance.kind == "CONCENTRATED":
            string_array.append(str(instance.nodes[0]))
        elif instance.kind == "DISTRIBUTED":
            string_array.append(str(instance.nodes[0]))
            string_array.append(str(instance.nodes[1]))
        for parameter in instance.parameters:
            string_array.append(str(parameter))

        f.write("    ".join(string_array) + "\n")

    f.write("#\n")
    f.write("#CONSTRAINS\n")

    for instance in constrain.instances:
        string_array = [instance.name,
                        str(instance.point.number),
                        str(instance.X),
                        str(instance.Y),
                        str(instance.RZ)]

        f.write("    ".join(string_array) + "\n")

    f.write("#\n")
    f.write("#ELEMENTS\n")

    for line in elements_matrix:
        string_array = [str(line[0]),
                        str(line[1]),
                        str(line[2]),
                        str(line[3])]

        f.write("    ".join(string_array) + "\n")

    f.write("#\n")
    f.write("#NODES\n")

    for i, line in enumerate(node_matrix):
        string_array = [str(i),
                        str(line[0]),
                        str(line[1]),
                        str(line[2])]
        f.write("    ".join(string_array) + "\n")

input_file.close()
