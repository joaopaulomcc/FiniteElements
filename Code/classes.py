# This file contains all the classes used by the software

import numpy as np


class Material:
    # stores the material name and physical properties
    # density, young's modulus and poisson's ratio

    instances = []

    def __init__(self, name, density, young, poisson):
        self.__class__.instances.append(self)
        self.name = name
        self.density = density
        self.young = young
        self.poisson = poisson


class Property:
    # A property has a kind (BAR, BEAM4 or BEAM6), a material
    # and some parameters, for the BAR the parameters
    # are only the area, for the BEAM the parameters
    # are the area and the moment of inertia

    instances = []

    def __init__(self, name, kind, material, area, inertia):
        self.__class__.instances.append(self)
        self.name = name
        self.kind = kind
        self.material = material
        self.area = area
        self.inertia = inertia


class Load:
    # A load have a name, a kind, CONCENTRATED, DISTRIBUTED or NO_LOAD
    # for the CONCENTRATED load the parameters are the FX, FY and MZ
    # components for the load, for the DISTRIBUTED load the parameters
    # are the fx, fy, mz components. For th NO_LOAD kind the parameter
    # is just a zero

    instances = []

    def __init__(self, name, kind, fx, fy, mz):
        self.__class__.instances.append(self)
        self.name = name
        self.kind = kind
        self.fx = fx
        self.fy = fy
        self.mz = mz


class AppliedLoad:

    # An AppliedLoad is created by applying a Load to a point, in the
    # case of a CONCENTRATED load, or a Line in the case of a
    # DISTRIBUTED load.

    instances = []

    def __init__(self, load, counter, nodes):
        self.__class__.instances.append(self)
        self.name = load.name + "_" + str(counter)
        self.kind = load.kind
        self.nodes = nodes
        self.fx = load.fx
        self.fy = load.fy
        self.mz = load.mz


class Constrain:
    # A constrain has a name, the node where it is applied
    # and it's values in the X axis, Y axis and rotation around
    # the Z axis, if a constrain does not exist, it's value is
    # the setting "FREE".

    instances = []

    def __init__(self, name, point, x, y, rz):
        self.__class__.instances.append(self)
        self.name = name
        self.point = point.number
        self.x = x
        self.y = y
        self.rz = rz


class Point:
    # A point is used to create Lines and apply Loads an constrains
    # to the structure. It stores it's number, the x, y coordinates
    # of the point and the Load applied to it

    instances = []

    def __init__(self, number, x, y, load):
        self.__class__.instances.append(self)
        self.number = number
        self.x = x
        self.y = y
        self.load = load


class Line:
    # A Line is a bar or a beam, it is created between two Points.
    # Every line has a number, a Property, a Load, the number (n)
    # of elements in which it should be divided

    instances = []

    def __init__(self, number, point_0, point_1, prop, load, n):
        self.__class__.instances.append(self)
        self.number = number
        self.point_0 = point_0
        self.point_1 = point_1
        self.prop = prop
        self.load = load
        self.n = n

    def mesh(self):
        # The mesh method returns an numpy array with th x,y coordinates of
        # the nodes inside the line

        x_0 = self.point_0.x
        y_0 = self.point_0.y
        x_1 = self.point_1.x
        y_1 = self.point_1.y
        n = self.n

        delta_x = (x_1 - x_0) / n
        delta_y = (y_1 - y_0) / n

        local_nodes = np.zeros((n + 1, 2))
        local_nodes[0][0] = x_0
        local_nodes[0][1] = y_0

        for i in range(1, n + 1):
            local_nodes[i][0] = local_nodes[i - 1][0] + delta_x
            local_nodes[i][1] = local_nodes[i - 1][1] + delta_y

        return local_nodes


class Node:
    # A Node

    instances = []

    def __init__(self, x, y, rz):
        self.__class__.instances.append(self)
        self.x = x
        self.y = y
        self.rz = rz


class Element:
    # An element is a finite element, it has a property, and the
    # nodes that it connects

    instances = []

    def __init__(self, prop, node_0, node_1):
        self.__class__.instances.append(self)
        self.prop = prop
        self.node_0 = node_0
        self.node_1 = node_1
