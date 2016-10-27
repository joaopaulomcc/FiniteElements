import time
import numpy as np
from scipy import linalg
import math

from classes import Material
from classes import Property
from classes import Load
from classes import AppliedLoad
from classes import Constrain
from classes import Point
from classes import Line
from classes import Node
from classes import Element

from Transient import transient

from functions import create_mesh
from functions import local2global
from functions import get_normalized_eigs
from functions import get_damp_matrix
from functions import damp_rayleigh
from functions import damp_bismarck
from element_functions import local_stiff_matrix
from element_functions import local_mass_matrix
from element_functions import local_geo_stiff_matrix
from element_functions import global_force_vector

from post_proc_functions import post_proc_static
from post_proc_functions import post_proc_modal
from post_proc_functions import post_proc_transdir
from post_proc_functions import post_proc_buckling

import matplotlib.animation as animation
import matplotlib.pyplot as plt

print("######################################################################")
print("##               FINITE ELEMENT ANALYSIS - João Paulo               ##")
print("######################################################################")

###############################################################################
###                                                                         ###
###                               PREPROCESSOR                              ###
###                                                                         ###
###############################################################################

###############################################################################
# Reading input file
while True:
    try:
        input_file_name = input("\nEnter the input file name: ")
        input_file = open(input_file_name, 'r')
        try:
            print("\nReading input file ...")
            star_time = time.clock()
            exec(input_file.read())
            input_file.close()
            print("Time used: " + str(round(time.clock() - star_time, 4)) +
                  "s")
            break
        except:
            print("ERROR: The input file was not valid, check for errors.")
    except OSError:
        print("ERROR: The file was not found")

###############################################################################
# Selecting Analysis type
while True:
    print("\nWhat kind of analysis do you want to perform?")
    print("[0] - Static Analysis")
    print("[1] - Modal Analysis")
    print("[2] - Transient Analysis - Direct Integration")
    print("[3] - Transient Analysis - Modal Analysis")
    print("[4] - Buckling Analysis")
    print("[5] - Modal Analysis - Reduced System")
    analysis_type = input("Choose one of the options above [#]: ")
    if analysis_type == '0':
        analysis_type = "STATIC"
        break
    elif analysis_type == '1':
        analysis_type = "MODAL"
        break
    elif analysis_type == '2':
        analysis_type = "TRANSDIR"
        break
    elif analysis_type == '3':
        analysis_type = "TRANSMOD"
        break
    elif analysis_type == '4':
        analysis_type = "BUCKLING"
        break
    elif analysis_type == '5':
        analysis_type = "MODALRED"
        break
    else:
        print("\nINVALID INPUT, try again.")

###############################################################################
# Creating Nodes, Elements and Applies Loads

star_time = time.clock()
print("\nCreating mesh ...")
lines_array = Line.instances
points_array = Point.instances
nodes_array, elements_array, app_loads_array = create_mesh(lines_array,
                                                           points_array)
print("Number of Nodes created: " + str(len(nodes_array)))
print("Number of Elements created: " + str(len(elements_array)))
print("Time used: " + str(round(time.clock() - star_time, 4)) + "s")

###############################################################################
# Creating Global Matrices

print("\nCreating System Matrices ...")
star_time = time.clock()

n_dof = 3 * len(nodes_array)  # Number of degrees of freedom
stiff_global = np.zeros((n_dof, n_dof))  # Global Stiffness Matrix

if analysis_type == "BUCKLING":
    geo_stiff_global = np.np.zeros((n_dof, n_dof))

elif analysis_type != "STATIC":
    mass_global = np.zeros((n_dof, n_dof))  # Global Mass Matrix

disp_global = np.zeros((n_dof, 1))  # Global Displacement vector
force_global = np.zeros((n_dof, 3))  # Global Force Vector
active_dofs = [False for x in range(n_dof)]  # Vector that controls what
#                                              are being used

###############################################################################
# Creating Local stiffness and mass matrices and build global matrices

for element in elements_array:
    # Creates global stiffness matrix
    stiff_local, map_vector = local_stiff_matrix(element,
                                                 nodes_array)
    local2global(stiff_local, stiff_global, map_vector)

    for position in map_vector:
        # Mark the active dofs in the active_dofs vector
        active_dofs[int(position)] = True

    if analysis_type == "BUCKLING":
        # Creates Global Geometric Stiffness Matrix
        geo_stiff_local = local_geo_stiff_matrix(element)
        local2global(geo_stiff_local, geo_stiff_global, map_vector)

    elif analysis_type != "STATIC":
        # Creates global mass matrix
        mass_local = local_mass_matrix(element)
        local2global(mass_local, mass_global, map_vector)

loads_array = AppliedLoad.instances

for app_load in app_loads_array:
    # Create the Global Force Vector
    global_force_vector(app_load, force_global, nodes_array)

###############################################################################
# Eliminate unused degrees of freedom, for example, in a node connected by a
# BAR element only the rotational degree of freedom is inactive

inactive_dofs = []

for i, value in enumerate(active_dofs):

    if not value:
        inactive_dofs.append(i)

stiff_global = np.delete(stiff_global, inactive_dofs, 0)
stiff_global = np.delete(stiff_global, inactive_dofs, 1)

if analysis_type == "BUCKLING":
    geo_stiff_global = np.delete(geo_stiff_global, inactive_dofs, 0)
    geo_stiff_global = np.delete(geo_stiff_global, inactive_dofs, 1)

elif analysis_type != "STATIC":
    mass_global = np.delete(mass_global, inactive_dofs, 0)
    mass_global = np.delete(mass_global, inactive_dofs, 1)

force_global = np.delete(force_global, inactive_dofs, 0)

###############################################################################
# Check what degrees of freedom are constrained

degrees_of_freedom = ["FREE" for x in range(n_dof)]

for constrain in Constrain.instances:

    node = constrain.node  # node is the number of the point
#                            where the constrain was applied

    degrees_of_freedom[node * 3] = constrain.x
    degrees_of_freedom[node * 3 + 1] = constrain.y
    degrees_of_freedom[node * 3 + 2] = constrain.rz

degrees_of_freedom = np.delete(degrees_of_freedom, inactive_dofs, 0)

constrained_dofs = []
unconstrained_dofs = []

for i, dof_value in enumerate(degrees_of_freedom):
    if dof_value == 'FREE':
        unconstrained_dofs.append(i)
    else:
        constrained_dofs.append(i)


print("Time used: " + str(round(time.clock() - star_time, 4)) + "s")

###############################################################################
###                                                                         ###
###                                  SOLVER                                 ###
###                                                                         ###
###############################################################################


###############################################################################
# Calculates Transient Response by direct integration

if analysis_type == "TRANSDIR":

    # Calculates eigenvalues and eingenvectors
    [eig_values, eig_vectors] = get_normalized_eigs(stiff_global,
                                                    mass_global)

    # Calculates Natural frequencies
    freq_vector = np.sqrt(abs(eig_values))  # abs is used to avoid numerical
    #                                         problems with eigenvalues that
    #                                         are zero but are calculated
    #                                         as very small negative numbers

    #Calculates Maximum time step
    max_t_step = (1 / (max(freq_vector) * math.pi))

    # Calculates Damping Matrix
    damp_matrix = get_damp_matrix(mass_global,
                                  stiff_global,
                                  eig_vectors,
                                  freq_vector)

    # Calculates initial state vectors

    disp_0 = np.zeros((len(degrees_of_freedom), 1))

    for i, dof_value in enumerate(degrees_of_freedom):

        if dof_value == "FREE":
            disp_0[i][0] = 0
        else:
            disp_0[i][0] = dof_value

    vel_0 = np.zeros((len(degrees_of_freedom), 1))

    # Ask user for simulation parameters

    while True:

        try:
            print("\nSimulation Parameters:")
            t_step = float(input("Enter time step (max = %.4e):" % max_t_step))
            t_0 = float(input("Enter initial time: "))
            t_f = float(input("Enter final time: "))
            break

        except ValueError:
            print("ERROR: The times mus be numbers")

    print("\nRunning Transient Analysis ...")
    star_time = time.clock()

    disp, vel, acc, force, time_arr = transient(mass_global,
                                                stiff_global,
                                                damp_matrix,
                                                force_global,
                                                constrained_dofs,
                                                disp_0,
                                                vel_0,
                                                t_step,
                                                t_0,
                                                t_f)

    print("Time used: " + str(round(time.clock() - star_time, 4)) + "s")

    # print(mass_global)
    # print(stiff_global)
    # print(damp_matrix)
    # print(force_global)
    # print(constrained_dofs)
    # print(disp_0)
    # print(vel_0)
    # print(t_step)
    # print(t_0)
    # print(t_f)

###############################################################################
# Calculates Reduced Matrices for Static and Modal Analysis

else:
    print("\nCalculating Reduced Matrices ...")
    star_time = time.clock()

    red_stiff_matrix = np.zeros((len(unconstrained_dofs),
                                 len(unconstrained_dofs)))
    red_force_vector = np.zeros((len(unconstrained_dofs), 1))

    if analysis_type == "BUCKLING":
        red_geo_stiff_matrix = np.zeros((len(unconstrained_dofs),
                                         len(unconstrained_dofs)))

    elif analysis_type != "STATIC":
        red_mass_matrix = np.zeros((len(unconstrained_dofs),
                                    len(unconstrained_dofs)))

    for i, i_g in enumerate(unconstrained_dofs):
        red_force_vector[i][0] = force_global[i_g][0]
        for j, j_g in enumerate(unconstrained_dofs):
            red_stiff_matrix[i][j] = stiff_global[i_g][j_g]

            if analysis_type == "BUCKLING":
                red_geo_stiff_matrix[i][j] = geo_stiff_global[i_g][j_g]

            elif analysis_type != "STATIC":
                red_mass_matrix[i][j] = mass_global[i_g][j_g]

    print("Time used: " + str(round(time.clock() - star_time, 4)) + "s")
###############################################################################
# Calculates Displacements for Static Analysis

if analysis_type == "STATIC":

    print("\nCalculating Displacements ...")
    star_time = time.clock()

    static_force_vector = red_force_vector[:, 0]
    displacements = linalg.solve(red_stiff_matrix, static_force_vector)

    print("Time used: " + str(round(time.clock() - star_time, 4)) + "s")

###############################################################################
# Calculates Natural Frequencies and Modes for Modal Analysis

elif analysis_type == "MODAL":
    print("\nCalculating System natural frequencies and modes ...")
    star_time = time.clock()

    [eig_values, eig_vectors] = get_normalized_eigs(red_stiff_matrix,
                                                    red_mass_matrix)
    freq_vector = np.sqrt(abs(eig_values))

    print("Time used: " + str(round(time.clock() - star_time, 4)) + "s")

    # Calculates Damping Matrix
    damp_matrix = get_damp_matrix(red_mass_matrix,
                                  red_stiff_matrix,
                                  eig_vectors,
                                  freq_vector)
###############################################################################
# Calculates Critical load for Buckling
elif analysis_type == "BUCKLING":
    print("\nCalculating SCritical load for Buckling ...")
    star_time = time.clock()

    [eig_values, eig_vectors] = linalg.eig(red_stiff_matrix,
                                           red_geo_stiff_matrix)

    print("Time used: " + str(round(time.clock() - star_time, 4)) + "s")

###############################################################################
###                                                                         ###
###                           POSTPROCESSOR                                 ###
###                                                                         ###
###############################################################################

###############################################################################
# Get coordinates from the original structure

nodes_orig_coord = np.zeros((len(nodes_array), 3))

for i, node in enumerate(nodes_array):
    nodes_orig_coord[i][0] = node.x
    nodes_orig_coord[i][1] = node.y
    nodes_orig_coord[i][2] = node.rz

###############################################################################
# Print and plot results
print()
print("######################################################################")
print("##                             RESULTS                              ##")
print("######################################################################")

if analysis_type == "STATIC":

    post_proc_static(title,
                     nodes_orig_coord,
                     displacements,
                     n_dof,
                     unconstrained_dofs,
                     active_dofs,
                     degrees_of_freedom,
                     nodes_array,
                     elements_array)

elif analysis_type == "MODAL":

    post_proc_modal(title,
                    nodes_orig_coord,
                    eig_vectors,
                    freq_vector,
                    n_dof,
                    unconstrained_dofs,
                    active_dofs,
                    degrees_of_freedom,
                    nodes_array,
                    elements_array)

elif analysis_type == "TRANSDIR":

    post_proc_transdir(title,
                       nodes_orig_coord,
                       disp,
                       time_arr,
                       nodes_array,
                       elements_array)


elif analysis_type == "BUCKLING":

    post_proc_buckling(title,
                       nodes_orig_coord,
                       eig_vectors,
                       freq_vector,
                       n_dof,
                       unconstrained_dofs,
                       active_dofs,
                       degrees_of_freedom,
                       nodes_array,
                       elements_array)




