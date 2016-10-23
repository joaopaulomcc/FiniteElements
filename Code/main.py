import time
import numpy as np

from classes import Material
from classes import Property
from classes import Load
from classes import AppliedLoad
from classes import Constrain
from classes import Point
from classes import Line
from classes import Node
from classes import Element

from functions import create_mesh
from functions import local2global
from element_functions import local_stiff_matrix
from element_functions import local_mass_matrix
from element_functions import global_force_vector

print("######################################################################")
print("##                    SOLVER - João Paulo                           ##")
print("######################################################################")

# Main loop

###############################################################################
# Reading input file
while True:
    try:
        input_file_name = input("Enter the input file name: ")
        input_file = open(input_file_name, 'r')
        try:
            print("\nReading input file ...")
            star_time = time.clock()
            exec(input_file.read())
            input_file.close()
            print("Time used: " + str(round(time.clock() - star_time, 4)) +
                  "s")
            break
        except SyntaxError:
            print("The input file was not valid, check for errors.")
    except OSError:
        print("ERROR: The file was not found\n")

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
print("Time used: " + str(round(time.clock() - star_time, 4)) + "s\n")

###############################################################################
# Creating Global Matrices

n_dof = 3 * len(nodes_array)  # Number of degrees of freedom
stiff_global = np.zeros((n_dof, n_dof))  # Global Stiffness Matrix

if analysis_type != "STATIC":
    mass_global = np.zeros((n_dof, n_dof))  # Global Mass Matrix
    # damp_global = np.zeros((n_dof, n_dof))  # Global Damping Matrix

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

    if analysis_type != "STATIC":
        # Creates global mas matrix
        mass_local = local_mass_matrix(element)
        local2global(mass_local, mass_global, map_vector)

loads_array = AppliedLoad.instances

for app_load in app_loads_array:
    # Create the Global Force Vecto
    global_force_vector(app_load, force_global, nodes_array)

###############################################################################
# Eliminate unused degrees of freedom, for example, in a node connected by a
# BAR element only the rotational degree of freedom is inactive

inactive_dofs = []

for i, value in enumerate(active_dofs):

    if not value:
        active_dofs.append(i)

stiff_global = np.delete(stiff_global, inactive_dofs, 0)
stiff_global = np.delete(stiff_global, inactive_dofs, 1)

if analysis_type != "STATIC":
    mass_global = np.delete(mass_global, inactive_dofs, 0)
    mass_global = np.delete(mass_global, inactive_dofs, 1)

force_global = np.delete(force_global, inactive_dofs, 0)

###############################################################################
# Check what degrees of freedom are constrained

# STOPPED AT LINE 756 OF SOLVER_2D
DOF = ["FREE" for x in range(N_DOF)]