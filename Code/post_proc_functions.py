import matplotlib.pyplot as plt
import numpy as np


def plot_results(title,
                 node_numx,
                 node_numy,
                 results_global,
                 nodes_array,
                 elements_array):

    scale_factor = float(input("\nEnter the Scale Factor for the displacements: "))

    if scale_factor == 0:
        print("ERROR: Scale Factor cannot be zero, set to 1")
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

    for element in elements_array:
        node_0 = element.node_0
        node_1 = element.node_1
        node_0_index = nodes_array.index(node_0)
        node_1_index = nodes_array.index(node_1)

        x = [node_numx[node_0_index], node_numx[node_1_index]]
        y = [node_numy[node_0_index], node_numy[node_1_index]]
        x_def = [node_res_x[node_0_index], node_res_x[node_1_index]]
        y_def = [node_res_y[node_0_index], node_res_y[node_1_index]]
        ax1.plot(x, y, 'ko-')
        ax1.plot(x_def, y_def, 'bs-')

    plt.axis('equal')
    plt.title(title)
    plt.show()

    return None


def post_proc_static(title,
                     nodes_orig_coord,
                     displacements,
                     n_dof,
                     unconstrained_dofs,
                     active_dofs,
                     degrees_of_freedom,
                     nodes_array,
                     elements_array):

    results_global_active = np.zeros(len(degrees_of_freedom))

    for i, i_g in enumerate(unconstrained_dofs):
        results_global_active[i_g] = displacements[i]

    results_global = np.zeros(n_dof)

    j = 0

    for i in range(n_dof):
        if active_dofs[i] == False:
            results_global[i] = 0
            j -= 1
        else:
            results_global[i] = results_global_active[j]
        j += 1

    print("\nDisplacements and Rotations:")
    print("#       Delta_X         Delta_Y         Delta_RZ")
    for i in range(0, len(results_global), 3):
        print('{0:4d}    {1:4e}    {2:4e}    {3:4e}'.format(i // 3, results_global[i], results_global[i + 1], results_global[i + 2]))

    plot_results(title,
                 nodes_orig_coord[:, 0],
                 nodes_orig_coord[:, 1],
                 results_global,
                 nodes_array,
                 elements_array)

    return None

def post_proc_modal():
    pass

def post_proc_transdir():
    pass