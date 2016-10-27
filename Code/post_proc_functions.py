import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
plt.style.use('ggplot')

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


def post_proc_modal(title,
                    nodes_orig_coord,
                    eig_vectors,
                    freq_vector,
                    n_dof,
                    unconstrained_dofs,
                    active_dofs,
                    degrees_of_freedom,
                    nodes_array,
                    elements_array):

    print("Modes:")
    print("   #    Frequencies (rad/s)")

    for j, freq in enumerate(freq_vector):
        print('{0:4}    {1:4e}'.format(j, freq))


    while True:

        n_mode = int(input("\nChoose a mode to plot (-1 to close): "))

        if n_mode == -1:
            break

        else:

            results_global_active = np.zeros(len(degrees_of_freedom))

            displacements = eig_vectors[:, n_mode]

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


            plot_results(title + "_Mode_" + str(n_mode),
                         nodes_orig_coord[:, 0],
                         nodes_orig_coord[:, 1],
                         results_global,
                         nodes_array,
                         elements_array)


def post_proc_transdir(title,
                       nodes_orig_coord,
                       displacements,
                       time_arr,
                       nodes_array,
                       elements_array):

    scale_factor = float(
        input("\nEnter the Scale Factor for the displacements: "))

    if scale_factor == 0:
        print("ERROR: Scale Factor cannot be zero, set to 1")
        scale_factor = 1

    node_res_x = np.zeros((len(nodes_orig_coord), len(time_arr)))
    node_res_y = np.zeros((len(nodes_orig_coord), len(time_arr)))
    node_res_rz = np.zeros((len(nodes_orig_coord), len(time_arr)))

    for i in range(len(nodes_orig_coord)):
        node_res_x[i, :] = nodes_orig_coord[i, 0] + scale_factor * \
                                                    displacements[i * 3, :]
        node_res_y[i, :] = nodes_orig_coord[i, 1] + scale_factor * \
                                                    displacements[i * 3 + 1, :]
        node_res_rz[i, :] = nodes_orig_coord[i, 2] + \
                            scale_factor * displacements[i * 3 + 2, :]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    def animate(i):

        ax1.clear()

        for element in elements_array:
            node_0 = element.node_0
            node_1 = element.node_1
            node_0_index = nodes_array.index(node_0)
            node_1_index = nodes_array.index(node_1)

            x = [nodes_orig_coord[node_0_index, 0],
                 nodes_orig_coord[node_1_index, 0]]
            y = [nodes_orig_coord[node_0_index, 1],
                 nodes_orig_coord[node_1_index, 1]]

            x_def = [node_res_x[node_0_index][i],
                     node_res_x[node_1_index][i]]
            y_def = [node_res_y[node_0_index][i],
                     node_res_y[node_1_index][i]]

            ax1.plot(x, y, 'ko-')
            ax1.plot(x_def, y_def, 'bs-')
            ax1.text(2, 0.65, "Time: %.4f" % time_arr[i])

        plt.axis('equal')
        plt.title(title)
        plt.ylim(-.4, 0.4)
        time_string = "time: " + str(round(time_arr[i], 4)) + " s"
        plt.text(0.02, 0.9, time_string, transform=ax1.transAxes,
                 fontsize=15)

    sample = [x for x in range(0, len(time_arr), int(len(time_arr)/ 1000))]
    ani = animation.FuncAnimation(fig, animate, sample)
    plt.show()

def post_proc_buckling(title,
                       nodes_orig_coord,
                       eig_vectors,
                       eig_values,
                       n_dof,
                       unconstrained_dofs,
                       active_dofs,
                       degrees_of_freedom,
                       nodes_array,
                       elements_array)

    print("Buckling Critical Load:")
    print("%.4e N" % max(eig_values))

    n_mode = int(eig_values.index(max(eig_values)))

    results_global_active = np.zeros(len(degrees_of_freedom))
    displacements = eig_vectors[:, n_mode]

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

    plot_results(title),
                 nodes_orig_coord[:, 0],
                 nodes_orig_coord[:, 1],
                 results_global,
                 nodes_array,
                 elements_array)