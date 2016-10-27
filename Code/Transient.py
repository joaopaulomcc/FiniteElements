import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt


def transient(mass_matrix,
              stiff_matrix,
              damp_matrix,
              force_vector,
              constrained_dofs,
              disp_0,
              vel_0,
              t_step,
              t_0,
              t_f):

    # This function uses a the central difference method to simulate the
    # the transient response of a structure

    m = mass_matrix
    k = stiff_matrix
    c = damp_matrix

    force_mag = force_vector[:, 0]
    force_freq = force_vector[:, 1]
    force_phase = force_vector[:, 2]

    n_t_steps = int(np.fix(((t_f - t_0) / t_step)))
    disp = np.zeros((len(disp_0), n_t_steps + 1))
    vel = np.zeros((len(disp_0), n_t_steps + 1))
    acc = np.zeros((len(disp_0), n_t_steps + 1))
    force = np.zeros((len(disp_0), n_t_steps + 1))

    disp[:, 0] = disp_0[:, 0]
    vel[:, 0] = vel_0[:, 0]
    force[:, 0] = force_mag * np.cos(force_freq * 0 + force_phase)

    #print(vel_fic)


    m_inv = linalg.inv(m)

    for i in range(n_t_steps):
        aux_matrix = force[:, i] - np.dot(c, vel[:, i]) - np.dot(k, disp[:, i])
        acc[:, i] = np.dot(m_inv, aux_matrix)

        # Sets constrained degrees of freedom acceleration to zero
        for j in range(len(acc[:, i])):
            if j in constrained_dofs:
                acc[j, i] = 0


        vel[:, i + 1] = vel[:, i] + t_step * acc[:, i]

        disp[:, i + 1] = disp[:, i] + t_step * vel[:, i + 1]
        force[:, i + 1] = force_mag * \
                          np.cos(force_freq * (i + i) * t_step + force_phase)

    aux_matrix = force[:, n_t_steps] - np.dot(c, vel[:, n_t_steps]) - \
                 np.dot(k, disp[:, n_t_steps])
    acc[:, n_t_steps] = np.dot(m_inv, aux_matrix)

    time = np.array([t_0 + i * t_step for i in range(n_t_steps + 1)])
    plt.plot(time, disp[10, :])
    plt.show()

#    print("M inv")
#    print(m_inv)
#    print("acc")
#    print(acc)
#    print("vel")
#    print(vel)
#    print("disp")
#    print(disp)
    return disp, vel, acc, force, time
