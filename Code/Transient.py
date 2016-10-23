import numpy as np
from scipy import linalg


def transient(mass_matrix,
              stiff_matrix,
              damp_matrix,
              force,
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

    force_mag = force[:, 0]
    force_freq = force[:, 1]
    force_phase = force[:, 2]

    n_t_steps = int(np.fix(((t_f - t_0) / t_step)))
    disp = np.zeros((len(disp_0), n_t_steps + 1))
    vel = np.zeros((len(disp_0), n_t_steps + 1))
    acc = np.zeros((len(disp_0), n_t_steps + 1))
    force = np.zeros((len(disp_0), n_t_steps + 1))

    disp[:, 0] = disp_0[:, 0]
    vel[:, 0] = vel_0[:, 0]
    force[:, 0] = force_mag * np.cos(force_freq * 0 + force_phase)
    vel_fic = disp[:, 0] - t_step * vel[:, 0] + 0.5 * (t_step**2) * acc[:, 0]
    vel_fic = (disp[:, 0] - vel_fic) / t_step

    m_inv = linalg.inv(m)

    for i in range(n_t_steps):
        aux_matrix = force[:, i] - np.dot(c, vel[:, i]) - np.dot(k, disp[:, i])
        acc[:, i] = np.dot(m_inv, aux_matrix)

        # NEEDS TO BE CHANGED TO THE GENERAL CASE
        acc[0, i] = 0
        acc[1, i] = 0
        if i == 0:
            vel[:, i + 1] = vel_fic + t_step * acc[:, i]
        else:
            vel[:, i + 1] = vel[:, i] + t_step * acc[:, i]

        disp[:, i + 1] = disp[:, i] + t_step * vel[:, i + 1]
        force[:, i + 1] = force_mag * \
                          np.cos(force_freq * (i + i) * t_step  + force_phase)

    aux_matrix = force[:, n_t_steps] - np.dot(c, vel[:, n_t_steps]) - \
                 np.dot(k, disp[:, n_t_steps])
    acc[:, n_t_steps] = np.dot(m_inv, aux_matrix)

    # time = np.array([t_0 + i * t_step for i in range(n_t_steps + 1)])
    # plt.plot(time, disp[8, :])
    # plt.show()

    return disp, vel, acc
