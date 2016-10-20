import numpy as np
from scipy import linalg

def transient(mass_matrix,
              stiff_matrix,
              damp_matrix,
              force_0,
              disp_0,
              vel_0,
              t_step,
              t_0,
              t_f):

    M = mass_matrix
    K = stiff_matrix
    C = damp_matrix

    print(M)
    print(K)
    print(C)
    n_t_steps = int((t_f - t_0) / t_step) + 1
    disp = np.zeros((len(disp_0), (n_t_steps + 1)))
    vel = np.zeros((len(disp_0), (n_t_steps + 1)))
    acc = np.zeros((len(disp_0), (n_t_steps + 1)))
    force = np.zeros((len(disp_0), (n_t_steps + 1)))

    disp[:, 0] = disp_0
    vel[:, 0] = vel_0

    vel_fic = disp[:, 0] - t_step * vel[:, 0] + 0.5 * (t_step**2) * acc[:, 0]
    vel_fic = (disp[:, 0] - vel_fic) / t_step

    for i in range(1, n_t_steps):
        v = vel[:, i].reshape(1,-1).T
        d = disp[:, 0].reshape(1,-1).T
        aux_matrix = force_0 - np.dot(C, v) - np.dot(K, d)
        acc[:, i - 1] = np.dot(linalg.inv(M), aux_matrix)[:,0]

        if i == 1:
            vel[:, i] = vel_fic + t_step * acc[:, i - 1]
        else:
            vel[:, i] = vel[:, i - 1] + t_step * acc[:, i - 1]

        disp[:, i] = disp[:, i - 1] + t_step * vel[:, i]

    aux_matrix = force[:, i] - np.dot(C, vel[:, i]) - np.dot(K, disp[:, 0])
    acc[:, n_t_steps - 1] = np.dot(linalg.inv(M), aux_matrix)

    return disp, vel, acc
