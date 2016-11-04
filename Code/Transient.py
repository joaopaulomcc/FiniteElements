import numpy as np
import sympy as sy
from scipy import linalg
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math


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
                          np.cos(force_freq * (i + 1) * t_step + force_phase)

    aux_matrix = force[:, n_t_steps] - np.dot(c, vel[:, n_t_steps]) - \
                 np.dot(k, disp[:, n_t_steps])
    acc[:, n_t_steps] = np.dot(m_inv, aux_matrix)

    time = np.array([t_0 + i * t_step for i in range(n_t_steps + 1)])

    while True:

        n = int(input("Choose a degree"
                      " of freedom to plot the response (-1 to close): "))

        if n == -1:
            break
        else:
            plt.plot(time, disp[n, :])
            plt.ylabel('Deslocamento (m)')
            plt.xlabel("tempo (s)")
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


def modal_transient(red_mass_matrix,
                    red_damp_matrix,
                    red_stiff_matrix,
                    red_force_vector,
                    eig_vectors,
                    t_0,
                    t_f,
                    t_step):

    # This function uses modal analysis to solve the response of a structure

    n_t_steps = int(np.fix(((t_f - t_0) / t_step)))

    gen_mass_matrix = np.dot(np.dot(eig_vectors.transpose(), red_mass_matrix),
                             eig_vectors)

    gen_damp_matrix = np.dot(np.dot(eig_vectors.transpose(), red_damp_matrix),
                             eig_vectors)

    gen_stiff_matrix = np.dot(np.dot(eig_vectors.transpose(), red_stiff_matrix),
                              eig_vectors)

    mag_force_vector = red_force_vector[:, 0]

    gen_force_vector = np.dot(eig_vectors.transpose(), mag_force_vector)

    disp = np.zeros((len(mag_force_vector), n_t_steps + 1))

    for i in range(len(mag_force_vector)):

        modal_mass = gen_mass_matrix[i][i]
        modal_damp = gen_damp_matrix[i][i]
        modal_stiff = gen_stiff_matrix[i][i]
        force_mag = gen_force_vector[i]
        force_freq = red_force_vector[i][1]
        force_phase = red_force_vector[i][2]

        mode_response = solve_ODE_num(modal_mass,
                                      modal_damp,
                                      modal_stiff,
                                      force_mag,
                                      force_freq,
                                      force_phase,
                                      t_0,
                                      n_t_steps,
                                      t_step)

        disp[i, :] = mode_response

    time = np.array([t_0 + x * t_step for x in range(n_t_steps + 1)])

    for i in range(len(time)):
        disp[:, i] = np.dot(eig_vectors, disp[:, i])

    while True:

        n = int(input("Choose a degree"
                      " of freedom to plot the response (-1 to close): "))

        if n == -1:
            break
        else:
            plt.plot(time, disp[n, :])
            plt.show()


def solve_ODE(modal_mass,
              modal_damp,
              modal_stiff,
              force_mag,
              force_freq,
              force_phase,
              t_0,
              n_t_steps,
              t_step):


    f = sy.Function('f')
    t = sy.Symbol("t", real=True)
    C1, C2 = sy.symbols("C1, C2")

    eq = sy.Eq(modal_mass * sy.diff(f(t), t, 2)
               + modal_damp * sy.diff(f(t), t)
               + modal_stiff * f(t)
               - force_mag * sy.sin(force_freq * t + force_phase))

    sol = sy.dsolve(eq, f(t),
                    hint="nth_linear_constant_coeff_undetermined_coefficients").rhs

    print(sol)

    cond_0 = sy.Eq(sol.subs(t, 0), 0)
    cond_1 = sy.Eq(sol.diff(t).subs(t, 0), 0)

    C1_C2 = sy.solve([cond_0, cond_1], (C1, C2))

    #print(C1_C2)

    y_sl1 = (sol.subs(C1_C2))
    # print("Solution with initial conditions:")
    # sy.pprint(y_sl1)

    f = sy.lambdify(t, y_sl1, "numpy")

    time = np.array([t_0 + i * t_step for i in range(n_t_steps + 1)])
    data = f(time)
    # print(f(0))
    # plt.plot(time, f(time))
    # plt.show()

    return data


def solve_ODE_num(modal_mass,
                  modal_damp,
                  modal_stiff,
                  force_mag,
                  force_freq,
                  force_phase,
                  t_0,
                  n_t_steps,
                  t_step):

    # This function solves numerically the modal ODE

    def deriv(x, t):
        z_p = x[1]
        x_p = (-modal_damp / modal_mass) * x[1] \
              + (-modal_stiff / modal_mass) * x[0] \
              + (force_mag / modal_mass) * math.cos(
            force_freq * t + force_phase)

        return [z_p, x_p]

    U_0 = [0, 0]

    time = np.array([t_0 + i * t_step for i in range(n_t_steps + 1)])
    Us = odeint(deriv, U_0, time)
    x = Us[:, 0]

    return x
