from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math

# x = np.array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
#               [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
#               [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
#               [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
#               [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]])
#
# y = np.random.rand(5, 11)
#
#
# fig = plt.figure()
# ax = fig.add_subplot(111, autoscale_on=False, xlim=(0, 11), ylim=(0, 1))
# ax.grid()
#
# line, = ax.plot([], [], 'o-', lw=2)
# time_template = 'time = %.1fs'
# time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
#
# def init():
#     line.set_data([], [])
#     time_text.set_text('')
#     return line, time_text
#
#
# def animate(i):
#     thisx = x[i]
#     thisy = y[i]
#
#     line.set_data(thisx, thisy)
#     time_text.set_text(time_template % (i))
#     return line, time_text
#
# ani = animation.FuncAnimation(fig, animate, np.arange(1, len(y)),
#                               interval=200, blit=True, init_func=init)
#
# plt.show()

t = np.array([[0.1],
              [0.1],
              [0.1]])

freq = np.array([[3],
                 [2],
                 [7]])

phase = np.array([[math.pi / 2],
                  [math.pi / 5],
                  [math.pi / 6]])

mag = np.array([[1000],
                [2000],
                [200]])

force = mag * np.sin(freq * t + phase)

force_mag = force[:, 0]
print(force_mag)

f = float("He")

print("\nReading input file ...")
star_time = time.clock()
print("Time used: " + str(round(time.clock() - star_time, 4)) + "s")