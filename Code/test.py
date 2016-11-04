from functions import damp_rayleigh
from functions import damp_bismarck
import numpy as np
import scipy as sc

m = np.array([[400, 0, 0],
              [0, 400, 0],
              [0, 0, 200]])

m = (1/386) * m

k = np.array([[2, -1, 0],
              [-1, 2, -1],
              [0, -1, 1]])

k = 610 * k

eigvalues, eigvectors = sc.linalg.eig(k, m)

freq = abs(np.sqrt(eigvalues))

print(freq)

damp_matrix, damp_vector =  damp_rayleigh(m, k, freq)

print(damp_matrix)
print(damp_vector)

damp_matrix, damp_vector =  damp_bismarck(m, k, freq)

print(damp_matrix)
print(damp_vector)