import numpy as np

a = np.array([[1, 2, 3],
              [1, 2, 3],
              [1, 2, 3]])

b = a[:,1].reshape(1,-1).T

print(b)