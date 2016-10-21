from classes import *
import numpy as np

def change_first_number(matrix):
    matrix[0][0] = 5
    matrix[1][1] = 6
    matrix[2][2] = 7


a = np.array([[1, 2, 3],
              [1, 2, 3],
              [1, 2, 3]])

print(a)
change_first_number(a)
print(a)

b = a[:,1].reshape(1,-1).T

print(b)

Node(1, 2, 3)
b = Node(4, 5, 6)
a = Node.instances
c = a.index(b)
print("index")
print(c)