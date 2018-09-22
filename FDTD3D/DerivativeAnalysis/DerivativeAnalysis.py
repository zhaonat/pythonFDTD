import numpy as np

A = list(range(100));
A = np.reshape(A, (10,10))

dy = np.roll(A, 1, axis = 1);
dx = np.roll(A,0, axis = 0);