import matplotlib.pyplot as plt;
import numpy as np
from numpy import matlib
from FDTD2D.DielectricGenerator2D import createShapes as cs

def PhotonicGrid(xCell, yCell, cellSize, epsilon):
    [Nx, Ny] = cellSize;
    r = int(Nx/3);
    cell = cs.createCircle(Nx, Ny, epsilon, r);
    grid = np.matlib.repmat(cell, xCell, yCell);

    return grid;

def PhotonicWaveguide():
    return None;
#
# xCell = 2; yCell = 2;
# cellSize = [50,50]; epsilon = 12;
#
# [Nx, Ny] = cellSize;
# r = int(Nx / 3);
# cell = cs.createCircle(Nx, Ny, epsilon, r);
# grid = np.matlib.repmat(cell, xCell, yCell);
# grid = PhotonicGrid(xCell, yCell, cellSize, 12);
# plt.imshow(grid)
# plt.show()
