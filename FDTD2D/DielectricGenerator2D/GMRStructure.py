import numpy as np
import matplotlib.pyplot as plt

def GMRWaveguide(Nx, Ny, structureLength, structureWidth, structureSpacing, numCells,diel):
    unitCellSpacing = structureLength + structureSpacing;

    structure = diel * np.ones((structureLength, structureWidth))
    spacing = np.ones((structureWidth, structureSpacing));
    unitcell = np.concatenate((structure, spacing), axis=1);
    GMR = unitcell
    for i in range(numCells):
        GMR = np.concatenate((GMR, unitcell), axis=1)

    ## concatenate vacuum on both sides
    print(GMR.shape)
    [Nx, ny] = GMR.shape;
    vacuum = np.ones((50, ny));
    print(vacuum.shape)
    simDomain = np.concatenate((vacuum, GMR, vacuum), axis=0)
    return simDomain


# Nx = 100;
# Ny  = 100;
#
# structureLength = 10;
# structureWidth = 10;
# structureSpacing = 10;
# diel = 12;
# eps = np.zeros((Nx,Ny))
#
# unitCellSpacing = structureLength+structureSpacing;
#
# structure = diel*np.ones((structureLength, structureWidth))
# spacing = np.ones((structureWidth, structureSpacing));
# unitcell = np.concatenate((structure, spacing), axis = 1);
# GMR = unitcell
# for i in range(int(Nx/unitCellSpacing)-1):
#     GMR = np.concatenate((GMR, unitcell), axis = 1)
#
# ## concatenate vacuum on both sides
# print(GMR.shape)
# [Nx, ny] = GMR.shape;
# vacuumedge = (Ny - structureWidth) / 2
# vacuum = np.ones((50,ny));
# print(vacuum.shape)
# simDomain = np.concatenate((vacuum, GMR, vacuum), axis = 0)
#
# plt.imshow(simDomain)
# plt.show()