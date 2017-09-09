import numpy as np;

def curl(Vx, Vy, Vz, roll =1):
    #dxVy - dyVx
    zcomp = (np.roll(Vy,roll, axis=2) - Vy) - (np.roll(Vx,roll, axis = 1) - Vx);
    #dzVx - dxVz
    ycomp = (np.roll(Vx,roll, axis=0) - Vx) - (np.roll(Vz,roll, axis=2) - Vz);
    #dxVy - dyVx
    xcomp = (np.roll(Vz,roll, axis=1) - Vz) - (np.roll(Vy,roll,axis=0) - Vy);
    return [xcomp, ycomp, zcomp];
