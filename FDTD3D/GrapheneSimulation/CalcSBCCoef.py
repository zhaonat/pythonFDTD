import numpy as np
## we want to encode a 2D grid defining where the SBC is and is not

def SBCCoefs3D():
    return None;

def squareHole(Nx, Ny):
    ## returns a grid containing the structure of the 2D sheet
    sheet = np.ones((Nx,Ny));
    center = [round(Nx/2), round(Ny/2)]
