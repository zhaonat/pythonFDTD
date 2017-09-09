import numpy as np

def createCircle(Nx, Ny, epsilon, r):
    eps = np.ones((Nx, Ny));
    # Create index arrays to z
    I, J = np.meshgrid(np.arange(eps.shape[0]),\
                       np.arange(eps.shape[1]))

    # calculate distance of all points to centre
    dist = np.sqrt((I - Nx/2 + 1) ** 2 + (J - Ny/2 + 1) ** 2)
    eps[np.where(dist < r**2)] = epsilon;
    return eps;