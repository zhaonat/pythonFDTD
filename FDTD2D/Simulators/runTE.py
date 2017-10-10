import numpy as np
import matplotlib.pyplot as plt
from FDTD2D.PML import UPML_coefs as umpl

def runTE(tsteps, dt, dx, dy, eps, Ez, Hx, Hy, J, Npml, sourcePos = [50,50]):
    ## courant constant
    plt.ion()
    fieldAmplitudes = list();
    plt.figure(figsize = (25, 20));
    storedFields = list();
    N = Ez.shape;
    [fij1, fij2, fij3, gij2, gij3] = umpl.UPML(N, Npml);

    Icex = 0; Icey = 0; Idz = 0;
    for t in range(tsteps):
        print('t= ' + str(t))

        ## now work on Dz
        # Hx[-1,:] = 0;
        # Hy[:,-1] = 0;
        derivy = (np.roll(Hx, -1, axis=1) - Hx) / dy;
        derivx = (np.roll(Hy, -1, axis=0) - Hy) / dx;
        Chz = derivx - derivy;

        Dz = eps * Ez;
        Idz += Dz;
        Dz = gij3 * Dz + gij2 * dt * Chz;


        Dz[sourcePos[0], sourcePos[1]] -= J[t];  # location of the source

        ## finally update Ez
        Ez = Dz / eps;  ## this is important for the next update of the curl # equations

        ## update H fields
        # Ez[:,0] = 0;
        # Ez[0,:] = 0;
        Cex = -(Ez - np.roll(Ez, 1, axis=1)) / dy;
        Cey = +(Ez - np.roll(Ez, 1, axis=0)) / dx;
        Icex += Cex;
        Icey += Cey;
        Hx = fij3 * Hx + (dt) * fij2 * Cex + fij1 * Icex;
        Hy = fij3 * Hy + (dt) * fij2 * Cey + fij1 * Icey;


        # insert point source
        if (t % 10 == 0):
            storedFields.append([Ez, Hx, Hy])
            imgplot = plt.imshow(Ez)
            imgplot.set_cmap('jet')
            plt.clim(-0.04, 0.04)
            plt.pause(0.001)
            plt.clf()

    plt.show()
    return storedFields
