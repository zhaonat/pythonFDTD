import numpy as np
import matplotlib.pyplot as plt

def runTE(tsteps, dt, dx, dy, eps, Ez, Hx, Hy, source):
    ## courant constant
    plt.ion()
    fieldAmplitudes = list();
    plt.figure;
    storedFields = list();
    for t in range(tsteps):
        print('t= ' + str(t))
        # update Hz, Ex, Ey
        # remember the yee grid and integer indexing
        ## current source needs to be scaled
        J = np.sin(2 * np.pi * t / 30) * (dt);

        Hx[-1, :] = 0;
        Hy[:, -1] = 0;
        deriv_y = (np.roll(Hx, -1, axis=1) - Hx) / dy;
        deriv_x = (np.roll(Hy, -1, axis=0) - Hy) / dx;
        Ez = Ez + np.multiply((dt / eps), (deriv_x - \
                                           (deriv_y)));
        Ez[source[0], source[1]] -= J / eps[
            source[0], source[1]];  # location of the source            # np.roll(Hx,-1)-Hx;
        ##vectorize this code
        Ez[:, 0] = 0;
        # deriv_y = dey/dy;
        deriv_y = (Ez - np.roll(Ez, 1, axis=1)) / dy;
        Ez[0, :] = 0;
        # derix_x = dex;
        deriv_x = (Ez - np.roll(Ez, 1, axis=0)) / dx;
        Hx = Hx - (dt) * deriv_y;
        Hy = Hy + (dt) * deriv_x;

        ## records field in probloc

        # insert point source
        if (t % 10 == 0):
            imgplot = plt.imshow(Hy)
            imgplot.set_cmap('jet')
            plt.pause(0.001)
            plt.clf()

    storedFields = np.array(storedFields);
    plt.show()
    return storedFields
