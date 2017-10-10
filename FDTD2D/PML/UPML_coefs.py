import numpy as np
import numpy.matlib

def UPML(N, Npml):
    Nx = N[0]; Ny = N[1];
    xn = np.zeros((Nx, 1))
    fi1 = np.zeros((Nx, 1));
    fj1 = np.zeros((Ny, 1));
    gi2 = np.ones((Nx, 1));
    gj2 = np.ones((Ny, 1));
    gi3 = np.ones((Nx, 1));
    gj3 = np.ones((Ny, 1));
    Npmlx = Npml[0];
    Npmly = Npml[1];

    for i in range(Npmlx):
        x_n = (1 / 3) * ((i + 1) / Npmlx) ** 3;
        xn[Nx - Npmlx + i] = x_n;
        xn[Npmlx - 1 - i] = x_n;

        fi1[Npmlx - 1 - i] = x_n;
        gi2[Npmlx - 1 - i] = 1 / (1 + x_n);
        gi3[Npmlx - 1 - i] = (1 - x_n) / (1 + x_n);

        fi1[Nx - Npmlx + i] = x_n;
        gi2[Nx - Npmlx + i] = 1 / (1 + x_n);
        gi3[Nx - Npmlx + i] = (1 - x_n) / (1 + x_n);

    for i in range(Npmly):
        x_n = (1 / 3) * ((i + 1) / Npml[1]) ** 3;
        fj1[Npmly - 1 - i] = x_n;
        gj2[Npmly - 1 - i] = 1 / (1 + x_n);
        gj3[Npmly - 1 - i] = (1 - x_n) / (1 + x_n);

        fj1[Ny - Npmly + i] = x_n;
        gj2[Ny - Npmly + i] = 1 / (1 + x_n);
        gj3[Ny - Npmly + i] = (1 - x_n) / (1 + x_n);

    gi3 = np.squeeze(gi3);
    gi2 = np.squeeze(gi2);
    gj3 = np.squeeze(gj3);
    gj2 = np.squeeze(gj2);

    ## creat the other f's from teh g's
    fi2 = gi2;
    fj2 = gj2;
    fj3 = gj3;
    fi3 = gi3;

    fi1 = np.matlib.repmat(fi1.T, Ny, 1);
    fj1 = np.matlib.repmat(fj1, 1, Nx);

    fi2 = np.matlib.repmat(fi2.T, Ny, 1);
    fi2 = np.rot90(fi2);
    fj2 = np.matlib.repmat(fj2.T, Nx, 1);

    fi3 = np.matlib.repmat(fi3, Ny, 1);
    fi3 = np.rot90(fi3)
    fj3 = np.matlib.repmat(fj3.T, Nx, 1);

    ## at this point all the fi's and fj's only have the coeffs on one side
    fij3 = np.multiply(fj3, fi3);
    fij2 = np.multiply(fj2, fi2);
    fij1 = np.multiply(fj1, fi1);
    fij1 = np.rot90(fij1)

    gi3x = np.matlib.repmat(gi3.T, Ny, 1);
    gi3x = np.rot90(gi3x);
    gj3x = np.matlib.repmat(gj3, Nx, 1);
    gij3 = np.multiply(gi3x, gj3x);

    gi2x = np.matlib.repmat(gi2.T, Ny, 1);
    gi2x = np.rot90(gi2x)
    gj2x = np.matlib.repmat(gj2, Nx, 1);
    gij2 = np.multiply(gi2x, gj2x);

    return [fij1, fij2, fij3, gij2, gij3];