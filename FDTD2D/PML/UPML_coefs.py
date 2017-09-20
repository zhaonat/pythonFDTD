import numpy as np

def UPML(N, Npml):
    Nx = N[0]; Ny = N[1];
    xn = np.zeros((Nx, 1))
    fi1 = np.zeros((Nx, 1))
    gi2 = np.ones((Nx, 1));
    gj2 = np.ones((Ny, 1));
    gi3 = np.ones((Nx, 1));
    gj3 = np.ones((Ny, 1));
    Npmlx = Npml[0];
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

    for i in range(Npml[1]):
        x_n = (1 / 3) * ((i + 1) / Npml[1]) ** 3;
        gj2[Npmlx - 1 - i] = 1 / (1 + x_n);
        gj3[Npmlx - 1 - i] = (1 - x_n) / (1 + x_n);

        gj2[Nx - Npmlx + i] = 1 / (1 + x_n);
        gj3[Nx - Npmlx + i] = (1 - x_n) / (1 + x_n);

    fj1 = fi1;

    gi3 = np.squeeze(gi3);
    gi2 = np.squeeze(gi2);
    gj3 = np.squeeze(gj3);
    gj2 = np.squeeze(gj2);

    fi2 = gi2;
    fj2 = gi2;
    fj3 = gj3;
    fi3 = gi3;

    fi1 = np.matlib.repmat(fi1.T, len(fi1), 1);
    fj1 = np.rot90(fi1);

    fi2 = np.matlib.repmat(fi2.T, len(fi2), 1);
    fj2 = np.rot90(fi2);

    fi3 = np.matlib.repmat(fi3.T, len(fi3), 1);
    fj3 = np.rot90(fi3);

    ## at this point all the fi's and fj's only have the coeffs on one side
    fij3 = fj3 * fi3;
    fij2 = fj2 * fi2;
    fij1 = fj1 * fi1;

    gi3x = np.matlib.repmat(gi3.T, len(gi3), 1);
    gj3x = np.matlib.repmat(gj3, len(gj3), 1);
    gj3x = np.rot90(gj3x);
    gij3 = np.multiply(gi3x, gj3x);

    gi2x = np.matlib.repmat(gi2.T, len(gi2), 1);
    gj2x = np.matlib.repmat(gj2, len(gj2), 1);
    gj2x = np.rot90(gj2x);
    gij2 = np.multiply(gi2x, gj2x);

    return [fij1, fij2, fij3, gij2, gij3];