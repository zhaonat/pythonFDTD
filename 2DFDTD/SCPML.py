## first, we need to be able to visualize a non-uniform grid
import numpy as np;
import matplotlib.pyplot as plt

## want to apply the pml as a matrix operation on a matrix HX of nxmxl
def sfactor(dPML, N, omega_0, m = 4):
    lnR = -16;  ## final desired reflection coefficient
    PmlBound = [dPML, N - dPML];
    sigma_max = (m + 1) * lnR / 2
    # need start and stop coordinates of the PML
    # abs(wrange - wpml);
    sfactor_domain = np.ones((N - 2 * dPML))
    d = np.linspace(0, dPML, dPML)
    pmlPolynomial = 1/(1 + (1j) * (sigma_max * d / omega_0) ** m);
    sfactor_domain = np.concatenate((np.flip(pmlPolynomial, 0),\
                                     sfactor_domain, pmlPolynomial), axis=0);

    sfactorMatrix = np.repeat(sfactor_domain, N);
    sfactorMatrix = np.reshape(sfactorMatrix, [N,N])
    return sfactorMatrix;

## EMLAB calculation of PML
Nx = 20; Ny = 20;

Npml = [5,5];
Nx2 = Nx+2*Npml[0];
Ny2 = Ny+2*Npml[1];
dt = 0.1;

## x-direction pml
## no 2x grid concept
sigx = np.zeros((Nx2, Ny2));
#left hand side (need negative values)
for nx in range(Npml[0]):
    nx1 = Npml[0]-nx-1;
    sigx[nx1, :] = (0.5/dt)*(nx/2/Npml[0])**3;
#right hand side
for nx in range(Npml[1]):
    nx1 = Nx2-Npml[1]+nx;
    sigx[nx1,:] = 0.5/dt*(nx/2/Npml[1])**3

## vertical y direction pml
sigy = np.zeros((Nx2, Ny2))
for ny in range(Npml[0]):
    ny1 = Npml[0] - ny -1 ;
    sigy[:, ny1] = (0.5/dt)*(ny/2/Npml[0])**3;
for ny in range(Npml[1]):
    ny1 = Ny2 - Npml[1] +ny;
    sigy[:, ny1] = (0.5/dt)*(ny/2/Npml[1])**3;

plt.imshow(sigx)
plt.hold(True)
plt.imshow(sigy)

import numpy as np
import matplotlib.pyplot as plt
import numpy as np;


def currSource(t, lambda_up, lambda_low, lambda_0, c, dx):
    omega_0 = 2*np.pi/lambda_0;
    sigma = (2/omega_0)*(lambda_0)/(lambda_up-lambda_low)
    #print(omega_0); print(sigma)
    return np.exp(-(t-4*sigma)**2/sigma**2)*np.sin(omega_0*t)*(10/dx)**2;


c = 1; # we are in units of nanometers;
lambda_0 = 800; omega_0 = 2*np.pi/lambda_0;
lambda_up = 1200;
lambda_low = 400;
dx = lambda_0/80; dy = dx;
dt = (1/2)**.5*dx/c;
xrange = [-400, 400];
yrange = [-400, 400];
xstart = np.array([-400,-400]); xend = np.array([400,400]);
source =  np.array([0,0])-xstart; source = [int(i/dx) for i in source];
probloc = np.array([200,0])-xstart; probloc = [int(i/dx) for i in probloc];

Nx = int((xrange[1]-xrange[0])/dx);
Ny = int((yrange[1]-yrange[0])/dy);
## create epsilon
eps = np.ones((Nx, Ny));
#eps[50:150, 50:150] = 4;

##initialize fields (consider a sparse initializiation)
Hx = np.zeros((Nx, Ny));
Hy = np.zeros((Nx, Ny));
Ez = np.zeros((Nx, Ny));
Npml = 5;
#spatial steps
#time steps
tsteps = 500;

## courant constant
plt.ion()
fieldAmplitudes = list();
plt.figure;
storedFields = list();
sfactormatrix = sfactor(Npml,Nx, omega_0, m = 4);
sx = sfactormatrix;
sy = np.flip(sfactormatrix, axis = 1);

for t in range(tsteps):
    print('t= '+str(t))
    #update Hz, Ex, Ey
    #remember the yee grid and integer indexing
    ## current source needs to be scaled
    J = np.sin(2*np.pi*t/20)*(1/dx)**2;
    #J = currSource((t+0.5)*dt, lambda_up, lambda_low, lambda_0, c, dx);
    #update E field components
        # TEz mode...only transverse x and y E field
    Hx[-1,:] = 0;
    Hy[:,-1] = 0;
    deriv_y = np.multiply(sy,(np.roll(Hx,-1, axis=1)-Hx)/dy);
    deriv_x = np.multiply(sx, (np.roll(Hy, -1, axis = 0)-Hy)/dx);
    Ez = Ez + np.multiply((dt/eps), (deriv_x - \
                                  (deriv_y)));
    ## enforce Dirichlet

    Ez[source[0],source[1]]-=J; #location of the source            # np.roll(Hx,-1)-Hx;
    ##vectorize this code
    Ez[:,0] = 0;
    deriv_y = np.multiply(sy,(Ez - np.roll(Ez,1,axis=1 ))/dy);
    Ez[0,:] = 0;
    deriv_x = np.multiply(sx,(Ez - np.roll(Ez,1, axis = 0))/dx);
    Hx = Hx - (dt) * deriv_y;
    Hy = Hy + (dt) * deriv_x ;

    ## records field in probloc
    storedFields.append(Ez[probloc[0], probloc[1]]);

    #insert point source
    if(t%10 == 0):
        imgplot = plt.imshow(np.real(Hy))
        imgplot.set_cmap('jet')
        plt.pause(0.001)
        plt.clf()

storedFields = np.array(storedFields);
FFT = np.fft.fft(storedFields)
plt.plot(np.abs(FFT));
plt.show()
