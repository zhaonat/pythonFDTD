import numpy as np
import matplotlib.pyplot as plt
from FDTD2D.DielectricGenerator2D import PhotonicCrystalGrid as pcg
from FDTD2D.Simulators import runTE

## fundamental parameter specificatoin (all 1)
c = 1; # we are in units of nanometers;
eps_0 = 1; mu_0 = 1;

lambda_0 = 800;
lambda_up = 1200;
lambda_low = 400;

dx = lambda_0/200; dy = dx;
dt = (1/2)**.5*dx/c;
xrange = [-400, 400];
yrange = [-400, 400];

xstart = np.array([-400,-400]); xend = np.array([400,400]);
source =  np.array([0, 0])-xstart; source = [int(i/dx) for i in source];
probloc = np.array([200,0])-xstart; probloc = [int(i/dx) for i in probloc];


Nx = int((xrange[1]-xrange[0])/dx);
Ny = int((yrange[1]-yrange[0])/dy);
## create epsilon

## structure specification

structureLength = 50;
structureWidth = 50;
structureSpacing = 0; numCells = 20; diel = 12;
vacuumSpace = int((Nx-10)/2);
xcells = 10; ycells = 10; cellSize = np.array([40,40])
eps = pcg.PhotonicGrid(xcells, ycells, cellSize, diel);
plt.imshow(eps)
plt.show()
[Nx, Ny] = eps.shape
print(eps.shape)
# plt.imshow(eps)
# plt.show()

##initialize fields (consider a sparse initializiation)
Npml = [10,10];
Hx = np.zeros((Nx, Ny));
Hy = np.zeros((Nx, Ny));
Ez = np.zeros((Nx, Ny));
#spatial steps
#time steps
tsteps = 1000;
timeRun = np.arange(1,1000+1)
J = np.sin(2*np.pi*timeRun/20);
def currSource(t, lambda_up, lambda_low, lambda_0, c, dx):
    omega_0 = 2*np.pi/lambda_0;
    sigma = (2/omega_0)*(lambda_0)/(lambda_up-lambda_low)
    #print(omega_0); print(sigma)
    return np.exp(-(t-4*sigma)**2/sigma**2)*np.sin(omega_0*t)*(10/dx)**2;

J = currSource(timeRun, lambda_up, lambda_low, 170, c, dx);
print(J)

storedFields = runTE.runTE(tsteps, dt, dx, dy, eps, Ez, Hx, Hy,J,Npml, sourcePos = [112,112])
