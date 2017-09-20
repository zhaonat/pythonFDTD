import numpy as np
import matplotlib.pyplot as plt
from FDTD2D.DielectricGenerator2D import GMRStructure as GMR
from FDTD2D.Simulators import runTE

c = 1; # we are in units of nanometers;
lambda_0 = 800;
lambda_up = 1200;
lambda_low = 400;
dx = lambda_0/80; dy = dx;
dt = (1/2)**.5*dx/c;
xrange = [-400, 400];
yrange = [-400, 400];
xstart = np.array([-400,-400]); xend = np.array([400,400]);
source =  np.array([-200, 0])-xstart; source = [int(i/dx) for i in source];
probloc = np.array([200,0])-xstart; probloc = [int(i/dx) for i in probloc];

Nx = int((xrange[1]-xrange[0])/dx);
Ny = int((yrange[1]-yrange[0])/dy);
## create epsilon
eps = GMR.GMRWaveguide(Nx, Ny,20,20, 10,10, 6);
[Nx, Ny] = eps.shape
print(eps.shape)
# plt.imshow(eps)
# plt.show()

##initialize fields (consider a sparse initializiation)
Hx = np.zeros((Nx, Ny));
Hy = np.zeros((Nx, Ny));
Ez = np.zeros((Nx, Ny));
#spatial steps
#time steps
tsteps = 1000;


storedFields = runTE.runTE(tsteps, dt, dx, dy, eps, Ez, Hx, Hy, source)
