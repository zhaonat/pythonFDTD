import numpy as np
## testing of the 3D FDTD equations

import matplotlib.pyplot as plt
import numpy as np;

def curl(Vx, Vy, Vz, roll =1):
    #dxVy - dyVx
    zcomp = (np.roll(Vy,roll, axis=2) - Vy) - (np.roll(Vx,roll, axis = 1) - Vx);
    #dzVx - dxVz
    ycomp = (np.roll(Vx,roll, axis=0) - Vx) - (np.roll(Vz,roll, axis=2) - Vz);
    #dxVy - dyVx
    xcomp = (np.roll(Vz,roll, axis=1) - Vz) - (np.roll(Vy,roll,axis=0) - Vy);
    return [xcomp, ycomp, zcomp];


## test the curl
def currSource(t, lambda_up, lambda_low, lambda_0, c, dx):
    omega_0 = 2*np.pi/lambda_0;
    sigma = (2/omega_0)*(lambda_0)/(lambda_up-lambda_low)
    #print(omega_0); print(sigma)
    return np.exp(-(t-4*sigma)**2/sigma**2)*np.sin(omega_0*t)*(10/dx)**2;


c = 1; # we are in units of nanometers;
lambda_0 = 800;
lambda_up = 1200;
lambda_low = 400;
dx = lambda_0/100; dy = dx; dz = dy;
dt = (1/3)**.5*dx/c;
xrange = [-400, 400];
yrange = [-400, 400];
zrange = [-400,400]
xstart = np.array([-400,-400]); xend = np.array([400,400]);
source =  np.array([-200, 0])-xstart; source = [int(i/dx) for i in source];
probloc = np.array([200,0])-xstart; probloc = [int(i/dx) for i in probloc];

Nx = int((xrange[1]-xrange[0])/dx);
Ny = int((yrange[1]-yrange[0])/dy);
Nz = int((yrange[1]-yrange[0])/dy);

## create epsilon
eps = np.ones((Nx, Ny, Nz));
#eps[50:150, 50:150, 50:150] = 4;


##initialize fields (consider a sparse initializiation)
Hx = np.zeros((Nx, Ny, Nz));
Hy = np.zeros((Nx, Ny, Nz));
Ez = np.zeros((Nx, Ny, Nz));
Ex = np.zeros((Nx, Ny, Nz));
Ey = np.zeros((Nx, Ny, Nz));
Hz = np.zeros((Nx, Ny, Nz));
#spatial steps
#time steps
tsteps = 500;

## courant constant
plt.ion()
fieldAmplitudes = list();
plt.figure;
storedFields = list();
jx = int(Nx/2); jy = int(Ny/2); jz = int(Nz/2);
x,y,z = np.mgrid[0:Nx, 0:Ny, 0:Nz];

dL = [dx, dy, dz];
print('dt: '+str(dt)+' dx: '+str(dx))
for t in range(tsteps):
    print('t= '+str(t))
    #update Hz, Ex, Ey
    #remember the yee grid and integer indexing
    ## current source needs to be scaled
    J = 1*np.sin(2*np.pi*t/30)
    #J = 0.5;
    ## we have six different components to update
    # 0 = dz, 1 = dy, 2 = dx

    # [chx, chy, chz] = curl(Hx, Hy, Hz, roll = -1);
    # #chx =  Chx(dL, Hz, Hy); chy =  Chy(dL, Hx, Hz);chz =  Chz(dL, Hx, Hy);
    #
    # # print(np.max(chx))
    # ## algorithm is unstable right now...somehow the source gets amplified as we do the different rolls
    # print('before: '+str(np.max(Ex)))
    # Ex -= (dt/dx)*(chx)
    # #print('after: '+str(np.max(Ex)))
    # Ey -= (dt/dy)*(chy)
    # Ez -= (dt/dz)*(chz)
    # Ez[jx,jy,jz] -= J;
    # #now at the 1/2 time step
    # [cex, cey, cez] = curl(Ex, Ey, Ez, roll = 1);
    #
    # #cex =  Cex(dL, Ez, Ey); cey =  Cey(dL, Ex, Ez);cez =  Cez(dL, Ex, Ey);
    #
    # Hx -= (dt/dx)*(cex)
    # Hy -= (dt/dy)*(cey)
    # Hz -= (dt/dz)*(cez)
    # print(np.max(Hz))

    #t =t part of the time step, get the curls from t-1/2 for the H field
    Chx =  ((dt/dy)*(np.roll(Hz, 1, axis = 1) - Hz) - (dt/dz)*(np.roll(Hy, 1, axis = 0)-Hy));
    Chy = -((dt/dx)*(np.roll(Hz, 1, axis = 2) - Hz) - (dt/dz)*(np.roll(Hx,1, axis = 0)-Hx));
    Chz = ((dt/dx)*(np.roll(Hy, 1, axis = 2) - Hy) - (dt/dy)*(np.roll(Hx, 1, axis = 1)-Hx));
    ## add dielectric
    Chx = np.multiply(1/eps, Chx);
    Chy = np.multiply(1/eps, Chy);
    Chz = np.multiply(1/eps, Chz);

    print(np.max(Chx))
    print(np.max(Ex))
    Ex -= Chx;
    print(np.max(Ex))
    Ey -=  Chy;
    Ez -=  Chz
    Ez[jx, jy, jz] -=J;
    # technically, we are t = t+1/2 here
    Hx +=  ((dt/dy)*(Ez - np.roll(Ez, -1, axis = 1) ) - (dt/dz)*(Ey- np.roll(Ey, -1, axis = 0)));
    Hy +=  -((dt/dx)*(Ez - np.roll(Ez, -1, axis = 2) ) - (dt/dz)*(Ex- np.roll(Ex, -1, axis = 0)));
    Hz +=  ((dt/dx)*(Ey - np.roll(Ey, -1, axis = 2) ) - (dt/dy)*(Ex- np.roll(Ex, -1, axis = 1)));


    #insert point source
    if(t%2 == 0):
        imgplot = plt.imshow(Ez[:,:,jz], cmap = 'jet')
        plt.clim(-0.01, 0.01)
        #plt.plot(Ex[:,:,jz])
        #imgplot = mlab.contour3d(x,y,z,Hz)
        plt.pause(0.001)
        #print(np.max(Ez))
        #print('current source power: '+str(J))
        #print(Ez[jx, jy, jz])
        plt.clf()

