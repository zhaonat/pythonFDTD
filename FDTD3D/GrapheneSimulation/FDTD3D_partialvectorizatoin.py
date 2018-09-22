import numpy as np
## testing of the 3D FDTD equations

import matplotlib.pyplot as plt
import numpy as np;

def curl(Vx, Vy, Vz, roll =1):
    zcomp = (np.roll(Vy,roll, axis=2) - Vy) - (np.roll(Vx,roll, axis = 1) - Vx);
    ycomp = (np.roll(Vx,roll, axis=0) - Vx) - (np.roll(Vz,roll, axis=2) - Vz);
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

## courant constant
plt.ion()
fieldAmplitudes = list();
plt.figure;
storedFields = list();
jx = int(Nx/2); jy = int(Ny/2); jz = int(Nz/2);
x,y,z = np.mgrid[0:Nx, 0:Ny, 0:Nz];

dL = [dx, dy, dz];
print('dt: '+str(dt)+' dx: '+str(dx))
tsteps = 500;
for t in range(tsteps):
    print('t= '+str(t))
    #update Hz, Ex, Ey
    #remember the yee grid and integer indexing
    ## current source needs to be scaled
    J = 1*np.sin(2*np.pi*t/30)

    CHX = np.zeros((Nx,Ny,Nz));
    CHY = np.zeros((Nx,Ny,Nz));
    CHZ = np.zeros((Nx,Ny,Nz));
    #t =t part of the time step, get the curls from t-1/2 for the H field
    for i in range(0,Nz):
        index = i+1;
        if(index >Nz-1):
            index = 0;
        Chx =  ((dt/dy)*(np.roll(Hz[i,:,:], 1, axis = 0) - Hz[i,:,:]) - (dt/dz)*(Hy[index, :,:]-Hy[i,:,:]));
        Chy = -((dt/dx)*(np.roll(Hz[i,:,:], 1, axis = 1) - Hz[i,:,:]) - (dt/dz)*(Hx[index, :,:]-Hx[i,:,:]));
        Chz = ((dt/dx)*(np.roll(Hy[i,:,:], 1, axis = 1) - Hy[i,:,:]) - (dt/dy)*(np.roll(Hx[i,:,:], 1, axis = 0)-Hx[i,:,:]));
        ## add dielectric
        # Chx = np.multiply(1/eps[i,:,:], Chx);
        # Chy = np.multiply(1/eps[i,:,:], Chy);
        # Chz = np.multiply(1/eps[i,:,:], Chz);
        CHX[i,:,:] = Chx;
        CHY[i,:,:] =  Chy;
        CHZ[i,:,:] =  Chz;

    Ex -= CHX;
    Ey -= CHY;
    Ez -= CHZ;
    Ez[jx, jy, jz] -=J;

    ## update H fields
    # # technically, we are t = t+1/2 here

    CEX = np.zeros((Nx,Ny,Nz));
    CEY = np.zeros((Nx,Ny,Nz));
    CEZ = np.zeros((Nx,Ny,Nz));
    for i in range(0,Nz):
        index = i-1;
        if(index <0):
            index = Nz-1;
        CEx =  ((dt/dy)*(Ez[i,:,:] - np.roll(Ez[i,:,:], -1, axis = 0)) - (dt/dz)*(Ey[i,:,:]- Ey[index,:,:]));
        CEy = -((dt/dx)*(Ez[i,:,:] - np.roll(Ez[i,:,:], -1, axis = 1)) - (dt/dz)*(Ex[i,:,:]- Ex[index,:,:]));
        CEz = ((dt/dx)*(Ey[i,:,:] - np.roll(Ey[i,:,:], -1, axis = 1)) - (dt/dy)*(Ex[i,:,:]- np.roll(Ex[i,:,:], -1, axis = 0)));

        CEX[i,:,:]= CEx;
        CEY[i,:,:] =  CEy;
        CEZ[i,:,:] =  CEz;

    Hx += CEX;
    Hy += CEY;
    Hz += CEZ;

    # Hx +=  ((dt/dy)*(Ez - np.roll(Ez, -1, axis = 1) ) - (dt/dz)*(Ey- np.roll(Ey, -1, axis = 0)));
    # Hy +=  -((dt/dx)*(Ez - np.roll(Ez, -1, axis = 2) ) - (dt/dz)*(Ex- np.roll(Ex, -1, axis = 0)));
    # Hz +=  ((dt/dx)*(Ey - np.roll(Ey, -1, axis = 2) ) - (dt/dy)*(Ex- np.roll(Ex, -1, axis = 1)));


    #insert point source
    if(t%2 == 0):
        imgplot = plt.imshow(((Ez[:,:,jz])), cmap = 'jet')
        print(np.max(Ez))
        #plt.plot(Ez[:,jy, jz])
        plt.clim(-0.01, 0.01)
        #plt.plot(Ex[:,:,jz])
        #imgplot = mlab.contour3d(x,y,z,Hz)
        plt.pause(0.001)
        #print(np.max(Ez))
        #print('current source power: '+str(J))
        #print(Ez[jx, jy, jz])
        plt.clf()

