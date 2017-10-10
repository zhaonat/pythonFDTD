import numpy as np

import numpy as np
from mayavi import mlab;
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
tsteps = 300;

tau = 0.1;
sigma_s = 1;
denom = dx*sigma_s+dt+2*tau;
## coefs of the sheet
c1 = (dt+2*tau)/(dx*sigma_s+dt+2*tau);
c2 = (dt+2*tau)/(dx*sigma_s+dt+2*tau);

fe1 = sigma_s*dt/(denom);
fe2 = sigma_s*dt/(denom);
fh11 = 1-2*dt/(denom); fh22 = 1-2*dt/(denom);
fh12 = (dt-2*tau)/denom; fh21 = (dt-2*tau)/denom;


## courant constant
plt.ion()
fieldAmplitudes = list();
plt.figure;
storedFields = list();
jx = int(Nx/2); jy = int(Ny/2); jz = int(Nz/2);
x,y,z = np.mgrid[0:Nx, 0:Ny, 0:Nz];

dL = [dx, dy, dz];
print('dt: '+str(dt)+' dx: '+str(dx))

Hx_sheet = np.zeros((Nx,Ny,2)); #third dimension denotes index, orwhich side of the sheet
Hy_sheet = np.zeros((Nx,Ny,2));
Ez_sheet = np.zeros((Nx,Ny,2));

sheet_loc_z = round(Nz/2);


for t in range(tsteps):
    print('t= '+str(t))
    #t =t part of the time step, get the curls from t-1/2 for the H field
    ## any Hy or Hx that appears will have a special update with the interface
    for i in range(Nz):
        if(i == sheet_loc_z):
            ## define F coefs for the Hx sheet
            # define derivatives of the sheet fields
            deriv_y_1 = np.roll(Ez_sheet[:, :, 0], 1, axis=1) - Ez_sheet[:, :, 0];
            deriv_y_2 = np.roll(Ez_sheet[:, :, 1], 1, axis=1) - Ez_sheet[:, :, 1];

            F1n_Hx = fh11 * (Hx_sheet[:, :, 0] + fh12 * Hx_sheet[:, :, 1]) - \
                     fe1 * (2 * Ey[:, :, sheet_loc_z]) - dz * deriv_y_1;
            F2n_Hx = fh21 * (Hx_sheet[:, :, 0] + fh22 * Hx_sheet[:, :, 1]) - \
                     fe2 * (2 * Ey[:, :, sheet_loc_z]) - dz * deriv_y_2;

            F1n_Hy = fh11 * (Hy_sheet[:, :, 0] + fh12 * Hy_sheet[:, :, 1]) - \
                     fe1 * (2 * Ey[:, :, sheet_loc_z]) - dz * deriv_y_1;
            F2n_Hy = fh21 * (Hy_sheet[:, :, 0] + fh22 * Hy_sheet[:, :, 1]) - \
                     fe2 * (2 * Ey[:, :, sheet_loc_z]) - dz * deriv_y_2 \
                ;
            Hx_sheet[:, :, 0] = (1 / (1 - c1 * c2)) * (F1n_Hx + c1 * F2n_Hx);
            Hx_sheet[:, :, 1] = (1 / (1 - c1 * c2)) * (F2n_Hx + c2 * F1n_Hx);  # use Fx_sheet to update the E fields

            Hy_sheet[:, :, 0] = (1 / (1 - c1 * c2)) * (F1n_Hy + c1 * F2n_Hy);
            Hy_sheet[:, :, 1] = (1 / (1 - c1 * c2)) * (F2n_Hy + c2 * F1n_Hy);  # use Fx_sheet to update the E fields
        else:
            Chx =  ((dt/dy)*(np.roll(Hz, 1, axis = 1) - Hz) - (dt/dz)*(np.roll(Hy, 1, axis = 0)-Hy));
            Chy = -((dt/dx)*(np.roll(Hz, 1, axis = 2) - Hz) - (dt/dz)*(np.roll(Hx,1, axis = 0)-Hx));
            Chz = ((dt/dx)*(np.roll(Hy, 1, axis = 2) - Hy) - (dt/dy)*(np.roll(Hx, 1, axis = 1)-Hx));
            ## add dielectric
            Chx = np.multiply(1/eps, Chx);
            Chy = np.multiply(1/eps, Chy);
            Chz = np.multiply(1/eps, Chz);

            Ex -= Chx;
            Ey -=  Chy;
            Ez -=  Chz
            J = 1*np.sin(2*np.pi*t/30)

        Ez[jx, jy, jz] -=J;

    ## update H fields with curl of Ez
    # technically, we are t = t+1/2 here
    Hx +=  ((dt/dy)*(Ez - np.roll(Ez, -1, axis = 1) ) - (dt/dz)*(Ey- np.roll(Ey, -1, axis = 0)));
    Hy +=  -((dt/dx)*(Ez - np.roll(Ez, -1, axis = 2) ) - (dt/dz)*(Ex- np.roll(Ex, -1, axis = 0)));
    Hz +=  ((dt/dx)*(Ey - np.roll(Ey, -1, axis = 2) ) - (dt/dy)*(Ex- np.roll(Ex, -1, axis = 1)));


    ## update fields connected to interface


    ## Time step the Ez sheet fields using the updated Hy and Hx sheet fields
    #calculate derivatives of the sheet fields for Ez update
    deriv_Hy_x = np.roll(Hy_sheet[:,:,0],1,axis = 0) - Hy_sheet[:,:,0];
    deriv_Hx_y = np.roll(Hx_sheet[:,:,0],1,axis = 1) - Hx_sheet[:,:,0];
    deriv_Hy_x2 = np.roll(Hx_sheet[:,:,1],1,axis = 1) - Hx_sheet[:,:,1];
    deriv_Hx_y2 = np.roll(Hx_sheet[:,:,1],1,axis = 1) - Hx_sheet[:,:,1];

    Ez_sheet_next_1 = (dt/eps)*(deriv_Hx_y-deriv_Hy_x)
    Ez_sheet_next_2 = (dt/eps)*(deriv_Hx_y2-deriv_Hy_x2)

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

