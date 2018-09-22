import numpy as np

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


tau = 0.1;
sigma_0 = 1;
## Hx sheet update coefficients
c1 = (dt) / (dz * sigma_0 + dt);
c2 = (dt) / (dz * sigma_0 + dt); ## if this is the courant number, there is instability introduced

fe1 = sigma_0 * c1;
fe2 = sigma_0 * c2;
fh11 = 1-2.0*c1;
fh22 = 1-2.0*c2
fh12 = c1;
fh21 = c2;


## Hy field coefficients on the sbc sheet
c1y = (dt) / (dz * sigma_0 - dt);
c2y = (dt) / (dz * sigma_0 - dt); ## if this is the courant number, there is instability introduced

fe1y = sigma_0 * c1y;
fe2y = sigma_0 * c2y;
fh11y = 1-2.0*c1y;
fh22y = 1-2.0*c2y
fh12y = c1y;
fh21y = c2y;



## courant constant
plt.ion()
fieldAmplitudes = list();
plt.figure(figsize = (20,20));
storedFields = list();
jx = int(Nx/4); jy = int(Ny/2); jz = int(Nz/4);
x,y,z = np.mgrid[0:Nx, 0:Ny, 0:Nz];

dL = [dx, dy, dz];
print('dt: '+str(dt)+' dx: '+str(dx))

## We have to define a special yee cell which contains the conducting surface
Hx_sheet = np.zeros((Nx,Ny,2)); #third dimension denotes index, orwhich side of the sheet
Hy_sheet = np.zeros((Nx,Ny,2));
Ez_sheet = np.zeros((Nx,Ny,2));
Ex_sheet = np.zeros((Nx,Ny));
Ey_sheet = np.zeros((Nx,Ny));
Hz_sheet = np.zeros((Nz,Ny));

sheet_loc_z = round(Nz/2);

#time steps
tsteps = 100;
StorageSheet = list();
for t in range(tsteps):
    print('t= '+str(t))
    #t =t part of the time step, get the curls from t-1/2 for the H field
    ## any Hy or Hx that appears will have a special update with the interface

    #############=========== E FIELD UPDATES =====================#################

    CHX = np.zeros((Nx,Ny,Nz));
    CHY = np.zeros((Nx,Ny,Nz));
    CHZ = np.zeros((Nx,Ny,Nz));
    for i in range(Nz):
        index = i + 1;
        if (index > Nz - 1):
            index = 0;
        if(i == sheet_loc_z): ## at the yee cell right before the surface cell
            Hy_int_bottom = Hy_sheet[:,:,0]; # these are updated in the H loop
            Hx_int_bottom = Hx_sheet[:,:,0];
            ## The border yee cell couples into the Hy, Hx, components of the cell above it in the z derivatives
            Chx = ((dt / dy) * (np.roll(Hz[i, :, :], 1, axis=0) - Hz[i, :, :]) - (dt / dz) * (Hy_int_bottom - Hy[i, :, :]));
            Chy = -((dt / dx) * (np.roll(Hz[i, :, :], 1, axis=1) - Hz[i, :, :]) - (dt / dz) * (Hx_int_bottom - Hx[i, :, :]));
            Chz = ((dt / dx) * (np.roll(Hy[i, :, :], 1, axis=1) - Hy[i, :, :]) - (dt / dy) * (np.roll(Hx[i, :, :], 1, axis=0) - Hx[i, :, :]));

            CHX[i, :, :] = Chx;
            CHY[i, :, :] = Chy;
            CHZ[i, :, :] = Chz;


            ####################UPDATE OF THE SURFACE YEE CELL########################################
            ## this is also where we do updates on the Ez sheet, something new from the 1D case
            ## Time step the Ez sheet fields using the updated Hy and Hx sheet fields
            # calculate derivatives of the sheet fields for Ez update
            deriv_Hy_x = np.roll(Hy_sheet[:, :, 0], 1, axis=1) - Hy_sheet[:, :, 0];
            deriv_Hx_y = np.roll(Hx_sheet[:, :, 0], 1, axis=0) - Hx_sheet[:, :, 0];
            deriv_Hy_x2 = np.roll(Hy_sheet[:, :, 1], 1, axis=1) - Hy_sheet[:, :, 1];
            deriv_Hx_y2 = np.roll(Hx_sheet[:, :, 1], 1, axis=0) - Hx_sheet[:, :, 1];



        elif(i == sheet_loc_z+1): ## yee cell right afte rthe sheet field update
            Hy_int_top = Hy_sheet[:,:,1];
            Hx_int_top = Hx_sheet[:,:,1];
            Chx = ((dt / dy) * (np.roll(Hz[i, :, :], 1, axis=0) - Hz[i, :, :]) - (dt / dz) * ( Hy[i, :, :]-Hy_int_top));
            Chy = -((dt / dx) * (np.roll(Hz[i, :, :], 1, axis=1) - Hz[i, :, :]) - (dt / dz) * (Hx[i, :, :]- Hx_int_top ));
            Chz = ((dt / dx) * (np.roll(Hy[i, :, :], 1, axis=1) - Hy[i, :, :]) - (dt / dy) * (np.roll(Hx[i, :, :], 1, axis=0) - Hx[i, :, :]));
            ## add dielectric
            # Chx = np.multiply(1/eps[i,:,:], Chx);
            # Chy = np.multiply(1/eps[i,:,:], Chy);
            # Chz = np.multiply(1/eps[i,:,:], Chz);
            CHX[i, :, :] = Chx;
            CHY[i, :, :] = Chy;
            CHZ[i, :, :] = Chz;


        else:

            Chx = ((dt / dy) * (np.roll(Hz[i, :, :], 1, axis=0) - Hz[i, :, :]) - (dt / dz) * (Hy[index, :, :] - Hy[i, :, :]));
            Chy = -((dt / dx) * (np.roll(Hz[i, :, :], 1, axis=1) - Hz[i, :, :]) - (dt / dz) * (Hx[index, :, :] - Hx[i, :, :]));
            Chz = ((dt / dx) * (np.roll(Hy[i, :, :], 1, axis=1) - Hy[i, :, :]) - (dt / dy) * (np.roll(Hx[i, :, :], 1, axis=0) - Hx[i, :, :]));
            ## add dielectric
            # Chx = np.multiply(1/eps[i,:,:], Chx);
            # Chy = np.multiply(1/eps[i,:,:], Chy);
            # Chz = np.multiply(1/eps[i,:,:], Chz);
            CHX[i, :, :] = Chx;
            CHY[i, :, :] = Chy;
            CHZ[i, :, :] = Chz;


    ## NOW UPDATE EVERYTHING
    # PLUS SIGN BECAUSE WE HAVE REVERSED THE CURL
    Ez_sheet[:, :, 0] += (dt / dx) * (deriv_Hx_y - deriv_Hy_x)
    Ez_sheet[:, :, 1] += (dt / dx) * (deriv_Hx_y2 - deriv_Hy_x2)

    ## UPDATE THE NORMAL FIELD COMPONENTS IN THE SURFACE CELL
    # there is no corresponding Hy or Hx for the Ex and Ey in the sbc cell, only the sheet fields.
    # derivative with the sheet field involved is not a full central difference, it's a backward diff, so a facto of 1/2
    Ex_sheet += (
    (dt / dy) * (np.roll(Hz_sheet, 1, axis=0) - Hz_sheet) - (dt / dz) *0.5* (Hy_int_bottom - Hy[sheet_loc_z, :, :]));
    Ey_sheet += -(
    (dt / dx) * (np.roll(Hz_sheet, 1, axis=1) - Hz_sheet) - (dt / dz) *0.5* (Hx_int_bottom - Hx[sheet_loc_z, :, :]));
    Ex -= CHX;
    Ey -= CHY;
    Ez -= CHZ;

    ## SOURECE UPDATE
    J = 1*np.sin(2*np.pi*t/30)
    Ez[jx, jy, jz] -=J;

    #############=========== H FIELD UPDATES =====================#################
    CEX = np.zeros((Nx,Ny,Nz));
    CEY = np.zeros((Nx,Ny,Nz));
    CEZ = np.zeros((Nx,Ny,Nz));
    for i in range(1,Nz):
        index = i-1;
        if(index <0):
            index = Nz-1;

        if( i== sheet_loc_z): # we are at the yee cell right before the surface conductivity cell
            ## The bordering yee cell coupled to Ey and Ex on the surface conductivity sheet
            Ez_top = Ez_sheet[:,:,0];
            Ez_bottom = Ez_sheet[:,:,1];

            ## the updates here do not touch any of the surface cell components...
            CEx =  ((dt/dy)*(Ez[i,:,:] - np.roll(Ez[i,:,:], -1, axis = 0)) - (dt/dz)*(Ey[i,:,:]- Ey[index, :, :]));
            CEy = -((dt/dx)*(Ez[i,:,:] - np.roll(Ez[i,:,:], -1, axis = 1)) - (dt/dz)*(Ex[i,:,:]- Ex[index, :, :]));
            CEz = ((dt/dx)*(Ey[i,:,:] - np.roll(Ey[i,:,:], -1, axis = 1)) - (dt/dy)*(Ex[i,:,:]- np.roll(Ex[i,:,:], -1, axis = 0)));
            CEX[i, :, :] = CEx;
            CEY[i, :, :] = CEy;

            CEZ[i, :, :] = CEz;


            ## define F coefs for the Hx sheet ## WITHIN HERE< WE ARE ACTUALLY GOING TO RUN THE UPDATE ON THE H FIELDS
            ## INSIDE THE CONDUCTIVITY CELL
            # define derivatives of the sheet fields
            deriv_y_1 = (np.roll(Ez_sheet[:, :, 0], 1, axis=1) - Ez_sheet[:, :, 0])/dz;
            deriv_y_2 = (np.roll(Ez_sheet[:, :, 1], 1, axis=1) - Ez_sheet[:, :, 1])/dz;

            F1n_Hx = fh11 * (Hx_sheet[:, :, 0] + fh12 * Hx_sheet[:, :, 1]) - \
                     fe1 * (2 * Ey[:, :, sheet_loc_z]) + dz * deriv_y_1;
            F2n_Hx = fh21 * (Hx_sheet[:, :, 0] + fh22 * Hx_sheet[:, :, 1]) + \
                     fe2 * (2 * Ey[:, :, sheet_loc_z+1] - dz * deriv_y_2);

            F1n_Hy = fh11y * (Hy_sheet[:, :, 0] + fh12y * Hy_sheet[:, :, 1]) - \
                     fe1y * (2 * Ey[:, :, sheet_loc_z]) + dz * deriv_y_1;
            F2n_Hy = fh21y * (Hy_sheet[:, :, 0] + fh22y * Hy_sheet[:, :, 1]) + \
                     fe2y * (2 * Ey[:, :, sheet_loc_z+1] - dz * deriv_y_2) ;
            F1n_Hy = -F1n_Hy; F2n_Hy = -F2n_Hy;


        elif(i == sheet_loc_z+1): ## NOW WE ARE AT THE CELL  ABOVE THE SURFACE CONDUCTIVITY YEE CELL
            ## Coupling only happesn on the z derivatives

            CEx =  ((dt/dy)*(Ez[i,:,:] - np.roll(Ez[i,:,:], -1, axis = 0)) - (dt/dz)*(Ey[i,:,:]- Ey_sheet));
            CEy = -((dt/dx)*(Ez[i,:,:] - np.roll(Ez[i,:,:], -1, axis = 1)) - (dt/dz)*(Ex[i,:,:]- Ex_sheet));
            CEz = ((dt/dx)*(Ey[i,:,:] - np.roll(Ey[i,:,:], -1, axis = 1)) - (dt/dy)*(Ex[i,:,:]- np.roll(Ex[i,:,:], -1, axis = 0)));

            CEX[i, :, :] = CEx;
            CEY[i, :, :] = CEy;
            CEZ[i, :, :] = CEz;


        else: ## execute normal update
            CEx =  ((dt/dy)*(Ez[i,:,:] - np.roll(Ez[i,:,:], -1, axis = 0)) - (dt/dz)*(Ey[i,:,:]- Ey[index,:,:]));
            CEy = -((dt/dx)*(Ez[i,:,:] - np.roll(Ez[i,:,:], -1, axis = 1)) - (dt/dz)*(Ex[i,:,:]- Ex[index,:,:]));
            CEz = ((dt/dx)*(Ey[i,:,:] - np.roll(Ey[i,:,:], -1, axis = 1)) - (dt/dy)*(Ex[i,:,:]- np.roll(Ex[i,:,:], -1, axis = 0)));
            CEX[i, :, :] = CEx;
            CEY[i, :, :] = CEy;
            CEZ[i, :, :] = CEz;
    ## once you calculate all the quantities, then do the updates
    ## We also need to update the Hz component of the surface cell
    Hz_sheet += ((dt / dx) * (Ey_sheet - np.roll(Ey_sheet, -1, axis=1)) - \
                 (dt / dy) * (Ex_sheet - np.roll(Ex_sheet, -1, axis=0)));

    Hx_sheet[:, :, 0] = (1 / (1 - c1 * c2)) * (F1n_Hx + c1 * F2n_Hx);
    Hx_sheet[:, :, 1] = (1 / (1 - c1 * c2)) * (F2n_Hx + c2 * F1n_Hx);  # use Fx_sheet to update the E fields

    StorageSheet.append([Hx_sheet, Hy_sheet, Hz_sheet])

    Hy_sheet[:, :, 0] = (1 / (1 - c1 * c2)) * (F1n_Hy + c1 * F2n_Hy);
    Hy_sheet[:, :, 1] = (1 / (1 - c1 * c2)) * (F2n_Hy + c2 * F1n_Hy);  # use Fx_sheet to update the E fields

    Hx += CEX;
    Hy += CEY;
    Hz += CEZ;


    if(t%2 == 0):
        plt.subplot(4,3,1)
        # print(np.max(Hx_sheet));
        # print(np.max(Ez))
        # print(np.max(Ex_sheet))
        plt.subplot(4,3,2)
        imgplot = plt.imshow(Hx_sheet[:,:,1], cmap = 'jet')
        plt.subplot(4,3,3)
        imgplot = plt.imshow(Ez_sheet[:,:,0], cmap = 'jet')
        plt.subplot(4,3,4)
        imgplot = plt.imshow(Ez_sheet[:,:,1], cmap = 'jet')
        plt.subplot(4,3,5)
        imgplot = plt.imshow(Hy_sheet[:,:,0], cmap = 'jet')
        plt.subplot(4,3,6)
        imgplot = plt.imshow(Hy_sheet[:,:,1], cmap = 'jet')
        plt.subplot(4,3,7)
        imgplot = plt.plot(Ez[:,jy,sheet_loc_z])
        plt.subplot(4,3,8)
        imgplot = plt.plot(Ey[:,jy,sheet_loc_z])
        plt.subplot(4,3,9)
        imgplot = plt.imshow(Ex[:,:,jz], cmap = 'jet')
        plt.subplot(400,3,10)
        plt.clim(-0.01, 0.01)
        imgplot = plt.imshow(Ey[:,:,jz], cmap = 'jet')
        plt.clim(-0.01, 0.01)
        plt.subplot(4,3,11)
        imgplot = plt.imshow(Ez[:,:,jz], cmap = 'jet')
        plt.clim(-0.01, 0.01)
        #plt.plot(Ex[:,:,0jz])
        #imgplot = mlab.co0ntour3d(x,y,z,Hz)
        plt.pause(0.001)
        #print(np.max(Ez))
        #print('current s00ource power: '+str(J))
        #print(Ez[jx, jy, jz])
        plt.clf()

