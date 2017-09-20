## first, we need to be able to visualize a non-uniform grid
import numpy as np;
import matplotlib.pyplot as plt
plt.close('All')
## want to apply the pml as a matrix operation on a matrix HX of nxmxl
## calculate the PML on a 2x grid, since the even values go for H
## and the odd values will go for E
## EMLAB calculation of PML

lambda_0 = 100; omega_0 = 2*np.pi/lambda_0;
lambda_up = 200;
lambda_low = 50;
dx = 800/100; dy = dx; c =1;
dt = (1/2)**.5*dx/c;
xrange = [-400, 400];
yrange = [-400, 400];
xstart = np.array([-400,-400]); xend = np.array([400,400]);
source =  np.array([0,0])-xstart; source = [int(i/dx) for i in source];
probloc = np.array([200,0])-xstart; probloc = [int(i/dx) for i in probloc];

Nx = int((xrange[1]-xrange[0])/dx);
Ny = int((yrange[1]-yrange[0])/dy);

c = 1; # we are in units of nanometers;

Npml = np.array([40,40]);
Nx2 = Nx+2*Npml[0];
Ny2 = Ny+2*Npml[1];
Nx2 = 2*Nx;
Ny2 = 2*Nx;

## x-direction pml
## no 2x grid concept
sigx = np.zeros((Nx2, Ny2));
#left hand side (need negative values)
for nx in range(2*Npml[0]):
    nx1 = 2*Npml[0]-nx-1;
    sigx[nx1, :] = (0.5/dt)*(nx/2/Npml[0])**4;
#right hand side
for nx in range(2*Npml[0]):
    nx1 = Nx2-2*Npml[1]+nx;
    sigx[nx1,:] =(0.5/dt)*(nx/2/Npml[1])**4

## vertical y direction pml
sigy = np.zeros((Nx2, Ny2))
for ny in range(2*Npml[1]):
    ny1 = 2*Npml[0] - ny-1 ;
    sigy[:, ny1] = (0.5/dt)*(ny/2/Npml[0])**4;
for ny in range(2*Npml[1]):
    ny1 = Ny2 - 2*Npml[1] +ny; #5 index in array for pml mapped value to go
    sigy[:, ny1] = (0.5/dt)*(ny/2/Npml[1])**4;

plt.imshow(sigx)
plt.figure()
plt.imshow(sigy)
plt.show()
## remember, no imaginary parts because we've transformed back into time domain


URxx =1; ## these are permability term
URyy = 1;

## CALCULATE UPDATE COEFFICIENTS with the PML incorporated
sigHx0 = sigx[1:Nx2:2, 0:Ny2:2]
sigHy0 = sigy[0:Nx2:2, 1:Ny2:2];
mHx0 = (1/dt)+sigHy0/2
mHx1 = np.divide(((1/dt)-sigHy0/2),mHx0);
mHx2 = -1/(URxx*mHx0) #update to curl term
mHx3 = -dt*np.divide(sigHx0,mHx0)

sigHx = sigx[1:Nx2:2, 0:Ny2:2]
sigHy = sigy[0:Nx2:2, 1:Ny2:2];
mHy0 = (1/dt)+sigHx/2; # allvalues greater than 1, #concave up like sigHy
mHy1 = np.divide(((1/dt)-sigHx/2),mHy0); #update to past Hy field, concave down
mHy2 = -np.squeeze(1/(mHy0)); #update to curl term
mHy3 = -(dt)*np.divide(sigHy,URyy*mHy0);

## Ez is in the center of the square yee cell
sigDx = sigx[1:Nx2:2, 1:Ny2:2];
sigDy = sigy[1:Nx2:2, 1:Ny2:2];
mDz0 = (1/dt) + (sigDx + sigDy)/2 + (sigDx*sigDy)*dt/2
mDz1 = (1/dt) - (sigDx + sigDy)/2 - (sigDx*sigDy)*dt/2;
mDz1 = np.divide(mDz1,mDz0);
mDz2 = 1/mDz0;
#mDz3 is 0
mDz4 = -(dt)*(sigDx*sigDy/mDz0);

plt.figure()
plt.plot(sigHy[1,:])
plt.plot(sigHx[:,1])
plt.figure()
plt.plot(sigHx[:,1])
plt.plot(sigDy[:,1])
# plt.subplot(221)
# plt.imshow(mDz0)
# plt.subplot(222)
# plt.imshow(mDz1)
# plt.subplot(223)
# plt.imshow(mDz2)
# plt.subplot(224)
# plt.imshow(mDz4)
plt.show()

#### TEST THE PML #######################################################
#########################################################################

def currSource(t, lambda_up, lambda_low, lambda_0, c, dx):
    omega_0 = 2*np.pi/lambda_0;
    sigma = (2/omega_0)*(lambda_0)/(lambda_up-lambda_low)
    #print(omega_0); print(sigma)
    return np.exp(-(t-4*sigma)**2/sigma**2)*np.sin(omega_0*t)*(10/dx)**2;



## create epsilon
eps = np.ones((Nx, Ny));
#eps[50:150, 50:150] = 4;

##initialize fields (consider a sparse initializiation)
Hx = np.zeros((Nx, Ny));
Hy = np.zeros((Nx, Ny));
Ez = np.zeros((Nx, Ny));
#spatial steps
#time steps

## courant constant
plt.ion()
fieldAmplitudes = list();
plt.figure;
storedFields = list();

Idz = 0; #integral collection for curl
Icey = 0; Icex = 0;

tsteps = 200;
source = [50, 50]
for t in range(tsteps):
    print('t= '+str(t))

    J = np.sin(2*np.pi*t/10)*dt;
    J = currSource((t+0.5)*dt, lambda_up, lambda_low, lambda_0, c, dx);
    ## code currently explodes at the origin

    ## update H fields
    Ez[:,0] = 0;
    Ez[0,:] = 0;
    Cex = (Ez - np.roll(Ez, 1, axis=1)) / dy;
    Cey = -(Ez - np.roll(Ez, 1, axis = 0)) / dx;
    Icex += Cex;
    Icey += Cey;
    Hx = mHx1*Hx +mHx2*Cex +mHx3*Icex; # add integral term
    Hy = mHy1*Hy + mHy2*Cey +mHy3*Icey;

    ## now work on Dz
    Hx[-1,:] = 0;
    Hy[:,-1] = 0;
    derivy = (np.roll(Hx, -1, axis=1) - Hx) / dy;
    derivx = (np.roll(Hy, -1, axis = 0) - Hy) / dx;
    Chz = derivx - derivy;

    Dz = np.squeeze(np.multiply(eps, Ez));
    Idz +=Dz;
    Dz = mDz1*Dz + mDz2 * Chz + mDz4*Idz;
    Dz = np.squeeze(Dz)
    Dz[source[0],source[1]]-=J; #location of the source

    ## finally update Ez
    Ez = (np.divide(Dz,eps)); ## this is important for the next update of the curl # equations

    ## records field in probloc
    storedFields.append(Ez[probloc[0], probloc[1]]);

    #insert point source
    if(t%10 == 0):
        print(np.max(Ez))
        plt.subplot(121)
        imgplot = plt.imshow(np.real(Hy))
        imgplot.set_cmap('jet')
        plt.subplot(122)
        plt.imshow(np.real(Idz))
        plt.pause(0.1)
        plt.clf()

# storedFields = np.array(storedFields);
# FFT = np.fft.fft(storedFields)
# plt.plot(np.abs(FFT));
plt.show()
