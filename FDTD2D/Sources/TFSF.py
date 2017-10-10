import numpy as np

## first, we need to be able to visualize a non-uniform grid
import numpy as np;
import matplotlib.pyplot as plt
import numpy.matlib
plt.close('All')
## want to apply the pml as a matrix operation on a matrix HX of nxmxl
## calculate the PML on a 2x grid, since the even values go for H
## and the odd values will go for E
## EMLAB calculation of PML

lambda_0 = 100; omega_0 = 2*np.pi/lambda_0;
lambda_up = 200;
lambda_low = 50;
dx = 800/200; dy = dx; c =1;
dt = (1/2)**.5*dx/c;
xrange = [-600, 600];
yrange = [-600, 600];
rstart = np.array([-400, -300]); rend = np.array([400, 300]);
source = np.array([0,0]) - rstart; source = [int(i / dx) for i in source];
source = [100, 75]
probloc = np.array([200,0]) - rstart; probloc = [int(i / dx) for i in probloc];
Nx = int((xrange[1]-xrange[0])/dx);
Ny = int((yrange[1]-yrange[0])/dy);
N = [Nx, Ny];
c = 1; # we are in units of nanometers;

Npml = 1*np.array([20,20]);
Nx2 = Nx+2*Npml[0];
Ny2 = Ny+2*Npml[1];


## x-direction pml
## no 2x grid concept
## need to make the code compatible for non-square grids

xn = np.zeros((Nx,1))
fi1 = np.zeros((Nx,1)); fj1 = np.zeros((Ny, 1));
gi2 = np.ones((Nx,1)); gj2 = np.ones((Ny, 1));
gi3 = np.ones((Nx,1)); gj3 = np.ones((Ny, 1));
Npmlx = Npml[0]; Npmly = Npml[1];

for i in range(Npmlx):
    x_n = (1/3)*((i+1)/Npmlx)**3;
    xn[Nx-Npmlx+i] = x_n;
    xn[Npmlx-1-i] = x_n;

    fi1[Npmlx-1-i] = x_n;
    gi2[Npmlx-1-i] = 1/(1+x_n);
    gi3[Npmlx-1-i] = (1-x_n)/(1+x_n);

    fi1[Nx-Npmlx+i] = x_n;
    gi2[Nx-Npmlx+i] = 1/(1+x_n);
    gi3[Nx-Npmlx+i] = (1-x_n)/(1+x_n);

for i in range(Npmly):
    x_n = (1/3)*((i+1)/Npml[1])**3;
    fj1[Npmly-1-i] = x_n;
    gj2[Npmly-1-i] = 1/(1+x_n);
    gj3[Npmly-1-i] = (1-x_n)/(1+x_n);

    fj1[Ny-Npmly+i] = x_n;
    gj2[Ny-Npmly+i] = 1/(1+x_n);
    gj3[Ny-Npmly+i] = (1-x_n)/(1+x_n);


gi3 = np.squeeze(gi3);
gi2 = np.squeeze(gi2);
gj3 = np.squeeze(gj3);
gj2 = np.squeeze(gj2);

## creat the other f's from teh g's
fi2 = gi2; fj2 = gj2;
fj3 = gj3; fi3 = gi3;

fi1 = np.matlib.repmat(fi1.T,Ny,1);
fj1 = np.matlib.repmat(fj1,1,Nx);

fi2 = np.matlib.repmat(fi2.T,Ny,1);
fi2 = np.rot90(fi2);
fj2 = np.matlib.repmat(fj2.T,Nx,1);

fi3 = np.matlib.repmat(fi3,Ny,1);
fi3 = np.rot90(fi3)
fj3 = np.matlib.repmat(fj3.T,Nx,1);

## at this point all the fi's and fj's only have the coeffs on one side
fij3 = np.multiply(fj3,fi3);
fij2 = np.multiply(fj2,fi2);
fij1 = np.multiply(fj1,fi1); fij1 = np.rot90(fij1)

gi3x = np.matlib.repmat(gi3.T,Ny,1); gi3x = np.rot90(gi3x);
gj3x = np.matlib.repmat(gj3,Nx,1);
gij3 = np.multiply(gi3x, gj3x);

gi2x = np.matlib.repmat(gi2.T,Ny,1); gi2x = np.rot90(gi2x)
gj2x = np.matlib.repmat(gj2,Nx,1);
gij2 = np.multiply(gi2x, gj2x);

###################INSERT THE TFSF######################################
########################################################################
#select a  single line of nodes
TFSF_x = 100;


#### TEST THE CODE #######################################################
#########################################################################

def currSource(t, lambda_up, lambda_low, lambda_0, c, dx):
    omega_0 = 2*np.pi/lambda_0;
    sigma = (2/omega_0)*(lambda_0)/(lambda_up-lambda_low)
    #print(omega_0); print(sigma)
    return np.exp(-(t-4*sigma)**2/sigma**2)*np.sin(omega_0*t)*(10/dx)**2;


## create epsilon
eps = np.ones((Nx, Ny));
eps[150:200,:] = 4;
##initialize fields (consider a sparse initializiation)
Hx = np.zeros((Nx, Ny));
Hy = np.zeros((Nx, Ny));
Ez = np.zeros((Nx, Ny));
#spatial steps
#time steps

## courant constant
plt.ion()
fieldAmplitudes = list();
plt.figure(figsize=(25, 20));
storedFields = list();

Idz = 0; #integral collection for curl
Icey = 0; Icex = 0;

tsteps = 1000;
source = [TFSF_x, 25]
time = np.arange(0,tsteps+1);
J = np.sin(2*np.pi*time/40);
J = np.exp(-(time-30)**2/100)
for t in range(tsteps):
    print('t= '+str(t))


    #J = currSource((t+0.5)*dt, lambda_up, lambda_low, lambda_0, c, dx);

    ## now work on Dz
    # Hx[-1,:] = 0;
    # Hy[:,-1] = 0;
    derivy = (np.roll(Hx, -1, axis=1) - Hx) / dy;
    derivx = (np.roll(Hy, -1, axis=0) - Hy) / dx;
    Chz = derivx - derivy;

    Dz = eps * Ez;
    Idz += Dz;
    Dz = gij3 * Dz + gij2 * dt * Chz;
    Dz[source[0], source[1]] -= J[t];  # location of the source

    ## finally update Ez
    Ez = Dz/eps;  ## this is important for the next update of the curl # equations
    Ez[TFSF_x, :] = J[t];

    ## update H fields
    #Ez[:,0] = 0;
    #Ez[0,:] = 0;
    Cex = -(Ez - np.roll(Ez, 1, axis=1)) / dy;
    Cey = +(Ez - np.roll(Ez, 1, axis = 0)) / dx;
    Icex += Cex; Icey += Cey;
    Hx = fij3*Hx + (dt)*fij2*Cex + fij1*Icex;
    Hy = fij3*Hy + (dt)*fij2*Cey + fij1*Icey;

    #update the Hx source for the TFSF
    delay = dy/2+dt/2
    #Hy[TFSF_x, :] += Ez[TFSF_x, :]

    ## records field in probloc
    storedFields.append(Ez[probloc[0], probloc[1]]);

    #insert point source
    if(t%15 == 0):
        print(np.max(Ez))
        plt.subplot(221)
        imgplot = plt.imshow(Ez)
        imgplot.set_cmap('jet')
        plt.clim(-1, 1)
        plt.subplot(222)
        plt.imshow(np.real(Idz))
        plt.clim(-0.04, 0.04)
        plt.subplot(223)
        plt.plot(Ez[40,:])
        plt.subplot(224)
        plt.plot(Hy[40,:])
        plt.pause(0.1)
        plt.clf()

# storedFields = np.array(storedFields);
# FFT = np.fft.fft(storedFields)
# plt.plot(np.abs(FFT));
plt.show()
