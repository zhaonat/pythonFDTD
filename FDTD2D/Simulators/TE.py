import numpy as np
import matplotlib.pyplot as plt
import numpy as np;


def currSource(t, lambda_up, lambda_low, lambda_0, c, dx):
    omega_0 = 2*np.pi/lambda_0;
    sigma = (2/omega_0)*(lambda_0)/(lambda_up-lambda_low)
    #print(omega_0); print(sigma)
    return np.exp(-(t-4*sigma)**2/sigma**2)*np.sin(omega_0*t)*(10/dx)**2;


c = 1; # we are in units of nanometers;
lambda_0 = 800;
lambda_up = 1200;
lambda_low = 400;
dx = lambda_0/160; dy = dx;
dt = (1/2)**.5*dx/c;
xrange = [-400, 400];
yrange = [-400, 400];
xstart = np.array([-400,-400]); xend = np.array([400,400]);
source =  np.array([200, 0])-xstart; source = [int(i/dx) for i in source];
probloc = np.array([200,0])-xstart; probloc = [int(i/dx) for i in probloc];
source = [80,80]
Nx = int((xrange[1]-xrange[0])/dx);
Ny = int((yrange[1]-yrange[0])/dy);
## create epsilon
eps = np.ones((Nx, Ny));
#eps[50:150, 50:150] = 4;


##initialize fields (consider a sparse initializiation)
Hx = np.zeros((Nx, Ny));
Hy = np.zeros((Nx, Ny));
Ez = np.zeros((Nx, Ny));
#spatial steps
#time steps
tsteps = 10000;

## courant constant
plt.ion()
fieldAmplitudes = list();
plt.figure;
storedFields = list();
for t in range(tsteps):
    print('t= '+str(t))
    #update Hz, Ex, Ey
    #remember the yee grid and integer indexing
    ## current source needs to be scaled
    J = np.sin(2*np.pi*t/40)*(dt)#*np.exp(-(t-20)**2/50);
    #J = currSource((t+0.5)*dt, lambda_up, lambda_low, lambda_0, c, dx);
    #update E field components
        # TEz mode...only transverse x and y E field
    Hx[-1,:] = 0;
    Hy[:,-1] = 0;
    deriv_y = (np.roll(Hx,-1, axis=1)-Hx)/dy;
    deriv_x = (np.roll(Hy, -1, axis = 0)-Hy)/dx;
    Dz = eps*Ez;
    Dz= Dz +((dt)* (deriv_x - \
                                  (deriv_y)));
    Ez = Dz/eps;
    ## enforce Dirichlet
    #if( t < 20):
    Ez[source[0],source[1]] -= J/eps[source[0], source[1]]; #location of the source            # np.roll(Hx,-1)-Hx;
    ##vectorize this code
    Ez[:,0] = 0;
    #deriv_y = dey/dy;
    deriv_y = (Ez - np.roll(Ez,1,axis=1 ))/dy;
    Ez[0,:] = 0;
    #derix_x = dex;
    deriv_x = (Ez - np.roll(Ez,1, axis = 0))/dx;
    Hx = Hx - (dt) * deriv_y;
    Hy = Hy + (dt) * deriv_x ;

    ## records field in probloc
    storedFields.append(Ez[probloc[0], probloc[1]]);

    #insert point source
    if(t%10 == 0):
        imgplot = plt.imshow(Ez)
        imgplot.set_cmap('jet')
        plt.pause(0.001)
        plt.clf()

storedFields = np.array(storedFields);
FFT = np.fft.fft(storedFields)
plt.plot(np.abs(FFT));
plt.show()
