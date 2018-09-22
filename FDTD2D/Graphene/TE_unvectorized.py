import numpy as np
import matplotlib.pyplot as plt

Nx = 200;
Ny = 200;
## create epsilon
eps = np.ones((Nx, Ny));

##initialize fields
Hx = np.zeros((Nx, Ny));
Hy = np.zeros((Nx, Ny));
Ez = np.zeros((Nx, Ny));
#spatial steps
dx = dy = 1;

#time steps
tsteps = 100;
dt = 1
c = 1;
## courant constant
cc = c*dt/min(dx, dy);
plt.ion()
fieldAmplitudes = list();
plt.figure(figsize=(10,10));
cc = 0.5;
for t in range(tsteps):
    print('t= '+str(t))
    #update Hz, Ex, Ey
    #remember the yee grid and integer indexing

    # update H field components
    deriv_y = np.zeros((Nx, Ny));
    deriv_x = np.zeros((Nx, Ny));
    for j in range(Nx):
        for k in range(Ny):
            indx = j - 1; indy = k - 1;
            if (indx < 0):
                indx = Nx - 1;
            if (indy < 0):
                indy = Ny - 1;
            # TEz mode...only transverse x and y E field
            deriv_y[k, j] = -(Hx[indy, j] - Hx[k, j]);
            deriv_x[k, j] = -(Hy[k, indx] - Hy[k, j]);

    Ez -= (deriv_x - deriv_y);
    Ez[100, 100] -= np.sin(2 * np.pi * t / 30) * (dt);


    #update E field components
    Curl_x = np.zeros((Nx, Ny));
    Curl_y = np.zeros((Nx, Ny));
    for j in range(Nx): #y axis
        for k in range(Ny): #x axis;
            indx = j + 1; indy = k + 1;
            if (indx > Nx - 1):
                indx = 0;
            if (indy > Ny - 1):
                indy = 0;

            deriv_y = -(Ez[indy, j] - Ez[k, j]);
            deriv_x = -(Ez[k, indx] - Ez[k, j]);
            Curl_x[k, j]= (cc) * (deriv_x); #curl Hz, dy(Hz)
            Curl_y[k, j] = (cc) * (deriv_y) ; #cul Hz dx(Hz)
    Hx -= Curl_y;
    Hy += Curl_x;



    #insert point source
    if(t%10 == 0):
        plt.subplot(121)
        imgplot = plt.imshow(Hy)
        imgplot.set_cmap('jet')
        plt.subplot(122)
        imgplot = plt.imshow(Ez)
        imgplot.set_cmap('jet')
        plt.pause(0.001)
        plt.clf()
plt.imshow(Ez)
plt.figure()
plt.plot(Ez[:, 50])
#get a 1D slice of the 2d field
plt.show()
