import numpy as np
import matplotlib.pyplot as plt

Nx = 200;
Ny = 200;
## create epsilon
eps = np.ones((Nx, Ny));

##initialize fields
Ex = np.zeros((Nx,Ny));
Ey = np.zeros((Nx, Ny));
Hz = np.zeros((Nx, Ny));
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
            indx = j + 1; indy = k + 1;
            if (indx > Nx - 1):
                indx = 0;
            if (indy > Ny - 1):
                indy = 0;
            # TEz mode...only transverse x and y E field
            deriv_y[k, j] = (Ex[indy, j] - Ex[k, j]);
            deriv_x[k, j] = (Ey[k, indx] - Ey[k, j]);

    Hz -= (deriv_x - deriv_y);
    Hz[100, 100] -= np.sin(2*np.pi*t/30)*(dt);


    #update E field components
    Curl_x = np.zeros((Nx, Ny));
    Curl_y = np.zeros((Nx, Ny));
    for j in range(Nx): #y axis
        for k in range(Ny): #x axis;
            indx = j - 1; indy = k - 1;
            if(indx < 0):
                indx = Nx - 1;
            if(indy < 0):
                indy = Ny - 1;
            deriv_y = -(Hz[indy, j] - Hz[k, j]);
            deriv_x = -(Hz[k, indx] - Hz[k, j]);
            Curl_x[k, j]= (cc) * (deriv_x); #curl Hz, dy(Hz)
            Curl_y[k, j] = (cc) * (deriv_y) ; #cul Hz dx(Hz)
    Ex += Curl_y;
    Ey -= Curl_x;



    #insert point source
    if(t%10 == 0):
        plt.subplot(121)
        imgplot = plt.imshow(Ey)
        imgplot.set_cmap('jet')
        plt.subplot(122)
        imgplot = plt.imshow(Hz)
        imgplot.set_cmap('jet')
        plt.pause(0.001)
        plt.clf()
plt.imshow(Hz)
plt.figure()
plt.plot(Hz[:,50])
#get a 1D slice of the 2d field
plt.show()
