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
tsteps = 200;
dt = 1
c = 1;
## courant constant
cc = c*dt/min(dx, dy);
imp0 = 377;
plt.ion()
fieldAmplitudes = list();
plt.figure;
for t in range(tsteps):
    print('t= '+str(t))
    #update Hz, Ex, Ey
    #remember the yee grid and integer indexing

    Hz[100, 100] += np.exp(-(t-30)*(t-30)/50)
    #update E field components
    for j in range(1,Nx):
        for k in range(1,Ny):
            Ex[j, k] = Ex[j, k] - (cc)*(Hz[j, k] - Hz[j, k - 1]) ;
            Ey[j, k] = Ey[j, k] + (cc)*(Hz[j ,k] - Hz[j - 1, k]) ;

    for j in range(Nx - 1):
        for k in range(Ny - 1):
            # TEz mode...only transverse x and y E field
            Hz[j, k] = Hz[j, k] + (cc) * (Ey[j + 1, k] - Ey[j, k] - \
                                          (Ex[j, k + 1] - Ex[j, k]));

    #insert point source
    plt.imshow(Ey)
    plt.pause(0.05)
    plt.clf()

plt.imshow(Hz)
plt.figure()
plt.plot(Hz[:,50])
#get a 1D slice of the 2d field
plt.show()
