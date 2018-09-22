import numpy as np
import matplotlib.pyplot as plt
## only Hx component sees the graphene if we put the sheet along x axis

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

## COEFFICIENTS FOR SHEET
sheet_loc_y= 85;
Hx_int_left = np.zeros((Nx,));
Hx_int_right = np.zeros((Nx,));
tau = 0;
#time steps
tsteps = 200;
dt = 1
c = 1;

## COEFFICIENTS FOR THE SHEET UPDATE


sigma_0 = 0.1;
denom = (dx * sigma_0) - dt;
c1 = (dt) / (dx * sigma_0 - dt);
c2 = (dt) / (dx * sigma_0 - dt); ## if this is the courant number, there is instability introduced

fe1 = sigma_0 * c1;
fe2 = sigma_0 * c2;
fh11 = 1-2.0*c1;
fh22 = 1-2.0*c2
fh12 = c1;
fh21 = c2;


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
            indx = j - 1;
            indy = k - 1;
            if (indx < 0):
                indx = Nx - 1;
            if (indy < 0):
                indy = Ny - 1;

            #Ez field yee updates couple into the sheet field
            if(k == sheet_loc_y):
                deriv_y[k, j] = +( Hx_int_left[j] - Hx[k-1,j]);
                deriv_x[k, j] = -(Hy[k, indx] - Hy[k,j]);
            elif(k == sheet_loc_y+1):
                deriv_y[k, j] = -(Hx_int_right[j]-Hx[k, j]);
                deriv_x[k, j] = -(Hy[k, indx] - Hy[k,j]);
            else:
                # TEz mode...only transverse x and y E field
                deriv_y[k, j] = -(Hx[indy, j] - Hx[k, j]);
                deriv_x[k, j] = -(Hy[k, indx] - Hy[k, j]);

    Ez -= (deriv_x - deriv_y);
    Ez[100, 100] -= np.sin(2 * np.pi * t / 30) * (dt);


    #update H field components
    Curl_x = np.zeros((Nx, Ny));
    Curl_y = np.zeros((Nx, Ny));
    F1n = np.zeros((Nx,));
    F2n = np.zeros((Nx,))
    for j in range(Nx): #y axis
        for k in range(Ny): #x axis;
            indx = j + 1;
            indy = k + 1;
            if (indx > Nx - 1):
                indx = 0;
            if (indy > Ny - 1):
                indy = 0;

            if(k == sheet_loc_y):
                deriv_x = -(Ez[k, indx] - Ez[k, j]);
                Curl_x[k, j] = (cc) * (deriv_x);  # curl Hz, dy(Hz) #Hy still exists at sheet loc
                 # DO NOT UPDATE THE FIELD AT THE SHEET LOC for Hx, IT DOES NOT EXIST
                ## consider the updates at k = sheet_loc_y, and k = sheet_loc_y+1
            elif(k == sheet_loc_y+1):
                deriv_y = -(Ez[indy, j] - Ez[k, j]);
                deriv_x = -(Ez[k, indx] - Ez[k, j]);
                Curl_x[k, j]= (cc) * (deriv_x); #curl Hz, dy(Hz)
                Curl_y[k, j] = (cc) * (deriv_y) ; #cul Hz dx(Hz)
            else:

                deriv_y = -(Ez[indy, j] - Ez[k, j]);
                deriv_x = -(Ez[k, indx] - Ez[k, j]);
                Curl_x[k, j]= (cc) * (deriv_x); #curl Hz, dy(Hz)
                Curl_y[k, j] = (cc) * (deriv_y) ; #cul Hz dx(Hz)

        F1n[j] = -2 * fe1 * Ez[sheet_loc_y, j] + (fh11 * Hx_int_left[j]) + (fh12 * (Hx_int_right[j]));
        F2n[j] = 2 * fe2 * Ez[sheet_loc_y+1 + 1, j] + (fh21 * Hx_int_left[j]) + (fh22 * (Hx_int_right[j]));
        # Ey is at time step n, but the interfaces are at n-1/2
        Hx_int_left[j] = (1 / (1 - c1 * c2)) * (F1n[j] + c1 * F2n[j]);  # time step to n+1/2 using F1n (at time step n)
        Hx_int_right[j] = (1 / (1 - c1 * c2)) * (F2n[j] + c2 * F1n[j]);
        print(np.max(F1n))


    Hx += -Curl_y;
    Hy += Curl_x;



    #insert point source
    if(t%2 == 0):
        # print(np.max(Hx_int_left));
        # print(np.max(Hx_int_right))
        # print(np.max(Hx))
        # print(np.max(Ez))
        plt.subplot(121)
        imgplot = plt.imshow(Hy)
        imgplot.set_cmap('jet')
        plt.subplot(122)
        imgplot = plt.imshow(Ez)
        imgplot.set_cmap('jet')
        plt.clim([-0.02, 0.02])
        plt.pause(0.001)
        plt.clf()
plt.imshow(Ez)
plt.figure()
plt.plot(Ez[:, 50])
#get a 1D slice of the 2d field
plt.show()
