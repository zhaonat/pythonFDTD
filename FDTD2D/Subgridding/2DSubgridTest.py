import numpy as np
import matplotlib.pyplot as plt
## test grid will be a coarse grid on left side, fine grid on right side

fullxrange = [-400, 400];
fullyrange = [-400, 400];
lambda_0 = 800;
dx_coarse = lambda_0/40;
dy_coarse = dx_coarse;
dx_fine = dx_coarse/2; ## refinement factor of 2;
dy_fine = dx_fine;
dt_coarse = (1/2)**.5*dx_coarse;
dt_fine = dt_coarse/2;

finexrange = [0, 400];
fineyrange = [-400, 400];
coarsexrange = [-400,0];
coarseyrange = [-400, 400]


Nxfine = int(800/dx_fine);
Nyfine = Nxfine;
Nxcoarse = int(800/dx_coarse);
Nycoarse = Nxcoarse;

Ez1 = np.zeros((Nxcoarse, Nycoarse)); Ez2 = np.zeros((Nxfine, Nyfine));
Hx1 = np.zeros((Nxcoarse, Nycoarse)); Hx2 = np.zeros((Nxfine, Nyfine));
Hy1 = np.zeros((Nxcoarse, Nycoarse)); Hy2 = np.zeros((Nxfine, Nyfine));

print(Ez1.shape);
print(Ez2.shape)
indx = int(Nxcoarse/2); indy = int(Nycoarse/2)
tsteps = 200;
fig, ax = plt.subplots(nrows=2, ncols=2)



## we need the dirichlet condition
for t in range(tsteps):
    print('timestep: '+str(t))
    ## step E components on coarse grid time step t
    Hx1[-1,:] = 0;
    deriv_x = dt_coarse*(np.roll(Hx1, -1, axis = 0) - Hx1)/dx_coarse;
    deriv_y = dt_coarse*(np.roll(Hy1, -1, axis = 1) - Hy1)/dy_coarse;
    Ez1 -= deriv_x - deriv_y
    ## insert simple point source in coarse grid
    J = np.sin(2*np.pi*t/20);
    Ez1[indx, 10] =J;

    ## step H components on coarse grid time step t+0.5
    Ez1[:,0] = 0;
    deriv_x = dt_coarse*(Ez1 - np.roll(Ez1, 1, axis = 0))/dx_coarse;
    Ez1[0,:] = 0;
    deriv_y = dt_coarse*(Ez1 - np.roll(Ez1, 1, axis = 1))/dy_coarse;
    Hx1 -= deriv_x;
    Hy1 += deriv_y;
    #insert point source
    # if(t%10 == 0):
    #     imgplot = plt.imshow(Ez1)
    #     imgplot.set_cmap('jet')
    #     plt.pause(0.01)
    #     plt.clf()

    ## INTERPOLATION at the quarter time step onto the fine grid
    # we can directly map elements onto the coarse grid to every other element on the fine grid
    ez1 = Ez1[:, -1]; hx1 = Hx1[:,-1]; hy1 = Hy1[:,-1];

    ## now insert 0's every other value
    insertIndices = np.arange(1,Nxcoarse+1);
    ez1 = np.insert(ez1,insertIndices ,0)
    hx1 = np.insert(hx1,insertIndices ,0)
    hy1 = np.insert(hy1,insertIndices ,0)


    # 0's represent values we need to interpolate
    # first step is to do a simple average
    ez = (np.roll(ez1, 2) + ez1) / 2
    ez = ez1 + np.roll(ez, 1)
    hx = (np.roll(hx1, -2) + hx1) / 2
    hx = hx1 + np.roll(hx, 1)
    hy = (np.roll(hy1, -2) + hy1) / 2
    hy = hx1 + np.roll(hy, 1)
    # plt.subplot(1,2,1)
    # plt.plot(hy)
    # plt.ylim([-0.2, 0.2])
    # plt.subplot(1,2,2)
    # plt.plot(hy1)
    # plt.ylim([-0.2, 0.2])
    # plt.pause(0.001)
    # plt.clf()

    ##Now we havr the subgrid field initial values...now we yee grid this
    Ez2[:,0] = ez;     Hx2[:,0] = hx;    Hy2[:,0] = hy;
    Hx2[-1,:] = 0;
    Hy2[:,-1] = 0;
    deriv_x = dt_fine*(np.roll(Hx2, -1, axis = 0) - Hx2)/dx_fine;
    deriv_y = dt_fine*(np.roll(Hy2, -1, axis = 1) - Hy2)/dy_fine;
    Ez2 -= deriv_x - deriv_y

    ## we have to be careful here
    # Ez1[:,0] = 0;

    deriv_x = dt_fine*(Ez2 - np.roll(Ez2, 1, axis = 0))/dx_fine;
    Ez2[0,:] = 0;
    deriv_y = dt_fine*(Ez2 - np.roll(Ez2, 1, axis = 1))/dy_fine;
    Hx2 -= deriv_x;
    Hy2 += deriv_y;

    ## Now we need to interpolate back onto the coarse grid
    ez11 = Ez2[0:-1:2,0];
    hx11 = Hx2[0:-1:2,0];
    hy11 = Hy2[0:-1:2,0];

    ## this reverse mapping is hard
    ez11 = (np.roll(ez11,1)+np.roll(ez11,2)+np.roll(ez11,3) +ez11)/9;
    hx11 = (np.roll(hx11,1) +hx11)/9;
    hy11 = (np.roll(hy11,1) +hy11)/9;

    Ez1[:,-1] = ez11;
    Hx1[:,-1] = hx11;
    Hy1[:,-1] = hy11;

    plt.subplot(1,2,1)
    plt.imshow(Ez1)
    plt.subplot(1,2,2)
    plt.imshow(Ez2)
    plt.pause(0.001)
    plt.clf()