## step 1, how do we visualize a nonuniform grid
import numpy as np;
import matplotlib.pyplot as plt
def interp(r1, r2, x):
    ans = r1[1] + (x-r1[0])*((r2[1]-r1[1])/(r2[0]-r1[0]));
    return ans;

fullxrange = [-400, 400];
lambda_0 = 800;
dx_coarse = lambda_0/100;
dx_fine = dx_coarse/2; ## refinement factor of 2;
c = 1;

#there are two time steps...
dt_coarse = dx_coarse;
dt_fine = dx_fine;

finexrange = [400, 800];

Nx = int(400/dx_coarse);
Nf = int(400/dx_fine);

xcoarse = np.linspace(-400, 400, Nx+1)
xfine = np.linspace(400,800, Nf+1);


Hz = np.zeros(Nx+1);
Ey = np.zeros(Nx+1);

## there is only one border points Ey[end] and ey[0]
ey = np.zeros(Nf+1);
hz = np.zeros(Nf+1);
## we need to be able to map the spatial positions of the coarse and fine grids


tsteps = 2000;
Hz_storage = list(); Hz_storage.append(Hz);
Ey_storage = list(); Ey_storage.append(Ey);

## we need to put in 3 sets of initial values...
hz_storage = list();
ey_storage = list();
for i in range(9):
    hz_storage.append(hz);
    ey_storage.append(ey);
Ey_prev = Ey;
for t in range(8,tsteps): #odd steps for main grid, even for subgrid
    # update Hz field on coarse grid
    print(t)
    Hz_next = Hz + (dt_coarse / dx_coarse) * (Ey - np.roll(Ey, 1))
    Hz = Hz_next;
    Hz_storage.append(Hz);

    # update Ey field on coarse grid
    Hz[0] = 0; #dirichlet
    Ey_next = Ey + (dt_coarse / dx_coarse) * (np.roll(Hz, -1) - Hz)
    Ey = Ey_next;

    # source injection
    #Ey[50] += np.exp(-(t - 30.) * (t - 30.) / 100);
    Ey_storage.append(Ey);
    ## Once we do the update, we need to analyze the border points between fine and coarse
    # now we transfer information to the fine grid
    hz[0] = (2/3)*(Hz[-1])+(1/3)*(hz[0]) #- fine grid values t-1
    ey[0] = (2/3)*(Ey[-1])+(1/3)*(ey[0]) #- fine grid interpolation at t -.5
    #hz[0] = hz_interp;
    #ey[0] = ey_interp;
    ##subgrid fields at t and t+0.5
    ##apparently, there's a time interpolation step?
    # es_tminushalf = (ey_storage[t-8][0]+ey_storage[t-4][0]+ey[0])/3;
    # hb_tminusquarter = (hz_storage[t-8][0]+hz_storage[t-4][0]+hz[0])/3;

    # ## use Yee update to get a different value of hb_tminusquarter
    # hb_tminusquarter2 = hz_storage[t-4][0] - (Ey[1]-(es_tminushalf));
    # hb_tminusquarter = hb_tminusquarter*0.65+0.35*hb_tminusquarter2

    # update ey field on the fine grid to time step t, which is the same time step as Ez on the coarse grid
    #hz[0] = hb_tminusquarter;
    hz[0] = 0
    ey_next = ey + (dt_fine / dx_fine) * (np.roll(hz, -1) - hz)
    ey = ey_next;
    ey[50] += np.exp(-(t - 30.) * (t - 30.) / 100);
    #ey[0] = es_tminushalf;

    #update hz field to t+0.5
    ey[-1] = 0;  # dirichlet
    hz_next = hz + (dt_fine / dx_fine) * (ey - np.roll(ey, 1))
    hz = hz_next;

    ## Now we the ghost values at the time level of t+1/2, which we do another interpolation for
    ey_tplushalf = ey[0] +0.5*(Ey[0]-Ey_prev[0])

    ## now we do another quartertime step to get e to t+1, the time step of HZ
    e_tplushalf = ey +(np.roll(hz,-1)-hz);


    ## update coarse E fields parallel and adjacent to the fine grid using newly updated H coarse fields


    ## store this field as the 'previous' field when we enter the next time step
    Ey_prev = Ey;

    ## at the end, everything has been stepped forward one coarse time step
    if(t%10 == 0):
        plt.plot(xcoarse, Hz)
        plt.hold(True)
        plt.plot(xfine, hz)
        plt.ylim([-0.5, 1])
        plt.pause(0.001)
        plt.clf()

    # plt.plot(fieldAmplitudes)
    # plt.pause(0.05)
