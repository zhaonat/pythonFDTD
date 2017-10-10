## step 1, how do we visualize a nonuniform grid
import numpy as np;
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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

hz_storage = list();
ey_storage = list();
for i in range(10):
    hz_storage.append(hz);
    ey_storage.append(ey);
Ey_prev = Ey;

#storage fields will contain fields for both coarse and fine grid
# at QUARTER STEPS
for t in range(8,tsteps): #odd steps for main grid, even for subgrid
    # update Hz field on coarse grid
    print(t)
    Ey[0] = 0;
    Hz_next = Hz + (dt_coarse / dx_coarse) * (Ey - np.roll(Ey, 1))
    Hz = Hz_next;
    Hz_storage.append(Hz);

    # update Ey field on coarse grid t--> t+dt
    Hz[0] = 0; #dirichlet
    Ey_next = Ey + (dt_coarse / dx_coarse) * (np.roll(Hz, -1) - Hz)
    Ey = Ey_next;

    # source injection
    Ey[50] += np.exp(-(t - 30.) * (t - 30.) / 100);
    Ey_storage.append(Ey);

    ## Once we do the update, we need to analyze the border points between fine and coarse
    # now we transfer information to the fine grid

    hz[0] = (2/3)*(Hz[-1])+(1/3)*(hz[1]) #- fine grid values t+ dt updated from t
    ey[0] = (2/3)*(Ey[-1])+(1/3)*(ey[1]) #- fine grid interpolation at t + dt/2 updated from t-dt/2

    ##time interpolate subgrid fields at t+dt/2 for hz and t+dt/4 for ey
    tsamples = [-4,0,8]; y = [ey_storage[t-4][0],ey_storage[t][0],ey_storage[t-8][0]];
    f = interp1d(tsamples,y, kind='quadratic')

    ## apply the interpolation in time to the border points
    es_tminushalf = f(-2);
    hb_tminusquarter =f(-2);

    # ## use Yee update to get a different value of the fine grids at t+dt/2 and t+dt/4
    hb_tminusquarter2 = hz_storage[t-4][0] - (Ey[1]-(es_tminushalf));
    #create weighted average
    hb_tminusquarter = hb_tminusquarter*0.65+0.35*hb_tminusquarter2

    ey[0] = es_tminushalf;
    hz[0] = hb_tminusquarter;

    ##### YEE UPDATES NO INTERPOLATION
    # update the fine grid values for the electric field
    hz[0] = 0; #ey
    ey_next = ey + (dt_fine / dx_fine) * (np.roll(hz, -1) - hz)
    ey = ey_next;
    ey_storage.append(ey) #ey at t

    #update hz field to h_int to t+0.25;
    ey[-1] = 0;  # dirichlet
    hz_next = hz + (dt_fine / dx_fine) * (ey - np.roll(ey, 1))
    hz = hz_next;
    hz_storage.append(hz) #at t+0.25;

    ## NOW THAT WE HAVE RUN THE UPDATE on the fine grid,
    ## we reverse interpolate back onto the coarse grid
    ## update all H coarse files at t+0.25;
    Hz[-1] = hz[0]*(1/3)+(2/3)*Hz[-2];

    ## now we do another quartertime step to get e to t+1, the time step of HZ
    ey = ey +(np.roll(hz,-1)-hz);
    hz = hz + (dt_fine / dx_fine) * (ey - np.roll(ey, 1))

    ## update coarse E fields and H fields via interpolation again
    # Hz[-1] = hz[0]*(1/3)+(2/3)*Hz[-1];
    # Ey[-1] = ey[0]*(1/3)+(2/3)*Ey[-2];

    ## why do we need the H fields at the quarter time steps???

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
