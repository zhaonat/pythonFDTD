## step 1, how do we visualize a nonuniform grid
import numpy as np;
import matplotlib.pyplot as plt

fullxrange = [-400, 400];
lambda_0 = 800;
dx_coarse = lambda_0/100;
dx_fine = dx_coarse/2; ## refinement factor of 2;
c = 1;

#there are two time steps...
dt_coarse = dx_coarse;
dt_fine = dx_fine;

finexrange = [-200, 200];

Nx = int(800/dx_coarse);
Nf = int(200/dx_fine);

xcoarse = np.linspace(-400, 400, Nx+1)
xfine = np.linspace(-200, 200, Nf+1);

borders = [-200, 200];
bordernodes = [25, 100];

Hz = np.zeros(Nx+1);
Ey = np.zeros(Nx+1);
ey = np.zeros(Nf+1);
hz = np.zeros(Nf+1);
## we need to be able to map the spatial positions of the coarse and fine grids


tsteps = 200;
for t in range(tsteps): #odd steps for main grid, even for subgrid
    # update Hz field
    if(t%2 == 1):
        print(t)
        Hz_next = np.zeros((Nx + 1,))
        deriv = Ey[1:] - np.roll(Ey, 1)[1:]
        deriv = Ey[0] + deriv;
        Hz_next = Hz + (dt_coarse / dx_coarse) * (Ey - np.roll(Ey, 1))

        Hz = Hz_next;

        # update Ey field
        Ey0 = Ey[-2]
        Ey_next = np.zeros((Nx + 1,));
        Ey_next = Ey + (dt_coarse / dx_coarse) * (np.roll(Hz, -1) - Hz)
        Ey = Ey_next;


        Ey[50] += np.exp(-(t - 30.) * (t - 30.) / 100);

        ## Once we do the update, we need to analyze the border points between fine and coarse
        # in 1D, every subgrid point is surrounded by two coarse points..
        for i in range(Nf-1):
            ey[i] = (Ey[24+i]+Ey[25+i])/2
            hz[i] = (Hz[24+i]+Hz[25+i])/2;


        ##apparently, there's a time interpolation step?
        #step ey and hz by a quarter steps
        hz_next = np.zeros((Nx + 1,))
        deriv = ey[1:] - np.roll(ey, 1)[1:]
        deriv = ey[0] + deriv;
        hz_next = hz + (dt_coarse / dx_coarse) * (ey - np.roll(ey, 1))
        hz = hz_next;

        # update Ey field
        ey0 = ey[-2]
        ey_next = np.zeros((Nf + 1,));
        ey_next = ey + (dt_fine / dx_fine) * (np.roll(hz, -1) - hz)
        ey = ey_next;


        #step Ey and Hz by a quarter step;
        Hz_next = np.zeros((Nx + 1,))
        deriv = Ey[1:] - np.roll(Ey, 1)[1:]
        deriv = Ey[0] + deriv;
        Hz_next = Hz + (dt_fine / dx_fine) * (Ey - np.roll(Ey, 1))

        Hz = Hz_next;

        # update Ey field
        Ey0 = Ey[-2]
        Ey_next = np.zeros((Nx + 1,));
        Ey_next = Ey + (dt_fine / dx_fine) * (np.roll(Hz, -1) - Hz)
        Ey = Ey_next;


        Ey[50] += np.exp(-(t - 30.) * (t - 30.) / 100);


        ## at the end, everything has been stepped forward one coarse time step
        plt.plot(xcoarse, Hz)
        plt.hold(True)
        plt.plot(xfine, hz)
        plt.ylim([-0.5, 1])
        plt.pause(0.001)
        plt.clf()
    else:
        None;
    # plt.plot(fieldAmplitudes)
    # plt.pause(0.05)
