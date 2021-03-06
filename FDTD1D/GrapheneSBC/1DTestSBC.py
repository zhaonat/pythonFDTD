import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as nl
import os
directory = 'D:\\Documents\\pythonFDTD\\FDTD1D\\GrapheneSBC\\Movies\\'

try:
    os.stat(directory)
except:
    os.mkdir(directory)

N = 100;

Hz = np.zeros((N+1, ));

x = np.linspace(-N/2,N/2,N+1);
Ey = np.zeros((N+1,))

plt.plot(x, Ey)
Ey = np.squeeze(Ey);
print(Ey.shape)
print(Hz.shape)
dx = 1;
dt = dx;
cc = 1;#dt/dx

## time step is SENSITIVE AT THE BOUNDARY???

plt.ion()
fieldAmplitudes = list();
plt.figure;
abc = (cc*dt/2 - dx)/(cc*dt/2 + dx);
sheet_loc= 50;
source_loc =70;
tau = 0;


## COEFFICIENTS FOR THE SHEET UPDATE
timeAlphaB = 1;

sigma_0 = 1;
denom = (dx * sigma_0) + dt/timeAlphaB;
c1 = (dt) / (dx * sigma_0 + dt/timeAlphaB);
c2 = (dt) / (dx * sigma_0 + dt/timeAlphaB); ## if this is the courant number, there is instability introduced

fe1 = sigma_0 * c1;
fe2 = sigma_0 * c2;
fh11 = 1-2.0*c1;
fh22 = 1-2.0*c2
fh12 = c1;
fh21 = c2;

## What is the courant number
print(str(c1))

Hz_int_right = 0;
Hz_int_left = 0;
Hinterface = list(); F1n = 0; F2n=0;
timeAlpha = 1.0;
tsteps = 200;
E_khalf = list();
for t in range(tsteps):
    print('time step: '+str(t))
    ########################## UPDATE THE HY FIELD ##############################
    Hz_next = np.zeros((N+1,))
    for i in range(1, N+1): #Hz is offset right of Ey
        index = i+1
        if(i == sheet_loc):
            ## Ey is at time n, get the update coefficients for the H interface before we replace Ey with the updated one
            # hz_int_right and hz_int_left are at time n-1/2
            F1n = -2 * fe1 * Ey[sheet_loc] + (fh11 * Hz_int_left) + (
            fh12 * (Hz_int_right));  # Ey is at time step n, but the interfaces are at n-1/2
            F2n = 2 * fe2 * Ey[sheet_loc + 1] + (fh21 * Hz_int_left) + (fh22 * (Hz_int_right));

            Hz_int_left = (1 / (1 - c1 * c2)) * (F1n + c1 * F2n);  # time step to n+1/2 using F1n (at time step n)
            Hz_int_right = (1 / (1 - c1 * c2)) * (F2n + c2 * F1n);
            continue;
        if(index > N):
            index = 0;
        Hz_next[i] = Hz[i] + (dt/dx)*(Ey[index] - (Ey[i]));


    #print(Hz_next)
    Hz_next[0] = Hz[1] + abc*(Hz_next[1] - Hz[0])
    Hz = Hz_next;

    ##################### update Ey field #######################
    #Hz[0] = 0;
    ## iterate through the updates since there's an interface
    Ey_next = np.zeros((N+1,));
    for i in range(0,N+1):
        if(i == sheet_loc+1):
            ## Ey field at sheet_loc and sheet_loc+1 have special updates
            Ey_next[i] = Ey[i] + (dt / dx/timeAlpha) * (Hz[i] - Hz_int_right);
        elif(i == sheet_loc):
            Ey_next[i] = Ey[i] + (dt / dx/timeAlpha) * (Hz_int_left - Hz[i-1]);
        else: # normal update
            index = i-1;
            if(index <0):
                index = N;
            Ey_next[i] = Ey[i] + (dt / dx/timeAlpha) * (Hz[i] - Hz[index]);

    # #E -field abc
    Ey_next[-1] = Ey[-1] + abc * (Ey_next[-2] - Ey[-1])

    #point source pulse
    Ey_next[source_loc] += np.exp(-(t - 20) * (t - 20) / (60*timeAlpha));


    print(str(Hz_int_left)+', '+str(Hz_int_right))
    Hinterface.append([Hz_int_left, Hz_int_right]);
    E_khalf.append((Hz_int_right-Hz_int_left)/sigma_0)
    ## calculate the discontinuity

    ##Hz_interfaces are now updated to n+1/2

    Ey = Ey_next; #update Ey to next timestep


    plt.figure()
    plt.plot(x,Hz)
    plt.plot(x,Ey)
    plt.ylim([-1, 1])
    plt.savefig(directory+'Fields_'+str(t)+'.jpg');
    plt.pause(0.000001)
    plt.clf()


Hinterface = np.array(Hinterface);
plt.figure()
# plt.plot(Hinterface);
plt.plot(E_khalf)
plt.show()
