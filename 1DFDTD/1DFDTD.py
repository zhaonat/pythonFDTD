import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as nl
N = 200;

Hz = np.zeros((N+1, ));
IC = lambda x: np.cos(np.pi*x/10)
x = np.linspace(-50,50,N+1);
Ey = np.zeros((N+1,))
for i in range(len(x)):
    if(abs(x[i])<= 5):
        Ey[i] = np.cos(x[i]*np.pi/10)**2;
    else:
        Ey[i] = 0;
plt.plot(x, Ey)
Ey = np.squeeze(Ey);
print(Ey.shape)
print(Hz.shape)
dt = 2;
dx = 2;
cc = 1;#dt/dx

tsteps = 60;
plt.ion()
fieldAmplitudes = list();
plt.figure;
abc = (cc*dt - dx)/(cc*dt + dx);

for t in range(tsteps):
    #update Hz field
    Hz_next = np.zeros((N+1,))
    deriv = Ey[1:]-np.roll(Ey,1)[1:]
    deriv = Ey[0]+deriv;
    Hz_next = Hz + (dt/dx)*(Ey - np.roll(Ey, 1))
    # for i in range(0, N):
    #     Hz_next[i] = Hz[i] + (dt/dx)*(Ey[i+1] - (Ey[i]));
    #H field abc, #Hz[0] is never updated during the loop
    #Hz_next[0] = Hz[1] + abc*(Hz_next[1] - Hz[0])
    Hz = Hz_next;

    #update Ey field
    Ey0 = Ey[-2]
    Ey_next = np.zeros((N+1,));
    Ey_next = Ey + (dt/dx)*(np.roll(Hz,-1) - Hz)

    # for i in range(1, N+1):
    #     Ey_next[i] = Ey[i] + (dt/dx)*(Hz[i] - Hz[i-1]);
    #E -field abc
    #Ey_next[N] = Ey0 +abc*(Ey_next[-2] - Ey[-1])
    #point source pulse
    #Ey[50] += np.exp(-(t - 30.) * (t - 30) / 100);
    Ey = Ey_next;

    fieldAmplitudes.append(nl.norm(Ey)**2+nl.norm(Hz)**2)

    plt.plot(x,Hz)
    plt.ylim([-0.5, 1])
    plt.pause(0.1)
    plt.clf()
    # plt.plot(fieldAmplitudes)
    # plt.pause(0.05)
    # plt.clf()