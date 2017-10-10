import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as nl
N = 200;

Hz = np.zeros((N+1, ));
IC = lambda x: np.cos(np.pi*x/10)
x = np.linspace(-50,50,N+1);
Ey = np.zeros((N+1,))
# for i in range(len(x)):
#     if(abs(x[i])<= 5):
#         Ey[i] = np.cos(x[i]*np.pi/10)**2;
#     else:
#         Ey[i] = 0;
# plt.plot(x, Ey)
Ey = np.squeeze(Ey);
print(Ey.shape)
print(Hz.shape)
dt = 2;
dx = 2;
cc = 1;#dt/dx

eps = np.ones((N+1,));
eps[125:150] = 4;

tsteps = 200;
plt.ion()
fieldAmplitudes = list();
plt.figure;
abc = (cc*dt - dx)/(cc*dt + dx);
mu_0 = 1;
for t in range(tsteps):
    print('time step: '+str(t))
    #update Hz field
    deriv = Ey[1:]-np.roll(Ey,1)[1:]
    deriv = Ey[0]+deriv;
    Hz_next = Hz + (dt/dx)*(Ey - np.roll(Ey, 1))

    Hz_next[0] = Hz[1] + abc*(Hz_next[1] - Hz[0])
    Hz = Hz_next;

    delay = dx/2+dt/2
    Hz[99] += -np.exp(-(t-30)**2/100)

    Dy = np.multiply(eps, Ey);
    #update Ey field
    Dy_next = Dy + (dt / dx) * (np.roll(Hz, -1) - Hz)

    #E -field abc
    Dy_next[-1] = Dy[-2] + abc * (Dy_next[-2] - Dy[-1])

    Ey = np.divide(Dy_next, eps);
    #point source pulse
    Ey[100] += np.exp(-(t - 30.+delay)**2/ 100);



    plt.plot(x,Hz)
    plt.plot(x,Ey)
    plt.ylim([-0.5, 1])
    plt.pause(0.0001)
    plt.clf()
    # plt.plot(fieldAmplitudes)
    # plt.pause(0.05)
    # plt.clf()