import numpy as np
import matplotlib.pyplot as plt
## 1D version
N = 209;
Hx = np.zeros((N + 1,));
IC = lambda x: np.cos(np.pi*x/20)
x = np.linspace(-40,40,N+1);
Ey = np.zeros((N+1,))

dx = 1;
cc = 1;#dt/dx
dt = dx;
Npmlx = 10;
Nx = len(Hx);

xn = np.zeros((Nx,1))
fi1 = np.zeros((Nx,1))
gi2 = np.ones((Nx,1))
gi3 = np.ones((Nx,1));
fi10 = np.zeros((Npmlx, 1))
gi20 = np.zeros((Npmlx, 1))
gi30 = np.zeros((Npmlx, 1))

for i in range(Npmlx):
    x_n = (1/3)*((i+1)/Npmlx)**3;
    xn[Nx-Npmlx+i] = x_n;
    xn[Npmlx-1-i] = x_n;

    fi1[Npmlx-1-i] = x_n;
    gi2[Npmlx-1-i] = 1/(1+x_n);
    gi3[Npmlx-1-i] = (1-x_n)/(1+x_n);

    fi10[i] = x_n;
    gi20[i] = 1/(1+x_n);
    gi30[i] = (1-x_n)/(1+x_n);

    fi1[Nx-Npmlx+i] = x_n;
    gi2[Nx-Npmlx+i] = 1/(1+x_n);
    gi3[Nx-Npmlx+i] = (1-x_n)/(1+x_n);

plt.plot(fi1)
plt.plot(gi2);
plt.plot(gi3)
plt.show()
gi3 = np.squeeze(gi3);
gi2 = np.squeeze(gi2);

## apply it to the fdtd


plt.ion()
fieldAmplitudes = list();
plt.figure;
abc = (cc*dt - dx)/(cc*dt + dx);
tsteps = 1000;
for t in range(tsteps):
    print(t)
    #J = np.exp(-(t-40)**2/100);
    J = 0.5*np.sin(2*np.pi*t/50)
    print(J)
    #update Hz field
    #Ey[-1] = 0;
    Hx_next = np.multiply(gi3,Hx) - (dt / dx) * \
                                    np.multiply(gi2, (Ey - np.roll(Ey, 1)))
    # for i in range(0, N):
    #     Hz_next[i] = Hz[i] + (dt/dx)*(Ey[i+1] - (Ey[i]));
    #H field abc, #Hz[0] is never updated during the loop
    #Hz_next[0] = Hz[1] + abc*(Hz_next[1] - Hz[0])
    Hx = Hx_next;

    #update Ey field
    #Hx[0] = 0;
    Ey_next = np.multiply(gi3,Ey) - \
              (dt/dx)*np.multiply(gi2,(np.roll(Hx, -1) - Hx))

    # for i in range(1, N+1):
    #     Ey_next[i] = Ey[i] + (dt/dx)*(Hz[i] - Hz[i-1]);
    #E -field abc
    #Ey_next[N] = Ey0 +abc*(Ey_next[-2] - Ey[-1])
    #point source pulse
    #Ey[50] += np.exp(-(t - 30.) * (t - 30) / 100);
    Ey = Ey_next;
    Ey[110] -= J;

    if(t%10 == 0):
        plt.plot(x,Ey)
        plt.ylim([-1, 1])
        plt.pause(0.00001)
        plt.clf()
    # plt.plot(fieldAmplitudes)
    # plt.pause(0.05)
    # plt.clf()