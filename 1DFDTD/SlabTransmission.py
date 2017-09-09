import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as nl


mu_0 = 4*np.pi*10**-7;
eps_0 = 8.85e-12;

def runFDTD(N, dt, dx, c, tsteps, eps, x):
    ## J is the current source
    abc = (c*dt-dx)/(c*dt+dx);
    Ey = np.zeros(N,);
    Hz = np.zeros(N,)
    ## initial condition on E
    #initial field profile has to be continuous!!!
    # for i in range(len(Ey)):
    #
    #     if ( i > 200 and i < 300):
    #         Ey[i] = np.cos((x[i]-250) * np.pi / 100) ** 2;
    #     else:
    #         Ey[i] = 0;
    # plt.figure;
    # plt.plot(Ey);
    # plt.show()
    plt.ion()
    fieldAmplitudes = list();
    plt.figure;
    omega = 2*np.pi/50;
    # update Ey field
    lambda_0 = 800;
    lambda_U = 1200;
    lambda_L = 600;
    sigma = (2 / omega) * (lambda_0 / (lambda_U + lambda_L))
    print(sigma)
    for t in range(tsteps):
        # update Hz field
        #D = np.multiply(eps, Ey);
        D = Ey;
        Hz_next = Hz + (dt*c / dx) * (D - np.roll(D, 1))
        Hz_next[0] = Hz[1] + abc*(Hz_next[1] - Hz[0])
        Hz = Hz_next;


        J = np.sin(omega*t)*np.exp(-(t-4*sigma)**2/sigma**2)
        print(J)
        Ey_next = Ey + ((dt *c/ dx)/eps) * (np.roll(Hz, -1) - Hz);
        Ey_next[100]-= (dt/1)*J;
        # E -field abc
        Ey_next[-1] = Ey[-2] +abc*(Ey_next[-2] - Ey[-1])

        Ey = Ey_next;
        fieldAmplitudes.append(nl.norm(Ey) ** 2 + nl.norm(Hz) ** 2)

        plt.plot(x, Ey)
        plt.ylim([-0.5, 1])
        plt.pause(0.1)
        plt.clf()

    return [Ey, Hz,J,eps]

def dielSlab(N, epsilon, xstart, xend, left, dx):
    eps = np.ones((N+1,));
    Nstart = locToGrid(0, dx, left);
    Nend = locToGrid(700, dx, left);
    print(str(Nstart)+' '+str(Nend))
    eps[Nstart:Nend] = epsilon;

    ## apply the staircase averaging
    for i in range(1, N-1):
        eps[i] = (eps[i-1]+eps[i+1])/2;
    return eps;

#map spatial location ot grid location
def locToGrid(loc, dx, left):
    x = int(abs(left-loc)/dx);
    print(int(abs(left-x)));
    return x;



## Run sim
xstart = -300;
xend = 1000;
dt = 1; dx = 1; c = 1;
N = int((xend-xstart)/dx);
x = np.linspace(xstart, xend, N)
tsteps = 300;
dielStart = 0; dielEnd = 700;

eps = dielSlab(N-1, 4, dielStart, dielEnd, xstart, dx);
# plt.figure()
# plt.plot(eps)
# plt.show()
[Ey, Hz,J, eps] = runFDTD(N, dt, dx, c, tsteps, eps, x);
plt.plot(Ey);
plt.plot(Hz);
plt.show()

