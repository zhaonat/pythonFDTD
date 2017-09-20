import numpy as np
import matplotlib.pyplot as plt
## linear interpolation of two points

def interp(r1, r2, x):
    ans = r1[1] + (x-r1[0])*((r2[1]-r1[1])/(r2[0]-r1[0]));
    return ans;

def bilinearinterp():
    #requires knowledge of four points on the plane
    return None;

def trilinearinterp():
    #requires knowledge of 9 points on the space
    return None;

r1 = [1,1];
r2 = [3,3];

test = list();
for i in range(20):
    test.append(interp(r1, r2, i))

plt.plot(test)
plt.show()