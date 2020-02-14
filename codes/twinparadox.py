import numpy as np
v0 = 0.3
a = 0.003
tr = v0 / a
nsteps = 10000
dt = tr  * 2.0 / nsteps
tprime = 0.0
t = 0.0
for i in range(nsteps):
    t = t + dt
    v = v0 - a * t
    x = v0 * t - a * t**2 / 2.0
    tprime += np.sqrt(1.0 - v** 2) * dt
    tdiff = tprime - (t-x)
    if(tdiff < 0):
        print(t - x)
        break
