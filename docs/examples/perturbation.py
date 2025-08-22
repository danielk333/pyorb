'''
Perturbation by modifying elements
===================================
'''

import pyorb

import numpy as np
import matplotlib.pyplot as plt

# We first create a standard orbit around the sun in SI units
orb = pyorb.Orbit(
    M0 = pyorb.M_sol,
    a = 1*pyorb.AU,
    e = 0.2,
    i = 0,
    omega = 0,
    Omega = 0,
    anom = 0,
    degrees = True,
)
print(orb)

# prepare a propagation
dt = 3600*24.0
num = 1000
r = np.empty((3, num))

for ti in range(num):
    orb.propagate(dt)

    # We need to calculate all the perturbations before perturbing it

    # Every step, perturb the inclination using the normalized velocity solar-component
    di = dt*1e-6*(orb.v[:, 0] @ orb.r[:, 0])/(orb.velocity*np.linalg.norm(orb.r))

    # Also modify the semi-major axis for a funky looking orbit
    da = (pyorb.AU - np.linalg.norm(orb.r))*dt*1e-7
    orb.i += di
    orb.a += da

    # now when we call orb.r the cartesian elements will be recalculated
    r[:, ti] = np.squeeze(orb.r)/pyorb.AU


fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(111, projection='3d')
ax.plot(r[0, :], r[1, :], r[2, :], '-b')
ax.set_title('Orbits', fontsize=22)
ax.set_xlabel('X-position [AU]', fontsize=20)
ax.set_ylabel('Y-position [AU]', fontsize=20)
ax.set_zlabel('Z-position [AU]', fontsize=20)
ax.view_init(14, -6)
plt.show()
