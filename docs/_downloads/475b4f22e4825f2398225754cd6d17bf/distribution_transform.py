'''
Distribution transformation
============================
'''

import pyorb

import numpy as np
import matplotlib.pyplot as plt

#turn on TeX interperter
plt.rc('text', usetex=True)

#for reproducibility
np.random.seed(12398748)

#We first create a standard orbit around the sun in SI units
orb = pyorb.Orbit(M0 = pyorb.M_sol)

#Create 1000 equal orbits
orb.add(num=1000, a=pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)
#calculate cartesian elements
orb.calculate_cartesian()

#Add a symmetric Gaussian distribution in the plane on the velocity
std = 3e3
orb.vx += np.random.randn(orb.num)*std
orb.vy += np.random.randn(orb.num)*std

#now when we call any Keplerian element, the distribution in kepler space will be calculated automatically

fig, axes = plt.subplots(1, 2, figsize=(10,6))
axes[0].plot(orb.vx*1e-3, orb.vy*1e-3,'.')
axes[0].set_title('Cartesian velocity', fontsize=22)
axes[0].set_xlabel('X-velocity $v_x$ [km/s]', fontsize=20)
axes[0].set_ylabel('Y-velocity $v_y$ [km/s]', fontsize=20)

axes[1].plot(orb.a/pyorb.AU, orb.e, '.')
axes[1].set_title('Keplerian orbital shape', fontsize=22)
axes[1].set_xlabel('Semi-major axis $a$ [AU]', fontsize=20)
axes[1].set_ylabel('Eccentricity $e$ [1]', fontsize=20)

plt.show()