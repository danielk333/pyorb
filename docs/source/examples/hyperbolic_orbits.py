'''
Hyperbolic orbits
==================

Hyperbolic orbits are also supported. However their anomalies function a bit differently then the anomalies for elliptical orbits.

'''
import matplotlib.pyplot as plt
import numpy as np

import pyorb

font_sz = 16

# Get the asymptotic true anomaly
theta_inf = pyorb.kepler.true_of_the_asymptote(e=1.2, degrees=True)

# dont go to the asymptotic value
theta_inf *= 0.9

num = 1000
orb = pyorb.Orbit(
    M0 = 1.0,
    G = pyorb.get_G(length='AU', mass='Msol', time='y'),
    a = 1.0, 
    e = 1.2, 
    i = 0, 
    omega = 0, 
    Omega = 0, 
    anom = np.linspace(-theta_inf, theta_inf, num=num),
    degrees = True,
    num = num,
    type = 'true',
)

print(orb)
print(orb[num//2])

r = orb.r

fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(121, projection='3d')

ax.plot(r[0, :], r[1, :], r[2, :], '-b')
ax.plot([0], [0], [0], 'or')

ax.set_title('Hyperbolic orbit', fontsize=font_sz + 2)
ax.set_xlabel('X-position [AU]', fontsize=font_sz)
ax.set_ylabel('Y-position [AU]', fontsize=font_sz)
ax.set_zlabel('Z-position [AU]', fontsize=font_sz)

ax = fig.add_subplot(122)
ax.plot(orb.anom, orb.velocity, '-b')
ax.set_title('Hyperbolic velocity', fontsize=font_sz + 2)
ax.set_xlabel('True anomaly [deg]', fontsize=font_sz)
ax.set_ylabel('Velocity [AU/y]', fontsize=font_sz)

plt.show()
