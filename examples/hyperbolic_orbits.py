'''
Hyperbolic orbits
==================

Hyperbolic orbits are also supported. However their anomalies function a bit differently then the anomalies for elliptical orbits.


'''
import matplotlib.pyplot as plt
import numpy as np

import pyorb

orb = pyorb.Orbit(
    M0 = 1.0,
    G = pyorb.get_G(length='AU', mass='Msol', time='y'),
    a = 1.0, 
    e = 1.2, 
    i = 0, 
    omega = 0, 
    Omega = 0, 
    anom = 0.0,
    degrees = True,
    type = 'true',
)

print(orb)
print(f'Orbit anomaly type: {orb.type}')

t_num = 1000
dt = 0.01 #y
r = np.empty((3,t_num))

for ti in range(t_num):
    r[:,ti] = np.squeeze(orb.r)
    orb.propagate(dt)


fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111, projection='3d')

ax.plot(r[0,:], r[1,:], r[2,:],  '-b')
ax.plot([0], [0], [0], 'or')

ax.set_title('Hyperbolic orbit', fontsize=22)
ax.set_xlabel('X-position [AU]', fontsize=20)
ax.set_ylabel('Y-position [AU]', fontsize=20)
ax.set_zlabel('Z-position [AU]', fontsize=20)
plt.show()