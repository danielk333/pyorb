'''
Create project logo
====================

This example shows how the project logo was created.

'''

import pyorb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

orb = pyorb.Orbit(
    M0=1,
    G=pyorb.get_G('AU', 'Msol', 'y'),
    num=300,
    degrees=True,
    a=1, 
    e=0.5, 
    i=0,
    omega=0, 
    Omega=0, 
    anom=np.linspace(0, 360, num=300),
)

steps = 10
incs = np.linspace(0, 69, steps)
omegas = np.linspace(0, 90, steps)

cmap = cm.get_cmap('rainbow')

fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(111, projection='3d')

for ind in range(steps):
    orb.i = incs[ind]
    orb.omega = omegas[ind]
    r = orb.r
    ax.plot(r[0, :], r[1, :], r[2, :], '-', color=cmap(float(ind)/steps))

# Hide grid lines
ax.grid(False)

# Hide axes ticks
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.set_axis_off()
ax.view_init(9, -148)

# fig.savefig('logo.svg', transparent=True, bbox_inches='tight')

plt.show()
