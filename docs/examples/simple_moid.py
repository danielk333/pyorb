"""
Simple moid
===============

"""

import matplotlib.pyplot as plt
import numpy as np
import pyorb
import time

t0 = time.time()

res = 100
orb_1 = pyorb.Orbit(
    M0=pyorb.M_sol,
    degrees=True,
    num=res,
    a=1 * pyorb.AU,
    e=0,
    i=0,
    omega=0,
    Omega=0,
    anom=np.linspace(0, 360, num=res),
    type="true",
)
orb_2 = pyorb.Orbit(
    M0=pyorb.M_sol,
    degrees=True,
    num=res,
    a=2 * pyorb.AU,
    e=0.4,
    i=20,
    omega=180,
    Omega=0,
    anom=np.linspace(0, 360, num=res),
    type="true",
)


dist = np.sum((orb_1.cartesian[:3, :, None] - orb_2.cartesian[:3, None, :]) ** 2, axis=0)
flat_dist = dist.flatten()

flat_ind = np.argmin(flat_dist)
final_ind = np.unravel_index(flat_ind, dist.shape)

dt = time.time() - t0

print(f"Execution time: {dt} s")
print(final_ind, np.sqrt(flat_dist[flat_ind]) / pyorb.AU)

fig, ax = plt.subplots()
ax.pcolormesh(dist)
plt.show()

fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(111, projection="3d")

ax.plot(
    orb_1.cartesian[0, :],
    orb_1.cartesian[1, :],
    orb_1.cartesian[2, :],
    "-b",
)
ax.plot(
    orb_2.cartesian[0, :],
    orb_2.cartesian[1, :],
    orb_2.cartesian[2, :],
    "-g",
)

ax.plot(
    orb_1.cartesian[0, final_ind[0]],
    orb_1.cartesian[1, final_ind[0]],
    orb_1.cartesian[2, final_ind[0]],
    "xr",
)
ax.plot(
    orb_2.cartesian[0, final_ind[1]],
    orb_2.cartesian[1, final_ind[1]],
    orb_2.cartesian[2, final_ind[1]],
    "xr",
)

plt.show()
