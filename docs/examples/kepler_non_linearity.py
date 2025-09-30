"""
The non-linearity of Kepler elements
====================================
"""

import pyorb

import numpy as np
import matplotlib.pyplot as plt

num = 50
vx, vy, vz = np.meshgrid(
    np.linspace(-1, 1, num=num) * 1e4,
    np.linspace(-1, 1, num=num) * 1e4,
    np.linspace(-1, 1, num=num) * 1e4,
)
orb = pyorb.Orbit(
    M0=pyorb.M_sol,
    num=num**3,
    x=pyorb.AU / np.sqrt(2),
    y=pyorb.AU / np.sqrt(2),
    z=0,
    vx=vx.flatten(),
    vy=vy.flatten(),
    vz=vz.flatten(),
    degrees=True,
)
v = np.sqrt(vx**2 + vy**2 + vz**2).flatten()


fig, axes = plt.subplots(1, 3, figsize=(10, 6))
axes[0].scatter(orb.vx, orb.vy, s=1, c=v)
axes[0].set_xlabel("Vx")
axes[0].set_ylabel("Vy")

axes[1].scatter(orb.vx, orb.vz, s=1, c=v)
axes[1].set_xlabel("Vx")
axes[1].set_ylabel("Vz")

axes[2].scatter(orb.vy, orb.vz, s=1, c=v)
axes[2].set_xlabel("Vy")
axes[2].set_ylabel("Vz")


fig, axes = plt.subplots(2, 3, figsize=(10, 6))
axes[0, 0].scatter(orb.a / pyorb.AU, orb.e, s=1, c=v)
axes[0, 0].set_xlabel("Semi-major axis [AU]")
axes[0, 0].set_ylabel("Eccentricity [1]")

axes[0, 1].scatter(orb.i, orb.omega, s=1, c=v)
axes[0, 1].set_xlabel("Inclination [deg]")
axes[0, 1].set_ylabel("Arg. of periapais [deg]")

axes[0, 2].scatter(orb.Omega, orb.anom, s=1, c=v)
axes[0, 2].set_xlabel("Long. of asc. node [deg]")
axes[0, 2].set_ylabel("True anomaly [deg]")

axes[1, 0].scatter(orb.e, orb.i, s=1, c=v)
axes[1, 0].set_xlabel("Eccentricity [1]")
axes[1, 0].set_ylabel("Inclination [deg]")

axes[1, 1].scatter(orb.a / pyorb.AU, orb.i, s=1, c=v)
axes[1, 1].set_xlabel("Semi-major axis [AU]")
axes[1, 1].set_ylabel("Inclination [deg]")

axes[1, 2].scatter(orb.e, orb.omega, s=1, c=v)
axes[1, 2].set_xlabel("Eccentricity [1]")
axes[1, 2].set_ylabel("Arg. of periapais [deg]")

plt.show()
