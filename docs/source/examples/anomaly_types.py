'''
Anomaly types
===============

To get another anomaly, either change anomaly type ``orb.type = 'mean'`` or use the property ``orb.mean_anomaly``.
'''
import matplotlib.pyplot as plt
import numpy as np

import pyorb

num = 500
[ecc, anom] = np.meshgrid(
    np.linspace(0,0.99,num=num), 
    np.linspace(0,360,num=num),
)

orb = pyorb.Orbit(
    M0 = pyorb.M_sol,
    num = num**2,
    a = 1*pyorb.AU, 
    e = ecc.reshape(num**2), 
    i = 0, 
    omega = 90, 
    Omega = 0, 
    anom = anom.reshape(num**2),
    degrees = True,
    type = 'true'
)

true = orb.true_anomaly.reshape(num,num)
mean = orb.mean_anomaly.reshape(num,num)

print(orb)
print(f'Orbit anomaly type: {orb.type}')

fig, ax = plt.subplots(1,1)

c = ax.pcolormesh(ecc, true, mean-true)
ax.set_xlabel('Eccentricity [1]')
ax.set_ylabel('True anomaly [deg]')
ax.set_title('Difference between True and mean anomaly [deg]')
fig.colorbar(c)

plt.show()