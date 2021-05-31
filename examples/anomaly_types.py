'''
Anomaly types
===============

To get another anomaly, either change anomaly type :code:`orb.type = 'mean'` or use the property
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
    auto_update = True, 
    direct_update = True,
    type = 'true'
)

true = orb.true_anomaly.reshape(num,num)
mean = orb.mean_anomaly.reshape(num,num)

print(orb)
print(f'Orbit anomaly type: {orb.type}')

fig, ax = plt.subplots(1,1)

c = ax.pcolormesh(ecc, true, mean-true)
fig.colorbar(c)

plt.show()