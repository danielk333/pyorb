'''
Array support
===============

Important to note
-------------------
- Properties act on ALL orbits in the class
- Only way to update individual orbits of a set is to use ``self.update`` with the ``inds`` keyword
- Iterations are passive, the objects are copies from the array so the array itself is NOT modified

'''

import pyorb

import numpy as np



orb = pyorb.Orbit(M0 = pyorb.M_sol, degrees=True)

orb.allocate(10)
orb.update(a=np.linspace(0.5,2,num=10)*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)

orb.add(a=4*pyorb.AU)
orb.add(num=4, a=1*pyorb.AU, e=np.linspace(0,0.2,num=4))
print(orb)
print('Orbit semi major axis [AU]:')
print(orb.a/pyorb.AU)

print('\nOrbit eccentricity:')
print(orb.e)
print('One of the orbits were not assigned an eccentricity, hence it is NaN')

print('\nOrbit index 14:')
print(orb[14])
print('All elements not assigned are NaN')

print('\nOrbit cartesian array, all x values:')
print(orb.cartesian[0,:])

print('Numpy warnings will trigger if calculation of a parameter is run and some items are NaN')
print('\nOrbit eccentric anomaly:')
print(orb.eccentric_anomaly)
orb.eccentric_anomaly = np.linspace(0,360,num=len(orb))
print(f'Orbit {orb.type} anomaly after setting:')
print(orb.anom)
print('One of the anomalies are still NaN since to convert from eccentric anomaly to True anomaly one needs an eccentricity')


orb.update(inds=[2,3,4], anom=[0,45,90])
print(f'Orbit {orb.type} anomaly after updating specific objects:')
print(orb.anom)