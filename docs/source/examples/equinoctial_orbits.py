'''
Equinoctial orbits
====================
'''

import pyorb

orb = pyorb.Orbit(M0 = pyorb.M_sol, degrees=True)

# equinoctial orbit
orb.update(a=1*pyorb.AU, e=0.01, i=0.9, omega=24, Omega=0, anom=22)
print("Original orbit object")
print(orb)

print("Equinoctial elements")
equi = orb.equinoctial
for ind, var in enumerate(orb.EQUINOCTIAL):
    print(f'{var:<2}: {equi[ind]}')

# change p variable

equi[3] += 0.1

# set values, kepler and cart automatically updated
orb.equinoctial = equi
print("New orbit object")
print(orb)

print("New equinoctial elements")
equi = orb.equinoctial
for ind, var in enumerate(orb.EQUINOCTIAL):
    print(f'{var:<2}: {equi[ind]}')

# Call transform manually

kep = orb.kepler[:, 0]
kep[pyorb.const.K_nu] = orb.mean_anomaly[0]

print("Manually new calculated elements")
equi_again = pyorb.kep_to_equi(kep, degrees=True)
for ind, var in enumerate(orb.EQUINOCTIAL):
    print(f'{var:<2}: {equi_again[ind]}')
