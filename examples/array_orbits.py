import pyorb

import numpy as np

#
# THE MOST IMPORTANT THING
# - properties act on ALL orbits in the class
# - Only way to update individual orbits of a set is to use .update
# - iterations are passive, the objects are copies from the array so the array itself is NOT modified


orb = pyorb.Orbit(M0 = pyorb.M_sol, degrees=True)

orb.allocate(10)
orb.update(a=np.linspace(0.5,2,num=10)*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)
# print(orb)

orb.add(a=4*pyorb.AU)
orb.add(num=4, a=1*pyorb.AU, e=np.linspace(0,0.2,num=4))
print(orb.a/pyorb.AU)
print(orb.e)
print(orb)
print(orb[14])

print(orb._cart[0,:])
print(orb.x)

print(orb.eccentric_anomaly)
orb.eccentric_anomaly = np.linspace(0,360,num=len(orb))
print(orb.anom)

# orb.update(inds=[2,3,4], anom=[0,45,90])
# for o in orb[2:5]:
#     print(o.r)

# for o in orb[1:2]:
#     print(o)

# print(orb.x)

# for o in orb[3:5]:
#     print(o)


# orb = pyorb.Orbit(M0 = pyorb.M_sol, degrees=True)
# orb.update(a=1*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=45)

# print(orb)