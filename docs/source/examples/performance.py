'''
Performance testing
=====================
'''
import timeit

str_format = '{} executions: average {:.3e} seconds per execution'

number = 1000
sizes = [1, 10, 100]

for sz in sizes:
    dt = timeit.timeit(
        'orb.calculate_cartesian()', 
        setup=
f'''
import pyorb
import numpy as np
orb = pyorb.Orbit(M0 = pyorb.M_sol, degrees=True)
orb.add(num={sz},a=np.linspace(0.5,2,num={sz})*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)
''', 
        number = number,
    )
    print('"orb.calculate_cartesian()" Performance')
    print(f'{sz} orbits')
    print(str_format.format(number, dt/number) + '\n')


for sz in sizes:
    dt = timeit.timeit(
        'orb.calculate_kepler()', 
        setup=
f'''
import pyorb
import numpy as np
orb = pyorb.Orbit(M0 = pyorb.M_sol, degrees=True)
orb.add(num={sz},x=np.linspace(0.5,2,num={sz})*pyorb.AU, y=0, z=0, vx=0, vy=2.9785e+04, vz=0)
''', 
        number = number,
    )
    print('"orb.calculate_kepler()" Performance')
    print(f'{sz} orbits')
    print(str_format.format(number, dt/number) + '\n')