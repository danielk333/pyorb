Getting started tutorial
===========================

Code

.. code-block:: python
    :linenos:

    import pyorb

    #We first create a standard orbit around the sun in SI units
    orb = pyorb.Orbit(M0 = pyorb.M_sol)

    #Lets switch to degrees for more human readable units, this can also be given 
    # at orbit creation as a keyword parameter
    orb.degrees = True

    #Currently the orbit has no values
    print(orb)

.. execute_code:: python
    :hide_code:

    import pyorb
    orb = pyorb.Orbit(M0 = pyorb.M_sol)
    orb.degrees = True
    print(orb)

Code

.. code-block:: python
    :linenos:

    #give it a circular orbit in the plane
    orb.update(a=1*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)
    print(orb)
    print(f'Orbital period: {orb.period/(3600.0*24)} days')

.. execute_code:: python
    :hide_code:

    import pyorb
    orb = pyorb.Orbit(M0 = pyorb.M_sol)
    orb.degrees = True
    orb.update(a=1*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)
    print(orb)
    print(f'Orbital period: {orb.period/(3600.0*24)} days')

Code

.. code-block:: python
    :linenos:

    #Now as soon as we try to look at any cartesian elements 
    # the orbit will be transformed to cartesian space and the 
    # Cartesian elements are stored
    print(f'Orbit X-position is {orb.x*1e-3} km')
    print(f'Orbit velocity vector is {orb.v*1e-3} km/s')
    print(orb)

.. execute_code:: python
    :hide_code:

    import pyorb
    orb = pyorb.Orbit(M0 = pyorb.M_sol)
    orb.degrees = True
    orb.update(a=1*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)
    print(f'Orbit X-position is {orb.x*1e-3} km')
    print(f'Orbit velocity vector is {orb.v*1e-3} km/s')
    print(orb)

Code

.. code-block:: python
    :linenos:

    #However, if we change one of the cartesian variables
    orb.x += 0.1*pyorb.AU

    #A flag will be raised in the class internally that 
    # the kepler elements needs recalculation

    #Converting a orbit instance to a string is 
    # intentionally **not** triggering a re-calculation
    print(orb)
    print('\n')

    #When we then try to get the kepler elements they are automatically recalculated 
    print(f'New kepler elements {orb.kepler}')
    print(orb)

.. execute_code:: python
    :hide_code:

    import pyorb
    orb = pyorb.Orbit(M0 = pyorb.M_sol)
    orb.degrees = True
    orb.update(a=1*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)
    a = orb.x
    orb.x += 0.1*pyorb.AU

    print(orb)
    print('\n')
    print(f'New kepler elements {orb.kepler}')
    print(orb)

As a standard, SI units are used but we can also create orbits with an arbitrary system of units

.. execute_code:: python
    :linenos:

    import pyorb

    #Some combinations are implement as standard, otherwise just pass a float 
    # that describes the conversion between SI and your unit of choice
    G_au = pyorb.get_G(length='AU', mass='kg', time='s')
    print(f'SI gravitation constant: {pyorb.G} m^3 kg^-1 s^-2')
    print(f'Alternative gravitation constant: {G_au} AU^3 kg^-1 s^-2')


    orb2 = pyorb.Orbit(M0 = pyorb.M_sol, G = G_au)
    orb2.update(a=1, e=0, i=0, omega=0, Omega=0, anom=0)

    #To calculate cartesian elements without trying to access any of them simply 
    # call the calculate_cartesian function
    orb2.calculate_cartesian()

    #Now we see that both the velocity and positions have changed to AU and AU/s
    print(orb2)

    #However, if we look at the orbital period, it is still given in seconds
    print(f'Orbital period: {orb2.period/(3600.0*24)} days\n')


    #We can also change this on the fly
    #A common system of units for dynamical astronomy is 
    # Astronomical units-Solar masses-years
    G_ast = pyorb.get_G(length='AU', mass='Msol', time='y')
    print(f'Astronomical gravitation constant: {G_ast} AU^3 Msol^-1 y^-2')
    orb2.G = G_ast

    #We also need to update the central mass
    orb2.M0 = 1.0

    #Since Kepler elements only have one variable with a physical quantity,
    # the semi-major-axis, this change only affects the cartesian elements.
    #Therefore we should recalculate the cartesian based on the current Keplerian
    orb2.calculate_cartesian()

    print(orb2)
    print(f'Orbital period: {orb2.period} years')

    #The orbital speed should be approximately 2pi as this is the 
    # circumference of a circle with radius 1 AU in units of AU
    print(f'Orbital speed: {orb2.speed} AU/y')

