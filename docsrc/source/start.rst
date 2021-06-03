Getting started tutorial
=========================

We first create a standard orbit around the sun in SI units

.. code-block:: python

    import pyorb
    orb = pyorb.Orbit(M0 = pyorb.M_sol)

Lets switch to degrees for more human readable units, this setting can also be given at orbit creation as a keyword parameter.

.. code-block:: python

    orb.degrees = True

The output from a ``print(orb)`` is

.. code-block:: 

    a    : nan   x : nan
    e    : nan   y : nan
    i    : nan   z : nan
    omega: nan   vx: nan
    Omega: nan   vy: nan
    anom : nan   vz: nan

Since we have not yet set any parameters for this object. Let us give it a circular orbit in the plane at 1 AU (i.e. at approximately the same distance as the Earth form the Sun).


.. code-block:: python

    orb.update(a=1*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)

The output from printing the following

.. code-block:: python

    print(orb)
    print(f'Orbital period: {orb.period/(3600.0*24)} days')

is now

.. code-block:: 

    a    : 1.4960e+11   x : 1.4960e+11
    e    : 0.0000e+00   y : 0.0000e+00
    i    : 0.0000e+00   z : 0.0000e+00
    omega: 0.0000e+00   vx: -0.0000e+00
    Omega: 0.0000e+00   vy: 2.9785e+04
    anom : 0.0000e+00   vz: 0.0000e+00
    Orbital period: [365.25137584] days

As we can see, the orbital period coincides very well with that of an Earth year and the object automatically calculated the cartesian coordinates as well. We can also look at a variety of implemented parameters available as attributes (see API documentation for complete list) such as ``print(f'Orbit velocity is {orb.velocity*1e-3} km/s')`` which would generate:

.. code-block:: 

    Orbit velocity is 29.785142169221025 km/s

Lets say we desire more control over when a conversion from Kepler to Cartesian is performed. We can then disable the direct conversion by

.. code-block:: python

    orb.direct_update = False

Now if we change the semi-major axis to be located at e.g. the orbit of Mars and print the object

.. code-block:: python

    orb.a = 1.5237*pyorb.AU
    print(orb)

We get a change in the kepler elements but not in the cartesian elements:

.. code-block:: 

    a    : 2.2794e+11   x : 1.4960e+11
    e    : 0.0000e+00   y : 0.0000e+00
    i    : 0.0000e+00   z : 0.0000e+00
    omega: 0.0000e+00   vx: -0.0000e+00
    Omega: 0.0000e+00   vy: 2.9785e+04
    anom : 0.0000e+00   vz: 0.0000e+00

These values are generated from the internal storage for the elements ``orb._kep`` and ``orb._cart``. However, as a safeguard to avoid inconsistent pairs of elements, the object knows that a change has been made to the kepler elements and hence the cartesian ones are out of date. So if e.g. we would print the x-axis location with ``print(f'Orbit X-axis: {orb.x/pyorb.AU}')`` the elements would be automatically updated

.. code-block:: 

    Orbit X-axis: [1.5237]

To disable this automatic conversion use following flag 

.. code-block:: python

    orb.auto_update = False

Then if we update the eccentricity, print the cartesian coordinates, change a cartesian coordinate and print the kepler coordinates:

.. code-block:: python

    orb.e = 0.5
    print(f'Cartesian: {orb.cartesian}')
    orb.vz = 30e3
    print(f'Kepler: {orb.kepler}')
    print('Both:')
    print(orb)

We see that any conversion has to be made manually and the pair can be inconsistent:

.. code-block:: 

    Cartesian: [[ 2.27942276e+11]
     [ 0.00000000e+00]
     [ 0.00000000e+00]
     [-0.00000000e+00]
     [ 2.41295901e+04]
     [ 0.00000000e+00]]
    Kepler: [[2.27942276e+11]
     [5.00000000e-01]
     [0.00000000e+00]
     [0.00000000e+00]
     [0.00000000e+00]
     [0.00000000e+00]]
    Both:
    a    : 2.2794e+11   x : 2.2794e+11
    e    : 5.0000e-01   y : 0.0000e+00
    i    : 0.0000e+00   z : 0.0000e+00
    omega: 0.0000e+00   vx: -0.0000e+00
    Omega: 0.0000e+00   vy: 2.4130e+04
    anom : 0.0000e+00   vz: 3.0000e+04

At this point, one would have to choose which set of coordinates is the one desired and use that as a basis for transformation. E.g. if we chose to use the kepler as base:

.. code-block:: python

    orb.calculate_cartesian()
    print(f'Orbit X-axis: {orb.x/pyorb.AU}')

We get a x-position consistent with a 0.5 eccentricity orbit:

.. code-block:: 

    Orbit X-axis: [0.76185]

This manual transformation should mainly be used if there are performance issues or if total control over the transformation is needed. Here we can also see another property of the orbit class: it is completely vectorized. Hence why the ``orb.x`` returns a numpy 1-length vector. Since the ``Orbit`` object can have multiple orbits in the same instance there are a few convenience functions to work with multiple orbits such as

.. code-block:: python

    import numpy as np

    orb.auto_update = True
    orb.direct_update = True

    orb.allocate(10)
    orb.update(
        a=np.linspace(0.5,2,num=10)*pyorb.AU, 
        e=0, 
        i=0, 
        omega=0, 
        Omega=0, 
        anom=0,
    )
    orb.add(num=2, a=4*pyorb.AU)

    print(orb)
    print('Orbit semi major axis [AU]:')
    for i,o in enumerate(orb):
        print(f'Item {i}: a={o.a/pyorb.AU} AU, e={o.e}')

Will generate:

.. code-block:: 

    12 Orbits
    Orbit semi major axis [AU]:
    Item 0: a=[0.5] AU, e=[0.]
    Item 1: a=[0.66666667] AU, e=[0.]
    Item 2: a=[0.83333333] AU, e=[0.]
    Item 3: a=[1.] AU, e=[0.]
    Item 4: a=[1.16666667] AU, e=[0.]
    Item 5: a=[1.33333333] AU, e=[0.]
    Item 6: a=[1.5] AU, e=[0.]
    Item 7: a=[1.66666667] AU, e=[0.]
    Item 8: a=[1.83333333] AU, e=[0.]
    Item 9: a=[2.] AU, e=[0.]
    Item 10: a=[4.] AU, e=[nan]
    Item 11: a=[4.] AU, e=[nan]

So there are quite a few steps to unpack there. 

Lets start with the ``orb.allocate(10)``: this command allocates space in the internal arrays used to store data and sets everything to ``nan``. This method remove all previous data and replaces it with arrays to support exactly 10 items.

Then we have the ``orb.update`` method. This method allows for smart assigning of parameters to internal items by use of the ``inds`` keyword argument. Only the items selected by ``inds`` are updated, and if no value is given, all items are updated. The ``inds`` parameter should be able to index a numpy array so any type that achieves that goal can be used. E.g. integers, logical arrays, integer arrays/lists and array slices can be used. Here we gave no inds, so all orbits are assigned the float values, e.g. all orbits will have an eccentricity of 0. The semi major axis however was given as an array of length 10, so here the individual items will get assigned to each element of this list. 

As opposed to the allocate function the ``orb.add`` method combines an insert at the end of the internal arrays with a call to ``orb.update``. Here we add two more orbits and give them both an semi major axis of 4 AU.

When printing the results we can see that each of the first 0-9 orbits have the linearly increasing semi major axis from the ``update`` method while the last two have the 4 AU values and ``nan`` as eccentricity values. Also the printing function reverts to only printing the number of orbits when more than 1 is used as not to clutter print statements.

Lastly we need to cover the loop. Iterating trough a ``Orbit`` object is slightly different than iterating trough a normal array. To save overhead internally the ``orb`` object is not a collection of 12 ``Orbit`` classes but the internal arrays are extended in the appropriate extra dimension, e.g. the internal ``orb._cart`` becomes a 6-by-12 matrix. Hence, to have access to all internal methods and properties when iterating trough an orbit a copy of the orbit with size 1 is created at the start of each iteration. This is the object ``o`` in the above loop.

**IMPORTANT NOTE**: This means that modifying a orbit instance in a loop using the ``Orbit`` iterator does NOT change the original instance. Also, as all properties return copies of the internal arrays one cannot iterate over the ``orb.cartesian`` numpy array and change the ``orb`` instance. To modify the internal variables in an external loop use the pointer to the internal array ``orb._cart``. This will modify the ``orb`` instance, but without triggering auto-update or direct-update. So a call to ``orb.calculate_kepler()`` will probably have to be performed after the iteration.

As a standard, SI units are used. But, we can also create orbits with an arbitrary system of units. For this there is an convenience function ``pyorb.get_G`` that generates the gravitational constant in the requested units. Some units are implement as standard and can be called by using its string representation (e.g. 'AU'), otherwise one can just pass a float to describe the conversion between the SI unit and this unit (e.g. 3600.0 for 'h'). For example:

.. code-block:: python
    :linenos:

    G_ast = pyorb.get_G(length='AU', mass='Msol', time='y')
    print(f'SI gravitation constant: {pyorb.G} m^3 kg^-1 s^-2')
    print(f'Astronomical gravitation constant: {G_ast} AU^3 Msol^-1 y^-2')

    orb2 = pyorb.Orbit(M0 = 1.0, G=G_ast)
    orb2.update(a=1, e=0, i=0, omega=0, Omega=0, anom=0)

    print(orb2)

Will generate:

.. code-block:: 

    SI gravitation constant: 6.6743e-11 m^3 kg^-1 s^-2
    Astronomical gravitation constant: 39.47812018693255 AU^3 Msol^-1 y^-2
    a    : 1.0000e+00   x : 1.0000e+00
    e    : 0.0000e+00   y : 0.0000e+00
    i    : 0.0000e+00   z : 0.0000e+00
    omega: 0.0000e+00   vx: -0.0000e+00
    Omega: 0.0000e+00   vy: 6.2832e+00
    anom : 0.0000e+00   vz: 0.0000e+00

Here we see all the units now conforming to the new gravitational constant. If we print the period in these units they should be approximently 1 year. Another interesting fact about these units is that the orbital speed should be approximately 2pi as this is the circumference of a circle with radius 1 AU in units of AU:

.. code-block:: python

    print(f'Orbital period: {orb2.period} years')
    print(f'Orbital speed: {orb2.speed} AU/y')

.. code-block:: 

    Orbital period: [1.00000377] years
    Orbital speed: [6.28316164] AU/y


This concludes the getting started guide. If there are any additions that should be added here: please open up an issue in the github repository. Have fun orbiting!