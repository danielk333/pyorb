'''Module for keeping physical constants without external dependencies.

Summary
-------

The module itself ( `const.py` ) is constructed as a Borg class instance [1]_
(ensuring that every instance has shared state) masquerading as a python module [2]_.
As such, this module ( `_const_docs.py` ) only documents its implementation,
the actual module ( `const.py` ) contains the implementation.


Mathematical constants
~~~~~~~~~~~~~~~~~~~~~~
pi
    You shouldn't have to ask
e
    Euler's number, not elementary charge
_2pi
    How it is usually used
_4pi
    How it is usually used

Fundamental physics constants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
    light speed in vacuum [m/s]
G
    Newtons gravitational constant [m^3 kg^-1 s^-2] from [3]_

WGS-84 Earth model
~~~~~~~~~~~~~~~~~~
wgs84_a
    Semi-major axis [m] (by definition)
wgs84_invf
    inverse flattening (by definition)
wgs84_f
    from whence flattening arises
wgs84_omega
    Nominal Mean Angular Velocity [rad/s]
wgs84_GM
    Geocentric Gravitational Constant [m3/s2]
wgs84_M_earth
    Earth mass [kg]

Astronomical quantities
~~~~~~~~~~~~~~~~~~~~~~~
AU
    Astronomical unit [m]
    mean distance between the sun and the earth (exact IAU definition)
M_sol
    Mass of the sun [kg] from [4]_
GM_sol
    Standard gravitational parameter of the sun m^3/s^2

Symbolic indices
~~~~~~~~~~~~~~~~
K_a, K_e, K_i, K_om, K_OM, K_nu
    Index :math:`0, 1, 2, 3, 4, 5` in that order for
    Semi-major axis, eccentricity, inclination, argument of periapsis,
    longitude of ascending node and anomaly variables
E_a, E_h, E_k, E_p, E_q, E_lam
    Index :math:`0, 1, 2, 3, 4, 5` in that order for
    Semi-major axis, h, k, p, q and longitude variables


Notes
-----
Standard gravitational parameter
    The gravitational parameter of the bodies is defined as
    :math:`GM = G M` in m^3 s^{-2} but the product GM
    is known to much higher precision than either
    :math:`G` or :math:`M`.


.. [1] Python Cookbook ch05s23:
    https://www.safaribooksonline.com/library/view/Python+Cookbook/0596001673/ch05s23.html
.. [2] Python Cookbook ch05s16:
    https://www.safaribooksonline.com/library/view/python-cookbook/0596001673/ch05s16.html
.. [3] NIST:
    http://physics.nist.gov/cgi-bin/cuu/Value?bg
.. [4] The Astronomical Almanac:
    http://asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf


'''
