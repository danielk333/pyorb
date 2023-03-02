#!/usr/bin/env python

'''Functions related to Keplerian orbital elements.

'''

# Python standard import


# Third party import
import numpy as np


# Local import


G = 6.6743e-11
'''float: Newtons gravitational constant [m^3 kg^-1 s^-2]

**Reference**: NIST - http://physics.nist.gov/cgi-bin/cuu/Value?bg
'''

AU = 1.495978707e11
'''float: Astronomical Unit [m]

The mean distance between the sun and the earth as defined in 
"International Astronomical Union, 31 August 2012"
'''


M_earth = 398600.5e9/G
'''float: Mass of the Earth using the WGS84 convention [kg]
'''

M_sol = 1.98847e30
"""float: The mass of the sun :math:`M_\\odot` given in kg

**Reference**: The Astronomical Almanac - 
http://asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf
"""

e_lim = 1e-9
"""float: The limit on eccentricity below witch an orbit is considered circular
"""

i_lim = np.pi*1e-9
"""float: The limit on inclination in radians below witch an orbit is 
considered not inclined
"""


def cart_to_equi(cart, mu=M_sol*G, degrees=False):
    '''TODO

    Converts set of Cartesian state vectors to set of Equinoctial orbital 
    elements.
    '''
    if not isinstance(cart, np.ndarray):
        raise TypeError('Input type {} not supported: must be {}'.format(
            type(cart), np.ndarray))
    if cart.shape[0] != 6:
        raise ValueError(
            f'Input data must have at least 6 variables along axis 0: \
            input shape is {cart.shape}')

    if len(cart.shape) < 2:
        input_is_vector = True
        try:
            cart.shape = (6, 1)
        except ValueError as e:
            print(f'Error {e} while trying to cast vector into single column.')
            print(f'Input array shape: {cart.shape}')
            raise
    else:
        input_is_vector = False

    if degrees:
        pass

    if input_is_vector:
        cart.shape = (6,)
        # o.shape = (6,)

    raise NotImplementedError()


def equi_to_cart(equi, mu=M_sol*G, degrees=False):
    '''TODO

    Converts set of Equinoctial orbital elements to set of Cartesian state 
    vectors.
    '''
    if not isinstance(equi, np.ndarray):
        raise TypeError('Input type {} not supported: must be {}'.format(
            type(equi), np.ndarray))
    if equi.shape[0] != 6:
        raise ValueError(
            f'Input data must have at least 6 variables along axis 0: \
            input shape is {equi.shape}')

    if len(equi.shape) < 2:
        input_is_vector = True
        try:
            equi.shape = (6, 1)
        except ValueError as e:
            print(f'Error {e} while trying to cast vector into single column.')
            print(f'Input array shape: {equi.shape}')
            raise
    else:
        input_is_vector = False

    if degrees:
        pass

    if input_is_vector:
        equi.shape = (6,)
        # o.shape = (6,)

    raise NotImplementedError()


def kep_to_equi(kep, degrees=False):
    '''TODO

    a, h, k, p, q, lambda0,

    **Reference:**
        Broucke, R.A., Cefola, P.J., 1972. On the equinoctial orbit elements. 
        Celestial Mechanics 5, 303–310. https://doi.org/10.1007/BF01228432

    '''

    om_bar = kep[3, ...] + kep[4, ...]
    om = kep[4, ...]
    hi = 0.5*kep[2, ...]
    lambda0 = om_bar + kep[5, ...]
    if degrees:
        om_bar = np.radians(om_bar)
        om = np.radians(om)
        hi = np.radians(hi)

    elems = np.empty(kep.shape, dtype=kep.dtype)

    elems[0, ...] = kep[0, ...]
    elems[1, ...] = kep[1, ...]*np.sin(om_bar)
    elems[2, ...] = kep[1, ...]*np.cos(om_bar)

    elems[3, ...] = np.sin(hi)*np.sin(om)
    elems[4, ...] = np.sin(hi)*np.cos(om)

    elems[5, ...] = lambda0

    return elems


def equi_to_kep(equi, degrees=False):
    '''TODO
    '''

    kep = np.empty(equi.shape, dtype=equi.dtype)

    kep[0, ...] = equi[0, ...]
    kep[1, ...] = np.sqrt(equi[1, ...]**2 + equi[2, ...]**2)
    hi = np.arcsin(np.sqrt(equi[3, ...]**2 + equi[4, ...]**2))
    kep[2, ...] = 2*hi
    kep[4, ...] = np.arctan2(equi[3, ...], equi[4, ...])

    lam = np.arctan2(equi[1, ...], equi[2, ...])
    if degrees:
        lam = np.degrees(lam)
        kep[4, ...] = np.degrees(kep[4, ...])
        kep[2, ...] = np.degrees(kep[2, ...])

    kep[3, ...] = lam - kep[4, ...]
    kep[5, ...] = equi[5, ...] - lam

    return kep


def cart_to_kep(cart, mu=M_sol*G, degrees=False):
    '''Converts set of Cartesian state vectors to set of Keplerian orbital 
    elements.

    Parameters
    ----------
    kep : numpy.ndarray kep
        (6, N) or (6, ) array of Cartesian state vectors where rows correspond 
        to :math:`x`, :math:`y`, :math:`z`, :math:`v_x`, :math:`v_y`, 
        :math:`v_z` and columns correspond to different objects.
    mu : float or numpy.ndarray
        Float or (N, ) array with the standard gravitational parameter of 
        objects. If :code:`mu` is a numpy vector, the element corresponding to 
        each column of :code:`cart` will be used for its element calculation, 
        Default value is in SI units a massless object orbiting the Sun.
    degrees : bool
        If :code:`True`, use degrees. Else all angles are given in radians.

    Returns
    -------
    numpy.ndarray
        (6, N) or (6, ) array of Keplerian orbital elements where rows 
        correspond to :math:`a`, :math:`e`, :math:`i`, :math:`\\omega`, 
        :math:`\\Omega`, :math:`\\nu` and columns correspond to different 
        objects.
    
    Notes
    -----
    Transform
        Basic functionality is largley based on standard literature like [1]_. 
        Some open source material is avalible like these `Orbital-mechanics notes on GitHub <https://orbital-mechanics.space/intro.html>`_ .

    Arbitrary constants used
       * :mod:`~pyorb.kepler.e_lim`: Used to determine circular orbits
       * :mod:`~pyorb.kepler.i_lim`: Used to determine non-inclined orbits

    Variables
        * :math:`a`: Semi-major axis
        * :math:`e`: Eccentricity
        * :math:`i`: Inclination
        * :math:`\\omega`: Argument of perihelion
        * :math:`\\Omega`: Longitude of ascending node
        * :math:`\\nu`: True anoamly
        * :math:`\\mu`: Standard gravitational parameter for the two body problem: :math:`G(M_1 + M_2)`.
    
    Units
       Using default standard gravitational parameter :code:`mu` 
       (:math:`\\mu`), all variables are in `SI Units 
       <https://www.nist.gov/pml/weights-and-measures/metric-si/si-units>`_ 
       If a :code:`mu` is given in other units, all other input variables 
       should also use the same unit system. Angles are by default given as 
       radians, all angles are radians internally in functions, input and 
       output angles can be both radians and degrees depending on the 
       :code:`degrees` boolean.


    Orientation of the ellipse in the coordinate system [2]_
        For 0 inclination :math:`i`: the ellipse is located in the x-y plane.

        The direction of motion as True anoamly :math:`\\nu`: increases for a 
        zero inclination :math:`i`: orbit is anti-coockwise, i.e. from +x 
        towards +y.

        If the eccentricity :math:`e`: is increased, the periapsis will lie in 
        +x direction.

        If the inclination :math:`i`: is increased, the ellipse will rotate 
        around the x-axis, so that +y is rotated toward +z.

        An increase in Longitude of ascending node :math:`\\Omega`: 
        corresponds to a rotation around the z-axis so that +x is rotated 
        toward +y.

        Changing argument of perihelion :math:`\\omega`: will not change the 
        plane of the orbit, it will rotate the orbit in the plane.

        The periapsis is shifted in the direction of motion.

    .. [1] A.E. Roy. "Orbital Motion"
    .. [2] D.A. Vallado. "Fundamentals of Astrodynamics and Applications"

    Examples
    --------
    
    Example of using the base conversion function to transform from
    kepler elements to cartesian coordinates

    >>> import pyorb
    >>> import numpy as np
    >>> G = pyorb.get_G(length='AU', mass='Msol', time='y')
    >>> G
    39.47812018693255
    >>> cart = np.array([
    ...     0.70710678,  0.70710678, 0.,
    ...     -4.4428662, 4.4428662, 0.,
    ... ])
    >>> kep = pyorb.cart_to_kep(
    ...     cart, 
    ...     mu=1*G, 
    ...     degrees=True,
    ... )

    See Also
    --------
    kep_to_cart : Convert from cartesian to kepler

    '''

    if not isinstance(cart, np.ndarray):
        raise TypeError('Input type {} not supported: must be {}'.format(
            type(cart), np.ndarray))
    if cart.shape[0] != 6:
        raise ValueError(
            f'Input data must have at least 6 variables along axis 0:\
             input shape is {cart.shape}')

    if len(cart.shape) < 2:
        input_is_vector = True
        try:
            cart.shape = (6, 1)
        except ValueError as e:
            print(f'Error {e} while trying to cast vector into single column.')
            print(f'Input array shape: {cart.shape}')
            raise
    else:
        input_is_vector = False

    ez = np.array([0, 0, 1], dtype=cart.dtype)

    o = np.empty(cart.shape, dtype=cart.dtype)

    r = cart[:3, :]
    v = cart[3:, :]
    rn = np.linalg.norm(r, axis=0)
    vn = np.linalg.norm(v, axis=0)

    vr = np.sum((r/rn)*v, axis=0)

    epsilon = vn**2*0.5 - mu/rn

    # ## ECCENTRICITY ###
    e = 1.0/mu*((vn**2 - mu/rn)*r - (rn*vr)*v)
    o[1, :] = np.linalg.norm(e, axis=0)

    # ## SEMI MAJOR AXIS ###
    # possible cases
    # e < 1
    # and
    # e >= 1
    e_hyp = o[1, :] >= 1
    o[0, :] = -mu/(2.0*epsilon)
    o[0, e_hyp] = -o[0, e_hyp]

    # ## ANUGLAR MOMENTUM ###
    h = np.cross(r, v, axisa=0, axisb=0, axisc=0)
    hn = np.linalg.norm(h, axis=0)
    o[2, :] = np.arccos(h[2, :]/hn)

    # possible cases
    eg = o[1, :] > e_lim  # e grater
    # i grater (include retrograde planar orbits)
    ig = np.logical_and(o[2, :] > i_lim, o[2, :] < np.pi - i_lim)
    el = np.logical_not(eg)  # e less equal
    il = np.logical_not(ig)  # i less equal

    # e > elim & i > ilim
    eg_ig = np.logical_and(eg, ig)

    # e > elim & i <= ilim
    eg_il = np.logical_and(eg, il)

    # e <= elim & i > ilim
    el_ig = np.logical_and(el, ig)

    # e <= elim & i <= ilim
    el_il = np.logical_and(el, il)

    # ## ASCENDING NODE ###
    # ascending node pointing vector
    n = np.empty_like(h)
    nn = np.empty_like(hn)
    n[:, ig] = np.cross(ez, h[:, ig], axisa=0, axisb=0, axisc=0)
    nn[ig] = np.linalg.norm(n[:, ig], axis=0)

    # ensure [0,2pi]
    ny_neg = np.logical_and(n[1, :] < 0.0, ig)
    o[4, ig] = np.arccos(n[0, ig]/nn[ig])
    o[4, ny_neg] = 2.0*np.pi - o[4, ny_neg]

    # non inclined: no ascending node
    o[4, il] = 0

    # ## ARGUMENT OF PERIAPSIS ###
    # circular orbits: no argument of periapsis
    o[3, el] = 0

    # elliptical and hyperbolic orbits
    # two cases
    cos_om = np.empty_like(hn)
    # first case: eg and ig (n-vector)
    # use vector angle between the two
    cos_om[eg_ig] = np.sum(
        n[:, eg_ig]*e[:, eg_ig],
        axis=0,
    )/(nn[eg_ig]*o[1, eg_ig])

    # second case: eg and il (no n-vector)
    # use e vector angle
    cos_om[eg_il] = e[0, eg_il]/o[1, eg_il]

    # remove unused array positions
    cos_om = cos_om[eg]

    # do not fail due to number precision fluctuation
    cos_om[cos_om > 1.0] = 1.0
    cos_om[cos_om < -1.0] = -1.0

    o[3, eg] = np.arccos(cos_om)

    # first case: e and n vector angle
    ez_neg = np.logical_and(e[2, :] < 0.0, eg_ig)
    o[3, ez_neg] = 2.0*np.pi - o[3, ez_neg]

    # second case: ex component
    # prograde
    ey_neg = np.logical_and(o[2, :] < np.pi*0.5, eg_il)
    ey_neg2 = np.logical_and(ey_neg, e[1, :] < 0.0)
    o[3, ey_neg2] = 2.0*np.pi - o[3, ey_neg2]

    # retrograde
    ey_neg = np.logical_and(o[2, :] > np.pi*0.5, eg_il)
    ey_neg2 = np.logical_and(ey_neg, e[1, :] >= 0.0)
    o[3, ey_neg2] = 2.0*np.pi - o[3, ey_neg2]

    # ## TRUE ANOMALY ###
    cos_nu = np.empty_like(hn)

    # three cases
    # elliptical and hyperbolic: (angle from periapsis using e and r)
    cos_nu[eg] = np.sum(e[:, eg]*r[:, eg], axis=0)/(o[1, eg]*rn[eg])

    # circular and inclined: (angle from periapsis using n and r)
    # if e=0 and omega := 0, with inclination +y -> +z perihelion is ascending 
    # node
    cos_nu[el_ig] = np.sum(
        (n[:, el_ig]/nn[el_ig])*(r[:, el_ig]/rn[el_ig]),
        axis=0,
    )

    # circular and planar: (use angle of position vector)
    cos_nu[el_il] = r[0, el_il]/rn[el_il]

    # do not fail due to number precision fluctuation
    cos_nu[cos_nu > 1.0] = 1.0
    cos_nu[cos_nu < -1.0] = -1.0

    o[5, :] = np.arccos(cos_nu)

    # ensure [0,2pi]
    # elliptical and hyperbolic
    tmp_ind_ = np.logical_and(vr < 0.0, eg)
    o[5, tmp_ind_] = 2.0*np.pi - o[5, tmp_ind_]

    # circular and inclined
    tmp_ind_ = np.logical_and(r[2, :] < 0.0, el_ig)
    o[5, tmp_ind_] = 2.0*np.pi - o[5, tmp_ind_]

    # circular and planar
    # prograde
    tmp_ind_ = np.logical_and(o[2, :] < np.pi*0.5, el_il)
    tmp_ind2_ = np.logical_and(tmp_ind_, r[1, :] < 0.0)
    o[5, tmp_ind2_] = 2.0*np.pi - o[5, tmp_ind2_]

    # if retrograde, its reversed
    tmp_ind_ = np.logical_and(o[2, :] > np.pi*0.5, el_il)
    tmp_ind2_ = np.logical_and(tmp_ind_, r[1, :] >= 0.0)
    o[5, tmp_ind2_] = 2.0*np.pi - o[5, tmp_ind2_]

    # # OUTPUT FORMATTING ##
    if degrees:
        o[2:, :] = np.degrees(o[2:, :])

    if input_is_vector:
        cart.shape = (6,)
        o.shape = (6,)

    return o


def orbit_total_angular_momentum(a, e, mu):
    '''Calculates the total angular momentum from orbital parameters.

    :param float/numpy.ndarray a: Semi-major axis of (>0) ellipse or (<0) 
        hyperbola.
    :param float/numpy.ndarray e: Eccentricity of ellipse (e<1), parabola 
        (e==1) or hyperbola (e>1).
    :param float/numpy.ndarray mu: Standard gravitation parameter 
        :math:`\\mu = G(m_1 + m_2)` of the orbit.

    :return: Total angular momentum.
    '''
    return np.sqrt(mu*a*(1 - e**2))


def parabolic_total_angular_momentum(q, mu):
    '''Calculates the total angular momentum from orbital parameters.

    :param float/numpy.ndarray q: Periapsis-distance of parabola.
    :param float/numpy.ndarray mu: Standard gravitation parameter 
        :math:`\\mu = G(m_1 + m_2)` of the orbit.

    :return: Total angular momentum.
    '''
    return np.sqrt(mu*2*q)


def true_to_eccentric(nu, e, degrees=False):
    '''Calculates the eccentric anomaly from the true anomaly.

    :param float/numpy.ndarray nu: True anomaly.
    :param float/numpy.ndarray e: Eccentricity of ellipse (e<1), parabola 
        (e==1) or hyperbola (e>1).
    :param bool degrees: If true degrees are used, else all angles are given in 
        radians

    :return: Eccentric anomaly.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _nu = np.radians(nu)
    else:
        _nu = nu

    if isinstance(_nu, np.ndarray) or isinstance(e, np.ndarray):
        if not isinstance(_nu, np.ndarray):
            _nu = np.ones(e.shape, dtype=e.dtype)*_nu
        if not isinstance(e, np.ndarray):
            e = np.ones(_nu.shape, dtype=_nu.dtype)*e
        if _nu.shape != e.shape:
            raise TypeError('Input dimensions does not agree')

        E = np.empty(_nu.shape, dtype=_nu.dtype)
        hyp = e > 1
        per = e == 1
        eli = e < 1

        E[hyp] = 2.0*np.arctanh(
            np.sqrt((e[hyp] - 1.0)/(e[hyp] + 1.0))*np.tan(_nu[hyp]*0.5)
        )
        E[per] = np.tan(_nu[per]*0.5)
        E[eli] = 2.0*np.arctan(
            np.sqrt((1.0 - e[eli])/(1.0 + e[eli]))*np.tan(_nu[eli]*0.5)
        )
    else:
        if e > 1:
            E = 2.0*np.arctanh(np.sqrt((e - 1.0)/(e + 1.0))*np.tan(_nu*0.5))
        elif e == 1:
            E = np.tan(_nu*0.5)
        else:
            E = 2.0*np.arctan(np.sqrt((1.0 - e)/(1.0 + e))*np.tan(_nu*0.5))

    E = np.mod(E + 2.*np.pi, 2.*np.pi)

    if degrees:
        E = np.degrees(E)

    return E


def eccentric_to_true(E, e, degrees=False):
    '''Calculates the true anomaly from the eccentric anomaly.

    :param float/numpy.ndarray E: elliptic, parabolic or hyperbolic eccentric 
        anomaly.
    :param float/numpy.ndarray e: Eccentricity of ellipse (e<1), 
        parabola (e==1) or hyperbola (e>1).
    :param bool degrees: If true degrees are used, else all angles are given in 
        radians

    :return: True anomaly.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _E = np.radians(E)
    else:
        _E = E

    if isinstance(_E, np.ndarray) or isinstance(e, np.ndarray):
        if not isinstance(_E, np.ndarray):
            _E = np.ones(e.shape, dtype=e.dtype)*_E
        if not isinstance(e, np.ndarray):
            e = np.ones(_E.shape, dtype=_E.dtype)*e
        if _E.shape != e.shape:
            raise TypeError('Input dimensions does not agree')

        nu = np.empty(_E.shape, dtype=_E.dtype)
        hyp = e > 1
        per = e == 1
        eli = e < 1

        nu[hyp] = 2.0*np.arctan(
            np.sqrt((e[hyp] + 1.0)/(e[hyp] - 1.0))*np.tanh(_E[hyp]*0.5)
        )
        nu[per] = 2.0*np.arctan(_E[per])
        nu[eli] = 2.0*np.arctan(
            np.sqrt((1.0 + e[eli])/(1.0 - e[eli]))*np.tan(_E[eli]*0.5)
        )
    else:
        if e > 1:
            nu = 2.0*np.arctan(np.sqrt((e + 1.0)/(e - 1.0))*np.tanh(_E*0.5))
        elif e == 1:
            nu = 2.0*np.arctan(_E)
        else:
            nu = 2.0*np.arctan(np.sqrt((1.0 + e)/(1.0 - e))*np.tan(_E*0.5))

    nu = np.mod(nu + 2.*np.pi, 2.*np.pi)
    if degrees:
        nu = np.degrees(nu)

    return nu


def eccentric_to_mean(E, e, degrees=False):
    '''Calculates the mean anomaly from the (elliptic, parabolic or hyperbolic) 
    eccentric anomaly using Kepler equation.

    :param float/numpy.ndarray E: elliptic, parabolic or hyperbolic eccentric 
        anomaly.
    :param float/numpy.ndarray e: Eccentricity of ellipse (e<1), 
        parabola (e==1) or hyperbola (e>1).
    :param bool degrees: If true degrees are used, else all angles are given in 
        radians

    :return: Mean anomaly.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _E = np.radians(E)
    else:
        _E = E

    if isinstance(_E, np.ndarray) or isinstance(e, np.ndarray):
        if not isinstance(_E, np.ndarray):
            _E = np.ones(e.shape, dtype=e.dtype)*_E
        if not isinstance(e, np.ndarray):
            e = np.ones(_E.shape, dtype=_E.dtype)*e
        if _E.shape != e.shape:
            raise TypeError('Input dimensions does not agree')

        M = np.empty(_E.shape, dtype=_E.dtype)
        hyp = e > 1
        per = e == 1
        eli = e < 1

        M[hyp] = e[hyp]*np.sinh(_E[hyp]) - _E[hyp]
        M[per] = _E[per] + _E[per]**3/3.0
        M[eli] = _E[eli] - e[eli]*np.sin(_E[eli])
    else:
        if e > 1:
            M = e*np.sinh(_E) - _E
        elif e == 1:
            M = _E + _E**3/3.0
        else:
            M = _E - e*np.sin(_E)

    if degrees:
        M = np.degrees(M)
    return M


def true_to_mean(nu, e, degrees=False):
    '''Transforms true anomaly to mean anomaly.

    :param float/numpy.ndarray nu: True anomaly.
    :param float/numpy.ndarray e: Eccentricity of ellipse (e<1), 
        parabola (e==1) or hyperbola (e>1).
    :param bool degrees: If true degrees are used, else all angles are given in 
        radians

    :return: Mean anomaly.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _nu = np.radians(nu)
    else:
        _nu = nu

    E = true_to_eccentric(_nu, e, degrees=False)
    M = eccentric_to_mean(E, e, degrees=False)

    if degrees:
        M = np.degrees(M)
    return M


def orbital_distance(h, mu, e, nu, degrees=False):
    '''Calculates the distance between the left focus point of an 
    ellipse (e<1), parabola (e==1) or hyperbola (e>1) and a
    point on the orbit defined by the true anomaly.

    :param float/numpy.ndarray h: Orbit total angular momentum.
    :param float/numpy.ndarray mu: Standard gravitation parameter.
    :param float/numpy.ndarray e: Eccentricity of ellipse (e<1), 
        parabola (e==1) or hyperbola (e>1).
    :param float/numpy.ndarray nu: True anomaly.
    :param bool degrees: If true degrees are used, else all angles are given in 
        radians

    :return: Radius from left focus point.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _nu = np.radians(nu)
    else:
        _nu = nu

    return h**2/(mu*(1 + e*np.cos(_nu)))


def elliptic_radius(E, a, e, degrees=False):
    '''Calculates the distance between the left focus point of an ellipse and a 
    point on the ellipse defined by the eccentric anomaly.

    :param float/numpy.ndarray E: Eccentric anomaly.
    :param float/numpy.ndarray a: Semi-major axis of ellipse.
    :param float/numpy.ndarray e: Eccentricity of ellipse.
    :param bool degrees: If true degrees are used, else all angles are given in 
        radians

    :return: Radius from left focus point.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _E = np.radians(E)
    else:
        _E = E

    return a*(1.0 - e*np.cos(_E))


def parabolic_radius(nu, q, degrees=False):
    '''Calculates the distance between the left focus point of an parabola and 
    a point on the parabola defined by the eccentric anomaly.

    :param float/numpy.ndarray nu: True anomaly.
    :param float/numpy.ndarray q: Periapsis-distance of parabola.
    :param bool degrees: If true degrees are used, else all angles are given in 
        radians

    :return: Radius from left focus point.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _nu = np.radians(nu)
    else:
        _nu = nu

    return 2.0*q/(1.0 + np.cos(_nu))


def hyperbolic_radius(nu, a, e, degrees=False):
    '''Calculates the distance between the left focus point of an hyperbola and 
        a point on the hyperbola defined by the hyperbolic anomaly.

    :param float/numpy.ndarray nu: True anomaly.
    :param float/numpy.ndarray a: Semi-transverse axis of hyperbola (positive).
    :param float/numpy.ndarray e: Eccentricity of hyperbola.
    :param bool degrees: If true degrees are used, else all angles are given in 
        radians

    :return: Radius from left focus point.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _nu = np.radians(nu)
    else:
        _nu = nu

    return a*(e**2 - 1.0)/(1.0 + e*np.cos(_nu))


def rot_mat_x(theta, dtype=np.float64):
    '''Generates the 3D transformation matrix for rotation around X-axis.

    :param float theta: Angle to rotate.
    :param numpy.dtype dtype: The data-type of the output matrix.

    :return: Rotation matrix
    :rtype: (3,3) numpy.ndarray
    '''
    R = np.zeros((3, 3), dtype=dtype)
    R[1, 1] = np.cos(theta)
    R[1, 2] = -np.sin(theta)
    R[2, 1] = np.sin(theta)
    R[2, 2] = np.cos(theta)
    R[0, 0] = 1.0
    return R


def rot_mat_y(theta, dtype=np.float64):
    '''Generates the 3D transformation matrix for rotation around Y-axis.

    :param float theta: Angle to rotate.
    :param numpy.dtype dtype: The data-type of the output matrix.

    :return: Rotation matrix
    :rtype: (3,3) numpy.ndarray
    '''
    R = np.zeros((3, 3), dtype=dtype)
    R[0, 0] = np.cos(theta)
    R[0, 2] = np.sin(theta)
    R[2, 0] = -np.sin(theta)
    R[2, 2] = np.cos(theta)
    R[1, 1] = 1.0
    return R


def rot_mat_z(theta, dtype=np.float64):
    '''Generates the 3D transformation matrix for rotation around Z-axis.

    :param float theta: Angle to rotate.
    :param numpy.dtype dtype: The data-type of the output matrix.

    :return: Rotation matrix
    :rtype: (3,3) numpy.ndarray
    '''
    R = np.zeros((3, 3), dtype=dtype)
    R[0, 0] = np.cos(theta)
    R[0, 1] = -np.sin(theta)
    R[1, 0] = np.sin(theta)
    R[1, 1] = np.cos(theta)
    R[2, 2] = 1.0
    return R


def laguerre_solve_kepler(E0, M, e, tol=1e-12, max_iter=5000, degree=5):
    '''Solve the Kepler equation using the The Laguerre Algorithm, a algorithm 
    that guarantees global convergence.

    Absolute numerical tolerance is defined as :math:`|f(E)| < tol` where 
    :math:`f(E) = M - E + e \\sin(E)` or where 
    :math:`f(E) = M + E - e \\sinh(E)`.

    # TODO: implement in C and bind using ctypes

    *Note:* Choice of polynomial degree does not matter significantly for 
        convergence rate.

    :param float E0: Initial guess for eccentric anomaly.
    :param float M: Mean anomaly.
    :param float e: Eccentricity of ellipse or hyperbola.
    :param float tol: Absolute numerical tolerance eccentric anomaly.
    :param int max_iter: Maximum number of iterations before solver is aborted.
    :param int degree: Polynomial degree in derivation of Laguerre Algorithm.
    :return: Eccentric anomaly and number of iterations.
    :rtype: tuple of (float, int)

    *Reference:* Conway, B. A. (1986). An improved algorithm due to Laguerre 
        for the solution of Kepler's equation. 
        Celestial mechanics, 39(2), 199-211.

    **Example:**

    .. code-block:: python

       import pyorb.kepler as pykep
       M = 3.14
       e = 0.8

       #Use mean anomaly as initial guess
       E, iterations = pykep.laguerre_solve_kepler(
          E0 = M,
          M = M,
          e = e,
          tol = 1e-12,
       )
    '''

    degree = float(degree)

    if e > 1:
        def _f(E): 
            return M + E - e*np.sinh(E)

        def _fp(E): 
            return 1.0 - e*np.cosh(E)

        def _fpp(E): 
            return -e*np.sinh(E)
    else:
        def _f(E): 
            return M - E + e*np.sin(E)

        def _fp(E): 
            return e*np.cos(E) - 1.0

        def _fpp(E): 
            return -e*np.sin(E)

    E = E0

    f_eval = _f(E)

    it_num = 0

    while np.abs(f_eval) >= tol:
        it_num += 1

        fp_eval = _fp(E)

        sqrt_term = np.sqrt(np.abs(
            (degree - 1.0)**2*fp_eval**2 - degree*(degree - 1.0)*f_eval*_fpp(E)
        ))

        denom_p = fp_eval + sqrt_term
        denom_m = fp_eval - sqrt_term

        if np.abs(denom_p) > np.abs(denom_m):
            delta = degree*f_eval/denom_p
        else:
            delta = degree*f_eval/denom_m

        E = E - delta

        f_eval = _f(E)

        if it_num >= max_iter:
            break

    return E, it_num


def _get_hyperbolic_kepler_guess(M, e):
    '''The different initial guesses for solving the Kepler equation based on 
    input mean anomaly

    :param float M: Mean anomaly in radians.
    :param float e: Eccentricity of hyperbola.
    :return: Guess for eccentric anomaly in radians.
    :rtype: float

    Reference: T. M. Burkardt and J. M. A. Danby, “The solutions of Kepler’s 
        equation. II,” Celestial Mechanics, vol. 31, pp. 317–328, Nov. 1983, 
        doi: 10.1007/BF01844230.
    '''

    E0 = np.log(2*M/e + 1.8)

    return E0


def _get_kepler_guess(M, e):
    '''The different initial guesses for solving the Kepler equation based on 
    input mean anomaly

    :param float M: Mean anomaly in radians.
    :param float e: Eccentricity of ellipse.
    :return: Guess for eccentric anomaly in radians.
    :rtype: float

    Reference: Esmaelzadeh, R., & Ghadiri, H. (2014). Appropriate starter for 
        solving the Kepler's equation. 
        International Journal of Computer Applications, 89(7).
    '''
    if M > np.pi:
        _M = 2.0*np.pi - M
    else:
        _M = M

    if _M < 0.25:
        E0 = _M + e*np.sin(_M)/(1.0 - np.sin(_M + e) + np.sin(_M))
    elif _M < 2.0:
        E0 = _M + e
    else:
        E0 = _M + (e*(np.pi - _M))/(1.0 + e)

    if M > np.pi:
        E0 = 2.0*np.pi - E0

    return E0


def kepler_guess(M, e):
    '''Guess the initial iteration point for newtons method.

    :param float/numpy.ndarray M: Mean anomaly in radians.
    :param float/numpy.ndarray e: Eccentricity of ellipse (e<1), 
        parabola (e==1) or hyperbola (e>1).
    :return: Guess for eccentric anomaly in radians.
    :rtype: numpy.ndarray or float

    *Reference:* Esmaelzadeh, R., & Ghadiri, H. (2014). Appropriate starter for 
        solving the Kepler's equation. 
        International Journal of Computer Applications, 89(7).
    *Reference:* T. M. Burkardt and J. M. A. Danby, “The solutions of Kepler’s 
        equation. II,” Celestial Mechanics, vol. 31, pp. 317–328, Nov. 1983, 
        doi: 10.1007/BF01844230.
    '''

    if isinstance(M, np.ndarray) or isinstance(e, np.ndarray):
        if not isinstance(M, np.ndarray):
            M = np.ones(e.shape, dtype=e.dtype)*M
        if not isinstance(e, np.ndarray):
            e = np.ones(M.shape, dtype=M.dtype)*e
        if M.shape != e.shape:
            raise TypeError('Input dimensions does not agree')

        E0 = np.empty(M.shape, dtype=M.dtype)

        out_it = E0.size
        Mit = np.nditer(M)
        eit = np.nditer(e)
        Eit = np.nditer(E0, op_flags=['readwrite'])

        for it in range(out_it):
            Mc = next(Mit)
            ec = next(eit)
            Ec = next(Eit)

            if ec > 1:
                E_calc = _get_hyperbolic_kepler_guess(Mc, ec)
            else:
                E_calc = _get_kepler_guess(Mc, ec)

            Ec[...] = E_calc

    else:
        if e > 1:
            E0 = _get_hyperbolic_kepler_guess(M, e)
        else:
            E0 = _get_kepler_guess(M, e)

    return E0


def mean_to_eccentric(M, e, solver_options=None, degrees=False):
    '''Calculates the eccentric anomaly from the mean anomaly by solving the 
    Kepler equation.

    :param float/numpy.ndarray M: Mean anomaly.
    :param float/numpy.ndarray e: Eccentricity of ellipse (e<1), 
        parabola (e==1) or hyperbola (e>1).
    :param dict solver_options: Options for the numerical solution of Kepler's 
        equation. See `pyorb.kepler.laguerre_solve_kepler` for information.
    :param bool degrees: If true degrees are used, else all angles are given in
        radians

    :return: True anomaly.
    :rtype: numpy.ndarray or float

    **Note:**
        * For parabolic orbits the equation is analytic and solvable, ref: 
            *Montenbruck, Oliver; Pfleger, Thomas (2009). Astronomy on the 
            Personal Computer. Springer-Verlag Berlin Heidelberg. 
            ISBN 978-3-540-67221-0. p 64*
        * For hyperbolic orbits the hyperbolic sine is used in the Kepler 
            equation.

    **Uses:**
       * :func:`~pyorb.kepler._get_kepler_guess`
       * :func:`~pyorb.kepler.laguerre_solve_kepler`
    '''
    if solver_options is None:
        solver_options = {}

    if degrees:
        _M = np.radians(M)
    else:
        _M = M

    if isinstance(_M, np.ndarray) or isinstance(e, np.ndarray):
        if not isinstance(_M, np.ndarray):
            _M = np.ones(e.shape, dtype=e.dtype)*_M
        if not isinstance(e, np.ndarray):
            e = np.ones(_M.shape, dtype=_M.dtype)*e
        if _M.shape != e.shape:
            raise TypeError(f'Input dimensions does not agree \
                M:{_M.shape} != e:{e.shape}')

        E = np.empty(_M.shape, dtype=_M.dtype)

        out_it = E.size
        Mit = np.nditer(_M)
        eit = np.nditer(e)
        Eit = np.nditer(E, op_flags=['readwrite'])

        for it in range(out_it):
            Mc = next(Mit)
            ec = next(eit)
            Ec = next(Eit)

            if ec > 1:
                E0 = _get_hyperbolic_kepler_guess(Mc, ec)
                E_calc, it_num = laguerre_solve_kepler(
                    E0, Mc, ec, **solver_options)
            elif ec == 1:
                A = 3.0/2.0*Mc
                B = np.cbrt(A + np.sqrt(A**2 + 1))
                E_calc = B - 1.0/B
            else:
                E0 = _get_kepler_guess(Mc, ec)
                E_calc, it_num = laguerre_solve_kepler(
                    E0, Mc, ec, **solver_options)

            Ec[...] = E_calc

    else:
        if e == 0:
            return M

        if e > 1:
            E0 = _get_hyperbolic_kepler_guess(_M, e)
            E, it_num = laguerre_solve_kepler(E0, _M, e, **solver_options)
        elif e == 1:
            A = 3.0/2.0*_M
            B = np.cbrt(A + np.sqrt(A**2 + 1))
            E = B - 1.0/B
        else:
            E0 = _get_kepler_guess(_M, e)
            E, it_num = laguerre_solve_kepler(E0, _M, e, **solver_options)

    if degrees:
        E = np.degrees(E)

    return E


def mean_to_true(M, e, solver_options=None, degrees=False):
    '''Transforms mean anomaly to true anomaly.

    :param float/numpy.ndarray M: Mean anomaly.
    :param float/numpy.ndarray e: Eccentricity of ellipse (e<1), 
        parabola (e==1) or hyperbola (e>1).
    :param dict solver_options: Options for the numerical solution of Kepler's 
        equation. See `pyorb.kepler.laguerre_solve_kepler` for information.
    :param bool degrees: If true degrees are used, else all angles are given 
        in radians

    :return: True anomaly.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _M = np.radians(M)
    else:
        _M = M

    E = mean_to_eccentric(_M, e, solver_options=solver_options, degrees=False)
    nu = eccentric_to_true(E, e, degrees=False)

    if degrees:
        nu = np.degrees(nu)

    return nu


def orbital_speed(r, a, mu):
    '''Calculates the orbital speed at a given radius for an Keplerian orbit 
    :math:`v = \\sqrt{\\mu \\left (\\frac{2}{r} - \\frac{1}{a} \\right )}`.

    :param float/numpy.ndarray r: Radius from the pericenter.
    :param float/numpy.ndarray a: Semi-major axis of (>0) ellipse or (<0) 
        hyperbola.
    :param float/numpy.ndarray mu: Standard gravitation parameter 
        :math:`\\mu = G(m_1 + m_2)` of the orbit.
    :return: Orbital speed.
    '''
    return np.sqrt(mu*(2.0/r - 1.0/a))


def orbital_period(a, mu):
    '''Calculates the orbital period of an Keplerian orbit based on the 
    semi-major axis :math:`P = 2\\pi\\sqrt{\\frac{a^3}{\\mu}}`.

    :param float/numpy.ndarray a: Semi-major axis of ellipse.
    :param float/numpy.ndarray mu: Standard gravitation parameter 
        :math:`\\mu = G(m_1 + m_2)` of the orbit.
    :return: Orbital period.
    '''
    return 2.0*np.pi*np.sqrt(a**3.0/mu)


def semi_major_axis(P, mu):
    '''Calculates the orbital semi-major axis of an Keplerian orbit based on 
    the orbital period :math:`a = \\mu^{\\frac{1}{3}}
    (\\frac{P}{2\\pi})^{\\frac{2}{3}}`.

    :param float/numpy.ndarray P: Orbital period
    :param float/numpy.ndarray mu: Standard gravitation parameter 
        :math:`\\mu = G(m_1 + m_2)` of the orbit.
    :return: semi-major axis.
    '''
    a = np.cbrt((P/(2.0*np.pi))**2*mu)
    return a


def true_of_the_asymptote(e, degrees=False):
    '''Calculate the True anomaly of the hyperbolic asymptotes.

    :param float/numpy.ndarray e: Eccentricity of hyperbola (e>1).
    :return: True anomaly of the hyperbolic asymptote.
    '''
    theta_inf = np.arccos(-1.0/e)

    if degrees:
        theta_inf = np.degrees(theta_inf)

    return theta_inf


def stumpff0(x):
    '''Calculates the Stumpff function number 0 value.

    **Reference**: Fundamentals of Celestial Mechanics 
        (Second Edition) (Hardback) [J.M.A. Danby - 1992]

    :param numpy.ndarray x: Stumpff input variable.
    :return: Stumpff function value.
    '''
    c0 = np.empty_like(x)
    inds = x > 0
    c0[inds] = np.cos(np.sqrt(x[inds]))

    c0[x == 0] = 1

    inds = x < 0
    c0[inds] = np.cosh(np.sqrt(-x[inds]))
    return c0


def stumpff1(x):
    '''Calculates the Stumpff function number 1 value.

    **Reference**: Fundamentals of Celestial Mechanics 
        (Second Edition) (Hardback) [J.M.A. Danby - 1992]

    :param numpy.ndarray x: Stumpff input variable.
    :return: Stumpff function value.
    '''
    c1 = np.empty_like(x)
    inds = x > 0
    c1[inds] = np.sin(np.sqrt(x[inds]))/np.sqrt(x[inds])

    c1[x == 0] = 1

    inds = x < 0
    c1[inds] = np.sinh(np.sqrt(-x[inds]))/np.sqrt(-x[inds])
    return c1


def stumpff2(x):
    '''Calculates the Stumpff function number 2 value.

    **Reference**: Fundamentals of Celestial Mechanics 
        (Second Edition) (Hardback) [J.M.A. Danby - 1992]

    :param numpy.ndarray x: Stumpff input variable.
    :return: Stumpff function value.
    '''
    c2 = np.empty_like(x)
    inds = x > 0
    c2[inds] = (1 - np.cos(np.sqrt(x[inds])))/x[inds]

    c2[x == 0] = 0.5

    inds = x < 0
    c2[inds] = (1 - np.cosh(np.sqrt(-x[inds])))/x[inds]
    return c2


def stumpff3(x):
    '''Calculates the Stumpff function number 3 value.

    **Reference**: Fundamentals of Celestial Mechanics 
        (Second Edition) (Hardback) [J.M.A. Danby - 1992]

    :param numpy.ndarray x: Stumpff input variable.
    :return: Stumpff function value.
    '''
    c3 = np.empty_like(x)
    inds = x > 0
    xsq = np.sqrt(x[inds])
    c3[inds] = (xsq - np.sin(xsq))/(xsq*x[inds])

    c3[x == 0] = 1/6

    inds = x < 0
    xsq = np.sqrt(-x[inds])
    c3[inds] = (xsq - np.sinh(xsq))/(xsq*x[inds])
    return c3


def stumpff(x):
    '''Calculate the 0-3 stumpff functions.

    These Stumpff function are used in the universal variable 
    formulation of the kepler orbit. 

    The `n`'th' Stumpff function is defined as:
    :math:`c_{n}(x) = \\sum^{\\infty}_{j=0} \\frac{(-1)^j x^j}{(n + 2j)!}`.

    However, it can be found that the first Stumpff functions can be expressed 
    in terms of trigonometric functions,
    as this is implemented here and only `c0, c1, c2, c3` are used in the 
    universal variable formulation of a kepler orbit.
    '''
    return stumpff0(x), stumpff1(x), stumpff2(x), stumpff3(x)


def euler_rotation_matrix(inc, omega, Omega, degrees=False):
    '''Generate the rotation matrix for the intrinsic rotation sequence Z-X-Z 
    used by keplerian orbital elements of (i, Omega, omega). 

    If any of the input angles are vectors the matrix will expand along a the 
    remaining axis.

    Appendix I (p. 483) of: Roithmayr, Carlos M.; Hodges, Dewey H. (2016), 
    Dynamics: Theory and Application of Kane's Method (1st ed.), 
    Cambridge University Press
    '''
    if isinstance(inc, np.ndarray):
        R = np.empty((3, 3, inc.size), dtype=np.float64)
    elif isinstance(omega, np.ndarray):
        R = np.empty((3, 3, omega.size), dtype=np.float64)
    elif isinstance(Omega, np.ndarray):
        R = np.empty((3, 3, Omega.size), dtype=np.float64)
    else:
        R = np.empty((3, 3), dtype=np.float64)

    if degrees:
        _omega = np.radians(omega)
        _inc = np.radians(inc)
        _Omega = np.radians(Omega)
    else:
        _omega = omega
        _inc = inc
        _Omega = Omega

    c1 = np.cos(_Omega)
    s1 = np.sin(_Omega)

    c2 = np.cos(_inc)
    s2 = np.sin(_inc)

    c3 = np.cos(_omega)
    s3 = np.sin(_omega)

    # col 0
    R[0, 0, ...] = c1*c3 - c2*s1*s3
    R[1, 0, ...] = c3*s1 + c1*c2*s3
    R[2, 0, ...] = s2*s3

    # col 1
    R[0, 1, ...] = -c1*s3 - c2*c3*s1
    R[1, 1, ...] = c1*c2*c3 - s1*s3
    R[2, 1, ...] = c3*s2

    # col 2
    R[0, 2, ...] = s1*s2
    R[1, 2, ...] = -c1*s2
    R[2, 2, ...] = c2

    return R


def kep_to_cart(kep, mu=M_sol*G, degrees=False):
    '''Converts set of Keplerian orbital elements to set of Cartesian state 
    vectors.

    Parameters
    ----------
    kep : numpy.ndarray kep
        (6, N) or (6, ) array of Keplerian orbital elements where rows 
        correspond to :math:`a`, :math:`e`, :math:`i`, :math:`\\omega`, 
        :math:`\\Omega`, :math:`\\nu` and columns correspond to different 
        objects.
    mu : float or numpy.ndarray
        Float or (N, ) array with the standard gravitational parameter of 
        objects. If :code:`mu` is a numpy vector, the element corresponding to 
        each column of :code:`cart` will be used for its element calculation, 
        Default value is in SI units a massless object orbiting the Sun.
    degrees : bool
        If :code:`True`, use degrees. Else all angles are given in radians.

    Returns
    -------
    numpy.ndarray
        (6, N) or (6, ) array of cartesian state vector(s) where rows 
        correspond to :math:`x`, :math:`y`, :math:`z`, :math:`v_x`, 
        :math:`v_y`, :math:`v_z` and columns correspond to different objects.
    
    Notes
    -----
    Transform
        Basic functionality is largley based on standard literature like [1]_. 
        Some open source material is avalible like these `Orbital-mechanics notes on GitHub <https://orbital-mechanics.space/intro.html>`_ .

    Variables
        * :math:`a`: Semi-major axis
        * :math:`e`: Eccentricity
        * :math:`i`: Inclination
        * :math:`\\omega`: Argument of perihelion
        * :math:`\\Omega`: Longitude of ascending node
        * :math:`\\nu`: True anoamly
        * :math:`\\mu`: Standard gravitational parameter for the two body problem: :math:`G(M_1 + M_2)`.
    
    Units
       Using default standard gravitational parameter :code:`mu` 
       (:math:`\\mu`), all variables are in `SI Units 
       <https://www.nist.gov/pml/weights-and-measures/metric-si/si-units>`_ 
       If a :code:`mu` is given in other units, all other input variables 
       should also use the same unit system. Angles are by default given as 
       radians, all angles are radians internally in functions, input and 
       output angles can be both radians and degrees depending on the 
       :code:`degrees` boolean.

    Orientation of the ellipse in the coordinate system [2]_
        For 0 inclination :math:`i`: the ellipse is located in the x-y plane.

        The direction of motion as True anoamly :math:`\\nu`: increases for a 
        zero inclination :math:`i`: orbit is anti-coockwise, i.e. from +x 
        towards +y.

        If the eccentricity :math:`e`: is increased, the periapsis will lie in 
        +x direction.

        If the inclination :math:`i`: is increased, the ellipse will rotate 
        around the x-axis, so that +y is rotated toward +z.

        An increase in Longitude of ascending node :math:`\\Omega`: 
        corresponds to a rotation around the z-axis so that +x is rotated 
        toward +y.

        Changing argument of perihelion :math:`\\omega`: will not change the 
        plane of the orbit, it will rotate the orbit in the plane.

        The periapsis is shifted in the direction of motion.

    .. [1] A.E. Roy. "Orbital Motion"
    .. [2] D.A. Vallado. "Fundamentals of Astrodynamics and Applications"

    Examples
    --------
    
    Example of using the base conversion function to transform from
    kepler elements to cartesian coordinates

    >>> import pyorb
    >>> import numpy as np
    >>> G = pyorb.get_G(length='AU', mass='Msol', time='y')
    >>> G
    39.47812018693255
    >>> kep = np.array([1, 0, 0, 0, 0, 45.0])
    >>> cart = pyorb.kep_to_cart(
    ...     kep, 
    ...     mu=1*G, 
    ...     degrees=True,
    ... )
    >>> cart
    array([ 0.70710678,  0.70710678,  0.        , -4.4428662 ,  4.4428662 ,
            0.        ])

    See Also
    --------
    cart_to_kep : Convert from cartesian to kepler

    '''
    if not isinstance(kep, np.ndarray):
        raise TypeError('Input type {} not supported: must be {}'.format(
            type(kep), np.ndarray))
    if kep.shape[0] != 6:
        raise ValueError(f'Input data must have at least 6 variables \
            along axis 0: input shape is {kep.shape}')

    if len(kep.shape) < 2:
        input_is_vector = True
        try:
            kep.shape = (6, 1)
        except ValueError as e:
            print(f'Error {e} while trying to cast vector into single column.')
            print(f'Input array shape: {kep.shape}')
            raise
    else:
        input_is_vector = False
    
    if not isinstance(mu, np.ndarray):
        mu = np.full((kep.shape[1],), mu, dtype=np.float64)

    x = np.empty(kep.shape, dtype=kep.dtype)

    if degrees:
        kep[2:, :] = np.radians(kep[2:, :])

    nu = kep[5, :]
    omega = kep[3, :]
    asc_node = kep[4, :]
    inc = kep[2, :]
    e = kep[1, :]
    a = kep[0, :]

    per = e == 1
    nper = np.logical_not(per)

    # In the consistent equations hyperbolas have negative semi-major axis
    a[e > 1] *= -1

    h = np.empty_like(a)
    h[nper] = orbit_total_angular_momentum(a[nper], e[nper], mu[nper])
    h[per] = parabolic_total_angular_momentum(a[per], mu[per])

    rn = orbital_distance(h, mu, e, nu, degrees=False)

    R = euler_rotation_matrix(inc, omega, asc_node, degrees=False)
    if len(R.shape) < 3:
        R.shape = R.shape + (1, )

    rx = rn*np.cos(nu)
    ry = rn*np.sin(nu)
    r = np.zeros((3, kep.shape[1]), dtype=kep.dtype)
    r[0, :] = R[0, 0, ...]*rx + R[0, 1, ...]*ry
    r[1, :] = R[1, 0, ...]*rx + R[1, 1, ...]*ry
    r[2, :] = R[2, 0, ...]*rx + R[2, 1, ...]*ry

    vn = mu/h
    vx = -vn*np.sin(nu)
    vy = vn*(e + np.cos(nu))

    v = np.zeros((3, kep.shape[1]), dtype=kep.dtype)
    v[0, :] = R[0, 0, ...]*vx + R[0, 1, ...]*vy
    v[1, :] = R[1, 0, ...]*vx + R[1, 1, ...]*vy
    v[2, :] = R[2, 0, ...]*vx + R[2, 1, ...]*vy

    x[:3, :] = r
    x[3:, :] = v

    if input_is_vector:
        x.shape = (6,)
        kep.shape = (6,)

    if degrees:
        kep[2:, ...] = np.degrees(kep[2:, ...])

    return x
