#!/usr/bin/env python

'''Functions related to Keplerian orbital elements.

'''

#Python standard import


#Third party import
import numpy as np


#Local import



G = 6.6743e-11
'''float: Newtons gravitational constant [m^3 kg^-1 s^-2]

**Reference**: NIST - http://physics.nist.gov/cgi-bin/cuu/Value?bg
'''

AU = 1.495978707e11
'''float: Astronomical Unit [m]

The mean distance between the sun and the earth as defined in "International Astronomical Union, 31 August 2012"
'''


M_earth = 398600.5e9/G
'''float: Mass of the Earth using the WGS84 convention [kg]
'''

M_sol = 1.98847e30
"""float: The mass of the sun :math:`M_\\odot` given in kg

**Reference**: The Astronomical Almanac - http://asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf
"""

e_lim = 1e-9
"""float: The limit on eccentricity below witch an orbit is considered circular
"""

i_lim = np.pi*1e-9
"""float: The limit on inclination in radians below witch an orbit is considered not inclined
"""

def cart_to_equi(cart, mu=M_sol*G, degrees=False):
    '''Converts set of Cartesian state vectors to set of Equinoctial orbital elements.
    '''
    if not isinstance(cart, np.ndarray):
        raise TypeError('Input type {} not supported: must be {}'.format( type(cart),np.ndarray ))
    if cart.shape[0] != 6:
        raise ValueError('Input data must have at least 6 variables along axis 0: input shape is {}'.format(cart.shape))
    
    if len(cart.shape) < 2:
        input_is_vector = True
        try:
            cart.shape=(6,1)
        except ValueError as e:
            print('Error {} while trying to cast vector into single column.'.format(e))
            print('Input array shape: {}'.format(cart.shape))
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
    '''Converts set of Equinoctial orbital elements to set of Cartesian state vectors.
    '''
    if not isinstance(equi, np.ndarray):
        raise TypeError('Input type {} not supported: must be {}'.format( type(equi),np.ndarray ))
    if equi.shape[0] != 6:
        raise ValueError('Input data must have at least 6 variables along axis 0: input shape is {}'.format(equi.shape))
    
    if len(equi.shape) < 2:
        input_is_vector = True
        try:
            equi.shape=(6,1)
        except ValueError as e:
            print('Error {} while trying to cast vector into single column.'.format(e))
            print('Input array shape: {}'.format(equi.shape))
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

    lam = kep[3,...] + kep[4,...]
    om = kep[4,...]
    hi = 0.5*kep[2,...]
    if degrees:
        lam = np.radians(lam)
        om = np.radians(om)
        hi = np.radians(hi)

    elems = np.empty(kep.shape, dtype=kep.dtype)

    elems[0,...] = kep[0,...]
    elems[1,...] = kep[1,...]*np.sin(lam)
    elems[2,...] = kep[1,...]*np.cos(lam)
    
    elems[3,...] = np.sin(hi)*np.sin(om)
    elems[4,...] = np.sin(hi)*np.cos(om)

    elems[5,...] = kep[5,...] + lam

    return elems

def equi_to_kep(equi, degrees=False):

    kep = np.empty(equi.shape, dtype=equi.dtype)

    kep[0,...] = equi[0,...]
    kep[1,...] = np.sqrt(equi[1,...]**2 + equi[2,...]**2)
    hi = np.arcsin(np.sqrt(equi[3,...]**2 + equi[4,...]**2))
    kep[2,...] = 2*hi
    kep[4,...] = np.arctan2(equi[3,...], equi[4,...])

    lam = np.arctan2(equi[1,...], equi[2,...])
    if degrees:
        lam = np.degrees(lam)
        kep[4,...] = np.degrees(kep[4,...])
        kep[2,...] = np.degrees(kep[2,...])
    
    kep[3,...] = lam - kep[4,...]
    kep[5,...] = equi[5,...] - lam

    return kep

def cart_to_kep(cart, mu=M_sol*G, degrees=False):
    '''Converts set of Cartesian state vectors to set of Keplerian orbital elements.

    **Units:**
       Using default :code:`mu`, all variables are in `SI Units <https://www.nist.gov/pml/weights-and-measures/metric-si/si-units>`_

       If a standard gravitational parameter :code:`mu` is given in other units, all other input variables should also use the same unit system.

       Angles are by default given as radians, all angles are radians internally in functions, input and output angles can be both radians and degrees depending on the :code:`degrees` boolean.


    **Variables:**
       * :math:`a`: Semi-major axis
       * :math:`e`: Eccentricity
       * :math:`i`: Inclination
       * :math:`\\omega`: Argument of perihelion
       * :math:`\\Omega`: Longitude of ascending node
       * :math:`\\nu`: True anomaly
       * :math:`\\mu`: Standard gravitational parameter for the two body problem: :math:`G(M_1 + M_2)`.


    **Orientation of the ellipse in the coordinate system:**
       * For zero inclination :math:`i`: the ellipse is located in the x-y plane.
       * The direction of motion as True anoamly :math:`\\nu`: increases for a zero inclination :math:`i`: orbit is anti-coockwise, i.e. from +x towards +y.
       * If the eccentricity :math:`e`: is increased, the periapsis will lie in +x direction.
       * If the inclination :math:`i`: is increased, the ellipse will rotate around the x-axis, so that +y is rotated toward +z.
       * An increase in Longitude of ascending node :math:`\\Omega`: corresponds to a rotation around the z-axis so that +x is rotated toward +y.
       * Changing argument of perihelion :math:`\\omega`: will not change the plane of the orbit, it will rotate the orbit in the plane.
       * The periapsis is shifted in the direction of motion.
       * True anomaly measures from the +x axis, i.e :math:`\\nu = 0` is located at periapsis and :math:`\\nu = \\pi` at apoapsis.
       * All anomalies and orientation angles reach between 0 and :math:`2\\pi`

       *Reference:* "Orbital Motion" by A.E. Roy.
    

    **Constants:**
       * :mod:`~pyorb.kepler.e_lim`: Used to determine circular orbits
       * :mod:`~pyorb.kepler.i_lim`: Used to determine non-inclined orbits


    :param numpy.ndarray cart: Cartesian state vectors where rows 1-6 correspond to :math:`x`, :math:`y`, :math:`z`, :math:`v_x`, :math:`v_y`, :math:`v_z` and columns correspond to different objects.
    :param float/numpy.ndarray mu: Standard gravitational parameter of objects. If `mu` is a numpy vector, the element corresponding to each column of `cart` will be used for its element calculation, Default value is in SI units a massless object orbiting the Sun.
    :param bool degrees: If `true`, degrees are used. Else all angles are given in radians.

    :return: Keplerian orbital elements where rows 1-6 correspond to :math:`a`, :math:`e`, :math:`i`, :math:`\\omega`, :math:`\\Omega`, :math:`\\nu` and columns correspond to different objects.
    :rtype: numpy.ndarray

    **Example:**
       
       Convert an orbit around the sun from Kepler elements to cartesian elements:

       .. code-block:: python

            import pyorb
            import numpy as np

            orb = np.array([1*pyorb.AU, 0.5, 30, 10, 55, 2], dtype=np.float64)

            x = pyorb.kep_to_cart(orb, mu=pyorb.M_sol*pyorb.G, degrees=True)

            print('Position in the solar-system for an orbit of')
            for j,var in enumerate(pyorb.Orbit.KEPLER):
                print(f'{var}: {orb[j]}')
            print('is:')
            for j,var in enumerate(pyorb.Orbit.CARTESIAN):
                print(f'{var}: {x[j]}')


    *Reference:* Daniel Kastinen Master Thesis: Meteors and Celestial Dynamics and references therein.
    '''

    if not isinstance(cart,np.ndarray):
        raise TypeError('Input type {} not supported: must be {}'.format( type(cart),np.ndarray ))
    if cart.shape[0] != 6:
        raise ValueError('Input data must have at least 6 variables along axis 0: input shape is {}'.format(x.shape))
    
    if len(cart.shape) < 2:
        input_is_vector = True
        try:
            cart.shape=(6,1)
        except ValueError as e:
            print('Error {} while trying to cast vector into single column.'.format(e))
            print('Input array shape: {}'.format(cart.shape))
            raise
    else:
        input_is_vector = False

    ez = np.array([0,0,1], dtype=cart.dtype)
    ex = np.array([1,0,0], dtype=cart.dtype)
    
    o = np.empty(cart.shape, dtype=cart.dtype)
    iter_n = cart.shape[1]

    r = cart[:3,:]
    v = cart[3:,:]
    rn = np.linalg.norm(r,axis=0)
    vn = np.linalg.norm(v,axis=0)

    vr = np.sum( (r/rn)*v , axis=0)

    epsilon = vn**2*0.5 - mu/rn


    ### ECCENTRICITY ###
    e = 1.0/mu*((vn**2 - mu/rn)*r - (rn*vr)*v)
    o[1,:] = np.linalg.norm(e,axis=0)

    ### SEMI MAJOR AXIS ###
    #possible cases
    # e < 1
    #and
    # e >= 1
    e_hyp = o[1,:] >= 1
    o[0,:] = -mu/(2.0*epsilon)
    o[0,e_hyp] = -o[0,e_hyp]

    ### ANUGLAR MOMENTUM ###
    h = np.cross(r, v, axisa=0, axisb=0, axisc=0)
    hn = np.linalg.norm(h, axis=0)
    o[2,:] = np.arccos(h[2,:]/hn)

    #possible cases
    eg = o[1,:] > e_lim #e grater
    ig = o[2,:] > i_lim #i grater
    el = np.logical_not(eg) #e less equal
    il = np.logical_not(ig) #i less equal

    #e > elim & i > ilim 
    eg_ig = np.logical_and(eg, ig)

    #e > elim & i <= ilim
    eg_il = np.logical_and(eg, il)

    #e <= elim & i > ilim
    el_ig = np.logical_and(el, ig)

    #e <= elim & i <= ilim
    el_il = np.logical_and(el, il)

    ### ASCENDING NODE ###
    #ascending node pointing vector
    n = np.empty_like(h)
    nn = np.empty_like(hn)
    n[:,ig] = np.cross(ez, h[:,ig], axisa=0, axisb=0, axisc=0)
    nn[ig] = np.linalg.norm(n[:,ig], axis=0)

    #ensure [0,2pi]
    ny_neg = np.logical_and(n[1,:] < 0.0, ig)
    o[4,ig] = np.arccos(n[0,ig]/nn[ig])
    o[4,ny_neg] = 2.0*np.pi - o[4,ny_neg]

    #non inclined: no ascending node
    o[4,il] = 0

    ### ARGUMENT OF PERIAPSIS ###
    # circular orbits: no argument of periapsis
    o[3,el] = 0

    # elliptical and hyperbolic orbits
    #two cases
    cos_om = np.empty_like(hn)
    #first case: eg and ig (n-vector)
    #use vector angle between the two
    cos_om[eg_ig] = np.sum(n[:,eg_ig]*e[:,eg_ig], axis=0)/(nn[eg_ig]*o[1, eg_ig])

    #second case: eg and il (no n-vector)
    #use e vector angle
    cos_om[eg_il] = e[0,eg_il]/o[1, eg_il]
    
    #remove unused array positions
    cos_om = cos_om[eg]
    
    #do not fail due to number precision fluctuation
    cos_om[cos_om > 1.0] = 1.0
    cos_om[cos_om < -1.0] = -1.0

    o[3,eg] = np.arccos(cos_om)

    #first case: e and n vector angle
    ez_neg = np.logical_and(e[2,:] < 0.0, eg_ig)
    o[3,ez_neg] = 2.0*np.pi - o[3,ez_neg]

    #second case: ex component
    ey_neg = np.logical_and(e[1,:] < 0.0, eg_il)
    o[3,ey_neg] = 2.0*np.pi - o[3,ey_neg]

    ### TRUE ANOMALY ###
    cos_nu = np.empty_like(hn)

    #three cases
    #elliptical and hyperbolic: (angle from periapsis using e and r)
    cos_nu[eg] = np.sum(e[:,eg]*r[:,eg], axis=0)/(o[1, eg]*rn[eg])

    #circular and inclined: (angle from periapsis using n and r)
    #if e=0 and omega := 0, with inclination +y -> +z perihelion is ascending node
    cos_nu[el_ig] = np.sum((n[:,el_ig]/nn[el_ig])*(r[:,el_ig]/rn[el_ig]), axis=0)

    #circular and planar: (use angle of position vector)
    cos_nu[el_il] = r[0,el_il]/rn[el_il]

    #do not fail due to number precision fluctuation
    cos_nu[cos_nu > 1.0] = 1.0
    cos_nu[cos_nu < -1.0] = -1.0

    o[5,:] = np.arccos(cos_nu)

    #ensure [0,2pi]
    #elliptical and hyperbolic
    tmp_ind_ = np.logical_and(vr < 0.0, eg)
    o[5,tmp_ind_] = 2.0*np.pi - o[5,tmp_ind_]

    #circular and inclined
    tmp_ind_ = np.logical_and(r[2,:] < 0.0, el_ig)
    o[5,tmp_ind_] = 2.0*np.pi - o[5,tmp_ind_]

    #circular and planar
    tmp_ind_ = np.logical_and(r[1,:] < 0.0, el_il)
    o[5,tmp_ind_] = 2.0*np.pi - o[5,tmp_ind_]


    ## OUTPUT FORMATTING ##
    if degrees:
        o[2:,:] = np.degrees(o[2:,:])

    if input_is_vector:
        cart.shape = (6,)
        o.shape = (6,)

    return o


def true_to_eccentric(nu, e, degrees=False):
    '''Calculates the eccentric anomaly from the true anomaly.

    :param float/numpy.ndarray nu: True anomaly.
    :param float/numpy.ndarray e: Eccentricity of ellipse.
    :param bool degrees: If true degrees are used, else all angles are given in radians

    :return: Eccentric anomaly.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _nu = np.radians(nu)
    else:
        _nu = nu

    E = 2.0*np.arctan( np.sqrt( (1.0 - e)/(1.0 + e) )*np.tan(_nu*0.5) )
    E = np.mod(E + 2.*np.pi,2.*np.pi)
    
    if degrees:
        E = np.degrees(E)
    
    return E


def eccentric_to_true(E, e, degrees=False):
    '''Calculates the true anomaly from the eccentric anomaly.

    :param float/numpy.ndarray E: Eccentric anomaly.
    :param float/numpy.ndarray e: Eccentricity of ellipse.
    :param bool degrees: If true degrees are used, else all angles are given in radians

    :return: True anomaly.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _E = np.radians(E)
    else:
        _E = E

    nu = 2.0*np.arctan( np.sqrt( (1.0 + e)/(1.0 - e) )*np.tan(_E*0.5) )
    nu = np.mod(nu + 2.*np.pi, 2.*np.pi)

    if degrees:
        nu = np.degrees(nu)
    
    return nu


def eccentric_to_mean(E, e, degrees=False):
    '''Calculates the mean anomaly from the eccentric anomaly using Kepler equation.

    :param float/numpy.ndarray E: Eccentric anomaly.
    :param float/numpy.ndarray e: Eccentricity of ellipse.
    :param bool degrees: If true degrees are used, else all angles are given in radians

    :return: Mean anomaly.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _E = np.radians(E)
    else:
        _E = E

    M = _E - e*np.sin(_E)

    if degrees:
        M = np.degrees(M)
    return M


def true_to_mean(nu, e, degrees=False):
    '''Transforms true anomaly to mean anomaly.
    
    **Uses:**
       * :func:`~pyorb.kepler.true_to_eccentric`
       * :func:`~pyorb.kepler.eccentric_to_mean`

    :param float/numpy.ndarray nu: True anomaly.
    :param float/numpy.ndarray e: Eccentricity of ellipse.
    :param bool degrees: If true degrees are used, else all angles are given in radians
    
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


def elliptic_radius(E, a, e, degrees=False):
    '''Calculates the distance between the left focus point of an ellipse and a point on the ellipse defined by the eccentric anomaly.

    :param float/numpy.ndarray E: Eccentric anomaly.
    :param float/numpy.ndarray a: Semi-major axis of ellipse.
    :param float/numpy.ndarray e: Eccentricity of ellipse.
    :param bool degrees: If true degrees are used, else all angles are given in radians
    
    :return: Radius from left focus point.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _E = np.radians(E)
    else:
        _E = E

    return a*(1.0 - e*np.cos( _E ))


def rot_mat_x(theta, dtype=np.float64):
    '''Generates the 3D transformation matrix for rotation around X-axis.
    
    :param float theta: Angle to rotate.
    :param numpy.dtype dtype: The data-type of the output matrix.

    :return: Rotation matrix
    :rtype: (3,3) numpy.ndarray
    '''
    R = np.zeros((3,3), dtype=dtype)
    R[1,1] = np.cos(theta)
    R[1,2] = -np.sin(theta)
    R[2,1] = np.sin(theta)
    R[2,2] = np.cos(theta)
    R[0,0] = 1.0
    return R


def rot_mat_y(theta, dtype=np.float64):
    '''Generates the 3D transformation matrix for rotation around Y-axis.
    
    :param float theta: Angle to rotate.
    :param numpy.dtype dtype: The data-type of the output matrix.

    :return: Rotation matrix
    :rtype: (3,3) numpy.ndarray
    '''
    R = np.zeros((3,3), dtype=dtype)
    R[0,0] = np.cos(theta)
    R[0,2] = np.sin(theta)
    R[2,0] = -np.sin(theta)
    R[2,2] = np.cos(theta)
    R[1,1] = 1.0
    return R


def rot_mat_z(theta, dtype=np.float64):
    '''Generates the 3D transformation matrix for rotation around Z-axis.
    
    :param float theta: Angle to rotate.
    :param numpy.dtype dtype: The data-type of the output matrix.

    :return: Rotation matrix
    :rtype: (3,3) numpy.ndarray
    '''
    R = np.zeros((3,3), dtype=dtype)
    R[0,0] = np.cos(theta)
    R[0,1] = -np.sin(theta)
    R[1,0] = np.sin(theta)
    R[1,1] = np.cos(theta)
    R[2,2] = 1.0
    return R


def laguerre_solve_kepler(E0, M, e, tol=1e-12, max_iter=5000, degree=5):
    '''Solve the Kepler equation using the The Laguerre Algorithm, a algorithm that guarantees global convergence.
    Adjusted for solving only real roots (non-hyperbolic orbits)
    
    Absolute numerical tolerance is defined as :math:`|f(E)| < tol` where :math:`f(E) = M - E + e \\sin(E)`.

    # TODO: implement hyperbolic solving.

    *Note:* Choice of polynomial degree does not matter significantly for convergence rate.

    :param float M: Initial guess for eccentric anomaly.
    :param float M: Mean anomaly.
    :param float e: Eccentricity of ellipse.
    :param float tol: Absolute numerical tolerance eccentric anomaly.
    :param int max_iter: Maximum number of iterations before solver is aborted.
    :param int degree: Polynomial degree in derivation of Laguerre Algorithm.
    :return: Eccentric anomaly and number of iterations.
    :rtype: tuple of (float, int)

    *Reference:* Conway, B. A. (1986). An improved algorithm due to Laguerre for the solution of Kepler's equation. Celestial mechanics, 39(2), 199-211.

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

    _f = lambda E: M - E + e*np.sin(E)
    _fp = lambda E: e*np.cos(E) - 1.0
    _fpp = lambda E: -e*np.sin(E)

    E = E0

    f_eval = _f(E)

    it_num = 0

    while np.abs(f_eval) >= tol:
        it_num += 1

        fp_eval = _fp(E)

        sqrt_term = np.sqrt(np.abs(
            (degree - 1.0)**2*fp_eval**2
            - degree*(degree - 1.0)*f_eval*_fpp(E)
        ))

        denom_p = fp_eval + sqrt_term
        denom_m = fp_eval - sqrt_term

        if np.abs(denom_p) > np.abs(denom_m):
            delta = degree*f_eval/denom_p
        else:
            delta = degree*f_eval/denom_m

        E = E - delta

        f_eval = _f(E)

        if it_num > max_iter:
            break

    return E, it_num


def _get_kepler_guess(M, e):
    '''The different initial guesses for solving the Kepler equation based on input mean anomaly
    
    :param float M: Mean anomaly in radians.
    :param float e: Eccentricity of ellipse.
    :return: Guess for eccentric anomaly in radians.
    :rtype: float

    Reference: Esmaelzadeh, R., & Ghadiri, H. (2014). Appropriate starter for solving the Kepler's equation. International Journal of Computer Applications, 89(7).
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
    :param float/numpy.ndarray e: Eccentricity of ellipse.
    :return: Guess for eccentric anomaly in radians.
    :rtype: numpy.ndarray or float

    *Reference:* Esmaelzadeh, R., & Ghadiri, H. (2014). Appropriate starter for solving the Kepler's equation. International Journal of Computer Applications, 89(7).
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

            E_calc = _get_kepler_guess(Mc, ec)

            Ec[...] = E_calc

    else:
        E0 = _get_kepler_guess(M, e)

    return E0


def mean_to_eccentric(M, e, tol=1e-12, degrees=False):
    '''Calculates the eccentric anomaly from the mean anomaly by solving the Kepler equation.

    :param float/numpy.ndarray M: Mean anomaly.
    :param float/numpy.ndarray e: Eccentricity of ellipse.
    :param float tol: Numerical tolerance for solving Keplers equation in units of radians.
    :param bool degrees: If true degrees are used, else all angles are given in radians
    
    :return: True anomaly.
    :rtype: numpy.ndarray or float

    **Uses:**
       * :func:`~pyorb.kepler._get_kepler_guess`
       * :func:`~pyorb.kepler.laguerre_solve_kepler`
    '''

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
            raise TypeError(f'Input dimensions does not agree M:{_M.shape} != e:{e.shape}')


        E = np.empty(_M.shape, dtype=_M.dtype)

        out_it = E.size
        Mit = np.nditer(_M)
        eit = np.nditer(e)
        Eit = np.nditer(E, op_flags=['readwrite'])

        for it in range(out_it):
            Mc = next(Mit)
            ec = next(eit)
            Ec = next(Eit)

            E0 = _get_kepler_guess(Mc, ec)
            E_calc, it_num = laguerre_solve_kepler(E0, Mc, ec, tol=tol)

            Ec[...] = E_calc

    else:
        if e == 0:
            return M

        E0 = _get_kepler_guess(_M, e)
        E, it_num = laguerre_solve_kepler(E0, _M, e, tol=tol)


    if degrees:
        E = np.degrees(E)

    return E


def mean_to_true(M, e, tol=1e-12, degrees=False):
    '''Transforms mean anomaly to true anomaly.
    
    **Uses:**
       * :func:`~pyorb.kepler.mean_to_eccentric`
       * :func:`~pyorb.kepler.eccentric_to_true`

    :param float/numpy.ndarray M: Mean anomaly.
    :param float/numpy.ndarray e: Eccentricity of ellipse.
    :param float tol: Numerical tolerance for solving Keplers equation in units of radians.
    :param bool degrees: If true degrees are used, else all angles are given in radians
    
    :return: True anomaly.
    :rtype: numpy.ndarray or float
    '''
    if degrees:
        _M = np.radians(M)
    else:
        _M = M

    E = mean_to_eccentric(_M, e, tol=tol, degrees=False)
    nu = eccentric_to_true(E, e, degrees=False)

    if degrees:
        nu = np.degrees(nu)

    return nu


def orbital_speed(r, a, mu):
    '''Calculates the orbital speed at a given radius for an Keplerian orbit :math:`v = \\sqrt{\\mu \\left (\\frac{2}{r} - \\frac{1}{a} \\right )}`.
    
    :param float/numpy.ndarray r: Radius from the pericenter.
    :param float/numpy.ndarray a: Semi-major axis of ellipse.
    :param float/numpy.ndarray mu: Standard gravitation parameter :math:`\\mu = G(m_1 + m_2)` of the orbit.
    :return: Orbital speed.
    '''
    return np.sqrt(mu*(2.0/r - 1.0/a))


def orbital_period(a, mu):
    '''Calculates the orbital period of an Keplerian orbit based on the semi-major axis :math:`P = 2\\pi\\sqrt{\\frac{a^3}{\\mu}}`.
    
    :param float/numpy.ndarray a: Semi-major axis of ellipse.
    :param float/numpy.ndarray mu: Standard gravitation parameter :math:`\\mu = G(m_1 + m_2)` of the orbit.
    :return: Orbital period.
    '''
    return 2.0*np.pi*np.sqrt(a**3.0/mu)


def semi_major_axis(P, mu):
    '''Calculates the orbital semi-major axis of an Keplerian orbit based on the orbital period :math:`a = \\mu^{\\frac{1}{3}}(\\frac{P}{2\\pi})^{\\frac{2}{3}}`.
    
    :param float/numpy.ndarray P: Orbital period
    :param float/numpy.ndarray mu: Standard gravitation parameter :math:`\\mu = G(m_1 + m_2)` of the orbit.
    :return: semi-major axis.
    '''
    np.cbrt((P/(2.0*np.pi))**2*mu)
    return a


def kep_to_cart(kep, mu=M_sol*G, degrees=False):
    '''Converts set of Keplerian orbital elements to set of Cartesian state vectors.
    
    **Units:**
       Using default :code:`mu`, all variables are in `SI Units <https://www.nist.gov/pml/weights-and-measures/metric-si/si-units>`_

       If a standard gravitational parameter :code:`mu` is given in other units, all other input variables should also use the same unit system.

       Angles are by default given as radians, all angles are radians internally in functions, input and output angles can be both radians and degrees depending on the :code:`degrees` boolean.

    **Orientation of the ellipse in the coordinate system:**
       * For zero inclination :math:`i`: the ellipse is located in the x-y plane.
       * The direction of motion as True anoamly :math:`\\nu`: increases for a zero inclination :math:`i`: orbit is anti-coockwise, i.e. from +x towards +y.
       * If the eccentricity :math:`e`: is increased, the periapsis will lie in +x direction.
       * If the inclination :math:`i`: is increased, the ellipse will rotate around the x-axis, so that +y is rotated toward +z.
       * An increase in Longitude of ascending node :math:`\\Omega`: corresponds to a rotation around the z-axis so that +x is rotated toward +y.
       * Changing argument of perihelion :math:`\\omega`: will not change the plane of the orbit, it will rotate the orbit in the plane.
       * The periapsis is shifted in the direction of motion.

       *Reference:* "Orbital Motion" by A.E. Roy.

    **Variables:**
       * :math:`a`: Semi-major axis
       * :math:`e`: Eccentricity
       * :math:`i`: Inclination
       * :math:`\\omega`: Argument of perihelion
       * :math:`\\Omega`: Longitude of ascending node
       * :math:`\\nu`: True anoamly

    **Uses:**
       * :func:`~pyorb.kepler.true_to_eccentric`
       * :func:`~pyorb.kepler.elliptic_radius`

    :param numpy.ndarray kep: Keplerian orbital elements where rows 1-6 correspond to :math:`a`, :math:`e`, :math:`i`, :math:`\\omega`, :math:`\\Omega`, :math:`\\nu` and columns correspond to different objects.
    :param float/numpy.ndarray mu: Standard gravitational parameter of objects. If `mu` is a numpy vector, the element corresponding to each column of `cart` will be used for its element calculation, Default value is in SI units a massless object orbiting the Sun.
    :param bool degrees: If `true`, degrees are used. Else all angles are given in radians.

    :return: Cartesian state vectors where rows 1-6 correspond to :math:`x`, :math:`y`, :math:`z`, :math:`v_x`, :math:`v_y`, :math:`v_z` and columns correspond to different objects.
    :rtype: numpy.ndarray

    **Example:**
       
       Convert ??

       .. code-block:: python

          import pyorb


    *Reference:* Daniel Kastinen Master Thesis: Meteors and Celestial Dynamics, "Orbital Motion" by A.E. Roy.
    '''
    if not isinstance(kep, np.ndarray):
        raise TypeError('Input type {} not supported: must be {}'.format( type(kep),np.ndarray ))
    if kep.shape[0] != 6:
        raise ValueError('Input data must have at least 6 variables along axis 0: input shape is {}'.format(kep.shape))
    
    if len(kep.shape) < 2:
        input_is_vector = True
        try:
            kep.shape=(6,1)
        except ValueError as e:
            print('Error {} while trying to cast vector into single column.'.format(e))
            print('Input array shape: {}'.format(kep.shape))
            raise
    else:
        input_is_vector = False

    x = np.empty(kep.shape, dtype=kep.dtype)

    if degrees:
        kep[2:,:] = np.radians(kep[2:,:])

    nu = kep[5,:]
    omega = kep[3,:]
    asc_node = kep[4,:]
    inc = kep[2,:]
    wf = nu + omega
    e = kep[1,:]
    a = kep[0,:]

    Ecc = true_to_eccentric(nu,e,degrees=False)

    rn = elliptic_radius(Ecc,a,e,degrees=False)

    r = np.zeros( (3, kep.shape[1]), dtype=kep.dtype )
    r[0,:] = np.cos(wf)
    r[1,:] = np.sin(wf)
    r = rn*r

    cos_Omega = np.cos(asc_node)
    sin_Omega = np.sin(asc_node)
    cos_i = np.cos(inc)
    sin_i = np.sin(inc)
    cos_w = np.cos(omega)
    sin_w = np.sin(omega)

    #order is important not to change variables before they are used
    x_tmp = r[0,:].copy()
    r[2,:] = r[1,:]*sin_i
    r[0,:] = cos_Omega*r[0,:] - sin_Omega*r[1,:]*cos_i
    r[1,:] = sin_Omega*x_tmp  + cos_Omega*r[1,:]*cos_i

    l1 = cos_Omega*cos_w - sin_Omega*sin_w*cos_i
    l2 = -cos_Omega*sin_w - sin_Omega*cos_w*cos_i
    m1 = sin_Omega*cos_w + cos_Omega*sin_w*cos_i
    m2 = -sin_Omega*sin_w + cos_Omega*cos_w*cos_i
    n1 = sin_w*sin_i
    n2 = cos_w*sin_i
    b = a*np.sqrt(1.0 - e**2)
    n = np.sqrt(mu/a**3)
    nar = n*a/rn
    
    v = np.zeros( (3, kep.shape[1]), dtype=kep.dtype )
    bcos_E = b*np.cos(Ecc)
    asin_E = a*np.sin(Ecc)

    v[0,:] = l2*bcos_E - l1*asin_E
    v[1,:] = m2*bcos_E - m1*asin_E
    v[2,:] = n2*bcos_E - n1*asin_E

    v = nar*v

    x[:3,:] = r
    x[3:,:] = v

    if input_is_vector:
        x.shape = (6,)
        kep.shape = (6,)
        if degrees:
            kep[2:] = np.degrees(kep[2:])
    else:
        if degrees:
            kep[2:,:] = np.degrees(kep[2:,:])


    return x

