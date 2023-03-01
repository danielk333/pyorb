#!/usr/bin/env python

'''Orbit class definition.
'''

# Python standard import
import copy

# Third party import
import numpy as np


# Local import
from .kepler import G as G_SI
from . import kepler as functions


class Orbit:
    '''Main encapsulating class for the concept of an orbit.


    :Pointers:
        All properties return copies of the internal variables rather then 
        return a pointer to the internal variable, e.g. `self.cartesian` 
        accesses and creates a copy of the internal storage `self._cart`, 
        while the pointer can technically be used this may cause unexpected 
        behavior in variables as internal parameters can be set without passing
        trough the setter functions and previously accessed variables can be 
        changed by the instance after access and should be avoided.

    :Copy:
        The `Orbit` class has a `self.copy` method that returns a completely 
        independent copy of the instance. The `copy` function from the `copy` 
        standard python package can also be used to the same effect.

    :Iteration:
        The `Orbit` instance has an iteration method that iterates trough each 
        orbit in the instance and upon each iteration returns an `Orbit` 
        instance with only the member of the current iteration step.


    :param numpy.dtype dtype: A numpy dtype that will be used to store all 
        orbit information. Defaults to `np.float64`. To save RAM footprint and
        if memory usage is a problem, consider if `np.float32` can support your
        required precision.
    :param bool direct_update: Toggles if the corresponding kepler/cartesian 
        elements are directly updated when a cartesian/kepler element is 
        changed. Consider disabling if performance is an issue or multiple 
        properties should be set prior to cartesian or kepler calculation. 
        To completely disable automatic conversion between kepler and cartesian
        elements, `direct_update` also needs to be disabled.
    :param bool auto_update: Advanced configuration, toggles automatic updating
        of the kepler/cartesian elements when accessed if any cartesian/kepler 
        element has changed trough a property setter since last access. Does 
        not affect behavior if `direct_update` is `True`.
    :param bool kepler_read_only: Toggles if the kepler variables can be only 
        read or the properties can also be set, e.g. if `kepler_read_only=True`
        an access of `self.a` is fine but `self.a = 1` would raise an 
        `AttributeError`.
    :param bool cartesian_read_only: Same as `kepler_read_only` but for the 
    cartesian variables.
    :param float G: Gravitational constant, should match the set of units used.
        Defaults to SI units. See :func:`pyorb.get_G` for more information.
    :param epoch: Storage for the epoch of the orbit. Not used internally for 
        any calculations.
    :param bool degrees: Toggles if angles are given in degrees instead of 
        radians. Radians are used by default, i.e `degrees=False`.
    :param str type: Specifies the type of anomaly is used to describe the 
        position on the orbit. The string values for the available options are 
        listed in `Orbit.ANOMALY`. The setter is not case sensitive.
    :param numpy.ndarray m: Masses of the objects for which to calculate 
        orbits. Used to calculate the gravitational parameter.
    :param int num: Integer describing the number of orbits in this instance. 
        As the orbit calculations are vectorized, for performance when 
        converting large amounts of orbits a single instance should be used 
        with many orbits.
    :param dict solver_options: Settings for the numerical solving of Kepler's 
        equation. See description below.


    :Kepler's equation:
        The Kepler equation is solved using the The Laguerre Algorithm, 
        an algorithm that guarantees global convergence. This algorithm has 
        been split into two, one branch adjusted for solving only real roots 
        (non-hyperbolic orbits) and for solving only imaginary roots 
        (hyperbolic orbits). The options for this algorithm are currently:

        - _float_ **tol**: Numerical tolerance for solving Kepler's equation 
            in units of radians.
        - _int_ **max_iter**: Maximum number of iterations to reach error 
            tolerance.
        - _int_ **degree**: Polynomial degree used in derivation of Laguerre 
            Algorithm.

    '''

    CARTESIAN = ['x', 'y', 'z', 'vx', 'vy', 'vz']
    KEPLER = ['a', 'e', 'i', 'omega', 'Omega', 'anom']
    EQUINOCTIAL = ['a', 'h', 'k', 'p', 'q', 'l']
    ANOMALY = ['true', 'eccentric', 'mean']

    UPDATE_KW = ['cartesian', 'kepler', 'm'] + KEPLER + CARTESIAN

    def __init__(self, M0, **kwargs):
        self.dtype = kwargs.pop('dtype', np.float64)

        self.auto_update = kwargs.pop('auto_update', True)
        self.direct_update = kwargs.pop('direct_update', True)

        self.G = kwargs.pop('G', G_SI)

        self.solver_options = dict(
            tol=1e-12,
            max_iter=5000,
            degree=5,
        )
        self.solver_options.update(kwargs.pop('solver_options', {}))

        self.M0 = M0
        self.epoch = kwargs.pop('epoch', None)
        self.degrees = kwargs.pop('degrees', False)
        self.__type = kwargs.pop('type', 'true')
        if self.__type not in Orbit.ANOMALY:
            raise ValueError(f'Anomaly type "{self.__type}" not recognized')

        self.kepler_read_only = False
        self.cartesian_read_only = False

        self.allocate(kwargs.pop('num', 1))
        if 'm' in kwargs:
            self.m[:] = kwargs['m']
            del kwargs['m']
        if self.num > 0:
            self.update(**kwargs)

        self.kepler_read_only = kwargs.pop('kepler_read_only', False)
        self.cartesian_read_only = kwargs.pop('cartesian_read_only', False)

    def __str__(self):
        if self.num == 1:
            str_ = '\n'.join([
                f'{kkey:<5}: {kval:<.4e}   {ckey:<2}: {cval:<.4e}'
                for kval, kkey, cval, ckey in
                zip(
                    self._kep[:, 0], Orbit.KEPLER,
                    self._cart[:, 0], Orbit.CARTESIAN
                )
            ])
        else:
            str_ = f'{self.num} Orbits'
        return str_

    def _cart_check(self):
        if self.auto_update and not self.direct_update:
            self.calculate_cartesian(
                np.logical_not(self.__cart_calculated),
            )

    def _kep_check(self):
        if self.auto_update and not self.direct_update:
            self.calculate_kepler(
                np.logical_not(self.__kep_calculated),
            )

    def __copy__(self):
        return self[:]

    def copy(self):
        return self[:]

    def __getitem__(self, inds):
        m = self.m[inds].copy()
        try:
            num = len(m)
        except TypeError:
            num = 1

        tmp_orb = Orbit(
            self.M0,
            G = self.G,
            num = num,
            m = m,
            degrees = self.degrees,
            solver_options = copy.deepcopy(self.solver_options),
            dtype = self.dtype,
            epoch = self.epoch,
            type = self.type,
        )
        tmp_orb.auto_update = self.auto_update
        tmp_orb.direct_update = self.direct_update
        tmp_orb.kepler_read_only = self.kepler_read_only
        tmp_orb.cartesian_read_only = self.cartesian_read_only

        tmp_orb._cart = self._cart[:, inds].copy()
        tmp_orb._kep = self._kep[:, inds].copy()

        if len(tmp_orb._cart.shape) == 1:
            tmp_orb._cart.shape = (6, 1)
            tmp_orb._kep.shape = (6, 1)

        return tmp_orb

    def __len__(self):
        return self.num

    def __iter__(self):
        self.__it = 0
        return self

    def __next__(self):
        if self.__it < self.num:
            self.__it += 1
            return self[self.__it-1]
        else:
            raise StopIteration

    @property
    def type(self):
        return self.__type

    @type.setter
    def type(self, val):
        val = val.lower().strip()
        if val not in Orbit.ANOMALY:
            raise ValueError(f'Anomaly type "{val}" not recognized')
        if self.__type != val:
            if val == 'true':
                self.anom = self.true_anomaly
            elif val == 'eccentric':
                self.anom = self.eccentric_anomaly
            elif val == 'mean':
                self.anom = self.mean_anomaly
            self.__type = val

    def propagate(self, dt):
        '''Propagate all orbits in time by modifying the anomaly of the orbit, 
        i.e. by updating the kepler elements
        '''
        if self.degrees:
            wrap_ = 360.0
        else:
            wrap_ = 2*np.pi
        self.mean_anomaly = np.mod(
            self.mean_anomaly + self.mean_motion*dt, wrap_)
        if self.direct_update:
            self.calculate_cartesian()

    def delete(self, inds):
        '''Deletes the selected orbits from this instance.

        :param int/list/numpy.ndarray/slice inds: Incidences of orbs to delete.
        :return: None
        '''
        self.m = np.delete(self.m, inds, axis=0)
        self._cart = np.delete(self._cart, inds, axis=1)
        self._kep = np.delete(self._kep, inds, axis=1)
        self._true_anomaly = np.delete(self._true_anomaly, inds, axis=0)
        self._eccentric_anomaly = np.delete(
            self._eccentric_anomaly, inds, axis=0)
        self._mean_anomaly = np.delete(self._mean_anomaly, inds, axis=0)
        self.__cart_calculated = np.delete(
            self.__cart_calculated, inds, axis=0)
        self.__kep_calculated = np.delete(self.__kep_calculated, inds, axis=0)

    def add(self, num=1, **kwargs):
        '''Adds orbits to this instance

        **Note**: Cannot add objects with BOTH cartesian and Keplerian elements

        :param int num: Number of orbits to add, defaults to 1
        :return: None

        :Keyword arguments:
            - Any valid Cartesian or Keplerian variable name can be given as a 
                key to update that parameter. 

        '''
        inds = slice(self.num, self.num + num, None)

        self.m = np.append(
            self.m,
            np.full((num,), 0, dtype=self.dtype),
            axis=0,
        )
        self._cart = np.append(
            self._cart,
            np.full((6, num), np.nan, dtype=self.dtype),
            axis=1,
        )
        self._kep = np.append(
            self._kep,
            np.full((6, num), np.nan, dtype=self.dtype),
            axis=1,
        )
        self._true_anomaly = np.append(
            self._true_anomaly,
            np.full((num,), np.nan, dtype=self.dtype),
            axis=0,
        )
        self._eccentric_anomaly = np.append(
            self._eccentric_anomaly,
            np.full((num,), np.nan, dtype=self.dtype),
            axis=0,
        )
        self._mean_anomaly = np.append(
            self._mean_anomaly,
            np.full((num,), np.nan, dtype=self.dtype),
            axis=0,
        )
        self.__cart_calculated = np.append(
            self.__cart_calculated,
            np.full((num,), False, dtype=np.bool_),
            axis=0,
        )
        self.__kep_calculated = np.append(
            self.__kep_calculated,
            np.full((num,), False, dtype=np.bool_),
            axis=0,
        )

        if 'm' in kwargs:
            self.m[inds] = kwargs['m']
            del kwargs['m']

        self.update(inds = inds, **kwargs)

    def allocate(self, num):
        '''Changes the current allocation of orbits to a fixed number.

        **Warning**: This method overrides all currently stored data to 
            allocate new space.

        :param int num: Number of orbits to change allocation to
        :return: None
        '''
        self.m = np.full((num,), 0, dtype=self.dtype)
        self._cart = np.full((6, num), np.nan, dtype=self.dtype)
        self._kep = np.full((6, num), np.nan, dtype=self.dtype)
        self._true_anomaly = np.full((num,), np.nan, dtype=self.dtype)
        self._eccentric_anomaly = np.full((num,), np.nan, dtype=self.dtype)
        self._mean_anomaly = np.full((num,), np.nan, dtype=self.dtype)
        self.__cart_calculated = np.full((num,), False, dtype=np.bool_)
        self.__kep_calculated = np.full((num,), False, dtype=np.bool_)

    def update(self, inds=slice(None, None, None), **kwargs):
        '''Calculates Keplerian elements based on Cartesian elements

        **Note**: Cannot update BOTH cartesian and Keplerian elements.

        :param int/list/numpy.ndarray/slice inds: Incidences of orbits to 
            update, defaults to all
        :return: None

        :Keyword arguments:
            - Any valid Cartesian or Keplerian variable name can be given as a 
                key to update that parameter. 
            - If `kepler` is given, it is interpreted as an input array the 
                same as using `orb.kepler = kepler` and no other updates are 
                done.
            - If `cartesian` is given, it is interpreted as an input array the 
                same as using `orb.cartesian = cartesian` and no other updates
                are done.
            - Mass can also be updated this way but this assumes that if kepler
                elements are updated, they stay constant due to the mass change
                and the cartesian are changed. Vice versa with updating mass 
                and cartesian elements. If only mass is updated 
                :code:`calculate_cartesian` or :code:`calculate_kepler` must be 
                called manually depending on which are to remain constant.
        '''

        cart_updated = False
        kep_updated = False

        if 'cartesian' in kwargs:
            if self.cartesian_read_only:
                raise AttributeError('Cannot update read only \
                    Cartesian elements')
            self._cart[:, inds] = kwargs['cartesian']
            if self.direct_update:
                self.calculate_kepler()
            else:
                self.__cart_calculated[inds] = True
                self.__kep_calculated[inds] = False
            return

        if 'kepler' in kwargs:
            if self.kepler_read_only:
                raise AttributeError('Cannot update read only Kepler elements')
            self._kep[:, inds] = kwargs['kepler']
            if self.direct_update:
                self.calculate_cartesian()
            else:
                self.__kep_calculated[inds] = True
                self.__cart_calculated[inds] = False
            return

        for ind, key in enumerate(Orbit.CARTESIAN):
            if key in kwargs:
                if self.cartesian_read_only:
                    raise AttributeError('Cannot update read only \
                        Cartesian elements')
                cart_updated = True

        for ind, key in enumerate(Orbit.KEPLER):
            if key in kwargs:
                if self.kepler_read_only:
                    raise AttributeError('Cannot update read only \
                        Cartesian elements')
                if cart_updated:
                    raise ValueError('Cannot update both cartesian and \
                        Keplerian elements simultaneously.')
                kep_updated = True

        if 'm' in kwargs:
            if self.direct_update:
                raise ValueError('Cannot set "m" and direct update, \
                    set "m" via property and manually call calculate.')
            self.m[inds] = kwargs['m']

        if cart_updated:
            for ind, key in enumerate(Orbit.CARTESIAN):
                if key in kwargs:
                    self._cart[ind, inds] = kwargs[key]

            if self.direct_update:
                self.calculate_kepler()
            else:
                self.__kep_calculated[inds] = False
                self.__cart_calculated[inds] = True
        if kep_updated:
            for ind, key in enumerate(Orbit.KEPLER):
                if key in kwargs:
                    self._kep[ind, inds] = kwargs[key]

            if self.direct_update:
                self.calculate_cartesian()
            else:
                self.__cart_calculated[inds] = False
                self.__kep_calculated[inds] = True

    def calculate_cartesian(self, inds=slice(None, None, None)):
        '''Calculates Cartesian elements based on Keplerian elements

        :param int/list/numpy.ndarray/slice inds: Incidences of orbits to 
            calculate, defaults to all
        :return: None
        '''
        do_inds = np.full((self.num,), False, dtype=np.bool_)
        do_inds[inds] = True
        do_inds[inds] = np.logical_not(
            np.any(np.isnan(self._kep[:, inds]), axis=0))

        if np.sum(do_inds) == 0:
            return

        self.__cart_calculated[do_inds] = True
        if self.type == 'true':
            kep_tmp = self._kep[:, do_inds]
        else:
            kep_tmp = self._kep[:, do_inds]
            self.calc_true_anomaly(do_inds)
            kep_tmp[5, :] = self._true_anomaly[do_inds]

        self._cart[:, do_inds] = functions.kep_to_cart(
            kep_tmp,
            mu=self.G*(self.M0 + self.m[do_inds]),
            degrees=self.degrees,
        )

    def calculate_kepler(self, inds=slice(None, None, None)):
        '''Calculates Keplerian elements based on Cartesian elements

        :param int/list/numpy.ndarray/slice inds: Incidences of orbits to 
            calculate, defaults to all
        :return: None
        '''
        do_inds = np.full((self.num,), False, dtype=np.bool_)
        do_inds[inds] = True
        do_inds[inds] = np.logical_not(
            np.any(np.isnan(self._cart[:, inds]), axis=0))

        if np.sum(do_inds) == 0:
            return

        self.__kep_calculated[do_inds] = True
        self._kep[:, do_inds] = functions.cart_to_kep(
            self._cart[:, do_inds],
            mu = self.G*(self.M0 + self.m[do_inds]),
            degrees = self.degrees,
        )

        if self.type == 'eccentric':
            self._kep[5, do_inds] = functions.true_to_eccentric(
                self._kep[5, do_inds],
                self.e[do_inds],
                degrees = self.degrees,
            )
        elif self.type == 'mean':
            self._kep[5, do_inds] = functions.true_to_mean(
                self._kep[5, do_inds],
                self.e[do_inds],
                degrees = self.degrees,
            )

    @property
    def num(self):
        '''Number of orbits
        '''
        nk = self._kep.shape[1]
        nc = self._cart.shape[1]
        assert nk == nc, f'Kepler "{nk}" and cartesian "{nc}" \
            sizes do not agree'
        return nk

    @property
    def r(self):
        '''Position vector
        '''
        self._cart_check()
        return self._cart[:3, :].copy()

    @r.setter
    def r(self, value):
        if self.cartesian_read_only:
            raise AttributeError('Cannot update read only Cartesian elements')
        self._cart[:3, :] = value

        if self.direct_update:
            self.calculate_kepler()
        else:
            self.__kep_calculated[:] = False
            self.__cart_calculated[:] = True

    @property
    def v(self):
        '''Velocity vector
        '''
        self._cart_check()
        return self._cart[3:, :].copy()

    @v.setter
    def v(self, value):
        if self.cartesian_read_only:
            raise AttributeError('Cannot update read only Cartesian elements')
        self._cart[3:, :] = value

        if self.direct_update:
            self.calculate_kepler()
        else:
            self.__kep_calculated[:] = False
            self.__cart_calculated[:] = True

    @property
    def cartesian(self):
        '''Cartesian state vector
        '''
        self._cart_check()
        return self._cart.copy()

    @cartesian.setter
    def cartesian(self, value):
        if self.cartesian_read_only:
            raise AttributeError('Cannot update read only Cartesian elements')
        if isinstance(value, np.ndarray):
            if len(value.shape) == 1:
                self._cart[:, :] = value[:, None]
            else:
                self._cart[:, :] = value[:, :]
        else:
            self._cart[:, :] = value

        if self.direct_update:
            self.calculate_kepler()
        else:
            self.__kep_calculated[:] = False
            self.__cart_calculated[:] = True

    @property
    def x(self):
        '''X Position
        '''
        self._cart_check()
        return self._cart[0, :].copy()

    @x.setter
    def x(self, value):
        if self.cartesian_read_only:
            raise AttributeError('Cannot update read only Cartesian elements')
        self._cart[0, :] = value

        if self.direct_update:
            self.calculate_kepler()
        else:
            self.__kep_calculated[:] = False

    @property
    def y(self):
        '''Y Position
        '''
        self._cart_check()
        return self._cart[1, :].copy()

    @y.setter
    def y(self, value):
        if self.cartesian_read_only:
            raise AttributeError('Cannot update read only Cartesian elements')
        self._cart[1, :] = value

        if self.direct_update:
            self.calculate_kepler()
        else:
            self.__kep_calculated[:] = False

    @property
    def z(self):
        '''Z Position
        '''
        self._cart_check()
        return self._cart[2, :].copy()

    @z.setter
    def z(self, value):
        if self.cartesian_read_only:
            raise AttributeError('Cannot update read only Cartesian elements')
        self._cart[2, :] = value

        if self.direct_update:
            self.calculate_kepler()
        else:
            self.__kep_calculated[:] = False

    @property
    def vx(self):
        '''Velocity X-component
        '''
        self._cart_check()
        return self._cart[3, :].copy()

    @vx.setter
    def vx(self, value):
        if self.cartesian_read_only:
            raise AttributeError('Cannot update read only Cartesian elements')
        self._cart[3, :] = value

        if self.direct_update:
            self.calculate_kepler()
        else:
            self.__kep_calculated[:] = False

    @property
    def vy(self):
        '''Velocity Y-component
        '''
        self._cart_check()
        return self._cart[4, :].copy()

    @vy.setter
    def vy(self, value):
        if self.cartesian_read_only:
            raise AttributeError('Cannot update read only Cartesian elements')
        self._cart[4, :] = value

        if self.direct_update:
            self.calculate_kepler()
        else:
            self.__kep_calculated[:] = False

    @property
    def vz(self):
        '''Velocity Z-component
        '''
        self._cart_check()
        return self._cart[5, :].copy()

    @vz.setter
    def vz(self, value):
        if self.cartesian_read_only:
            raise AttributeError('Cannot update read only Cartesian elements')
        self._cart[5, :] = value

        if self.direct_update:
            self.calculate_kepler()
        else:
            self.__kep_calculated[:] = False

    @property
    def kepler(self):
        '''Keplerian state vector
        '''
        self._kep_check()
        return self._kep.copy()

    @kepler.setter
    def kepler(self, value):
        if self.kepler_read_only:
            raise AttributeError('Cannot update read only Kepler elements')
        if isinstance(value, np.ndarray):
            if len(value.shape) == 1:
                self._kep[:, :] = value[:, None]
            else:
                self._kep[:, :] = value[:, :]
        else:
            self._kep[:, :] = value

        if self.direct_update:
            self.calculate_cartesian()
        else:
            self.__cart_calculated[:] = False
            self.__kep_calculated[:] = True

    @property
    def a(self):
        '''Semi-major axis
        '''
        self._kep_check()
        return self._kep[0, :].copy()

    @a.setter
    def a(self, value):
        if self.kepler_read_only:
            raise AttributeError('Cannot update read only Kepler elements')
        self._kep[0, :] = value

        if self.direct_update:
            self.calculate_cartesian()
        else:
            self.__cart_calculated[:] = False

    @property
    def e(self):
        '''Eccentricity
        '''
        self._kep_check()
        return self._kep[1, :].copy()

    @e.setter
    def e(self, value):
        if self.kepler_read_only:
            raise AttributeError('Cannot update read only Kepler elements')
        self._kep[1, :] = value

        if self.direct_update:
            self.calculate_cartesian()
        else:
            self.__cart_calculated[:] = False

    @property
    def i(self):
        '''Inclination
        '''
        self._kep_check()
        return self._kep[2, :].copy()

    @i.setter
    def i(self, value):
        if self.kepler_read_only:
            raise AttributeError('Cannot update read only Kepler elements')
        self._kep[2, :] = value

        if self.direct_update:
            self.calculate_cartesian()
        else:
            self.__cart_calculated[:] = False

    @property
    def omega(self):
        '''Argument of perihelion
        '''
        self._kep_check()
        return self._kep[3, :].copy()

    @omega.setter
    def omega(self, value):
        if self.kepler_read_only:
            raise AttributeError('Cannot update read only Kepler elements')
        self._kep[3, :] = value

        if self.direct_update:
            self.calculate_cartesian()
        else:
            self.__cart_calculated[:] = False

    @property
    def Omega(self):
        '''Longitude of the ascending node
        '''
        self._kep_check()
        return self._kep[4, :].copy()

    @Omega.setter
    def Omega(self, value):
        if self.kepler_read_only:
            raise AttributeError('Cannot update read only Kepler elements')
        self._kep[4, :] = value

        if self.direct_update:
            self.calculate_cartesian()
        else:
            self.__cart_calculated[:] = False

    @property
    def anom(self):
        '''The orbital anomaly, depending on :code:`self.type` 
        it is either the True, Eccentric or Mean anomaly
        '''
        self._kep_check()
        return self._kep[5, :].copy()

    @anom.setter
    def anom(self, value):
        if self.kepler_read_only:
            raise AttributeError('Cannot update read only Kepler elements')
        self._kep[5, :] = value

        if self.direct_update:
            self.calculate_cartesian()
        else:
            self.__cart_calculated[:] = False

    @property
    def mean_motion(self):
        '''Mean motion
        '''
        if self.degrees:
            norm_ = 360.0
        else:
            norm_ = np.pi*2.0

        return norm_/self.period

    @mean_motion.setter
    def mean_motion(self, value):
        if self.kepler_read_only:
            raise AttributeError('Cannot update read only Kepler elements')
        if self.degrees:
            norm_ = 360.0
        else:
            norm_ = np.pi*2.0
        self.period = norm_/value

        if self.direct_update:
            self.calculate_cartesian()
        else:
            self.__cart_calculated[:] = False

    def calc_true_anomaly(self, inds=slice(None, None, None)):
        '''Calculates the true anomaly and stores it in 
        :code:`self._true_anomaly`.
        '''
        if self.type == 'eccentric':
            self._true_anomaly[inds] = functions.eccentric_to_true(
                self.anom[inds],
                self.e[inds],
                degrees = self.degrees,
            )
        elif self.type == 'mean':
            self._true_anomaly[inds] = functions.mean_to_true(
                self.anom[inds],
                self.e[inds],
                solver_options = self.solver_options,
                degrees = self.degrees,
            )
        else:
            self._true_anomaly[inds] = self.anom[inds]

    def calc_mean_anomaly(self, inds=slice(None, None, None)):
        '''Calculates the mean anomaly and stores it in 
        :code:`self._mean_anomaly`.
        '''
        if self.type == 'eccentric':
            self._mean_anomaly[inds] = functions.eccentric_to_mean(
                self.anom[inds],
                self.e[inds],
                degrees = self.degrees,
            )
        elif self.type == 'true':
            self._mean_anomaly[inds] = functions.true_to_mean(
                self.anom[inds],
                self.e[inds],
                degrees = self.degrees,
            )
        else:
            self._mean_anomaly[inds] = self.anom[inds]

    def calc_eccentric_anomaly(self, inds=slice(None, None, None)):
        '''Calculates the eccentric anomaly and stores it in 
        :code:`self._eccentric_anomaly`.
        '''
        if self.type == 'mean':
            self._eccentric_anomaly = functions.mean_to_eccentric(
                self.anom[inds],
                self.e[inds],
                solver_options = self.solver_options,
                degrees = self.degrees,
            )
        elif self.type == 'true':
            self._eccentric_anomaly[inds] = functions.true_to_eccentric(
                self.anom[inds],
                self.e[inds],
                degrees = self.degrees,
            )
        else:
            self._eccentric_anomaly[inds] = self.anom[inds]

    @property
    def eccentric_anomaly(self):
        '''Eccentric Anomaly
        '''
        if self.type == 'eccentric':
            return self.anom

        # save it as a local variable to allow speedy get in case someone
        # wants to repeatedly access the value without recalculation
        #
        # But then its on the users head if something changes since this does
        # not get updated until someone calls calc_eccentric_anomaly
        self.calc_eccentric_anomaly()
        return self._eccentric_anomaly.copy()

    @eccentric_anomaly.setter
    def eccentric_anomaly(self, value):
        if self.kepler_read_only:
            raise AttributeError('Cannot update read only Kepler elements')
        if self.type == 'eccentric':
            self.anom = value
        elif self.type == 'true':
            self.anom = functions.eccentric_to_true(
                value,
                self.e,
                degrees = self.degrees,
            )
        elif self.type == 'mean':
            self.anom = functions.eccentric_to_mean(
                value,
                self.e,
                degrees = self.degrees,
            )

        if self.direct_update:
            self.calculate_cartesian()
        else:
            self.__cart_calculated[:] = False

    @property
    def true_anomaly(self):
        '''true Anomaly
        '''
        if self.type == 'true':
            return self.anom

        self.calc_true_anomaly()
        return self._true_anomaly.copy()

    @true_anomaly.setter
    def true_anomaly(self, value):
        if self.kepler_read_only:
            raise AttributeError('Cannot update read only Kepler elements')
        if self.type == 'true':
            self.anom = value
        elif self.type == 'eccentric':
            self.anom = functions.true_to_eccentric(
                value,
                self.e,
                degrees = self.degrees,
            )
        elif self.type == 'mean':
            self.anom = functions.true_to_mean(
                value,
                self.e,
                degrees = self.degrees,
            )

        if self.direct_update:
            self.calculate_cartesian()
        else:
            self.__cart_calculated[:] = False

    @property
    def mean_anomaly(self):
        '''Mean Anomaly
        '''
        if self.type == 'mean':
            return self.anom

        self.calc_mean_anomaly()
        return self._mean_anomaly.copy()

    @mean_anomaly.setter
    def mean_anomaly(self, value):
        if self.kepler_read_only:
            raise AttributeError('Cannot update read only Kepler elements')
        if self.type == 'mean':
            self.anom = value
        elif self.type == 'true':
            self.anom = functions.mean_to_true(
                value,
                self.e,
                solver_options = self.solver_options,
                degrees = self.degrees,
            )
        elif self.type == 'eccentric':
            self.anom = functions.mean_to_eccentric(
                value,
                self.e,
                solver_options = self.solver_options,
                degrees = self.degrees,
            )

        if self.direct_update:
            self.calculate_cartesian()
        else:
            self.__cart_calculated[:] = False

    @property
    def apoapsis(self):
        return (1.0 + self.e)*self.a

    @property
    def periapsis(self):
        return (1.0 - self.e)*self.a

    @property
    def mean_longitude(self):
        return self.mean_anomaly + self.omega + self.Omega

    @property
    def equinoctial(self):
        kep = self.kepler
        kep[5, ...] = self.mean_anomaly
        elems = functions.kep_to_equi(kep, degrees=self.degrees)
        if self.num == 1:
            elems.shape = (6,)
        return elems

    @equinoctial.setter
    def equinoctial(self, vals):
        elems = functions.equi_to_kep(vals, degrees=self.degrees)
        tmp_elms = elems[:5, ...]
        if len(tmp_elms.shape) == 1:
            tmp_elms.shape = (5, 1)
        self._kep[:5, ...] = tmp_elms
        self.mean_anomaly = elems[5, ...]

    @property
    def period(self):
        '''Orbital period
        '''
        self._kep_check()
        return functions.orbital_period(self.a, self.G*(self.M0 + self.m))

    @period.setter
    def period(self):
        if self.kepler_read_only:
            raise AttributeError('Cannot update read only Kepler elements')
        self.a = functions.semi_major_axis(
            self.period,
            self.G*(self.M0 + self.m),
        )
        if self.direct_update:
            self.calculate_cartesian()
        else:
            self.__cart_calculated[:] = False

    @property
    def velocity(self):
        '''Orbital velocity (from cartesian)
        '''
        self._cart_check()
        return np.linalg.norm(self.v, axis=0)

    @velocity.setter
    def velocity(self, value):
        if self.cartesian_read_only:
            raise AttributeError('Cannot update read only Cartesian elements')
        self.v *= value/np.linalg.norm(self.v, axis=0)

        if self.direct_update:
            self.calculate_kepler()
        else:
            self.__kep_calculated[:] = False

    @property
    def speed(self):
        '''Orbital speed (from kepler)
        '''
        self._kep_check()
        return functions.orbital_speed(
            functions.elliptic_radius(
                self.eccentric_anomaly,
                self.a,
                self.e,
                degrees=self.degrees,
            ),
            self.a,
            self.G*(self.M0 + self.m),
        )
