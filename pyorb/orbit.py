#!/usr/bin/env python

'''Orbit class definition

'''

#Python standard import


#Third party import
import numpy as np


#Local import
from .kepler import G as G_SI
from . import kepler as functions



class Orbit:

    CARTESIAN = ['x', 'y', 'z', 'vx', 'vy', 'vz']
    KEPLER = ['a', 'e', 'i', 'omega', 'Omega', 'anom']
    ANOMALY = ['true', 'eccentric', 'mean']
    
    def __init__(self, M0, **kwargs):
        self.dtype = kwargs.pop('dtype', np.float64)

        self.G = kwargs.pop('G', G_SI)
        self.tol = kwargs.pop('tol', 1e-12)
        self.M0 = M0
        self.epoch = kwargs.pop('epoch', None)
        self.degrees = kwargs.pop('degrees', False)
        self.type = kwargs.pop('type', 'true')
        if self.type not in Orbit.ANOMALY:
            raise ValueError(f'Anomaly type "{self.type}" not recognized')

        self.allocate(kwargs.pop('num',1))
        if 'm' in kwargs:
            self.m[:] = kwargs['m']
            del kwargs['m']
        if self.num > 0:
            self.update(**kwargs)


    def __str__(self):
        if self.num == 1:
            str_  = '\n'.join([
                f'{kkey:<5}: {kval:<.4e}   {ckey:<2}: {cval:<.4e}' 
                for kval, kkey, cval, ckey in 
                zip(self._kep[:,0], Orbit.KEPLER, self._cart[:,0], Orbit.CARTESIAN)
            ])
        else:
            str_ = f'{self.num} Orbits'
        return str_

    def _cart_check(self):
        self.calculate_cartesian(
            np.logical_not(self.__cart_calculated),
        )

    def _kep_check(self):
        self.calculate_kepler(
            np.logical_not(self.__kep_calculated),
        )


    def __getitem__(self, inds):
        tmp_orb = Orbit(
            self.M0,
            G = self.G,
            m = self.m[inds],
            degrees = self.degrees,
            tol = self.tol, 
            dtype = self.dtype,
            epoch = self.epoch, 
            type = self.type,
        )
        tmp_orb._cart = self._cart[:,inds]
        tmp_orb._kep = self._kep[:,inds]

        if len(tmp_orb._cart.shape) == 1:
            tmp_orb._cart.shape = (6,1)
            tmp_orb._kep.shape = (6,1)
            
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


    def propagate(self, dt):
        '''Propagate all orbits in time by modifying the anomaly of the orbit, i.e. by updating the kepler elements
        '''
        if self.degrees:
            wrap_ = 360.0
        else:
            wrap_ = 2*np.pi
        self.mean_anomaly = np.mod(self.mean_anomaly + self.mean_motion*dt, wrap_)


    def add(self, num=1, **kwargs):
        '''Adds orbits to this instance

        **Note**: Cannot add objects with BOTH cartesian and Keplerian elements.

        :param int num: Number of orbits to add, defaults to 1
        :return: None

        :Keyword arguments:
            - Any valid Cartesian or Keplerian variable name can be given as a key to update that parameter. 

        '''
        inds = slice(self.num, self.num + num, None)

        self.m = np.append(
            self.m, 
            np.full((num,), 0, dtype=self.dtype), 
            axis=0,
        )
        self._cart = np.append(
            self._cart, 
            np.full((6,num), np.nan, dtype=self.dtype),
            axis=1,
        )
        self._kep = np.append(
            self._kep,
            np.full((6,num), np.nan, dtype=self.dtype),
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
            np.full((num,), False, dtype=np.bool),
            axis=0,
        )
        self.__kep_calculated = np.append(
            self.__kep_calculated,
            np.full((num,), False, dtype=np.bool),
            axis=0,
        )

        if 'm' in kwargs:
            self.m[inds] = kwargs['m']
            del kwargs['m']

        self.update(inds = inds, **kwargs)


    def allocate(self, num):
        '''Changes the current allocation of orbits to a fixed number.

        **Warning**: This method overrides all currently stored data to allocate new space.

        :param int num: Number of orbits to change allocation to
        :return: None
        '''
        self.m = np.full((num,), 0, dtype=self.dtype)
        self._cart = np.full((6,num), np.nan, dtype=self.dtype)
        self._kep = np.full((6,num), np.nan, dtype=self.dtype)
        self._true_anomaly = np.full((num,), np.nan, dtype=self.dtype)
        self._eccentric_anomaly = np.full((num,), np.nan, dtype=self.dtype)
        self._mean_anomaly = np.full((num,), np.nan, dtype=self.dtype)
        self.__cart_calculated = np.full((num,), False, dtype=np.bool)
        self.__kep_calculated = np.full((num,), False, dtype=np.bool)


    def update(self, inds=slice(None, None, None), **kwargs):
        '''Calculates Keplerian elements based on Cartesian elements

        **Note**: Cannot update BOTH cartesian and Keplerian elements.

        :param int/list/numpy.ndarray/slice inds: Incidences of orbits to update, defaults to all
        :return: None

        :Keyword arguments:
            - Any valid Cartesian or Keplerian variable name can be given as a key to update that parameter. 
            - Mass can also be updated this way but this requires manual execution of either :code:`calculate_cartesian` or :code:`calculate_kepler` depending on which are to remain constant.
        '''

        cart_updated = False
        kep_updated = False

        if 'cartesian' in kwargs:
            self._cart[:, inds] = kwargs['cartesian']
            self.__cart_calculated[inds] = True
            self.__kep_calculated[inds] = False
            return

        if 'kepler' in kwargs:
            self._kep[:, inds] = kwargs['kepler']
            self.__kep_calculated[inds] = True
            self.__cart_calculated[inds] = False
            return

        for ind, key in enumerate(Orbit.CARTESIAN):
            if key in kwargs:
                self._cart[ind, inds] = kwargs[key]
                cart_updated = True

        for ind, key in enumerate(Orbit.KEPLER):
            if key in kwargs:
                if cart_updated:
                    raise ValueError('Cannot update both cartesian and Keplerian elements simultaneously.')
                self._kep[ind, inds] = kwargs[key]
                kep_updated = True

        if 'm' in kwargs:
            self.m[inds] = kwargs['m']
            self.__kep_calculated[inds] = False
            self.__cart_calculated[inds] = False
        else:
            if cart_updated:
                self.__kep_calculated[inds] = False
                self.__cart_calculated[inds] = True
            if kep_updated:
                self.__cart_calculated[inds] = False
                self.__kep_calculated[inds] = True


    def calculate_cartesian(self, inds=slice(None,None,None)):
        '''Calculates Cartesian elements based on Keplerian elements

        :param int/list/numpy.ndarray/slice inds: Incidences of orbits to calculate, defaults to all
        :return: None
        '''
        do_inds = np.full((self.num,), False, dtype=np.bool)
        do_inds[inds] = True
        do_inds[inds] = np.logical_not(np.any(np.isnan(self._kep[:,inds]), axis=0))

        if np.sum(do_inds) == 0:
            return

        self.__cart_calculated[do_inds] = True
        if self.type == 'true':
            kep_tmp = self._kep[:,do_inds]
        else:
            kep_tmp = self._kep[:,do_inds]
            self.calc_true_anomaly(do_inds)
            kep_tmp[5,:] = self._true_anomaly[:,do_inds]

        self._cart[:,do_inds] = functions.kep_to_cart(
            kep_tmp,
            mu=self.G*(self.M0 + self.m[do_inds]),
            degrees=self.degrees,
        )


    def calculate_kepler(self, inds=slice(None,None,None)):
        '''Calculates Keplerian elements based on Cartesian elements

        :param int/list/numpy.ndarray/slice inds: Incidences of orbits to calculate, defaults to all
        :return: None
        '''
        do_inds = np.full((self.num,), False, dtype=np.bool)
        do_inds[inds] = True
        do_inds[inds] = np.logical_not(np.any(np.isnan(self._cart[:,inds]), axis=0))

        if np.sum(do_inds) == 0:
            return

        self.__kep_calculated[do_inds] = True
        self._kep[:,do_inds] = functions.cart_to_kep(
            self._cart[:,do_inds],
            mu = self.G*(self.M0 + self.m[do_inds]),
            degrees = self.degrees,
        )

        if self.type == 'eccentric':
            self._kep[5,do_inds] = functions.true_to_eccentric(
                self._kep[5,do_inds], 
                self.e[do_inds],
                degrees = self.degrees, 
            )
        elif self.type == 'mean':
            self._kep[5,do_inds] = functions.true_to_mean(
                self._kep[5,do_inds], 
                self.e[do_inds], 
                degrees = self.degrees, 
            )


    @property
    def num(self):
        '''Number of orbits
        '''
        nk = self._kep.shape[1]
        nc = self._cart.shape[1]
        assert nk == nc, f'Kepler "{nk}" and cartesian "{nc}" sizes do not agree'
        return nk


    @property
    def r(self):
        '''Position vector
        '''
        self._cart_check()
        return self._cart[:3,:]
    @r.setter
    def r(self, value):
        self.__kep_calculated[:] = False
        self.__cart_calculated[:] = True
        self._cart[:3,:] = value

    @property
    def v(self):
        '''Velocity vector
        '''
        self._cart_check()
        return self._cart[3:,:]
    @v.setter
    def v(self, value):
        self.__kep_calculated[:] = False
        self.__cart_calculated[:] = True
        self._cart[3:,:] = value


    @property
    def cartesian(self):
        '''Cartesian state vector
        '''
        self._cart_check()
        return self._cart.copy()
    @cartesian.setter
    def cartesian(self, value):
        self.__kep_calculated[:] = False
        self.__cart_calculated[:] = True
        self._cart[:,:] = value


    @property
    def x(self):
        '''X Position
        '''
        self._cart_check()
        return self._cart[0,:]
    @x.setter
    def x(self, value):
        self.__kep_calculated[:] = False
        self._cart[0,:] = value


    @property
    def y(self):
        '''Y Position
        '''
        self._cart_check()
        return self._cart[1,:]
    @y.setter
    def y(self, value):
        self.__kep_calculated[:] = False
        self._cart[1,:] = value


    @property
    def z(self):
        '''Z Position
        '''
        self._cart_check()
        return self._cart[2,:]
    @z.setter
    def z(self, value):
        self.__kep_calculated[:] = False
        self._cart[2,:] = value


    @property
    def vx(self):
        '''Velocity X-component
        '''
        self._cart_check()
        return self._cart[3,:]
    @vx.setter
    def vx(self, value):
        self.__kep_calculated[:] = False
        self._cart[3,:] = value


    @property
    def vy(self):
        '''Velocity Y-component
        '''
        self._cart_check()
        return self._cart[4,:]
    @vy.setter
    def vy(self, value):
        self.__kep_calculated[:] = False
        self._cart[4,:] = value


    @property
    def vz(self):
        '''Velocity Z-component
        '''
        self._cart_check()
        return self._cart[5,:]
    @vz.setter
    def vz(self, value):
        self.__kep_calculated[:] = False
        self._cart[5,:] = value


    @property
    def kepler(self):
        '''Keplerian state vector
        '''
        self._kep_check()
        return self._kep.copy()
    @kepler.setter
    def kepler(self, value):
        self.__cart_calculated[:] = False
        self.__kep_calculated[:] = True
        self._kep[:,:] = value


    @property
    def a(self):
        '''Semi-major axis
        '''
        self._kep_check()
        return self._kep[0,:]
    @a.setter
    def a(self, value):
        self.__cart_calculated[:] = False
        self._kep[0,:] = value


    @property
    def e(self):
        '''Eccentricity
        '''
        self._kep_check()
        return self._kep[1,:]
    @e.setter
    def e(self, value):
        self.__cart_calculated[:] = False
        self._kep[1,:] = value


    @property
    def i(self):
        '''Inclination
        '''
        self._kep_check()
        return self._kep[2,:]
    @i.setter
    def i(self, value):
        self.__cart_calculated[:] = False
        self._kep[2,:] = value


    @property
    def omega(self):
        '''Argument of perihelion
        '''
        self._kep_check()
        return self._kep[3,:]
    @omega.setter
    def omega(self, value):
        self.__cart_calculated[:] = False
        self._kep[3,:] = value


    @property
    def Omega(self):
        '''Longitude of the ascending node
        '''
        self._kep_check()
        return self._kep[4,:]
    @Omega.setter
    def Omega(self, value):
        self.__cart_calculated[:] = False
        self._kep[4,:] = value


    @property
    def anom(self):
        '''The orbital anomaly, depending on :code:`self.type` it is either the True, Eccentric or Mean anomaly
        '''
        self._kep_check()
        return self._kep[5,:]
    @anom.setter
    def anom(self, value):
        self.__cart_calculated[:] = False
        self._kep[5,:] = value


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
        self.__cart_calculated[:] = False
        if self.degrees:
            norm_ = 360.0
        else:
            norm_ = np.pi*2.0
        self.period = norm_/value



    def calc_true_anomaly(self, inds=slice(None,None,None)):
        '''Calculates the true anomaly and stores it in :code:`self._true_anomaly`.
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
                tol = self.tol,
                degrees = self.degrees,
            )
        else:
            self._true_anomaly[inds] = self.anom[inds]


    def calc_mean_anomaly(self, inds=slice(None,None,None)):
        '''Calculates the mean anomaly and stores it in :code:`self._mean_anomaly`.
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



    def calc_eccentric_anomaly(self, inds=slice(None,None,None)):
        '''Calculates the eccentric anomaly and stores it in :code:`self._eccentric_anomaly`.
        '''
        if self.type == 'mean':
            self._eccentric_anomaly = functions.mean_to_eccentric(
                self.anom[inds], 
                self.e[inds],
                tol = self.tol,
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

        #save it as a local variable to allow speedy get in case someone 
        #wants to repeatedly access the value without recalculation
        #
        #But then its on the users head if something changes since this does 
        #not get updated until someone calls calc_eccentric_anomaly
        self.calc_eccentric_anomaly()
        return self._eccentric_anomaly

    @eccentric_anomaly.setter
    def eccentric_anomaly(self, value):
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
        self.__cart_calculated[:] = False


    @property
    def true_anomaly(self):
        '''true Anomaly
        '''
        if self.type == 'true':
            return self.anom
        
        self.calc_true_anomaly()
        return self._true_anomaly

    @true_anomaly.setter
    def true_anomaly(self, value):
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
        self.__cart_calculated[:] = False


    @property
    def mean_anomaly(self):
        '''Mean Anomaly
        '''
        if self.type == 'mean':
            return self.anom

        self.calc_mean_anomaly()
        return self._mean_anomaly


    @mean_anomaly.setter
    def mean_anomaly(self, value):
        if self.type == 'mean':
            self.anom = value
        elif self.type == 'true':
            self.anom = functions.mean_to_true(
                value, 
                self.e,
                tol = self.tol,
                degrees = self.degrees,
            )
        elif self.type == 'eccentric':
            self.anom = functions.mean_to_eccentric(
                value, 
                self.e,
                tol = self.tol,
                degrees = self.degrees,
            )
        self.__cart_calculated[:] = False


    @property
    def period(self):
        '''Orbital period
        '''
        self._kep_check()
        return functions.orbital_period(self.a, self.G*(self.M0 + self.m))

    @period.setter
    def period(self):
        self.a = functions.semi_major_axis(self.period, self.G*(self.M0 + self.m))
        self.__cart_calculated[:] = False


    @property
    def velocity(self):
        '''Orbital velocity
        '''
        if not self.__cart_calculated:
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
        else:
            return np.linalg.norm(self.v)

    @velocity.setter
    def velocity(self, value):
        self.v *= value/np.linalg.norm(self.v, axis=0)


Orbit.speed = Orbit.velocity