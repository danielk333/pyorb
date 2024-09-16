#!/usr/bin/python
# :encoding=utf-8

# Documentation for this module is in `_const_docs.py` due to
# difficulties of extracting docstrings when masquerading instances as modules

# Recipes from Python Cookbook
# const class masquerading as a module:
#   https://www.safaribooksonline.com/library/view/python-cookbook/0596001673/ch05s16.html
# Borg class ensuring that every instance has shared state:
#   https://www.safaribooksonline.com/library/view/Python+Cookbook/0596001673/ch05s23.html
import numpy
import sys


class _const:

    class ConstError(TypeError):
        pass
    _shared_state = {}

    def __init__(self):
        self.__dict__ = self._shared_state

    def __setattr__(self, name, value):
        if name in self.__dict__:
            if self.__dict__[name] == value:
                return
            raise self.ConstError("Can't rebind const({})".format(name))
        self.__dict__[name] = value

    def __delattr__(self, name):
        if name in self.__dict__:
            raise self.ConstError("Can't unbind const({})".format(name))
        raise NameError(name)

_self = _const()


# Mathematical constants
_self.pi = numpy.pi                     # You shouldn't have to ask
_self.e = numpy.e                       # Euler's number, not elementary charge
_self._2pi = 2*numpy.pi                 # How it is usually used
_self._4pi = 4*numpy.pi                 # How it is usually used

# Fundamental physics constants
_self.c = 299792458.0                   # light speed in vacuum [m/s]
_self.G = 6.6743e-11                    # Newtons gravitational constant [m^3 kg^-1 s^-2]
# **Reference**: NIST - http://physics.nist.gov/cgi-bin/cuu/Value?bg

# WGS-84 Earth model
_self.wgs84_a = 6378137.0               # Semi-major axis [m] (by definition)
_self.wgs84_invf = 298.257223563        # inverse flattening (by definition)
_self.wgs84_f = 1/_self.wgs84_invf      # from whence flattening arises
_self.wgs84_omega = 7292115e-11         # Nominal Mean Angular Velocity [rad/s]
# Geocentric Gravitational Constant [m3/s2]
_self.wgs84_GM = 3986004.418e8
_self.wgs84_M_earth = _self.wgs84_GM / _self.G

# Astronomical quantities
_self.AU = 1.495978707e11               # Astronomical unit [m]
                                        # mean distance between the sun and the earth
_self.M_sol  = 1.98847e30               # [kg]
_self.GM_sol = 1.32712440018e20         # m^3/s^2


# Symbolic names for the indices that indicate keplerian and equinoctial elements
_self.K_a, _self.K_e, _self.K_i, _self.K_om, _self.K_OM, _self.K_nu = 0, 1, 2, 3, 4, 5
_self.E_a, _self.E_h, _self.E_k, _self.E_p, _self.E_q, _self.E_lam  = 0, 1, 2, 3, 4, 5


sys.modules[__name__] = _self
