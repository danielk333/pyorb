#!/usr/bin/env python

'''Unit handling

'''

# Python standard import


# Third party import
import numpy as np


# Local import
from .kepler import G as G_SI
from .kepler import M_sol, AU


def get_G(length, mass, time):

    if isinstance(mass, str):
        mass = mass.lower().replace('_', '').replace(' ', '').replace('-', '')
        if mass == 'kg':
            _mass = 1.0
        elif mass == 'g':
            _mass = 1e-3
        elif mass == 'msun' or mass == 'msol':
            _mass = M_sol
        else:
            raise TypeError(f'Unit "{mass}" not recognized')
    elif isinstance(mass, float):
        _mass = mass
    else:
        raise TypeError(f'Type "{type(mass)}" for mass not supported')

    if isinstance(time, str):
        time = time.lower()
        if time == 's':
            _time = 1.0
        elif time == 'h':
            _time = 3600.0
        elif time == 'd':
            _time = 3600.0*24
        elif time == 'y':
            _time = 3600.0*24*365.25
        else:
            raise TypeError(f'Unit "{time}" not recognized')
    elif isinstance(time, float):
        _time = time
    else:
        raise TypeError(f'Type "{type(time)}" for time not supported')

    if isinstance(length, str):
        length = length.lower()
        if length == 'm':
            _length = 1.0
        elif length == 'cm':
            _length = 1e-2
        elif length == 'km':
            _length = 1e3
        elif length == 'au':
            _length = AU
        elif length == 'pc':
            _length = 3.08567758149137e16  # IAU 2012 exact SI def
        else:
            raise TypeError(f'Unit "{length}" not recognized')
    elif isinstance(length, float):
        _length = length
    else:
        raise TypeError(f'Type "{type(length)}" for length not supported')

    return G_SI*(_mass*_time**2/_length**3)


def angle_units(in_arg_inds, in_arg_keys, out_arg_inds):
    '''Wrapper to automatically convert input arguemnts from degrees to radians 
    and back to degrees if the keyword argument `degrees = True`.
    '''

    def angle_converter_warpper(func):
        def wrapped_func(*args, **kwargs):
            if not kwargs.pop('degrees', False):
                return func(*args, **kwargs)

            args = list(args)
            if in_arg_inds is not None:
                for ind in in_arg_inds:
                    args[ind] = np.radians(args[ind])
            if in_arg_keys is not None:
                for key in in_arg_keys:
                    if key in kwargs:
                        kwargs[key] = np.radians(kwargs[key])

            ret = func(*args, **kwargs)

            ret = list(ret)
            if out_arg_inds is not None:
                for ind in out_arg_inds:
                    ret[ind] = np.degrees(ret[ind])

                return ret
        return wrapped_func
    return angle_converter_warpper
