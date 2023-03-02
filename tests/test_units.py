#!/usr/bin/env python

'''
Test Units
'''
import unittest
import numpy as np

import pyorb


class TestUnits(unittest.TestCase):

    def test_variables(self):
        var = ['G', 'AU', 'M_earth', 'M_sol']
        for key in var:
            assert hasattr(pyorb, key)

    def test_get_G(self):
        G_SI = pyorb.G

        with self.assertRaises(TypeError):
            _ = pyorb.get_G(1, 1.0, 1.0)
        with self.assertRaises(TypeError):
            _ = pyorb.get_G(2.0, 'test', 1.0)
        with self.assertRaises(TypeError):
            class A:
                pass
            _ = pyorb.get_G(2.0, A, 1.0)

        assert G_SI == pyorb.get_G(1.0, 1.0, 1.0)
        assert G_SI == 0.5*pyorb.get_G(1.0, 2.0, 1.0)
        assert G_SI == pyorb.get_G('m', 'kg', 's')
        assert G_SI == pyorb.get_G(length='m', mass='kg', time='s')

        pc = 3.08567758149137e16  # IAU 2012 exact SI def [meter]
        G_ast = pyorb.get_G(length='km', mass='Msol', time='s')
        # Change one unit of "km" to "pc"
        # to get (km/s)^2 pc / M_sol
        G_ast /= pc*1e-3

        assert np.abs(G_ast - 4.30091e-3) < 1e-6
