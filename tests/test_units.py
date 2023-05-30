#!/usr/bin/env python

'''
Test Units
'''
import unittest
import numpy as np
import numpy.testing as nt

import pyorb
import pyorb.unit as unit


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

class TestWrappers(unittest.TestCase):

    def test_angle_units_pos_arg(self):

        @unit.angle_units([1], None, None, degrees=False)
        def test_func(text, th):
            nt.assert_almost_equal(th, np.pi)

        test_func('hello', 180, degrees=True)
        test_func('hello', np.pi, degrees=False)
        test_func('hello', np.pi)

    def test_angle_units_multi_arg(self):

        @unit.angle_units([0, 1], None, None, degrees=False)
        def test_func(ph, th):
            nt.assert_almost_equal(th, np.pi)
            nt.assert_almost_equal(ph, 2*np.pi)

        test_func(360, 180, degrees=True)
        test_func(2*np.pi, np.pi, degrees=False)
        test_func(2*np.pi, np.pi)

    def test_angle_units_kw_arg(self):

        @unit.angle_units(None, ['th'], None, degrees=False)
        def test_func(text='foo', th=0):
            nt.assert_almost_equal(th, np.pi)

        test_func(th=180, degrees=True)
        test_func(th=np.pi, degrees=False)
        test_func(th=np.pi)

    def test_angle_units_out(self):

        @unit.angle_units(None, None, True, degrees=False)
        def test_func():
            return np.pi

        nt.assert_almost_equal(test_func(degrees=True), 180)
        nt.assert_almost_equal(test_func(degrees=False), np.pi)
        nt.assert_almost_equal(test_func(), np.pi)

        @unit.angle_units(None, None, [0], degrees=False)
        def test_func():
            return np.pi, 'hellos'

        nt.assert_almost_equal(test_func(degrees=True)[0], 180)
        nt.assert_almost_equal(test_func(degrees=False)[0], np.pi)
        nt.assert_almost_equal(test_func()[0], np.pi)

    def test_angle_units_all(self):

        @unit.angle_units([0], ['th'], True, degrees=False)
        def test_func(ph, th=0):
            nt.assert_almost_equal(ph, np.pi)
            nt.assert_almost_equal(th, 2*np.pi)
            return 3*np.pi

        ret = test_func(180, th=360, degrees=True)
        nt.assert_almost_equal(ret, 180*3)
        ret = test_func(np.pi, th=2*np.pi, degrees=False)
        nt.assert_almost_equal(ret, 3*np.pi)
        ret = test_func(np.pi, th=2*np.pi)
        nt.assert_almost_equal(ret, 3*np.pi)
