#!/usr/bin/env python

'''
Test Orbit class
'''

import time

import unittest
import numpy as np
import numpy.testing as nt

import pyorb
from pyorb import Orbit

class TestOrbit(unittest.TestCase):

    def setUp(self):
        self.cart_orb = dict(x=8e-01,y=0,z=0,vx=0,vy=7.7,vz=0)
        self.kep_orb = dict(a=1,e=0.2,i=0,omega=0,Omega=0,anom=0)
        self.G = pyorb.get_G(length='AU', mass='Msol', time='y')
        self.M = 1

    def test_init(self):
        orb = Orbit(0)

    def test_init_exception(self):
        with self.assertRaises(TypeError):
            orb = Orbit()
            
    def test_dtype(self):
        orb = Orbit(1)
        assert orb.cartesian.dtype==np.float64
        assert orb.dtype==np.float64
        orb = Orbit(1, dtype = np.float32)
        assert orb.cartesian.dtype==np.float32
        assert orb.kepler.dtype==np.float32
        assert orb.dtype==np.float32

        orb.add(num=10)
        orb.x = 2.0

        assert orb.cartesian.dtype==np.float32
        assert orb.kepler.dtype==np.float32
        assert orb.x.dtype==np.float32
        assert orb.dtype==np.float32


    def test_disable_update_anoms(self):
        orb = Orbit(M0=self.M, G=self.G, direct_update=False, auto_update=False, **self.cart_orb)
        orb.calculate_kepler()

        mu0 = orb.mean_anomaly
        orb.x += 0.1
        mu = orb.mean_anomaly
        nt.assert_array_equal(mu0, mu)

        mu0 = orb.mean_anomaly
        nu0 = orb.anom

        orb.vx -= 3
        orb.calculate_kepler()

        mu = orb.mean_anomaly
        nu = orb.anom

        with self.assertRaises(AssertionError):
            nt.assert_array_equal(mu0, mu)
        with self.assertRaises(AssertionError):
            nt.assert_array_equal(nu0, nu)

        nu0 = orb.anom
        cart0 = orb.cartesian
        orb.mean_anomaly = 0
        nu = orb.anom
        cart = orb.cartesian
        with self.assertRaises(AssertionError):
            nt.assert_array_equal(nu0, nu)
        nt.assert_array_equal(cart0, cart)


    def test_disable_update(self):
        orb = Orbit(M0=self.M, G=self.G, direct_update=False, auto_update=False, **self.cart_orb)
        assert np.all(np.isnan(orb.kepler))
        orb.x += 1
        assert np.all(np.isnan(orb.kepler))

        orb.calculate_kepler()

        kep0 = orb.kepler
        print('Before cart change')
        print(orb)

        orb.x += 1

        kep = orb.kepler
        print('After cart change')
        print(orb)

        nt.assert_array_equal(kep0, kep)

        orb = Orbit(M0=self.M, G=self.G, direct_update=False, auto_update=False, **self.kep_orb)
        assert np.all(np.isnan(orb.cartesian))
        orb.a += 1
        assert np.all(np.isnan(orb.cartesian))

        orb.a = 1
        orb.calculate_cartesian()

        cart0 = orb.cartesian
        print('Before kep change')
        print(orb)

        orb.a += 1

        cart = orb.cartesian
        print('After kep change')
        print(orb)

        nt.assert_array_equal(cart0, cart)


    def test_direct_update_cart(self):
        orb = Orbit(M0=self.M, G=self.G, direct_update=True, **self.cart_orb)
        assert orb.direct_update is True
        assert np.all(np.logical_not(np.isnan(orb.kepler)))

        kep0 = orb.kepler
        print('Before cart change')
        print(orb)

        orb.x += 1

        kep = orb.kepler
        print('After cart change')
        print(orb)

        with self.assertRaises(AssertionError):
            nt.assert_array_equal(kep0, kep)

        orb = Orbit(M0=self.M, G=self.G, direct_update=False, **self.cart_orb)
        assert orb.direct_update is False
        assert np.all(np.isnan(orb._kep)), 'internally before access should be nan'
        assert np.all(np.logical_not(np.isnan(orb.kepler))), 'With auto_update should be calculated'

        kep0 = orb.kepler
        print('Before cart change')
        print(orb)

        orb.x = 2

        kep = orb._kep.copy()
        print('After cart change')
        print(orb)

        kep2 = orb.kepler
        
        nt.assert_array_equal(kep0, kep)
        with self.assertRaises(AssertionError):
            nt.assert_array_equal(kep0, kep2)


    def test_direct_update_cart(self):
        orb = Orbit(M0=self.M, G=self.G, direct_update=True, **self.kep_orb)
        assert orb.direct_update is True
        assert np.all(np.logical_not(np.isnan(orb.cartesian)))

        kep0 = orb.cartesian
        print('Before cart change')
        print(orb)

        orb.a += 1

        kep = orb.cartesian
        print('After cart change')
        print(orb)

        with self.assertRaises(AssertionError):
            nt.assert_array_equal(kep0, kep)

        orb = Orbit(M0=self.M, G=self.G, direct_update=False, **self.kep_orb)
        assert orb.direct_update is False
        assert np.all(np.isnan(orb._cart)), 'internally before access should be nan'
        assert np.all(np.logical_not(np.isnan(orb.cartesian))), 'With auto_update should be calculated'

        kep0 = orb.cartesian
        print('Before cart change')
        print(orb)

        orb.a = 2

        kep = orb._cart.copy()
        print('After cart change')
        print(orb)

        kep2 = orb.cartesian
        
        nt.assert_array_equal(kep0, kep)
        with self.assertRaises(AssertionError):
            nt.assert_array_equal(kep0, kep2)



    def test_kepler_read_only(self):
        orb = Orbit(M0=self.M, G=self.G, kepler_read_only=True, **self.cart_orb)

        a = orb.a

        with self.assertRaises(AttributeError):
            orb.a = 1

        orb.x += 0.1

        with self.assertRaises(AssertionError):
            nt.assert_array_equal(a, orb.a)


    def test_kepler_read_only(self):
        orb = Orbit(M0=self.M, G=self.G, cartesian_read_only=True, **self.cart_orb)

        x = orb.x

        with self.assertRaises(AttributeError):
            orb.x = 2

        orb.a += 0.1

        with self.assertRaises(AssertionError):
            nt.assert_array_equal(x, orb.x)

    def test_solver_options(self):
        orb1 = Orbit(M0=self.M, G=self.G, solver_options=dict(max_iter=1), **self.kep_orb)
        orb2 = Orbit(M0=self.M, G=self.G, solver_options=dict(max_iter=500), **self.kep_orb)

        orb1.mean_anomaly += 0.545657871
        orb2.mean_anomaly += 0.545657871

        assert orb1.true_anomaly != orb2.true_anomaly

        orb1 = Orbit(M0=self.M, G=self.G, solver_options=dict(tol=1e-12), **self.kep_orb)
        orb2 = Orbit(M0=self.M, G=self.G, solver_options=dict(tol=1e-1), **self.kep_orb)

        orb1.mean_anomaly += 0.545657871
        orb2.mean_anomaly += 0.545657871

        assert orb1.true_anomaly != orb2.true_anomaly

        orb1 = Orbit(M0=self.M, G=self.G, solver_options=dict(tol=1e-12), **self.kep_orb)
        orb2 = Orbit(M0=self.M, G=self.G, solver_options=dict(tol=1e-9), **self.kep_orb)

        orb1.mean_anomaly += 0.545657871
        orb2.mean_anomaly += 0.545657871

        nt.assert_almost_equal(orb1.eccentric_anomaly, orb2.eccentric_anomaly, decimal=8)

    def test_mass(self):
        orb1 = Orbit(M0=self.M, G=self.G, m=0.1, **self.kep_orb)
        orb2 = Orbit(M0=self.M, G=self.G, m=0.0, **self.kep_orb)

        with self.assertRaises(AssertionError):
            nt.assert_almost_equal(orb1.cartesian, orb2.cartesian)

        orb1.m[:] = 0
        orb1.calculate_cartesian()
        nt.assert_almost_equal(orb1.cartesian, orb2.cartesian)
