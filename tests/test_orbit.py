#!/usr/bin/env python

'''
Test Orbit class
'''
import unittest
import numpy as np
import numpy.testing as nt

import pyorb
from pyorb import Orbit


class TestOrbit(unittest.TestCase):

    def setUp(self):
        self.cart_orb = dict(x=8e-01, y=0, z=0, vx=0, vy=7.7, vz=0)
        self.kep_orb = dict(a=1, e=0.2, i=0, omega=0, Omega=0, anom=0)
        self.G = pyorb.get_G(length='AU', mass='Msol', time='y')
        self.M = 1

    def test_init(self):
        orb = Orbit(0)

    def test_init_exception(self):
        with self.assertRaises(TypeError):
            orb = Orbit()

    def test_dtype(self):
        orb = Orbit(1)
        assert orb.cartesian.dtype == np.float64
        assert orb.dtype == np.float64
        orb = Orbit(1, dtype = np.float32)
        assert orb.cartesian.dtype == np.float32
        assert orb.kepler.dtype == np.float32
        assert orb.dtype == np.float32

        orb.add(num=10)
        orb.x = 2.0

        assert orb.cartesian.dtype == np.float32
        assert orb.kepler.dtype == np.float32
        assert orb.x.dtype == np.float32
        assert orb.dtype == np.float32

    def test_disable_update_anoms(self):
        orb = Orbit(M0=self.M, G=self.G, direct_update=False,
                    auto_update=False, **self.cart_orb)
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
        orb = Orbit(M0=self.M, G=self.G, direct_update=False,
                    auto_update=False, **self.cart_orb)
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

        orb = Orbit(M0=self.M, G=self.G, direct_update=False,
                    auto_update=False, **self.kep_orb)
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

    def test_direct_update_kep(self):
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
        assert np.all(np.isnan(orb._kep)
                      ), 'internally before access should be nan'
        assert np.all(np.logical_not(np.isnan(orb.kepler))
                      ), 'With auto_update should be calculated'

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
        assert np.all(np.isnan(orb._cart)
                      ), 'internally before access should be nan'
        assert np.all(np.logical_not(np.isnan(orb.cartesian))
                      ), 'With auto_update should be calculated'

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

    def test_cart_read_only(self):
        orb = Orbit(M0=self.M, G=self.G,
                    kepler_read_only=True, **self.cart_orb)

        a = orb.a

        with self.assertRaises(AttributeError):
            orb.a = 1

        orb.x += 0.1

        with self.assertRaises(AssertionError):
            nt.assert_array_equal(a, orb.a)

    def test_kepler_read_only(self):
        orb = Orbit(M0=self.M, G=self.G,
                    cartesian_read_only=True, **self.cart_orb)

        x = orb.x

        with self.assertRaises(AttributeError):
            orb.x = 2

        orb.a += 0.1

        with self.assertRaises(AssertionError):
            nt.assert_array_equal(x, orb.x)

    def test_solver_options(self):
        orb1 = Orbit(M0=self.M, G=self.G, solver_options=dict(
            max_iter=1), **self.kep_orb)
        orb2 = Orbit(M0=self.M, G=self.G, solver_options=dict(
            max_iter=500), **self.kep_orb)

        orb1.mean_anomaly += 0.545657871
        orb2.mean_anomaly += 0.545657871

        assert orb1.true_anomaly != orb2.true_anomaly

        orb1 = Orbit(M0=self.M, G=self.G, solver_options=dict(
            tol=1e-12), **self.kep_orb)
        orb2 = Orbit(M0=self.M, G=self.G, solver_options=dict(
            tol=1e-1), **self.kep_orb)

        orb1.mean_anomaly += 0.545657871
        orb2.mean_anomaly += 0.545657871

        assert orb1.true_anomaly != orb2.true_anomaly

        orb1 = Orbit(M0=self.M, G=self.G, solver_options=dict(
            tol=1e-12), **self.kep_orb)
        orb2 = Orbit(M0=self.M, G=self.G, solver_options=dict(
            tol=1e-9), **self.kep_orb)

        orb1.mean_anomaly += 0.545657871
        orb2.mean_anomaly += 0.545657871

        nt.assert_almost_equal(orb1.eccentric_anomaly,
                               orb2.eccentric_anomaly, decimal=8)

    def test_mass(self):
        orb1 = Orbit(M0=self.M, G=self.G, m=0.1, **self.kep_orb)
        orb2 = Orbit(M0=self.M, G=self.G, m=0.0, **self.kep_orb)

        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(orb1.cartesian, orb2.cartesian)

        orb1.m[:] = 0
        orb1.calculate_cartesian()
        nt.assert_array_almost_equal(orb1.cartesian, orb2.cartesian)

    def test_num(self):
        orb = Orbit(M0=self.M, G=self.G, num=10)
        assert orb.num == 10
        assert orb._kep.shape[1] == 10
        assert orb._cart.shape[1] == 10
        i = 0
        for o in orb:
            i += 1
        assert i == 10

    def test_update(self):
        orb = Orbit(M0=self.M, G=self.G, m=0.1, **self.kep_orb)

        with self.assertRaises(ValueError):
            orb.update(a=2, m=3)

        with self.assertRaises(ValueError):
            orb.update(a=2, x=3)

        orb.kepler_read_only = True
        with self.assertRaises(AttributeError):
            orb.update(a=2)
        orb.kepler_read_only = False

        orb.cartesian_read_only = True
        with self.assertRaises(AttributeError):
            orb.update(x=3)
        orb.cartesian_read_only = False

        orb.direct_update = False
        kep0 = orb.kepler
        cart0 = orb.cartesian
        orb.update(a=2, m=3)
        kep0[0] = 2
        nt.assert_array_almost_equal(orb.kepler, kep0)
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(orb.cartesian, cart0)

    def test_allocate(self):
        orb = Orbit(M0=self.M, G=self.G, m=0.1, **self.kep_orb)
        assert orb.num == 1
        orb.allocate(10)
        assert orb.num == 10
        orb.add(num=5)
        assert orb.num == 15
        orb.allocate(10)
        assert orb.num == 10
        assert np.all(np.isnan(orb.cartesian))
        assert np.all(np.isnan(orb.kepler))
        nt.assert_array_almost_equal(orb.m, np.zeros_like(orb.m))

    def test_propagate(self):
        orb = Orbit(M0=self.M, G=self.G, m=0.1, **self.kep_orb)
        kep0 = orb.kepler
        cart0 = orb.cartesian
        orb.propagate(orb.period)
        nt.assert_array_almost_equal(orb.kepler, kep0)
        nt.assert_array_almost_equal(orb.cartesian, cart0)

        orb = Orbit(M0=self.M, G=self.G, a=1, e=0,
                    i=0, omega=0, Omega=0, anom=0)
        cart0 = orb.cartesian
        orb.propagate(orb.period*0.5)
        nt.assert_almost_equal(-orb.cartesian[0], cart0[0])
        nt.assert_array_almost_equal(orb.cartesian[1:3], cart0[1:3])
        nt.assert_almost_equal(-orb.cartesian[4], cart0[4])
        nt.assert_almost_equal(orb.cartesian[3], cart0[3])
        nt.assert_almost_equal(orb.cartesian[5], cart0[5])

    def test_delete(self):
        orb = Orbit(M0=self.M, G=self.G, m=0.1, **self.kep_orb)

        a0 = self.kep_orb['a']
        inds = [2, 4, 6]

        orb.add(num=9, **self.kep_orb)
        a = orb.a
        a[inds] = a0 + 2
        orb.a = a

        orb.delete(inds)

        assert orb.num == 10-3
        assert np.all(orb.a < a0 + 1)

    def test_type(self):
        orb = Orbit(M0=self.M, G=self.G, m=0.1, degrees=True, **self.kep_orb)
        orb.anom = 90
        anom = orb.anom
        m0 = orb.mean_anomaly
        orb.type = 'mean'

        nt.assert_almost_equal(m0, orb.mean_anomaly)
        nt.assert_almost_equal(m0, orb.anom)
        nt.assert_almost_equal(anom, orb.true_anomaly)


class TestOrbitProperties(unittest.TestCase):

    def setUp(self):
        self.kep_orb = dict(
            a=1, e=0.2, i=3, 
            omega=10, Omega=20, anom=90, 
            degrees=True,
        )
        self.G = pyorb.get_G(length='AU', mass='Msol', time='y')
        self.M0 = 1
        self.m = 0.1
        self.orb = Orbit(M0=self.M0, G=self.G, m=self.m, **self.kep_orb)
        self.cart = self.orb.cartesian
        self.kep = self.orb.kepler

    def test_r(self):
        nt.assert_array_almost_equal(self.orb.r, self.orb.cartesian[:3])

    def test_v(self):
        nt.assert_array_almost_equal(self.orb.v, self.orb.cartesian[3:])

    def test_velocity(self):
        vel0 = self.orb.velocity
        nt.assert_almost_equal(vel0, np.linalg.norm(self.orb.cartesian[3:]))
        self.orb.velocity = vel0*1.1
        nt.assert_almost_equal(self.orb.velocity, vel0*1.1)
        nt.assert_almost_equal(np.linalg.norm(
            self.orb.cartesian[3:]), vel0*1.1)

    def test_speed_vs_velocity(self):
        nt.assert_almost_equal(self.orb.speed, self.orb.velocity, decimal=2)

    def test_speed(self):
        # circular orbits have same speed
        self.orb.e = 0
        s0 = self.orb.speed
        self.orb.anom = 0
        nt.assert_almost_equal(s0, self.orb.speed)

    def test_period(self):
        p0 = self.orb.period
        self.orb.a += 1
        with self.assertRaises(AssertionError):
            nt.assert_almost_equal(p0, self.orb.period)

        nt.assert_almost_equal(self.orb.period*self.orb.mean_motion, 360.0)
        self.orb.degrees = False
        nt.assert_almost_equal(self.orb.period*self.orb.mean_motion, 2*np.pi)

    def test_a(self):
        assert self.orb.a == self.kep_orb['a']

    def test_set_a(self):
        self.orb.a = 2
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(self.orb.cartesian, self.cart)

    def test_e(self):
        assert self.orb.e == self.kep_orb['e']

    def test_set_e(self):
        self.orb.e = 0.1
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(self.orb.cartesian, self.cart)

    def test_i(self):
        assert self.orb.i == self.kep_orb['i']

    def test_set_i(self):
        self.orb.i = 20
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(self.orb.cartesian, self.cart)

    def test_omega(self):
        assert self.orb.omega == self.kep_orb['omega']

    def test_set_omega(self):
        self.orb.omega = 0
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(self.orb.cartesian, self.cart)

    def test_Omega(self):
        assert self.orb.Omega == self.kep_orb['Omega']

    def test_set_Omega(self):
        self.orb.Omega = 0
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(self.orb.cartesian, self.cart)

    def test_anom(self):
        assert self.orb.anom == self.kep_orb['anom']

    def test_set_anom(self):
        self.orb.anom = 0
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(self.orb.cartesian, self.cart)

    def test_x(self):
        assert self.orb.x == self.cart[0]

    def test_set_x(self):
        self.orb.x = 2
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(self.orb.kepler, self.kep)

    def test_y(self):
        assert self.orb.y == self.cart[1]

    def test_set_y(self):
        self.orb.y = 2
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(self.orb.kepler, self.kep)

    def test_z(self):
        assert self.orb.z == self.cart[2]

    def test_set_z(self):
        self.orb.z = 2
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(self.orb.kepler, self.kep)

    def test_vx(self):
        assert self.orb.vx == self.cart[3]

    def test_set_vx(self):
        self.orb.vx = 2
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(self.orb.kepler, self.kep)

    def test_vy(self):
        assert self.orb.vy == self.cart[4]

    def test_set_vy(self):
        self.orb.vy = 2
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(self.orb.kepler, self.kep)

    def test_vz(self):
        assert self.orb.vz == self.cart[5]

    def test_set_vz(self):
        self.orb.vz = 2
        with self.assertRaises(AssertionError):
            nt.assert_array_almost_equal(self.orb.kepler, self.kep)
