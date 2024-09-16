#!/usr/bin/env python

'''
Test basic kepler functions
'''

import unittest
import pytest
import numpy as np
import numpy.testing as nt

import pyorb.kepler as kep
import pyorb


class TestAlternativeParameters(unittest.TestCase):

    def setUp(self):
        self.a = 7500e3         # around 800 km above Earth
        self.orb_init = np.array(
            [self.a, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
        self.m = 1.0
        self.mu = pyorb.GM_earth

        s2 = np.sqrt(2)/2
        s3 = np.sqrt(3)/2
        t3 = 1/np.sqrt(3)
        a = self.a
        cases = np.array([
            # kepl a e i o O n equi a h k p q l
            [[a, 0, 0,   0,  0, 0], [a,  0,  0,  0, 0,   0]],
            [[a, 0, 0,  30,  0, 0], [a,  0,  0,  0, 0,  30]],
            [[a, 0, 90,  0,  0, 0], [a,  0,  0,  0, 1,   0]],
            [[a, 0, 90,  0, 60, 0], [a,  0,  0, s3, .5, 60]],
            [[a, 0, 90, 45, 30, 0], [a,  0,  0, .5, s3, 75]],
            [[a, .5, 0, 45,  0, 0], [a, s2/2, s2/2,  0, 0,  45]],
            [[a, .5, 60, 30, 30, 0],[a, s3/2, .25, t3/2, t3*s3, 60]]
        ])
        cart = kep.kep_to_cart(cases[:, 0].T, mu=self.mu, degrees=True).T
        self.testcases = dict(kepl=cases[:, 0], equi=cases[:, 1], cart=cart)

    def test_cart_to_equi(self):
        total_tests = len(self.testcases['kepl'])
        for ix in range(len(self.testcases['kepl'])):
            _cart = self.testcases['cart'][ix]
            _equi = kep.cart_to_equi(_cart, mu=self.mu, degrees=True)
            nt.assert_array_almost_equal(
                _equi,
                self.testcases['equi'][ix],
                decimal = 6,
                err_msg = f"from cart to equi case {ix}/{total_tests}",
            )

    def test_equi_to_cart(self):
        total_tests = len(self.testcases['equi'])
        for ix in range(len(self.testcases['equi'])):
            _equi = self.testcases['equi'][ix]
            _cart = kep.equi_to_cart(_equi, mu=self.mu, degrees=True)
            nt.assert_array_almost_equal(
                _cart,
                self.testcases['cart'][ix],
                decimal = 6,
                err_msg = f"from equi to cart case {ix}/{total_tests}",
            )

    def test_equi_cart_consistency(self):
        total_tests = len(self.testcases['equi'])*2
        test_num = 0

        for ix in range(len(self.testcases['equi'])):
            _equi0 = self.testcases['equi'][ix].copy()
            _cart = kep.equi_to_cart(_equi0, degrees=True)
            _equi = kep.cart_to_equi(_cart, degrees=True)

            nt.assert_array_almost_equal(
                _equi,
                _equi0,
                decimal = 6,
                err_msg = (
                    f"from equi to cart to equi case {ix} [{test_num=}/{total_tests}]\n"
                    f"Input kepler {self.testcases['kepl'][ix]}\n"
                    f"Input cart {self.testcases['cart'][ix]}\n"
                ),
            )
            test_num += 1

        for ix in range(len(self.testcases['cart'])):
            _cart0 = self.testcases['cart'][ix]
            _equi = kep.cart_to_equi(_cart0, degrees=True)
            _cart = kep.equi_to_cart(_equi, degrees=True)

            # we are comparing absolute values to decimals in SI so we can relax to account
            # for numerical unstabillity of some of the trig functions
            nt.assert_array_almost_equal(
                _cart,
                _cart0,
                decimal = 1,
                err_msg = (
                    f"from cart to equi to cart case {ix} [{test_num=}/{total_tests}]\n"
                    f"Input kepler {self.testcases['kepl'][ix]}\n"
                    f"Input equi {self.testcases['equi'][ix]}\n"
                ),
            )
            test_num += 1

    def test_kep_to_equi(self):
        orb_init = self.orb_init.copy()

        x = kep.kep_to_equi(orb_init, degrees=True)
        try:
            nt.assert_array_almost_equal(x, orb_init)
        except AssertionError:
            print(x)
            raise

        # one at a time
        for ix in range(len(self.testcases['kepl'])):
            try:
                assert np.allclose(
                    kep.kep_to_equi(self.testcases['kepl'][ix], degrees=True),
                    self.testcases['equi'][ix],
                )
            except AssertionError:
                print(f"kep_to_equi, case {ix}")
                print(f"input {self.testcases['kepl'][ix]}")
                print(f"     =>  {kep.kep_to_equi(self.testcases['kepl'][ix], degrees=True)}")
                print(f" expected {self.testcases['equi'][ix]}")
                raise

        # Then all at once
        assert np.allclose(
            kep.kep_to_equi(self.testcases['kepl'].T, degrees=True).T,
            self.testcases['equi'],
        )

    def test_equi_to_kep(self):
        orb_init = self.orb_init.copy()

        x = kep.equi_to_kep(orb_init, degrees=True)
        try:
            # This test case looks the same in Kepler and Equinoctial
            nt.assert_array_almost_equal(x, orb_init)
        except AssertionError:
            print('Problem!')
            raise

        # one at a time
        for ix in range(len(self.testcases['equi'])):
            try:
                nt.assert_(
                    kep.kep_equivalent(
                        kep.equi_to_kep(self.testcases['equi'][ix], degrees=True),
                        self.testcases['kepl'][ix]
                    ),
                    f'Error in test case [{ix}]',
                )
            except AssertionError:
                print(f"equi_to_kep, case {ix}")
                print(f"input {self.testcases['equi'][ix]}")
                print(f"      =>  {kep.equi_to_kep(self.testcases['equi'][ix], degrees=True)}")
                print(f" expected {self.testcases['kepl'][ix]}")
                raise

        # Then all at once
        # assert np.allclose(kep.equi_to_kep(self.testcases['equi'].T, degrees=True).T,
        #                 self.testcases['kepl'])

    def test_equi_kep_consistency(self):
        total_tests = len(self.testcases['equi'])*2
        test_num = 0

        for ix in range(len(self.testcases['equi'])):
            _equi0 = self.testcases['equi'][ix].copy()
            _kep = kep.equi_to_kep(_equi0, degrees=True)
            _equi = kep.kep_to_equi(_kep, degrees=True)

            nt.assert_array_almost_equal(
                _equi,
                _equi0,
                decimal = 6,
                err_msg = f"from equi to kep to equi case {ix} [{test_num=}/{total_tests}]",
            )
            test_num += 1

        for ix in range(len(self.testcases['kepl'])):
            _kep0 = self.testcases['kepl'][ix]
            _equi = kep.kep_to_equi(_kep0, degrees=True)
            _kep = kep.equi_to_kep(_equi, degrees=True)

            if _kep0[kep.K_e] <= kep.e_lim:
                # This means all aop information will be lost and put into anomaly/longitude
                _kep0[kep.K_nu] += _kep0[kep.K_om]
                _kep0[kep.K_om] = 0

            nt.assert_array_almost_equal(
                _kep,
                _kep0,
                decimal = 6,
                err_msg = f"from kep to equi to kep case {ix} [{test_num=}/{total_tests}]",
            )
            test_num += 1


class TestKepCart(unittest.TestCase):

    def setUp(self):
        self.M_e = 5.972e24
        self.R_e = 6371e3
        self.a = self.R_e*1.5
        self.orb_init = np.array(
            [self.a, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
        self.m = 1.0
        self.mu = (pyorb.M_earth + self.m)*pyorb.G

    def test_cart_kep_loop_degradation(self):

        orb_init = np.array(
            [self.a, 0.2, 23.0, 136.0, 44.0, 10.0], dtype=np.float64)

        x = kep.kep_to_cart(orb_init, mu=self.mu, degrees=True)
        x_ref = x.copy()
        for ind in range(200):
            o = kep.cart_to_kep(x, mu=self.mu, degrees=True)
            x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        nt.assert_array_almost_equal(np.ones_like(x), x_ref/x, decimal = 1)

    def test_orientation_convention(self):

        e = 0.5
        orb0 = np.array([self.a, e, 0, 0, 0, 0], dtype=np.float64)
        q = self.a*(1 - e)
        Q = self.a*(1 + e)
        L = self.a*(1 - e**2)

        orb = orb0.copy()
        x = kep.kep_to_cart(orb, mu=self.mu, degrees=True)
        nt.assert_almost_equal(x[0], q, decimal = 6)

        orb = orb0.copy()
        orb[5] = 90
        x = kep.kep_to_cart(orb, mu=self.mu, degrees=True)
        nt.assert_almost_equal(x[1], L, decimal = 6)

        orb = orb0.copy()
        orb[3] = 180
        x = kep.kep_to_cart(orb, mu=self.mu, degrees=True)
        nt.assert_almost_equal(x[0], -q, decimal = 6)

        orb[5] = 180
        x = kep.kep_to_cart(orb, mu=self.mu, degrees=True)
        nt.assert_almost_equal(x[0], Q, decimal = 6)

        orb = orb0.copy()
        orb[3] = 90
        orb[4] = 90
        orb[5] = 180
        x = kep.kep_to_cart(orb, mu=self.mu, degrees=True)
        nt.assert_almost_equal(x[0], Q, decimal = 6)

        orb = orb0.copy()
        orb[2] = 90
        x = kep.kep_to_cart(orb, mu=self.mu, degrees=True)
        nt.assert_almost_equal(x[0], q, decimal = 6)
        nt.assert_almost_equal(x[1], 0, decimal = 6)

        orb[3] = 90
        x = kep.kep_to_cart(orb, mu=self.mu, degrees=True)
        nt.assert_almost_equal(x[0], 0, decimal = 6)
        nt.assert_almost_equal(x[1], 0, decimal = 6)
        nt.assert_almost_equal(x[2], q, decimal = 6)

        orb[4] = 180
        x = kep.kep_to_cart(orb, mu=self.mu, degrees=True)
        nt.assert_almost_equal(x[0], 0, decimal = 6)
        nt.assert_almost_equal(x[1], 0, decimal = 6)
        nt.assert_almost_equal(x[2], q, decimal = 6)

        orb = orb0.copy()
        orb[2] = 45
        orb[3] = 90
        orb[4] = 90
        x = kep.kep_to_cart(orb, mu=self.mu, degrees=True)
        nt.assert_almost_equal(x[0], -q*np.sqrt(2)*0.5, decimal = 6)
        nt.assert_almost_equal(x[1], 0, decimal = 6)
        nt.assert_almost_equal(x[2], q*np.sqrt(2)*0.5, decimal = 6)

    def test_planar_detection(self):
        '''Longitude of ascending node gets transfered to argument of periapsis
        in planar orbits. Test if this is true for retro-grade and pro-grade
        planar orbits.
        '''

        orb_in = np.array([self.a, 0.5, 0, 0, 45, 0], dtype=np.float64)
        orb_ref = np.array([self.a, 0.5, 0, 45, 0, 0], dtype=np.float64)

        x0 = kep.kep_to_cart(orb_in, mu=self.mu, degrees=True)
        orb_out = kep.cart_to_kep(x0, mu=self.mu, degrees=True)

        nt.assert_almost_equal(orb_ref, orb_out)

        orb_in = np.array([self.a, 0.5, 180, 0, 45, 0], dtype=np.float64)
        orb_ref = np.array([self.a, 0.5, 180, 360-45, 0, 0], dtype=np.float64)

        x1 = kep.kep_to_cart(orb_in, mu=self.mu, degrees=True)
        orb_out = kep.cart_to_kep(x1, mu=self.mu, degrees=True)

        # should be same but reversed velocity
        x1[3:] = -x1[3:]
        nt.assert_almost_equal(x0, x1)

        nt.assert_almost_equal(orb_ref, orb_out)

    @pytest.mark.slow
    def test_cart_kep_inverse(self):
        a = np.linspace(self.a, self.a + self.R_e, num=2, dtype=np.float64)
        e = np.linspace(0.0, 0.99, num=10, dtype=np.float64)
        inc = np.linspace(0.0, 180.0, num=10, dtype=np.float64)
        omega = np.linspace(0.0, 360.0, num=10, dtype=np.float64)
        Omega = omega.copy()
        nu = np.array([0, 45, 180], dtype=np.float64)
        av, ev, incv, omegav, Omegav, nuv = np.meshgrid(
            a, e, inc, omega, Omega, nu)

        ait = np.nditer(av)
        eit = np.nditer(ev)
        incit = np.nditer(incv)
        omegait = np.nditer(omegav)
        Omegait = np.nditer(Omegav)
        nuit = np.nditer(nuv)

        o = np.empty((6, av.size), dtype=np.float64)
        for ind in range(av.size):
            o[0, ind] = next(ait)
            o[1, ind] = next(eit)
            o[2, ind] = next(incit)
            o[3, ind] = next(omegait)
            o[4, ind] = next(Omegait)
            o[5, ind] = next(nuit)

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)
        o_calc = kep.cart_to_kep(x, mu=self.mu, degrees=True)
        x_calc = kep.kep_to_cart(o_calc, mu=self.mu, degrees=True)

        for ind in range(av.size):
            try:
                nt.assert_array_almost_equal(
                    x[:3, ind]/o[0, ind], x_calc[:3, ind]/o[0, ind],
                    decimal = 6,
                )
                nt.assert_array_almost_equal(
                    x[3:, ind]/o[0, ind], x_calc[3:, ind]/o[0, ind],
                    decimal = 6,
                )
            except AssertionError:
                print(ind)
                print(x[:3, ind] - x_calc[:3, ind])
                print(o[:, ind])
                print(o_calc[:, ind])
                print(o_calc[:, ind] - o[:, ind])
                raise

        el = o[1, :] < kep.e_lim

        il = o[2, :] < np.degrees(kep.i_lim)
        il_ret = o[2, :] > 180 - np.degrees(kep.i_lim)
        o_orig = o.copy()

        o[3, el] = 0.0
        o[4, il] = 0.0
        o[4, il_ret] = 0.0
        o[5, el] += o_orig[3, el]

        inds = np.logical_and(il, el)
        o[5, np.logical_and(il, el)] += o_orig[4, inds]

        inds = np.logical_and(il, np.logical_not(el))
        o[3, np.logical_and(il, np.logical_not(el))] += o_orig[4, inds]

        inds = np.logical_and(il_ret, el)
        o[5, np.logical_and(il_ret, el)] -= o_orig[4, inds]

        inds = np.logical_and(il_ret, np.logical_not(el))
        o[3, np.logical_and(il_ret, np.logical_not(el))] -= o_orig[4, inds]

        o[5, :] = np.mod(o[5, :] + 360.0, 360.0)

        for ind in range(av.size):
            test = np.isclose(o_calc[0, ind]/self.R_e,
                              o[0, ind]/self.R_e, atol=1e-3)
            test = np.logical_and(test, np.isclose(
                o_calc[1, ind], o[1, ind], atol=1e-4))
            test_ang = np.logical_or(
                np.isclose(o_calc[2:, ind], o[2:, ind], atol=1e-7),
                np.isclose(np.mod(o_calc[2:, ind]+180.0, 360.0),
                           np.mod(o[2:, ind]+180.0, 360.0), atol=1e-7),
            )
            test = np.logical_and(test, test_ang)
            if not np.all(test):

                # x_ = kep.kep_to_cart(o_orig[:,ind], mu=self.mu, degrees=True)
                # o_calc_ = kep.cart_to_kep(x_, mu=self.mu, degrees=True)
                # x_calc_ = kep.kep_to_cart(o_calc_, mu=self.mu, degrees=True)

                print(ind)
                print(o_orig[:, ind])
                print(o[:, ind])
                print(o_calc[:, ind])
                print(o_calc[:, ind] - o[:, ind])
            assert np.all(test)

    def test_cart_to_kep_circ(self):
        res = 100
        nu = np.linspace(0, 360.0, num=res, dtype=np.float64)
        x = np.empty((6, res), dtype=np.float64)

        v = kep.orbital_speed(self.a, self.a, self.mu)
        r = self.a

        o_ref = np.empty((6, res), dtype=np.float64)

        for i in range(res):
            o_ref[:5, i] = self.orb_init
            o_ref[5, i] = nu[i]

            x[0, i] = r*np.cos(np.radians(nu[i]))
            x[1, i] = r*np.sin(np.radians(nu[i]))
            x[2, i] = 0.0
            x[3, i] = -v*np.sin(np.radians(nu[i]))
            x[4, i] = v*np.cos(np.radians(nu[i]))
            x[5, i] = 0.0

        o = kep.cart_to_kep(x, mu=self.mu, degrees=True)

        nt.assert_array_almost_equal(o_ref[2:, :], o[2:, :], decimal=7)
        nt.assert_array_almost_equal(o_ref[1, :], o[1, :], decimal=4)
        nt.assert_array_almost_equal(
            o_ref[0, :]/self.R_e, o[0, :]/self.R_e, decimal=3)

    def test_kep_to_cart_circ(self):
        res = 100
        nu = np.linspace(0, 360.0, num=res, dtype=np.float64)

        o = np.empty((6, res), dtype=np.float64)
        for i in range(res):
            o[:5, i] = self.orb_init
            o[5, i] = nu[i]

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        r0 = np.ones((res,), dtype=np.float64)
        r0 *= self.a/self.R_e
        r = np.sqrt(np.sum(x[:3, :]**2, axis=0))/self.R_e

        nt.assert_array_almost_equal(r0, r, decimal=3)

        v0 = np.ones((res,), dtype=np.float64)
        v0 *= kep.orbital_speed(self.a, self.a, self.mu)*1e-3
        v = np.sqrt(np.sum(x[3:, :]**2, axis=0))*1e-3

        nt.assert_array_almost_equal(v0, v, decimal=3)

        # check perpendicular vel to rad

        dot = np.sum(x[3:, :]*x[:3, :], axis=0)
        nt.assert_array_almost_equal(dot, np.zeros(
            dot.shape, dtype=np.float64), decimal=3)

        # check 0 inc in xy
        nt.assert_array_almost_equal(x[2, :], np.zeros(
            (res,), dtype=np.float64), decimal=3)
        nt.assert_array_almost_equal(x[5, :], np.zeros(
            (res,), dtype=np.float64), decimal=3)

    def test_kep_to_cart_ecc(self):
        self.orb_init[1] = 0.8

        nu = np.degrees(
            np.array([0, np.pi*0.5, 1.5*np.pi, np.pi], dtype=np.float64))
        res = nu.size

        o = np.empty((6, res), dtype=np.float64)
        for i in range(res):
            o[:5, i] = self.orb_init
            o[5, i] = nu[i]

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        r = np.sqrt(np.sum(x[:3, :]**2, axis=0))
        v = np.sqrt(np.sum(x[3:, :]**2, axis=0))

        # test true anomaly = 0 is periapsis
        assert np.all(r[0] < r[1:]), 'nu=0 is not periapsis'
        assert np.all(r[3] > r[:3]), 'nu=pi is not apoapsis'

        assert x[0, 0] > 0.0, 'periapsis NOT to +x'
        assert x[0, 3] < 0.0, 'apoapsis NOT to +x'

        # test velocity, fast at peri slow at apo
        assert np.all(v[0] > v[1:]), 'periapsis vel not fastest'
        assert np.all(v[3] < v[:3]), 'apoapsis vel not slowest'

        # test velocity is perp at peri and apo
        nt.assert_almost_equal(x[3, 0]*1e-3, 0.0, decimal=3)
        nt.assert_almost_equal(x[3, 3]*1e-3, 0.0, decimal=3)

        # test ellipse oriented along x
        nt.assert_almost_equal(x[1, 0]/self.R_e, 0.0, decimal=3)
        nt.assert_almost_equal(x[1, 3]/self.R_e, 0.0, decimal=3)

        # test orbital motion counter-clockwise
        assert x[1, 1] > 0.0, 'orbital motion is not counter-clockwise'
        assert x[1, 2] < 0.0, 'orbital motion is not counter-clockwise'

        # test periapsis distance and apoapsis distance
        nt.assert_almost_equal(
            r[0]/self.R_e,
            (1.0 - self.orb_init[1])*self.orb_init[0]/self.R_e,
            decimal = 3,
        )
        nt.assert_almost_equal(
            r[3]/self.R_e,
            (1.0 + self.orb_init[1])*self.orb_init[0]/self.R_e,
            decimal = 3,
        )

    def test_kep_to_cart_inc(self):
        self.orb_init[2] = 45.0

        nu = np.degrees(np.array([np.pi*0.5], dtype=np.float64))
        res = nu.size

        o = np.empty((6, res), dtype=np.float64)
        for i in range(res):
            o[:5, i] = self.orb_init
            o[5, i] = nu[i]

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        assert x[2, 0] > 0.0 and x[1, 0] > 0.0, '+y inclines towards +z'

    def test_kep_to_cart_omega(self):
        self.orb_init[1] = 0.8
        self.orb_init[3] = 90.0

        nu = np.degrees(np.array([0.0], dtype=np.float64))
        res = nu.size

        o = np.empty((6, res), dtype=np.float64)
        for i in range(res):
            o[:5, i] = self.orb_init
            o[5, i] = nu[i]

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        # rotates counterclockwise in orbital plane
        r = np.sqrt(np.sum(x[:3, :]**2, axis=0))

        nt.assert_almost_equal(
            r[0]/self.R_e,
            (1.0 - self.orb_init[1])*self.orb_init[0]/self.R_e,
            decimal = 3,
        )

    def test_kep_to_cart_omega_inc(self):
        self.orb_init[1] = 0.8
        self.orb_init[2] = 90.0
        self.orb_init[3] = 90.0

        nu = np.degrees(np.array([0.0], dtype=np.float64))
        res = nu.size

        o = np.empty((6, res), dtype=np.float64)
        for i in range(res):
            o[:5, i] = self.orb_init
            o[5, i] = nu[i]

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        # rotates counterclockwise in orbital plane
        r = np.sqrt(np.sum(x[:3, :]**2, axis=0))

        nt.assert_almost_equal(
            r[0]/self.R_e,
            (1.0 - self.orb_init[1])*self.orb_init[0]/self.R_e,
            decimal = 3,
        )

        # check periapsis lies on +z axis
        nt.assert_almost_equal(
            x[2, 0]/self.R_e,
            np.linalg.norm(x[:, 0])/self.R_e,
            decimal=3,
        )

    def test_kep_to_cart_Omega_inc(self):
        self.orb_init[1] = 0.8
        self.orb_init[2] = 45.0
        self.orb_init[3] = 90.0
        self.orb_init[4] = 90.0

        nu = np.degrees(np.array([0.0], dtype=np.float64))
        res = nu.size

        o = np.empty((6, res), dtype=np.float64)
        for i in range(res):
            o[:5, i] = self.orb_init
            o[5, i] = nu[i]

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        # rotates counterclockwise in orbital plane
        r = np.sqrt(np.sum(x[:3, :]**2, axis=0))

        nt.assert_almost_equal(
            r[0]/self.R_e,
            (1.0 - self.orb_init[1])*self.orb_init[0]/self.R_e,
            decimal = 3,
        )

        # check quadrant of periapsis
        assert x[2, 0] > 0.0 and x[0, 0] < 0.0 and np.abs(
            x[1, 0]/self.R_e) < 1e-3


if __name__ == '__main__':
    unittest.main(verbosity=2)
