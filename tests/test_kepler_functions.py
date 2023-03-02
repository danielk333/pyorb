#!/usr/bin/env python

'''
Test basic kepler functions
'''

import unittest
import numpy as np
import numpy.testing as nt

import pyorb.kepler as kep


class TestAnomaliesHyperbolic(unittest.TestCase):

    def test_true_to_hyperbolic(self):
        nu = kep.true_to_eccentric(0.0, np.array([0.1, 1.2]))
        nt.assert_array_almost_equal(nu, np.zeros_like(nu))

        nu = kep.true_to_eccentric(0.0, 1.2)
        nt.assert_almost_equal(nu, 0.0)

    def test_true_to_parabolic(self):
        nu = kep.true_to_eccentric(0.0, np.array([0.1, 1]))
        nt.assert_array_almost_equal(nu, np.zeros_like(nu))

        nu = kep.true_to_eccentric(0.0, 1)
        nt.assert_almost_equal(nu, 0.0)

    def test_hyperbolic_to_true(self):
        nu = kep.eccentric_to_true(0.0, np.array([0.1, 1.2]))
        nt.assert_array_almost_equal(nu, np.zeros_like(nu))

        nu = kep.eccentric_to_true(0.0, 1.2)
        nt.assert_almost_equal(nu, 0.0)

    def test_parabolic_to_true(self):
        nu = kep.eccentric_to_true(0.0, np.array([0.1, 1]))
        nt.assert_array_almost_equal(nu, np.zeros_like(nu))

        nu = kep.eccentric_to_true(0.0, 1)
        nt.assert_almost_equal(nu, 0.0)

    def test_mean_to_parabolic(self):
        E = kep.mean_to_eccentric(0.0, np.array([0.1, 1.0]))
        nt.assert_array_almost_equal(E, np.zeros_like(E))

        E = kep.mean_to_eccentric(0.0, 1.2)
        nt.assert_almost_equal(E, 0.0)

    def test_mean_to_hyperbolic(self):
        E = kep.mean_to_eccentric(0.0, np.array([0.1, 1.2]))
        nt.assert_array_almost_equal(E, np.zeros_like(E))

        E = kep.mean_to_eccentric(0.0, 1.2)
        nt.assert_almost_equal(E, 0.0)

    def test_parabolic_to_mean(self):
        M = kep.eccentric_to_mean(0.0, np.array([0.1, 1.0]))
        nt.assert_array_almost_equal(M, np.zeros_like(M))

        M = kep.eccentric_to_mean(0.0, 1)
        nt.assert_almost_equal(M, 0.0)

    def test_hyperbolic_to_mean(self):
        M = kep.eccentric_to_mean(0.0, np.array([0.1, 1.2]))
        nt.assert_array_almost_equal(M, np.zeros_like(M))

        M = kep.eccentric_to_mean(0.0, 1.2)
        nt.assert_almost_equal(M, 0.0)

    def test_parabolic_radius(self):
        nu = np.linspace(np.pi/2, 0, num=100)
        # Semi-latus rectum p = 2q

        r = kep.parabolic_radius(nu, 1.0)
        nt.assert_almost_equal(r[0], 2)
        nt.assert_almost_equal(r[-1], 1)

        r = kep.parabolic_radius(-nu, 1.0)
        nt.assert_almost_equal(r[0], 2)
        nt.assert_almost_equal(r[-1], 1)
        with np.errstate(divide='ignore'):
            r = kep.parabolic_radius(np.pi, 1.0)
            assert np.isinf(r)
            r = kep.parabolic_radius(-np.pi, 1.0)
            assert np.isinf(r)

    def test_hyperbolic_to_true_inverse(self):
        E0 = np.linspace(0, np.pi/2, num=100)
        e = np.ones_like(E0)*1.2

        nu0 = kep.eccentric_to_true(E0, e)
        E = kep.true_to_eccentric(nu0, e)
        nu = kep.eccentric_to_true(E, e)

        nt.assert_array_almost_equal(E0, E)
        nt.assert_array_almost_equal(nu0, nu)

    def test_parabolic_to_true_inverse(self):
        E0 = np.linspace(0, np.pi/2, num=100)
        e = np.ones_like(E0)

        nu0 = kep.eccentric_to_true(E0, e)
        E = kep.true_to_eccentric(nu0, e)
        nu = kep.eccentric_to_true(E, e)

        nt.assert_array_almost_equal(E0, E)
        nt.assert_array_almost_equal(nu0, nu)

    def test_parabolic_to_mean_inverse(self):
        E0 = np.linspace(0, np.pi/2, num=100)
        e = np.ones_like(E0)

        M0 = kep.eccentric_to_mean(E0, e)
        E = kep.mean_to_eccentric(M0, e)
        M = kep.eccentric_to_mean(E, e)

        nt.assert_array_almost_equal(E0, E)
        nt.assert_array_almost_equal(M0, M)

    def test_hyperbolic_to_mean_inverse(self):
        E0 = np.linspace(0, np.pi/2, num=100)
        e = np.ones_like(E0)*1.2

        M0 = kep.eccentric_to_mean(E0, e)
        E = kep.mean_to_eccentric(M0, e)
        M = kep.eccentric_to_mean(E, e)

        nt.assert_array_almost_equal(E0, E)
        nt.assert_array_almost_equal(M0, M)


class TestAnomalies(unittest.TestCase):

    def test_true_to_eccentric_coencides(self):

        e_test = np.linspace(0, 0.9, num=100)
        for e in e_test:
            E = kep.true_to_eccentric(0.0, e)
            self.assertAlmostEqual(E, 0.0)

            E = kep.true_to_eccentric(np.pi, e)
            self.assertAlmostEqual(E, np.pi)

    def test_eccentric_to_true_hand_calc(self):
        e = np.linspace(0.1, 0.9, num=100)

        # For E = pi/2 hand calc and Pythagoras gives 
        # :math:`\\nu = \pi - tan^{-1}(\sqrt{e^{-2} - 1})`.
        E = np.ones(e.shape)*np.pi*0.5
        nu0 = np.pi - np.arctan(np.sqrt(1.0/e**2 - 1.0))
        nu = kep.eccentric_to_true(E, e)
        nt.assert_array_almost_equal(nu, nu0, decimal=7)

    def test_true_to_eccentric_hand_calc(self):
        e = np.linspace(0.1, 0.9, num=100)

        # For E = pi/2 hand calc and Pythagoras gives 
        # :math:`\\nu = \pi - tan^{-1}(\sqrt{e^{-2} - 1})`.
        nu = np.pi - np.arctan(np.sqrt(1.0/e**2 - 1.0))
        E0 = np.ones(e.shape)*np.pi*0.5
        E = kep.true_to_eccentric(nu, e)
        nt.assert_array_almost_equal(E, E0, decimal=7)

    def test_eccentric_to_true_coencides(self):
        e_test = np.linspace(0, 0.9, num=100)
        for e in e_test:
            nu = kep.eccentric_to_true(0.0, e)
            self.assertAlmostEqual(nu, 0.0)

            nu = kep.eccentric_to_true(np.pi, e)
            self.assertAlmostEqual(nu, np.pi)

    def test_eccentric_true_inverse(self):
        nu0 = 1.2345
        e = 0.5
        nu = kep.eccentric_to_true(kep.true_to_eccentric(nu0, e), e)
        self.assertAlmostEqual(nu0, nu)

    def test_eccentric_to_mean(self):
        E = np.linspace(0, 2.0*np.pi, num=100, dtype=np.float64)
        M = kep.eccentric_to_mean(E, 0.0)
        nt.assert_array_almost_equal(E, M, decimal=7)

        e = np.linspace(0, 0.9, num=100, dtype=np.float64)
        M = kep.eccentric_to_mean(np.pi, e)
        M0 = np.ones(e.shape, dtype=np.float64)*np.pi
        nt.assert_array_almost_equal(M0, M, decimal=7)

    def test_mean_to_eccentric(self):
        M = np.linspace(0, 2.0*np.pi, num=100, dtype=np.float64)
        E = kep.mean_to_eccentric(M, 0.0)
        nt.assert_array_almost_equal(E, M, decimal=7)

        e = np.linspace(0, 0.9, num=100, dtype=np.float64)
        E = kep.mean_to_eccentric(np.pi, e)
        E0 = np.ones(e.shape, dtype=np.float64)*np.pi
        nt.assert_array_almost_equal(E0, E, decimal=7)

    def test_mean_eccentric_inverse_array(self):
        M = np.linspace(0.0, 2.0*np.pi, num=100, dtype=np.float64)
        e = np.linspace(0, 0.99, num=100, dtype=np.float64)

        Mv, ev = np.meshgrid(M, e)

        E_test = kep.mean_to_eccentric(Mv, ev)
        M_test = kep.eccentric_to_mean(E_test, ev)
        nt.assert_array_almost_equal(M_test, Mv, decimal=7)

    def test_mean_eccentric_inverse_float(self):
        M = np.linspace(0.0, 2.0*np.pi, num=100, dtype=np.float64)
        e = np.linspace(0, 0.99, num=100, dtype=np.float64)

        for eit in e:
            for Mit in M:
                E_test = kep.mean_to_eccentric(Mit, eit)
                M_test = kep.eccentric_to_mean(E_test, eit)
                nt.assert_almost_equal(M_test, Mit, decimal=7)

    def test_true_to_mean(self):
        nu = np.linspace(0, 1.5*np.pi, num=100, dtype=np.float64)
        M = kep.true_to_mean(nu, 0.0)
        nt.assert_array_almost_equal(nu, M, decimal=7)

        e = np.linspace(0, 0.9, num=100, dtype=np.float64)
        M = kep.true_to_mean(np.pi, e)
        M0 = np.ones(e.shape, dtype=np.float64)*np.pi
        nt.assert_array_almost_equal(M0, M, decimal=7)

    def test_mean_to_true(self):
        M = np.linspace(0, 1.5*np.pi, num=100, dtype=np.float64)
        nu = kep.mean_to_true(M, 0.0)
        nt.assert_array_almost_equal(nu, M, decimal=7)

        e = np.linspace(0, 0.9, num=100, dtype=np.float64)
        nu = kep.mean_to_true(np.pi, e)
        nu0 = np.ones(e.shape, dtype=np.float64)*np.pi
        nt.assert_array_almost_equal(nu, nu0, decimal=7)

    def test_mean_true_inverse_float(self):
        M = np.linspace(0.0, 1.5*np.pi, num=100, dtype=np.float64)
        e = np.linspace(0, 0.99, num=100, dtype=np.float64)

        for eit in e:
            for Mit in M:
                nu_test = kep.mean_to_true(Mit, eit)
                M_test = kep.true_to_mean(nu_test, eit)
                nt.assert_almost_equal(M_test, Mit, decimal=7)

    def test_mean_true_inverse_array(self):
        M = np.linspace(0.0, 1.5*np.pi, num=100, dtype=np.float64)
        e = np.linspace(0, 0.99, num=100, dtype=np.float64)

        Mv, ev = np.meshgrid(M, e)

        nu_test = kep.mean_to_true(Mv, ev)
        M_test = kep.true_to_mean(nu_test, ev)
        nt.assert_array_almost_equal(M_test, Mv, decimal=7)


class TestEllipticOrbit(unittest.TestCase):

    def test_elliptic_radius(self):
        E = np.linspace(0, 2.0*np.pi, num=100, dtype=np.float64)

        r = kep.elliptic_radius(E, 1.0, 0.0)
        r0 = np.ones(E.shape, dtype=np.float64)

        nt.assert_array_almost_equal(r, r0, decimal=7)

        # test periapsis and apoapsis
        e = np.linspace(0, 0.9, num=100, dtype=np.float64)
        r = kep.elliptic_radius(0.0, 1.0, e)
        nt.assert_array_almost_equal(r, 1.0 - e, decimal=7)

        r = kep.elliptic_radius(np.pi, 1.0, e)
        nt.assert_array_almost_equal(r, 1.0 + e, decimal=7)

    def test_speed_hand_calc(self):
        v = kep.orbital_speed(r=0.5, a=1.0, mu=1.0)
        self.assertAlmostEqual(v**2.0, 3.0)

        a_n = np.array([1.0, 1.0], dtype=np.float64)
        r_n = np.array([0.5, 0.5], dtype=np.float64)
        v_n = kep.orbital_speed(r=r_n, a=a_n, mu=1.0)
        ref_n = np.array([3.0, 3.0], dtype=np.float64)

        nt.assert_almost_equal(v_n**2.0, ref_n, decimal=7)

    def test_speed_numpy(self):
        x = np.linspace(1, 3, num=100, dtype=np.float64)
        y = np.ones((100, 10), dtype=np.float64)

        v = kep.orbital_speed(r=x, a=3.0, mu=1.0)
        self.assertEqual(v.shape, x.shape)

        v = kep.orbital_speed(r=0.5, a=x, mu=1.0)
        self.assertEqual(v.shape, x.shape)

        v = kep.orbital_speed(r=0.5, a=1.0, mu=x)
        self.assertEqual(v.shape, x.shape)

        v = kep.orbital_speed(r=x, a=x, mu=1.0)
        self.assertEqual(v.shape, x.shape)

        v = kep.orbital_speed(r=y, a=y, mu=1.0)
        self.assertEqual(v.shape, y.shape)

        v = kep.orbital_speed(r=y, a=1.0, mu=1.0)
        self.assertEqual(v.shape, y.shape)

    def test_period_hand_calc(self):
        t = kep.orbital_period(a=1.0, mu=1.0)
        self.assertAlmostEqual(t, 2.0*np.pi)

        a_n = np.array([1.0, 1.0], dtype=np.float64)
        t_n = kep.orbital_period(a=a_n, mu=1.0)
        ref_n = np.array([2.0*np.pi, 2.0*np.pi], dtype=np.float64)

        nt.assert_almost_equal(t_n, ref_n, decimal=7)

    def test_period_numpy(self):
        x = np.linspace(1, 3, num=100, dtype=np.float64)
        y = np.ones((100, 10), dtype=np.float64)

        t = kep.orbital_period(a=x, mu=1.0)
        self.assertEqual(t.shape, x.shape)

        t = kep.orbital_period(a=1.0, mu=x)
        self.assertEqual(t.shape, x.shape)

        t = kep.orbital_period(a=x, mu=x)
        self.assertEqual(t.shape, x.shape)

        t = kep.orbital_period(a=y, mu=1.0)
        self.assertEqual(t.shape, y.shape)

        t = kep.orbital_period(a=y, mu=y)
        self.assertEqual(t.shape, y.shape)


if __name__ == '__main__':
    unittest.main(verbosity=2)
