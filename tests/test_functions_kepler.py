#!/usr/bin/env python

'''

'''

import sys
import os
sys.path.insert(0, os.path.abspath('.'))
import time

import unittest
import numpy as np
import numpy.testing as nt

import dasst.functions.kepler as kep


class TestKeplerSolver(unittest.TestCase):

    def test_kepler_guess(self):
        E = n.linspace(0.0, 2.0*n.pi, num=100, dtype=n.float)
        e = n.linspace(0, 0.99, num=100, dtype=n.float)

        Ev, ev = n.meshgrid(E, e)

        Mv = dpt.eccentric2mean(Ev, ev)

        E0 = dpt.kepler_guess(Mv, ev)
        
        #the initial guess SHOULD be at least 30 degrees to true
        test_accuracy = n.abs(E0 - Ev)*180.0/n.pi < 30.0
        nt.assert_array_equal(
            test_accuracy,
            n.full(test_accuracy.shape, True)
        )

    def test_laguerre_solve_kepler(self):
        E = n.linspace(0.0, 2.0*n.pi, num=300, dtype=n.float)
        e = n.linspace(0, 0.99, num=500, dtype=n.float)

        for I, eit in enumerate(e):
            for J, Eit in enumerate(E):
                M = dpt.eccentric2mean(Eit, eit)

                E0 = dpt.kepler_guess(M, eit)
                E_calc, it = dpt.laguerre_solve_kepler(E0, M, eit, tol=1e-12)
                fun_err = n.abs(M - E_calc + eit*n.sin(E_calc))

                nt.assert_almost_equal(Eit, E_calc, decimal = 1e-9)
                assert fun_err < 1e-11





class TestAnomalies(unittest.TestCase):

    def test_true2eccentric_coencides(self):

        e_test = n.linspace(0,0.9,num=100)
        for e in e_test:
            E = dpt.true2eccentric(0.0, e)
            self.assertAlmostEqual(E, 0.0)

            E = dpt.true2eccentric(n.pi, e)
            self.assertAlmostEqual(E, n.pi)

    def test_eccentric2true_hand_calc(self):
        e = n.linspace(0.1,0.9,num=100)

        #For E = pi/2 hand calc and Pythagoras gives :math:`\\nu = \pi - tan^{-1}(\sqrt{e^{-2} - 1})`.
        E = n.ones(e.shape)*n.pi*0.5
        nu0 = n.pi - n.arctan(n.sqrt(1.0/e**2 - 1.0))
        nu = dpt.eccentric2true(E, e)
        nt.assert_array_almost_equal(nu, nu0, decimal=7)

    def test_true2eccentric_hand_calc(self):
        e = n.linspace(0.1,0.9,num=100)

        #For E = pi/2 hand calc and Pythagoras gives :math:`\\nu = \pi - tan^{-1}(\sqrt{e^{-2} - 1})`.
        nu = n.pi - n.arctan(n.sqrt(1.0/e**2 - 1.0))
        E0 = n.ones(e.shape)*n.pi*0.5
        E = dpt.true2eccentric(nu, e)
        nt.assert_array_almost_equal(E, E0, decimal=7)


    def test_eccentric2true_coencides(self):
        e_test = n.linspace(0,0.9,num=100)
        for e in e_test:
            nu = dpt.eccentric2true(0.0, e)
            self.assertAlmostEqual(nu, 0.0)

            nu = dpt.eccentric2true(n.pi, e)
            self.assertAlmostEqual(nu, n.pi)

    def test_eccentric_true_inverse(self):
        nu0 = 1.2345
        e = 0.5
        nu = dpt.eccentric2true(dpt.true2eccentric(nu0, e), e)
        self.assertAlmostEqual(nu0,nu)

    def test_eccentric2mean(self):
        E = n.linspace(0, 2.0*n.pi, num=100, dtype=n.float)
        M = dpt.eccentric2mean(E, 0.0)
        nt.assert_array_almost_equal(E, M, decimal=7)

        e = n.linspace(0, 0.9, num=100, dtype=n.float)
        M = dpt.eccentric2mean(n.pi, e)
        M0 = n.ones(e.shape, dtype=n.float)*n.pi
        nt.assert_array_almost_equal(M0, M, decimal=7)


    def test_mean2eccentric(self):
        M = n.linspace(0, 2.0*n.pi, num=100, dtype=n.float)
        E = dpt.mean2eccentric(M, 0.0)
        nt.assert_array_almost_equal(E, M, decimal=7)

        e = n.linspace(0, 0.9, num=100, dtype=n.float)
        E = dpt.mean2eccentric(n.pi, e)
        E0 = n.ones(e.shape, dtype=n.float)*n.pi
        nt.assert_array_almost_equal(E0, E, decimal=7)

    def test_mean_eccentric_inverse_array(self):
        M = n.linspace(0.0, 2.0*n.pi, num=100, dtype=n.float)
        e = n.linspace(0, 0.99, num=100, dtype=n.float)

        Mv, ev = n.meshgrid(M, e)

        E_test = dpt.mean2eccentric(Mv, ev)
        M_test = dpt.eccentric2mean(E_test, ev)
        nt.assert_array_almost_equal(M_test, Mv, decimal=7)

    def test_mean_eccentric_inverse_float(self):
        M = n.linspace(0.0, 2.0*n.pi, num=100, dtype=n.float)
        e = n.linspace(0, 0.99, num=100, dtype=n.float)

        for eit in e:
            for Mit in M:
                E_test = dpt.mean2eccentric(Mit, eit)
                M_test = dpt.eccentric2mean(E_test, eit)
                nt.assert_almost_equal(M_test, Mit, decimal=7)

    def test_true2mean(self):
        nu = n.linspace(0, 1.5*n.pi, num=100, dtype=n.float)
        M = dpt.true2mean(nu, 0.0)
        nt.assert_array_almost_equal(nu, M, decimal=7)

        e = n.linspace(0, 0.9, num=100, dtype=n.float)
        M = dpt.true2mean(n.pi, e)
        M0 = n.ones(e.shape, dtype=n.float)*n.pi
        nt.assert_array_almost_equal(M0, M, decimal=7)

    def test_mean2true(self):
        M = n.linspace(0, 1.5*n.pi, num=100, dtype=n.float)
        nu = dpt.mean2true(M, 0.0)
        nt.assert_array_almost_equal(nu, M, decimal=7)

        e = n.linspace(0, 0.9, num=100, dtype=n.float)
        nu = dpt.mean2true(n.pi, e)
        nu0 = n.ones(e.shape, dtype=n.float)*n.pi
        nt.assert_array_almost_equal(nu, nu0, decimal=7)

    def test_mean_true_inverse_float(self):
        M = n.linspace(0.0, 1.5*n.pi, num=100, dtype=n.float)
        e = n.linspace(0, 0.99, num=100, dtype=n.float)

        for eit in e:
            for Mit in M:
                nu_test = dpt.mean2true(Mit, eit)
                M_test = dpt.true2mean(nu_test, eit)
                nt.assert_almost_equal(M_test, Mit, decimal=7)

    def test_mean_true_inverse_array(self):
        M = n.linspace(0.0, 1.5*n.pi, num=100, dtype=n.float)
        e = n.linspace(0, 0.99, num=100, dtype=n.float)

        Mv, ev = n.meshgrid(M, e)

        nu_test = dpt.mean2true(Mv, ev)
        M_test = dpt.true2mean(nu_test, ev)
        nt.assert_array_almost_equal(M_test, Mv, decimal=7)




class TestOrbits(unittest.TestCase):

    def test_elliptic_radius(self):
        E = n.linspace(0, 2.0*n.pi, num=100, dtype=n.float)

        r = dpt.elliptic_radius(E, 1.0, 0.0)
        r0 = n.ones(E.shape, dtype=n.float)

        nt.assert_array_almost_equal(r, r0, decimal=7)

        #test periapsis and apoapsis
        e = n.linspace(0, 0.9, num=100, dtype=n.float)
        r = dpt.elliptic_radius(0.0, 1.0, e)
        nt.assert_array_almost_equal(r, 1.0 - e, decimal=7)

        r = dpt.elliptic_radius(n.pi, 1.0, e)
        nt.assert_array_almost_equal(r, 1.0 + e, decimal=7)

    def test_speed_hand_calc(self):
        v = dpt.orbital_speed(r=0.5, a=1.0, mu=1.0)
        self.assertAlmostEqual(v**2.0, 3.0)

        a_n = n.array([1.0, 1.0], dtype=n.float)
        r_n = n.array([0.5, 0.5], dtype=n.float)
        v_n = dpt.orbital_speed(r=r_n, a=a_n, mu=1.0)
        ref_n = n.array([3.0, 3.0], dtype=n.float)

        nt.assert_almost_equal(v_n**2.0, ref_n, decimal=7)

    def test_speed_numpy(self):
        x = n.linspace(1, 3, num=100, dtype=n.float)
        y = n.ones((100,10), dtype=n.float)

        v = dpt.orbital_speed(r=x, a=3.0, mu=1.0)
        self.assertEqual(v.shape, x.shape)

        v = dpt.orbital_speed(r=0.5, a=x, mu=1.0)
        self.assertEqual(v.shape, x.shape)

        v = dpt.orbital_speed(r=0.5, a=1.0, mu=x)
        self.assertEqual(v.shape, x.shape)

        v = dpt.orbital_speed(r=x, a=x, mu=1.0)
        self.assertEqual(v.shape, x.shape)

        v = dpt.orbital_speed(r=y, a=y, mu=1.0)
        self.assertEqual(v.shape, y.shape)

        v = dpt.orbital_speed(r=y, a=1.0, mu=1.0)
        self.assertEqual(v.shape, y.shape)

    def test_period_hand_calc(self):
        t = dpt.orbital_period(a=1.0, mu=1.0)
        self.assertAlmostEqual(t, 2.0*n.pi)

        a_n = n.array([1.0, 1.0], dtype=n.float)
        t_n = dpt.orbital_period(a=a_n, mu=1.0)
        ref_n = n.array([2.0*n.pi, 2.0*n.pi], dtype=n.float)

        nt.assert_almost_equal(t_n, ref_n, decimal=7)

    def test_period_numpy(self):
        x = n.linspace(1, 3, num=100, dtype=n.float)
        y = n.ones((100,10), dtype=n.float)

        t = dpt.orbital_period(a=x, mu=1.0)
        self.assertEqual(t.shape, x.shape)

        t = dpt.orbital_period(a=1.0, mu=x)
        self.assertEqual(t.shape, x.shape)

        t = dpt.orbital_period(a=x, mu=x)
        self.assertEqual(t.shape, x.shape)

        t = dpt.orbital_period(a=y, mu=1.0)
        self.assertEqual(t.shape, y.shape)

        t = dpt.orbital_period(a=y, mu=y)
        self.assertEqual(t.shape, y.shape)


class TestKepCart(unittest.TestCase):

    def setUp(self):
        self.M_e = 5.972e24
        self.R_e = 6371e3
        self.a = self.R_e*1.5
        self.orb_init = n.array([self.a, 0.0, 0.0, 0.0, 0.0], dtype=n.float)
        self.m = n.array([1.0])

    def test_cart_kep_loop(self):

        orb_init = n.array([self.a, 0.2, 23.0, 136.0, 44.0, 10.0], dtype=n.float)

        x = dpt.kep2cart(orb_init, m=self.m, M_cent=self.M_e, radians=False)
        x_ref = x.copy()
        for ind in range(200):
            o = dpt.cart2kep(x, m=self.m, M_cent=self.M_e, radians=False)
            x = dpt.kep2cart(o, m=self.m, M_cent=self.M_e, radians=False)

        nt.assert_array_almost_equal(x, x_ref, decimal = 6)


    def test_cart_kep_inverse(self):
        a = n.linspace(self.a, self.a + self.R_e, num=5, dtype=n.float64)
        e = n.linspace(0.0, 0.99, num=10, dtype=n.float64)
        inc = n.linspace(0.0, n.pi, num=10, dtype=n.float64)
        omega = n.linspace(0.0, 2.0*n.pi, num=10, dtype=n.float64)
        Omega = omega.copy()
        nu = omega.copy()
        av, ev, incv, omegav, Omegav, nuv = n.meshgrid(a, e, inc, omega, Omega, nu)

        ait = n.nditer(av)
        eit = n.nditer(ev)
        incit = n.nditer(incv)
        omegait = n.nditer(omegav)
        Omegait = n.nditer(Omegav)
        nuit = n.nditer(nuv)

        o = n.empty((6,av.size), dtype=n.float64)
        for ind in range(av.size):
            o[0,ind] = next(ait)
            o[1,ind] = next(eit)
            o[2,ind] = next(incit)
            o[3,ind] = next(omegait)
            o[4,ind] = next(Omegait)
            o[5,ind] = next(nuit)

        x = dpt.kep2cart(o, m=self.m, M_cent=self.M_e, radians=False)
        o_calc = dpt.cart2kep(x, m=self.m, M_cent=self.M_e, radians=False)
        x_calc = dpt.kep2cart(o_calc, m=self.m, M_cent=self.M_e, radians=False)

        for ind in range(av.size):
            try:
                nt.assert_array_almost_equal(x[:3,ind]/o[0,ind], x_calc[:3,ind]/o[0,ind], decimal = 6)
                nt.assert_array_almost_equal(x[3:,ind]/o[0,ind], x_calc[3:,ind]/o[0,ind], decimal = 6)
            except AssertionError:
                print(ind)
                print(x[:3,ind] - x_calc[:3,ind])
                print(o[:,ind])
                print(o_calc[:,ind])
                print(o_calc[:,ind] - o[:,ind])
                raise

        x = dpt.kep2cart(o, m=self.m, M_cent=self.M_e, radians=False)
        o_calc = dpt.cart2kep(x, m=self.m, M_cent=self.M_e, radians=False)

        o[3, o[1,:] < 1e-3] *= 0.0
        o[4, o[2,:] < 180e-3] *= 0.0

        for ind in range(av.size):
            test = n.isclose(o_calc[0,ind]/self.R_e, o[0,ind]/self.R_e, atol=1e-3)
            test = n.logical_and(test, n.isclose(o_calc[1,ind], o[1,ind], atol=1e-4))
            test_ang = n.logical_or(
                n.isclose(o_calc[2:,ind], o[2:,ind], atol=1e-7),
                n.isclose(n.mod(o_calc[2:,ind]+180.0, 360.0), n.mod(o[2:,ind]+180.0, 360.0), atol=1e-7),
            )
            test = n.logical_and(test, test_ang)
            if not n.all(test):
                print(ind)
                print(o[:,ind])
                print(o_calc[:,ind])
                print(o_calc[:,ind] - o[:,ind])
            assert n.all(test)
        



    def test_cart2kep_circ(self):
        res = 100
        nu = n.linspace(0, 360.0, num=res, dtype=n.float)
        x = n.empty((6,res),dtype=n.float)
        
        v = dpt.orbital_speed(self.a, self.a, scipy.constants.G*(self.m[0] + self.M_e))
        r = self.a

        o_ref = n.empty((6,res),dtype=n.float)
        
        for i in range(res):
            o_ref[:5,i] = self.orb_init
            o_ref[5,i] = nu[i]

            x[0,i] = r*n.cos(n.radians(nu[i]))
            x[1,i] = r*n.sin(n.radians(nu[i]))
            x[2,i] = 0.0
            x[3,i] = -v*n.sin(n.radians(nu[i]))
            x[4,i] = v*n.cos(n.radians(nu[i]))
            x[5,i] = 0.0

        o = dpt.cart2kep(x, m=self.m, M_cent=self.M_e, radians=False)

        nt.assert_array_almost_equal(o_ref[2:,:], o[2:,:], decimal=7)
        nt.assert_array_almost_equal(o_ref[1,:], o[1,:], decimal=4)
        nt.assert_array_almost_equal(o_ref[0,:]/self.R_e, o[0,:]/self.R_e, decimal=3)


    def test_kep2cart_circ(self):
        res = 100
        nu = n.linspace(0, 360.0, num=res, dtype=n.float)

        o = n.empty((6,res),dtype=n.float)
        for i in range(res):
            o[:5,i] = self.orb_init
            o[5,i] = nu[i]

        x = dpt.kep2cart(o, m=self.m, M_cent=self.M_e, radians=False)

        r0 = n.ones((res,), dtype=n.float)
        r0 *= self.a/self.R_e
        r = n.sqrt(n.sum(x[:3,:]**2, axis=0))/self.R_e

        nt.assert_array_almost_equal(r0, r, decimal=3)

        v0 = n.ones((res,), dtype=n.float)
        v0 *= dpt.orbital_speed(self.a, self.a, scipy.constants.G*(self.m[0] + self.M_e))*1e-3
        v = n.sqrt(n.sum(x[3:,:]**2, axis=0))*1e-3

        nt.assert_array_almost_equal(v0, v, decimal=3)

        #check perpendicular vel to rad

        dot = n.sum(x[3:,:]*x[:3,:], axis=0)
        nt.assert_array_almost_equal(dot, n.zeros(dot.shape, dtype=n.float), decimal=3)

        #check 0 inc in xy
        nt.assert_array_almost_equal(x[2,:], n.zeros((res,), dtype=n.float), decimal=3)
        nt.assert_array_almost_equal(x[5,:], n.zeros((res,), dtype=n.float), decimal=3)


    def test_kep2cart_ecc(self):
        self.orb_init[1] = 0.8
        
        nu = n.degrees(n.array([0, n.pi*0.5, 1.5*n.pi, n.pi], dtype=n.float))
        res = nu.size

        o = n.empty((6,res),dtype=n.float)
        for i in range(res):
            o[:5,i] = self.orb_init
            o[5,i] = nu[i]

        x = dpt.kep2cart(o, m=self.m, M_cent=self.M_e, radians=False)

        r = n.sqrt(n.sum(x[:3,:]**2, axis=0))
        v = n.sqrt(n.sum(x[3:,:]**2, axis=0))

        #test true anomaly = 0 is periapsis
        assert n.all(r[0] < r[1:]), 'nu=0 is not periapsis'
        assert n.all(r[3] > r[:3]), 'nu=pi is not apoapsis'

        assert x[0,0] > 0.0, 'periapsis NOT to +x'
        assert x[0,3] < 0.0, 'apoapsis NOT to +x'

        #test velocity, fast at peri slow at apo
        assert n.all(v[0] > v[1:]), 'periapsis vel not fastest'
        assert n.all(v[3] < v[:3]), 'apoapsis vel not slowest'
        
        #test velocity is perp at peri and apo
        nt.assert_almost_equal(x[3,0]*1e-3, 0.0, decimal=3)
        nt.assert_almost_equal(x[3,3]*1e-3, 0.0, decimal=3)

        #test ellipse oriented along x
        nt.assert_almost_equal(x[1,0]/self.R_e, 0.0, decimal=3)
        nt.assert_almost_equal(x[1,3]/self.R_e, 0.0, decimal=3)
        
        #test orbital motion counter-clockwise
        assert x[1,1] > 0.0, 'orbital motion is not counter-clockwise'
        assert x[1,2] < 0.0, 'orbital motion is not counter-clockwise'

        #test periapsis distance and apoapsis distance
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


    def test_kep2cart_inc(self):
        self.orb_init[2] = 45.0
        
        nu = n.degrees(n.array([n.pi*0.5], dtype=n.float))
        res = nu.size

        o = n.empty((6,res),dtype=n.float)
        for i in range(res):
            o[:5,i] = self.orb_init
            o[5,i] = nu[i]

        x = dpt.kep2cart(o, m=self.m, M_cent=self.M_e, radians=False)

        assert x[2,0] > 0.0 and x[1,0] > 0.0, '+y inclines towards +z'


    def test_kep2cart_omega(self):
        self.orb_init[1] = 0.8
        self.orb_init[3] = 90.0
        
        nu = n.degrees(n.array([0.0], dtype=n.float))
        res = nu.size

        o = n.empty((6,res),dtype=n.float)
        for i in range(res):
            o[:5,i] = self.orb_init
            o[5,i] = nu[i]

        x = dpt.kep2cart(o, m=self.m, M_cent=self.M_e, radians=False)

        #rotates counterclockwise in orbital plane
        r = n.sqrt(n.sum(x[:3,:]**2, axis=0))

        nt.assert_almost_equal(
            r[0]/self.R_e,
            (1.0 - self.orb_init[1])*self.orb_init[0]/self.R_e,
            decimal = 3,
        )


    def test_kep2cart_omega_inc(self):
        self.orb_init[1] = 0.8
        self.orb_init[2] = 90.0
        self.orb_init[3] = 90.0
        
        nu = n.degrees(n.array([0.0], dtype=n.float))
        res = nu.size

        o = n.empty((6,res),dtype=n.float)
        for i in range(res):
            o[:5,i] = self.orb_init
            o[5,i] = nu[i]

        x = dpt.kep2cart(o, m=self.m, M_cent=self.M_e, radians=False)

        #rotates counterclockwise in orbital plane
        r = n.sqrt(n.sum(x[:3,:]**2, axis=0))

        nt.assert_almost_equal(
            r[0]/self.R_e,
            (1.0 - self.orb_init[1])*self.orb_init[0]/self.R_e,
            decimal = 3,
        )

        #check periapsis lies on +z axis
        nt.assert_almost_equal(
            x[2,0]/self.R_e,
            n.linalg.norm(x[:,0])/self.R_e,
            decimal=3,
        )


    def test_kep2cart_Omega_inc(self):
        self.orb_init[1] = 0.8
        self.orb_init[2] = 45.0
        self.orb_init[3] = 90.0
        self.orb_init[4] = 90.0
        
        nu = n.degrees(n.array([0.0], dtype=n.float))
        res = nu.size

        o = n.empty((6,res),dtype=n.float)
        for i in range(res):
            o[:5,i] = self.orb_init
            o[5,i] = nu[i]

        x = dpt.kep2cart(o, m=self.m, M_cent=self.M_e, radians=False)

        #rotates counterclockwise in orbital plane
        r = n.sqrt(n.sum(x[:3,:]**2, axis=0))

        nt.assert_almost_equal(
            r[0]/self.R_e,
            (1.0 - self.orb_init[1])*self.orb_init[0]/self.R_e,
            decimal = 3,
        )

        #check quadrant of periapsis
        assert x[2,0] > 0.0 and x[0,0] < 0.0 and n.abs(x[1,0]/self.R_e) < 1e-3



if __name__ == '__main__':
    unittest.main(verbosity=2)