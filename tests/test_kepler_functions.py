#!/usr/bin/env python

'''
Test basic kepler functions
'''

import time

import unittest
import numpy as np
import numpy.testing as nt

import pyorb.kepler as kep
import pyorb

@unittest.skip("Not yet implemented")
class TestAlternativeParameters(unittest.TestCase):

    def setUp(self):
        pass

    def test_cart_to_equi(self):
        assert False

    def test_equi_to_cart(self):
        assert False

    def test_equi_cart_consistency(self):
        assert False

    def test_kep_to_equi(self):
        assert False

    def test_equi_to_kep(self):
        assert False

    def test_equi_kep_consistency(self):
        assert False


class TestRotations(unittest.TestCase):

    def setUp(self):
        mat = np.eye(3, dtype=np.float64)
        self.x = mat[:,0]
        self.y = mat[:,1]
        self.z = mat[:,2]
        self.theta = np.pi/2


    def test_rot_mat_x(self):
        R = kep.rot_mat_x(self.theta, dtype=np.float32)
        assert R.dtype==np.float32
        
        R = kep.rot_mat_x(self.theta)
        assert R.dtype==np.float64

        nt.assert_almost_equal(self.z, R @ self.y)


    def test_rot_mat_y(self):
        R = kep.rot_mat_y(self.theta, dtype=np.float32)
        assert R.dtype==np.float32
        
        R = kep.rot_mat_y(self.theta)
        assert R.dtype==np.float64

        nt.assert_almost_equal(-self.z, R @ self.x)


    def test_rot_mat_z(self):
        R = kep.rot_mat_z(self.theta, dtype=np.float32)
        assert R.dtype==np.float32
        
        R = kep.rot_mat_z(self.theta)
        assert R.dtype==np.float64

        nt.assert_almost_equal(self.y, R @ self.x)



class TestKeplerSolver(unittest.TestCase):

    def test_kepler_guess(self):
        E = np.linspace(0.0, 2.0*np.pi, num=100, dtype=np.float64)
        e = np.linspace(0, 0.99, num=100, dtype=np.float64)

        Ev, ev = np.meshgrid(E, e)

        Mv = kep.eccentric_to_mean(Ev, ev)

        E0 = kep.kepler_guess(Mv, ev)
        
        #the initial guess SHOULD be at least 30 degrees to true
        test_accuracy = np.abs(E0 - Ev)*180.0/np.pi < 30.0
        nt.assert_array_equal(
            test_accuracy,
            np.full(test_accuracy.shape, True)
        )

    def test_laguerre_solve_kepler(self):
        E = np.linspace(0.0, 2.0*np.pi, num=300, dtype=np.float64)
        e = np.linspace(0, 0.99, num=500, dtype=np.float64)

        for I, eit in enumerate(e):
            for J, Eit in enumerate(E):
                M = kep.eccentric_to_mean(Eit, eit)

                E0 = kep.kepler_guess(M, eit)
                E_calc, it = kep.laguerre_solve_kepler(E0, M, eit, tol=1e-12)
                fun_err = np.abs(M - E_calc + eit*np.sin(E_calc))

                nt.assert_almost_equal(Eit, E_calc, decimal = 1e-9)
                assert fun_err < 1e-11

    @unittest.skip("Not yet implemented")
    def test_laguerre_solve_hyperbolic_kepler(self):
        assert False

    @unittest.skip("Not yet implemented")
    def test_hyperbolic_kepler_guess(self):
        assert False

@unittest.skip("Not yet implemented")
class TestAnomaliesHyperbolic(unittest.TestCase):

    def test_true_to_hyperbolic(self):
        assert False

    def test_hyperbolic_to_true(self):
        assert False

    def test_mean_to_hyperbolic(self):
        assert False

    def test_hyperbolic_to_mean(self):
        assert False

    def test_true_to_hyperbolic_coencides(self):
        assert False
        e_test = [1, 2]
        for e in e_test:
            E = kep.true_to_eccentric(0.0, e)
            self.assertAlmostEqual(E, 0.0)


class TestAnomalies(unittest.TestCase):

    def test_true_to_eccentric_coencides(self):

        e_test = np.linspace(0,0.9,num=100)
        for e in e_test:
            E = kep.true_to_eccentric(0.0, e)
            self.assertAlmostEqual(E, 0.0)

            E = kep.true_to_eccentric(np.pi, e)
            self.assertAlmostEqual(E, np.pi)

    def test_eccentric_to_true_hand_calc(self):
        e = np.linspace(0.1,0.9,num=100)

        #For E = pi/2 hand calc and Pythagoras gives :math:`\\nu = \pi - tan^{-1}(\sqrt{e^{-2} - 1})`.
        E = np.ones(e.shape)*np.pi*0.5
        nu0 = np.pi - np.arctan(np.sqrt(1.0/e**2 - 1.0))
        nu = kep.eccentric_to_true(E, e)
        nt.assert_array_almost_equal(nu, nu0, decimal=7)

    def test_true_to_eccentric_hand_calc(self):
        e = np.linspace(0.1,0.9,num=100)

        #For E = pi/2 hand calc and Pythagoras gives :math:`\\nu = \pi - tan^{-1}(\sqrt{e^{-2} - 1})`.
        nu = np.pi - np.arctan(np.sqrt(1.0/e**2 - 1.0))
        E0 = np.ones(e.shape)*np.pi*0.5
        E = kep.true_to_eccentric(nu, e)
        nt.assert_array_almost_equal(E, E0, decimal=7)


    def test_eccentric_to_true_coencides(self):
        e_test = np.linspace(0,0.9,num=100)
        for e in e_test:
            nu = kep.eccentric_to_true(0.0, e)
            self.assertAlmostEqual(nu, 0.0)

            nu = kep.eccentric_to_true(np.pi, e)
            self.assertAlmostEqual(nu, np.pi)

    def test_eccentric_true_inverse(self):
        nu0 = 1.2345
        e = 0.5
        nu = kep.eccentric_to_true(kep.true_to_eccentric(nu0, e), e)
        self.assertAlmostEqual(nu0,nu)

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




class TestOrbits(unittest.TestCase):

    def test_elliptic_radius(self):
        E = np.linspace(0, 2.0*np.pi, num=100, dtype=np.float64)

        r = kep.elliptic_radius(E, 1.0, 0.0)
        r0 = np.ones(E.shape, dtype=np.float64)

        nt.assert_array_almost_equal(r, r0, decimal=7)

        #test periapsis and apoapsis
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
        y = np.ones((100,10), dtype=np.float64)

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
        y = np.ones((100,10), dtype=np.float64)

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


class TestKepCart(unittest.TestCase):

    def setUp(self):
        self.M_e = 5.972e24
        self.R_e = 6371e3
        self.a = self.R_e*1.5
        self.orb_init = np.array([self.a, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
        self.m = 1.0
        self.mu = (pyorb.M_earth + self.m)*pyorb.G

    def test_cart_kep_loop_degradation(self):

        orb_init = np.array([self.a, 0.2, 23.0, 136.0, 44.0, 10.0], dtype=np.float64)

        x = kep.kep_to_cart(orb_init, mu=self.mu, degrees=True)
        x_ref = x.copy()
        for ind in range(200):
            o = kep.cart_to_kep(x, mu=self.mu, degrees=True)
            x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        nt.assert_array_almost_equal(x, x_ref, decimal = 6)



    def test_orientation_convention(self):

        e = 0.5
        orb0 = np.array([self.a, e, 0, 0, 0, 0], dtype=np.float64)
        q = self.a*(1 - e)
        Q = self.a*(1 + e)
        l = self.a*(1 - e**2)

        orb = orb0.copy()
        x = kep.kep_to_cart(orb, mu=self.mu, degrees=True)
        nt.assert_almost_equal(x[0], q, decimal = 6)

        orb = orb0.copy()
        orb[5] = 90
        x = kep.kep_to_cart(orb, mu=self.mu, degrees=True)
        nt.assert_almost_equal(x[1], l, decimal = 6)

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
        '''Longitude of ascending node gets transfered to argument of periapsis in planar orbits. Test if this is true for retro-grade and pro-grade planar orbits.
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

        #should be same but reversed velocity
        x1[3:] = -x1[3:]
        nt.assert_almost_equal(x0, x1)

        nt.assert_almost_equal(orb_ref, orb_out)


    def test_cart_kep_inverse(self):
        a = np.linspace(self.a, self.a + self.R_e, num=2, dtype=np.float64)
        e = np.linspace(0.0, 0.99, num=10, dtype=np.float64)
        inc = np.linspace(0.0, 180.0, num=10, dtype=np.float64)
        omega = np.linspace(0.0, 360.0, num=10, dtype=np.float64)
        Omega = omega.copy()
        nu = np.array([0,45,180], dtype=np.float64)
        av, ev, incv, omegav, Omegav, nuv = np.meshgrid(a, e, inc, omega, Omega, nu)

        ait = np.nditer(av)
        eit = np.nditer(ev)
        incit = np.nditer(incv)
        omegait = np.nditer(omegav)
        Omegait = np.nditer(Omegav)
        nuit = np.nditer(nuv)

        o = np.empty((6,av.size), dtype=np.float64)
        for ind in range(av.size):
            o[0,ind] = next(ait)
            o[1,ind] = next(eit)
            o[2,ind] = next(incit)
            o[3,ind] = next(omegait)
            o[4,ind] = next(Omegait)
            o[5,ind] = next(nuit)

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)
        o_calc = kep.cart_to_kep(x, mu=self.mu, degrees=True)
        x_calc = kep.kep_to_cart(o_calc, mu=self.mu, degrees=True)

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

        el = o[1,:] < kep.e_lim
        
        il = o[2,:] < np.degrees(kep.i_lim)
        il_ret = o[2,:] > 180 - np.degrees(kep.i_lim)
        o_orig = o.copy()

        o[3, el] = 0.0
        o[4, il] = 0.0
        o[4, il_ret] = 0.0
        o[5, el] += o_orig[3, el]

        o[5, np.logical_and(il,el)] += o_orig[4, np.logical_and(il,el)]
        o[3, np.logical_and(il,np.logical_not(el))] += o_orig[4, np.logical_and(il,np.logical_not(el))]

        o[5, np.logical_and(il_ret,el)] -= o_orig[4, np.logical_and(il_ret,el)]
        o[3, np.logical_and(il_ret,np.logical_not(el))] -= o_orig[4, np.logical_and(il_ret,np.logical_not(el))]

        o[5,:] = np.mod(o[5,:] + 360.0, 360.0)

        for ind in range(av.size):
            test = np.isclose(o_calc[0,ind]/self.R_e, o[0,ind]/self.R_e, atol=1e-3)
            test = np.logical_and(test, np.isclose(o_calc[1,ind], o[1,ind], atol=1e-4))
            test_ang = np.logical_or(
                np.isclose(o_calc[2:,ind], o[2:,ind], atol=1e-7),
                np.isclose(np.mod(o_calc[2:,ind]+180.0, 360.0), np.mod(o[2:,ind]+180.0, 360.0), atol=1e-7),
            )
            test = np.logical_and(test, test_ang)
            if not np.all(test):

                # x_ = kep.kep_to_cart(o_orig[:,ind], mu=self.mu, degrees=True)
                # o_calc_ = kep.cart_to_kep(x_, mu=self.mu, degrees=True)
                # x_calc_ = kep.kep_to_cart(o_calc_, mu=self.mu, degrees=True)

                print(ind)
                print(o_orig[:,ind])
                print(o[:,ind])
                print(o_calc[:,ind])
                print(o_calc[:,ind] - o[:,ind])
            assert np.all(test)
        


    def test_cart_to_kep_circ(self):
        res = 100
        nu = np.linspace(0, 360.0, num=res, dtype=np.float64)
        x = np.empty((6,res),dtype=np.float64)
        
        v = kep.orbital_speed(self.a, self.a, self.mu)
        r = self.a

        o_ref = np.empty((6,res),dtype=np.float64)
        
        for i in range(res):
            o_ref[:5,i] = self.orb_init
            o_ref[5,i] = nu[i]

            x[0,i] = r*np.cos(np.radians(nu[i]))
            x[1,i] = r*np.sin(np.radians(nu[i]))
            x[2,i] = 0.0
            x[3,i] = -v*np.sin(np.radians(nu[i]))
            x[4,i] = v*np.cos(np.radians(nu[i]))
            x[5,i] = 0.0

        o = kep.cart_to_kep(x, mu=self.mu, degrees=True)

        nt.assert_array_almost_equal(o_ref[2:,:], o[2:,:], decimal=7)
        nt.assert_array_almost_equal(o_ref[1,:], o[1,:], decimal=4)
        nt.assert_array_almost_equal(o_ref[0,:]/self.R_e, o[0,:]/self.R_e, decimal=3)


    def test_kep_to_cart_circ(self):
        res = 100
        nu = np.linspace(0, 360.0, num=res, dtype=np.float64)

        o = np.empty((6,res),dtype=np.float64)
        for i in range(res):
            o[:5,i] = self.orb_init
            o[5,i] = nu[i]

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        r0 = np.ones((res,), dtype=np.float64)
        r0 *= self.a/self.R_e
        r = np.sqrt(np.sum(x[:3,:]**2, axis=0))/self.R_e

        nt.assert_array_almost_equal(r0, r, decimal=3)

        v0 = np.ones((res,), dtype=np.float64)
        v0 *= kep.orbital_speed(self.a, self.a, self.mu)*1e-3
        v = np.sqrt(np.sum(x[3:,:]**2, axis=0))*1e-3

        nt.assert_array_almost_equal(v0, v, decimal=3)

        #check perpendicular vel to rad

        dot = np.sum(x[3:,:]*x[:3,:], axis=0)
        nt.assert_array_almost_equal(dot, np.zeros(dot.shape, dtype=np.float64), decimal=3)

        #check 0 inc in xy
        nt.assert_array_almost_equal(x[2,:], np.zeros((res,), dtype=np.float64), decimal=3)
        nt.assert_array_almost_equal(x[5,:], np.zeros((res,), dtype=np.float64), decimal=3)


    def test_kep_to_cart_ecc(self):
        self.orb_init[1] = 0.8
        
        nu = np.degrees(np.array([0, np.pi*0.5, 1.5*np.pi, np.pi], dtype=np.float64))
        res = nu.size

        o = np.empty((6,res),dtype=np.float64)
        for i in range(res):
            o[:5,i] = self.orb_init
            o[5,i] = nu[i]

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        r = np.sqrt(np.sum(x[:3,:]**2, axis=0))
        v = np.sqrt(np.sum(x[3:,:]**2, axis=0))

        #test true anomaly = 0 is periapsis
        assert np.all(r[0] < r[1:]), 'nu=0 is not periapsis'
        assert np.all(r[3] > r[:3]), 'nu=pi is not apoapsis'

        assert x[0,0] > 0.0, 'periapsis NOT to +x'
        assert x[0,3] < 0.0, 'apoapsis NOT to +x'

        #test velocity, fast at peri slow at apo
        assert np.all(v[0] > v[1:]), 'periapsis vel not fastest'
        assert np.all(v[3] < v[:3]), 'apoapsis vel not slowest'
        
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


    def test_kep_to_cart_inc(self):
        self.orb_init[2] = 45.0
        
        nu = np.degrees(np.array([np.pi*0.5], dtype=np.float64))
        res = nu.size

        o = np.empty((6,res),dtype=np.float64)
        for i in range(res):
            o[:5,i] = self.orb_init
            o[5,i] = nu[i]

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        assert x[2,0] > 0.0 and x[1,0] > 0.0, '+y inclines towards +z'


    def test_kep_to_cart_omega(self):
        self.orb_init[1] = 0.8
        self.orb_init[3] = 90.0
        
        nu = np.degrees(np.array([0.0], dtype=np.float64))
        res = nu.size

        o = np.empty((6,res),dtype=np.float64)
        for i in range(res):
            o[:5,i] = self.orb_init
            o[5,i] = nu[i]

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        #rotates counterclockwise in orbital plane
        r = np.sqrt(np.sum(x[:3,:]**2, axis=0))

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

        o = np.empty((6,res),dtype=np.float64)
        for i in range(res):
            o[:5,i] = self.orb_init
            o[5,i] = nu[i]

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        #rotates counterclockwise in orbital plane
        r = np.sqrt(np.sum(x[:3,:]**2, axis=0))

        nt.assert_almost_equal(
            r[0]/self.R_e,
            (1.0 - self.orb_init[1])*self.orb_init[0]/self.R_e,
            decimal = 3,
        )

        #check periapsis lies on +z axis
        nt.assert_almost_equal(
            x[2,0]/self.R_e,
            np.linalg.norm(x[:,0])/self.R_e,
            decimal=3,
        )


    def test_kep_to_cart_Omega_inc(self):
        self.orb_init[1] = 0.8
        self.orb_init[2] = 45.0
        self.orb_init[3] = 90.0
        self.orb_init[4] = 90.0
        
        nu = np.degrees(np.array([0.0], dtype=np.float64))
        res = nu.size

        o = np.empty((6,res),dtype=np.float64)
        for i in range(res):
            o[:5,i] = self.orb_init
            o[5,i] = nu[i]

        x = kep.kep_to_cart(o, mu=self.mu, degrees=True)

        #rotates counterclockwise in orbital plane
        r = np.sqrt(np.sum(x[:3,:]**2, axis=0))

        nt.assert_almost_equal(
            r[0]/self.R_e,
            (1.0 - self.orb_init[1])*self.orb_init[0]/self.R_e,
            decimal = 3,
        )

        #check quadrant of periapsis
        assert x[2,0] > 0.0 and x[0,0] < 0.0 and np.abs(x[1,0]/self.R_e) < 1e-3



if __name__ == '__main__':
    unittest.main(verbosity=2)