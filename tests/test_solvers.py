#!/usr/bin/env python

'''
Test basic kepler functions
'''

import unittest
import numpy as np
import numpy.testing as nt

import pyorb.kepler as kep


class TestKeplerSolver(unittest.TestCase):

    def test_hyperbolic_kepler_guess(self):
        E = np.linspace(0.0, np.pi, num=300, dtype=np.float64)
        e = np.linspace(1.001, 10, num=500, dtype=np.float64)

        Ev, ev = np.meshgrid(E, e)

        Mv = kep.eccentric_to_mean(Ev, ev)

        E0 = kep.kepler_guess(Mv, ev)

        # the initial guess SHOULD be at least 35 degrees to true
        test_accuracy = np.abs(E0 - Ev)*180.0/np.pi < 35.0
        nt.assert_array_equal(
            test_accuracy,
            np.full(test_accuracy.shape, True)
        )

    def test_laguerre_solve_hyperbolic_kepler(self):
        E = np.linspace(0.0, np.pi, num=300, dtype=np.float64)
        e = np.linspace(1.001, 10, num=500, dtype=np.float64)

        for I, eit in enumerate(e):
            for J, Eit in enumerate(E):
                M = kep.eccentric_to_mean(Eit, eit)

                E0 = kep.kepler_guess(M, eit)
                E_calc, it = kep.laguerre_solve_kepler(E0, M, eit, tol=1e-12)
                M_calc = kep.eccentric_to_mean(E_calc, eit)
                fun_err = np.abs(M - M_calc)

                nt.assert_almost_equal(Eit, E_calc, decimal = 1e-9)
                assert fun_err < 1e-11

    def test_kepler_guess(self):
        E = np.linspace(0.0, 2.0*np.pi, num=100, dtype=np.float64)
        e = np.linspace(0, 0.99, num=100, dtype=np.float64)

        Ev, ev = np.meshgrid(E, e)

        Mv = kep.eccentric_to_mean(Ev, ev)

        E0 = kep.kepler_guess(Mv, ev)

        # the initial guess SHOULD be at least 30 degrees to true
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
