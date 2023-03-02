#!/usr/bin/env python

'''
Test basic kepler functions
'''

import unittest
import numpy as np
import numpy.testing as nt

import pyorb.kepler as kep


class TestRotations(unittest.TestCase):

    def setUp(self):
        mat = np.eye(3, dtype=np.float64)
        self.x = mat[:, 0]
        self.y = mat[:, 1]
        self.z = mat[:, 2]
        self.theta = np.pi/2

    def test_rot_mat_x(self):
        R = kep.rot_mat_x(self.theta, dtype=np.float32)
        assert R.dtype == np.float32

        R = kep.rot_mat_x(self.theta)
        assert R.dtype == np.float64

        nt.assert_almost_equal(self.z, R @ self.y)

    def test_rot_mat_y(self):
        R = kep.rot_mat_y(self.theta, dtype=np.float32)
        assert R.dtype == np.float32

        R = kep.rot_mat_y(self.theta)
        assert R.dtype == np.float64

        nt.assert_almost_equal(-self.z, R @ self.x)

    def test_rot_mat_z(self):
        R = kep.rot_mat_z(self.theta, dtype=np.float32)
        assert R.dtype == np.float32

        R = kep.rot_mat_z(self.theta)
        assert R.dtype == np.float64

        nt.assert_almost_equal(self.y, R @ self.x)
