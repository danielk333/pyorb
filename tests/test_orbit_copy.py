"""Copies keep the latest representation of each lazily updated orbit."""
import copy

import numpy as np
import numpy.testing as nt
import pytest

from pyorb import Orbit


def mixed_orbits(anomaly_type):
    orbit = Orbit(1., G=1., num=3, direct_update=False, type=anomaly_type,
                  a=[2., 3., 4.], e=.2, i=.3, omega=.4, Omega=.5, anom=.6)
    orbit.calculate_cartesian()
    orbit.update(inds=[0, 2], a=[5., 6.])
    orbit.update(inds=[1], vx=.2)
    return orbit


@pytest.mark.parametrize("anomaly_type", ["true", "eccentric", "mean"])
@pytest.mark.parametrize("selection", [0, slice(None), slice(None, None, -1),
                                       [2, 0], np.array([False, True, True])])
@pytest.mark.parametrize("first", ["kepler", "cartesian"])
def test_slice_preserves_mixed_pending_updates(anomaly_type, selection, first):
    original = mixed_orbits(anomaly_type)
    reference = mixed_orbits(anomaly_type)
    expected = {"cartesian": reference.cartesian[:, selection].reshape(6, -1),
                "kepler": reference.kepler[:, selection].reshape(6, -1)}
    result = original[selection]
    nt.assert_allclose(getattr(result, first), expected[first])
    nt.assert_allclose(result.cartesian, expected["cartesian"])
    nt.assert_allclose(result.kepler, expected["kepler"])
    original_before = original.cartesian
    result.a *= 1.01
    assert np.isfinite(result.cartesian).all()
    nt.assert_array_equal(original.cartesian, original_before)
    result.add(a=2., e=.1, i=.2, omega=.3, Omega=.4, anom=.5)
    result.delete(-1)
    assert np.isfinite(result.cartesian).all()


@pytest.mark.parametrize("clone", [copy.copy, lambda orbit: orbit.copy(),
                                   lambda orbit: next(iter(orbit))])
def test_copy_and_iteration_preserve_latest_kepler_state(clone):
    original = mixed_orbits("true")
    result = clone(original)
    nt.assert_array_equal(result.a, original.a[:result.num])
