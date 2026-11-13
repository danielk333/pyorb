"""Retune orbital size through period and angular speed."""
import numpy as np
import numpy.testing as nt
import pytest

from pyorb import Orbit


@pytest.mark.parametrize("attribute", ["period", "mean_motion"])
@pytest.mark.parametrize("degrees", [False, True])
@pytest.mark.parametrize("direct", [False, True])
@pytest.mark.parametrize("num", [1, 3])
@pytest.mark.parametrize("scalar", [False, True])
def test_retuned_circular_orbits_follow_the_requested_period(attribute, degrees, direct, num, scalar):
    orbit = Orbit(1., G=1., num=num, m=np.arange(num)*.1,
                  a=1., e=0., i=0., omega=0., Omega=0., anom=0.,
                  degrees=degrees, direct_update=direct)
    period = 3. if scalar else np.arange(num) + 3.
    norm = 360. if degrees else 2*np.pi
    value = period if attribute == "period" else norm/period
    setattr(orbit, attribute, value)
    mu = 1. + orbit.m
    expected_a = (mu*(period/(2*np.pi))**2)**(1./3.)
    nt.assert_allclose(orbit.a, expected_a, rtol=1e-14)
    nt.assert_allclose(orbit.period, np.broadcast_to(period, (num,)), rtol=1e-14)
    nt.assert_allclose(orbit.mean_motion, norm/period, rtol=1e-14)
    orbit.propagate(period/4)
    expected = np.zeros((6, num))
    expected[1] = expected_a
    expected[3] = -np.sqrt(mu/expected_a)
    nt.assert_allclose(orbit.cartesian, expected, rtol=1e-13, atol=1e-13)


@pytest.mark.parametrize("attribute", ["period", "mean_motion"])
def test_retuning_respects_read_only_elements(attribute):
    orbit = Orbit(1., G=1., a=1., e=.1, i=0., omega=0., Omega=0., anom=0.,
                  kepler_read_only=True)
    with pytest.raises(AttributeError, match="read only"):
        setattr(orbit, attribute, 3.)
