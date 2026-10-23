"""Element conversion can reuse callers' orbital element buffers."""
import numpy as np
import numpy.testing as nt
import pytest

from pyorb import kep_to_cart


@pytest.mark.parametrize("eccentricity", [0., .2, 1., 1.2])
@pytest.mark.parametrize("vector", [False, True])
@pytest.mark.parametrize("degrees", [False, True])
def test_repeated_conversion_keeps_periapsis_elements(eccentricity, vector, degrees):
    a, mu = 1e7, 3.986004418e14
    elements = np.array([a, eccentricity, 0., 0., 0., 0.])
    if not vector:
        elements = np.column_stack([elements, elements])
    original = elements.copy()
    radius = a if eccentricity == 1 else a * abs(1 - eccentricity)
    expected = np.array([radius, 0., 0., 0., np.sqrt(mu * (1 + eccentricity) / radius), 0.])
    if not vector:
        expected = np.column_stack([expected, expected])
    for _ in range(2):
        nt.assert_allclose(kep_to_cart(elements, mu=mu, degrees=degrees), expected, rtol=1e-14)
        nt.assert_array_equal(elements, original)


@pytest.mark.parametrize("degrees", [False, True])
@pytest.mark.parametrize("readonly", [False, True])
def test_conversion_accepts_shared_and_readonly_element_buffers(degrees, readonly):
    backing = np.array([[1e7, 2e7, 3e7, 4e7], [.2, 1.2, 1., .3],
                        [27.11, 33.2, 51.7, 72.8], [35.73, 61.3, 23.1, 17.9],
                        [42.22, 18.4, 94.6, 62.3], [17.33, 24.1, 35.8, 42.7]])
    if not degrees:
        backing[2:] = np.radians(backing[2:])
    original = backing.copy()
    elements = backing[:, ::2]
    elements.flags.writeable = not readonly
    result = kep_to_cart(elements, degrees=degrees)
    assert np.isfinite(result).all()
    nt.assert_array_equal(backing, original)
