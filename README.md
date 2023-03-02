# PyOrb

PyOrb is a lightweight package designed to convert back and forth between cartesian and kepler coordinates seamlessly and in a physically consistent manner, following defined rules. It provides a convenience class for handling orbits and is tested for special cases such as planar and circular orbits.

See full documentation [here](https://danielk.developer.irf.se/pyorb/).

![Example interactive orbit](/docs/source/static/example.gif)

## Feature list

Current features:
- Clear definition of an orbit, consistent throughout the code, including planar and circular orbits
- Kepler to Cartesian conversion
- Cartesian to Kepler conversion
- Can handle hyperbolic orbits
- All function handles all special cases (e.g. planar, circular and parabolic orbits)
- Convenient ``Orbit`` class or storing orbits and seamlessly convert between Kepler and Cartesian elements
- Access to all types of orbit anomalies
- Vectorized function for increased performance
- Access to alternative parameterizations such as Equinoctial elements


On the upcoming feature list:
- C-implementation of conversion function for performance
- Converting of orbits to a byte-stream
- Saving orbits to file (binary or HDFS 5)


## To install


```bash
pip install pyorb
```
or to do the "nightly" build:

```bash
git clone https://github.com/danielk333/pyorb
cd pyorb
git checkout develop
pip install .
```

Alternatively, if you are following updates closely you can install using ``pip install -e .`` so that in the future a ``git pull`` will update the library.


## Example

```python
import pyorb

orb = pyorb.Orbit(M0 = pyorb.M_sol, degrees=True)
orb.update(a=1*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)

# Convert and get cartesian elements
print(orb.cartesian)

# Make eccentric and place at aphelion
orb.e = 0.2
orb.anom = 180

# print cartesian position in AU at aphelion after the above changes
print(orb.r/pyorb.AU)
```


## Ellipse and angle definitions

Variables:
 - **a**: Semi-major axis
 - **e**: Eccentricity
 - **i**: Inclination
 - **omega**: Argument of perihelion
 - **Omega**: Longitude of the ascending node
 - **nu**: True anomaly
 - **E**: Elliptic, parabolic or hyperbolic eccentric anomaly
 - **M**: Mean anomaly


Orientation of the ellipse in the coordinate system and angle definitions:
 - For zero inclination: the ellipse is located in the x-y plane.
 - The direction of motion as True anomaly increases for a zero inclination orbit is anti-coockwise, i.e. from +x towards +y.
 - If the eccentricity is increased for an unrotated orbit, the periapsis will lie in +x direction.
 - If the inclination is increased, the ellipse will rotate around the x-axis, so that +y is rotated toward +z.
 - An increase in Longitude of the ascending node corresponds to a rotation around the z-axis so that +x is rotated toward +y.
 - Changing argument of perihelion will not change the plane of the orbit, it will rotate the orbit in the plane.
 - Changing argument of perihelion will rotate the periapsis in the direction of motion.
 - True anomaly measures from the +x axis, i.e **nu = 0** is located at periapsis and **nu = pi** at apoapsis.
 - All anomalies and orientation angles reach between **0** and **2pi**
 - If the inclination is **0** or **pi** the longitude of the ascending node is always zero (the rotation is described by only argument of perihelion).
 - If the eccentricity is zero, the argument of perihelion is always zero (the rotation is described by only the longitude of the ascending node).
 - If both **e=0** and **i=0** or **i=pi**: the position on the circle is only described by the anomaly.
 - The eccentric anomaly is used for elliptic, parabolic and hyperbolic cases but the kepler equation changes accordingly.
 - For parabolic and hyperbolic orbits the true and eccentric anomaly wraps at **pi** to the same trajectory, not the mirror version.

 Shape definitions:
 - The Semi-major axis is always positive. 
 - In the case of a parabolic orbit, as the Semi-major axis is undefined it is used as the periapsis distance instead.

## Notes

### Disabling direct conversion

There are two toggle flags in the ``pyorb.Orbit`` class for changing the conversion behavior: ``direct_update`` and ``auto_update`` that are ``True`` by default. 

Disabling ``direct_update`` will stop automatic conversion between elements if any element is changed. This would allow for e.g. 
```python
orb.a = 1
orb.omega = 0
```
without any conversion to be done. However, as the kepler elements changed, the class has internally tracked this change and if ``auto_update=True`` once an access to a cartesian property is performed, e.g. ``print(orb.x)``, the conversion is performed so that the pair of cartesian-kepler elements are never contradictory.

If also ``auto_update`` is disabled, the update between kepler and cartesian needs to be manually by calling 
```python
orb.calculate_cartesian()
```
or 
```python
orb.calculate_kepler()
```

### Using conversion functions directly

The ``pyorb.cart_to_kep`` or ``pyorb.kep_to_cart`` uses True anomaly and takes only numpy arrays ordered as per the function documentation.

### Frames

Remember that an Keplerian orbit only makes sense in an inertial frame if gravitation dominated physics is your concern.

### Array orbits

- Properties act on ALL orbits in the class
- Only way to update individual orbits of a set is to use ``self.update`` with the ``inds`` keyword
- Iterations are passive, the objects are copies from the array so the array itself is NOT modified
