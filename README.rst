PyOrb
=========

Feature list
-------------

Current features:
* Clear definition of an orbit, consistent throughout the code, including planar and circular orbits
* Kepler to Cartesian conversion
* Cartesian to Kepler conversion
* All function handles all special cases (e.g. planar and circular orbits)
* Convenient ``Orbit`` class or storing orbits and seamlessly convert between Kepler and Cartesian elements
* Access to all types of orbit anomalies
* Vectorized function for increased performance
* Access to alternative parameterizations such as Equinoctial elements

On the upcoming feature list:
* Can handle hyperbolic orbits
* C-implementation of conversion function for performance
* Converting of orbits to a byte-stream
* Saving orbits to file (binary or HDFS 5)



To install
-----------------

.. code-block:: bash

    pip install pyorb

or to do the "nightly" build:

.. code-block:: bash

    git clone https://github.com/danielk333/pyorb
    cd pyorb
    git checkout develop
    pip install .

Alternatively, if you are following updates closely you can install using ``pip install -e .`` so that in the future a ``git pull`` will update the library.


Example
---------

.. code-block:: python
    :linenos:

    import pyorb

    orb = pyorb.Orbit(M0 = pyorb.M_sol, degrees=True)
    orb.update(a=1*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)

    #Convert and get cartesian elements
    print(orb.cartesian)

    #Make eccentric and place at aphelion
    orb.e = 0.2
    orb.anom = 180

    #print cartesian position in AU at aphelion after the above changes
    print(orb.r/pyorb.AU)
