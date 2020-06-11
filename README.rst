PyOrb
=========

Feature list
-------------

* 


To install
-----------------

.. code-block:: bash

   pip install pyorb


Example
---------

.. execute_code:: python
    :linenos:

    import pyorb

    orb = pyorb.Orbit(M0 = pyorb.M_sol)
    orb.update(a=1*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)

    #Convert and get cartesian elements

    print(orb.cartesian)