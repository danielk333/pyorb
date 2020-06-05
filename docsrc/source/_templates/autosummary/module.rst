{{ name | escape | underline }}

.. automodule:: {{ fullname }}

.. currentmodule:: {{ fullname }}

Module summary
--------------

{% if classes %}
.. rubric:: Classes

.. autosummary::
    :toctree: .
    {% for class in classes %}
    {{ class }}
    {% endfor %}

{% endif %}

{% if functions %}
.. rubric:: Functions

.. autosummary::
    :toctree: .
    {% for function in functions %}
    {{ function }}
    {% endfor %}

{% endif %}

Contents
----------

{% if classes %}
{% for class in classes %}

.. rubric:: {{ class }}

.. autoclass:: {{ class }}
   :show-inheritance:
   :noindex:
   :members:
   :inherited-members:

{% endfor %}
{% endif %}

{% if functions %}
.. rubric:: Functions
{% for function in functions %}

.. autofunction:: {{ function }}
   :noindex:

{% endfor %}
{% endif %}
