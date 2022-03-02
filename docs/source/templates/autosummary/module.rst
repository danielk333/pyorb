{{ fullname | escape }}
{{ (fullname | escape | length)*"=" }}

.. automodule:: {{ fullname }}

.. currentmodule:: {{ fullname }}

Summary
-------

{% if modules %}
Modules
^^^^^^^

.. autosummary::
    :toctree: .
    {% for module in modules %}
    {{ module }}
    {% endfor %}

{% endif %}

{% if classes %}
Classes
^^^^^^^

.. autosummary::
    :toctree: .
    {% for class in classes %}
    {{ class }}
    {% endfor %}

{% endif %}

{% if functions %}
Functions
^^^^^^^^^

.. autosummary::
    :toctree: .
    {% for function in functions %}
    {{ function }}
    {% endfor %}

{% endif %}


{% block classes %}
{% if classes %}
Classes
-------

{% for item in classes %}
.. autoclass:: {{ item }}
    :noindex:
{%- endfor %}
{% endif %}
{% endblock %}


{% block functions %}
{% if functions %}
Functions
---------

{% for item in functions %}
.. autofunction:: {{ item }}
    :noindex:
{%- endfor %}
{% endif %}
{% endblock %}