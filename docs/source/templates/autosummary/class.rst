{{ fullname | escape }}
{{ (fullname | escape | length)*"=" }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

Summary
-------

{% if attributes %}
Attributes
^^^^^^^^^^

.. autosummary::
  :toctree: .
{% for item in attributes %}
  ~{{ name }}.{{ item }}
{%- endfor %}
{% endif %}


{% if methods %}
Methods
^^^^^^^

.. autosummary::
  :toctree: .
{% for item in methods %}
  ~{{ name }}.{{ item }}
{%- endfor %}
{% endif %}


{% block attributes %}
{% if attributes %}
Attributes
----------

{% for item in attributes %}
.. autoattribute:: {{ name }}.{{ item }}
    :noindex:
{%- endfor %}
{% endif %}
{% endblock %}


{% block methods %}
{% if methods %}
Methods
-------

{% for item in methods %}
.. automethod:: {{ name }}.{{ item }}
    :noindex:
{%- endfor %}
{% endif %}
{% endblock %}