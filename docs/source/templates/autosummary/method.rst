{{ fullname | escape }}
{{ (fullname | escape | length)*"=" }}

.. currentmodule:: {{ module }}

{{ objname }}
{{ (objname | escape | length)*"-" }}

.. automethod:: {{ objname }}
