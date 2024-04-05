########
Examples
########

See `/LocVib/example/` directory for some examples of typical runs.


.. code-block:: console

   /LocVib/example/
   ├── 1_composition.py
   ├── 2_locmodes.py
   ├── 3_couplings.py
   ├── Ala10
   │   ├── coord
   │   ├── restart
   │   └── snf.out

As a typical example of a vibrational spectra calculation in terms of localized modes, we take the results of an **SNF** calculation (https://reiher.ethz.ch/software/snf.html) of an *Ala10* molecule (harmonic approximation).

Composition
===========

Here we see the script for the first part of the example `1_composition.py`:

.. dropdown:: Open Python Script - 1_composition.py

   .. literalinclude:: ../example/1_composition.py
       :language: python

When we run the script we get:

   >>> /LocVib/example/ % python 1_composition.py

.. dropdown:: Open Output: python 1_composition.py

   .. literalinclude:: examples_output/1_output.txt
       :language: text

locmodes
========

Here we see the script for the first part of the example `2_locmodes.py`:

.. dropdown:: Open Python Script - 2_locmodes.py

   .. literalinclude:: ../example/2_locmodes.py
       :language: python

When we run the script we get:

   >>> /LocVib/example/% python 2_locmodes.py

.. dropdown:: Open Output: python 2_locmodes.py

   .. literalinclude:: examples_output/2_output.txt
       :language: text


couplings
=========

Here we see the script for the first part of the example `3_couplings.py`:

.. dropdown:: Open Python Script - 3_couplings.py

   .. literalinclude:: ../example/3_couplings.py
       :language: python

When we run the script we get:

   >>> /LocVib/example/% python 3_couplings.py

.. dropdown:: Open Output: python 3_couplings.py

   .. literalinclude:: examples_output/3_output.txt
       :language: text
