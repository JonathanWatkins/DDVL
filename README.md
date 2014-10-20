DDVL
====

Densiy Driven Vortex Lattice Simulation

This code requires the intel compiler and boost libraries.


Changes in version 2
--------------------

Version 2 now reads the system geometry from file. It supports particle types that can be integrated (on not) separately. Particles can also be initialised with a starting velocity.



Installation
------------

*On Windows*

To compile with cmake. Run

cmake . -G "NMake Makefiles"

Then 

nmake


Visualising
-----------

To view the simulation use the vis.nb mathematica notebook.
