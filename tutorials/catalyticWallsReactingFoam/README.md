All tutorials in this folder use the 
<span style="font-family:Courier;">catalyticZonesReactingFoam</span> 
solver. This solver is an extension of rhoReactingFoam with the option to include 
several catalytic zones, which are modeled using a pseudo-homogeneous approach.

## Instructions

Each tutorial can be run by executing `./Allrun` from the case directory. The case directory is cleaned via `./Allclean`.

## Case description

Inspiration for these tutorials was found in the following paper:

Matteo Maestri, Alberto Cuoci (2013), 
"Coupling CFD with detailed microkinetic modeling in heterogeneous catalysis",
*Chemical Engineering Science 96*, pp. 106-117. 
[https://doi.org/10.1016/j.ces.2013.03.048]

However, since the kinetic model is not identical to the one used in above paper,
no direct comparison is possible.

**annulus**
Isothermal H2 fuel rich combustion in an annular reactor. Part of the wall 
consists of the Rh catalyst.

**twoSpheres**
Adiabatic H2 fuel rich combustion in a channel with two catalytic Rh spheres.
