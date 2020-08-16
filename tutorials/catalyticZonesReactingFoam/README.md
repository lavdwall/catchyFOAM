All tutorials in this folder use the 
<span style="font-family:Courier;">catalyticZonesReactingFoam</span> 
solver. This solver is an extension of rhoReactingFoam with the option to include 
several catalytic zones, which are modeled using a pseudo-homogeneous approach.

## Instructions

Each tutorial can be run by executing `./Allrun` from the case directory. If 
Cantera is installed, the simulations can be compared with a 1D plug flow reactor 
simulation by executing `python3 validate.py`. 
The case directory is cleaned via `./Allclean`.

To clean, run and validate all 
<span style="font-family:Courier;">catalyticZonesReactingFoam</span> 
tutorials in this folder:
```
./Allclean
./Allrun
./Allvalidate
```

## Case description

**isoOxidativeCouplingMethane**
Isothermal 1D packed bed simulation for oxidative coupling of methane. This case
is used to test the implementation of the catalytic chemistry. All cells are 
selected as catalytic zone.

**adiOxidativeCouplingMethane**
Adiabatic 1D packed bed simulation for oxidative coupling of methane. This case
is used to test the implementation of the catalytic chemistry in combination
with reaction heat. Only the first part of the geometry is the bed,
downstream of it, gas phase chemistry is responsible for further reactions.

**porousMediaIsoOxidativeCouplingMethane**
This case is identical to the isoOxidativeCouplingMethane case, but  includes 
the use of the porousMedia function object to calculate the pressure drop in the
catalytic zone.

**methanePartialOxidation**
This case is similar to the isoOxidativeCouplingMethane case, but with chemistry
for methane partial oxidation. The case is primarily used to test the 
implementation of coverage-dependent surface reaction rates.


*NOTE: The agreement between these simulations and the validation cases from Cantera can*
*be further improved by disabling diffusion (put As to 0 in mechanism/transportProperties),* 
*refining the mesh, optimizing discretization settings etc.*