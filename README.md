# Molecular Modeling Experiments

Mostly experiments using [SAMSON](https://www.samson-connect.net/) to model simple machines.

Based on Tom Moore's [Diamond Machine Parts](https://github.com/mooreth/Diamond_Machine_Parts) projects.

## Next steps

* create molecules in classes, instead of structural models
    * structural models should be created at a higher level, and molecules added to them

## Simulating with LAMMPS

You'll need LAMMPS input scripts for this step. They can be adapted from the examples here, but the main simulation input especially may need adjustment for whatever motion you want to induce.

These simulations also rely on a CH.airebo file (can remain unchanged).

. [Install LAMMPS](https://docs.lammps.org/Install.html)
. Before exporting your SAMSON model, minimize it pretty well. The SAMSON minimizer can deal with larger divergences than LAMMPS can.
. Export your SAMSON model using a SAMSON2LAMMPS python script.
  . You can modify the examples in this repo for different atom groups as needed.
. Run LAMMPS to minimize energy, bring to temperature, and then simulate your structure
  . `lmps -in ./Energy_min.input`
  . `lmps -in ./Bring_to_temp.input`
  . `lmps -in ./simulation.input`