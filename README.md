# Molecular Modeling Experiments

Mostly experiments using [SAMSON](https://www.samson-connect.net/) to model simple machines.

Based on Tom Moore's [Diamond Machine Parts](https://github.com/mooreth/Diamond_Machine_Parts) projects.

## Scripts

Use the python files in the scripts directory to autogenerate larger structures. These can be triggered from within the python interpeter of SAMSON.

Note that SAMSON's python interpreter will cache imported modules. If you want to iteratively make changes to the scripts, run these lines first to make sure local modules get re-loaded:

- This doesn't seem to work actually. Gotta investigate more.
```
get_ipython().run_line_magic("load_ext", "autoreload")
get_ipython().run_line_magic("autoreload", "2")  # Reload all modules before executing
```

## Next steps

* linear gear teeth
* utilities script with calculate_atomic_distance, translate_molecule, rotate_molecule, remove bond, etc.

## Future Experiments

* rack and pinion
* atomCAD
* [planar pseudogears](http://apm.bplaced.net/w/index.php?title=Linear_reciprocative_pseudogears)

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