# The "coarse-grained amyloids MD simulation" project

This project uses the lammps\_multistate\_rods library to make coarse-grained multi-scale hybrid MD-MC simulations of amyloidogenesis.

The project is written in python3 and contains scipts and programs to set up everything for the simulations, submit them as jobs on remote clusters and exctract, analyse and visualise the data from those simulations.

## Installation / setup

Study the scripts in the *setup_scripts* folder, that should be enough.

### Requirements

The requirements are:
1. **Python 3** (min 3.5, it's **NOT Python2 compatible**) with *mpi4py* and *matplotlib* (for analysis tools), along with libraries needed for the *lammps_multistate_rods* library:
    ```
    pip install mpi4py matplotlib ...
    ```

2. **Jupyter** (or equivalent to run IPython notebooks) is necessary to run the analysis notebooks
	* The following terminal line can be of use to enable Jupyter to use a *venv kernel* that has everything installed (presumably the system python won't be used):
	```
    python -m ipykernel install --user --name=<kernel-name>
    ```
    Then that kernel can be chosen in Jupyter/IPython to run the code (menu under "Kernel").

