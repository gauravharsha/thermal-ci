# thermal-ci
#### Configuration Interaction Theory for Finite Temperature Quantum Chemistry

A repository to generate finite temperature energy expectation values for many-electron systems using a CISD like approach within the Thermofield Formalism published here.
The code is largely written in Python 3.6, however, the essential modules that drives the imaginary time and chemical potential evolution are generated from drudge in Fortran.
These Fortran modules are then converted to python modules with the help of the utility `f2py`

#### Installation
To convert the Fortran files to python modules, run the `makemodule` file in bash.


