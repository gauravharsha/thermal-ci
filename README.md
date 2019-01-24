# thermal-ci
### Configuration Interaction Theory for Finite Temperature Quantum Chemistry

A repository to generate finite temperature energy expectation values for many-electron systems using a CISD like approach within the Thermofield Formalism (see our [paper](https://arxiv.org/abs/1901.06753) for more details).
The code is largely written in Python 3.6, however, the essential modules that drives the imaginary time and chemical potential evolution are generated from drudge in Fortran.
These Fortran modules are then converted to python modules with the help of the utility `f2py`.

**_Note:_** The current version of this project is not really optimized and parallelized.

### Installation
As such, this package is a bunch of python functions and does not require any installation. However, the Fortran subroutines need to be compiled and converted into python modules.

To convert the Fortran files to python modules, run the `makemodule` csh script.
```bashscript
./makemodule
```

### Using the thermal-ci package
1. To run thermal-ci calculations for any hamiltonian, one needs as an input the one- and two-electron integrals in the spin-orbital basis, stored as an HDF (.h5) file with the keys `h1` and `eri` respectively.

   `h1`   :   A 2-dimensional array of size (N,N)
   
   `eri`  :   A 4-dimensional array of size (N,N,N,N)
   
    Such integrals can be generated using standard Quantum Chemistry packages such as [pyscf](https://github.com/pyscf/pyscf). For example, for 2-site Hubbard model with U/t = 1, we can have a data file named `hub_2s_u1_data.h5`.

2. Once the input data is prepared, open the `Input` file and edit / set the the appropriate parameters such as the name of the input data file, average number of electrons, final value of \beta up to which you need to perform the evolution, number of \beta grid points, etc. (see the Input file).

3. Finally, use the commands `./Cov_Run_CISD` or `./FixRef_Run_CISD` to perform the thermal-ci calculations in the covariant or the fixed-reference approach respectively. Again, details about these two approaches can be found in our [paper](https://arxiv.org/abs/1901.06753). However, the recommended approach, or rather the method which gives sane results for all the cases, is the covariant approach.
