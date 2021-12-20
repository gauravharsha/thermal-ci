# tfd-cisd
### Thermofield dynamics based configuration interaction with singles and doubles

[![Python application](https://github.com/gauravharsha/thermal-ci/actions/workflows/python-app.yml/badge.svg)](https://github.com/gauravharsha/thermal-ci/actions/workflows/python-app.yml)

Python package to computing grand-canonical ensemble properties of many-body electronic systems at finite-temperatures using configuration interaction theory with singles and doubles. The theoretical framework is presented our paper

> Thermofield Theory for Finite-Temperature Quantum Chemistry, G. Harsha, T. M. Henderson, and G. E. Scuseria, [J. Chem. Phys. 150, 154109 (2019)](https://aip.scitation.org/doi/10.1063/1.5089560) ([arxiv link](https://arxiv.org/abs/1901.06753)).

The code requires Python 3.7 or higher. Computationally intensive modules and functions that drives the imaginary time and chemical potential evolution are written in Fortran. The equations were generated using a modified [drudge](https://github.com/tschijnmo/drudge) class called `ThermofieldDrudge` which, along with the drudge scripts for generating equations, is included under `drudge` directory.

## Installation
If you already have the packages listed in `requirements.txt`, you can simply install this package by
```bashscript
python setup.py install
```

## Usage
1. Prepare the one- and two-electron integrals in the basis in which mean-field Hamiltonian is diagonal. We will also need these mean-field eigenvalues. In the current implementation, the thermal-ci package expects these integrals in the spin-orbital basis, stored as an HDF (.h5) file with the keys `h1`, `eri` and `eigs` respectively.

   `h1`   :   2D array of size (N,N)
   
   `eri`  :   4D array of size (N,N,N,N)
   
   `eigs` :   1D array of size (N)
   
    The integrals can be generated using standard Quantum Chemistry packages such as [pyscf](https://github.com/pyscf/pyscf). For example, for 2-site Hubbard model with U/t = 1, we can have a data file named `hub_2s_u1_data.h5`.

    Each of these objects in the HDF5 integrals file should have the same list of attributes.

2. With the integrals' file ready, the details about the system and the thermal evolution needs to be specified in the `Input` file. As an example, see the the Input file in the `scripts` subdirectory.

3. Scripts to perform the imaginary time evolution while maintaining a fixed number of particles, and compute the internal energies within this thermal CISD approximation are included in the  `scripts` subdirectory.

    `main_cisd_fixedN.py`           :   Python script for imaginary time evolution using covariant CISD.

    `main_cisd_fixedN_fixref.py`    :   Python script for imaginary time evolution using fixed-reference CISD.

(*NOTE*: For details on covariant and fixed-reference theory, see the article linked above.)
    
(*NOTE*: Fixed-reference CISD code has not been tested rigorously.)
