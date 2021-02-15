# thermal-ci
### Configuration Interaction Theory for Finite Temperature Quantum Chemistry

Python library for thermofield based finite temperature configuration interaction theory in electronic many-body systems. This library uses the configuration interaction singles and doubles (CISD) theory. The theoretical framework is presented our paper, Thermofield theory for finite temperature quantum chemistry ([arxiv link](https://arxiv.org/abs/1901.06753)).

The code is largely written in Python 3.6, however, the essential modules that drives the imaginary time and chemical potential evolution are generated from drudge in Fortran.
These Fortran modules are then converted to python modules with the help of the utility `f2py`.

### Installation
The installation process is straightforward. First, install the required python packages, listed in `requirements.txt` and then execute
```bashscript
python3 setup.py install
```

### Using the thermal-ci package
1. To run the thermal CISD calculations, we need to provide the one- and two-electron integrals as inputs. In the current implementation, the thermal-ci package expects these integrals in the spin-orbital basis, stored as an HDF (.h5) file with the keys `h1` and `eri` respectively.

   `h1`   :   A 2-dimensional array of size (N,N)
   
   `eri`  :   A 4-dimensional array of size (N,N,N,N)
   
    Such integrals can be generated using standard Quantum Chemistry packages such as [pyscf](https://github.com/pyscf/pyscf). For example, for 2-site Hubbard model with U/t = 1, we can have a data file named `hub_2s_u1_data.h5`.

2. With the integrals' file ready, the details about the system and the thermal evolution needs to be specified in the `Input` file. As an example, see the the Input file in the scripts folder.

3. A script to compute the internal energy of electronic many-body systems using the thermal CISD approximation can be found in the `scripts` directory.
