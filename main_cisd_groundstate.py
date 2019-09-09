import sys
sys.path.append('./fort_src/')
sys.path.append('./src/')
sys.path.append('./input/')

import numpy as np, h5py
from scipy.optimize import fixed_point, broyden1, minimize
from scipy.special import comb
import time

from iofuncs import *
from inttran import *
from rcisd import *

# from rCISDEqns impo

#
# GLOBAL VARIABLES
#

mu_step_0 = +5e-2
len_t1 = 0
len_t2 = 0




#
# MAIN
# 

def main():
    start_time = time.time()
    
    #################################################################
    #                       READ INPUT                              #
    #       Get the data first - including the attributes           #
    #################################################################

    print('\n==============================================================')
    print('\tReading Input')
    print('==============================================================')

    # Initialize the Evolution module
    evol = IOps(inp_file='InputRCI')

    # Extract the parameters
    n_elec = evol.n_elec
    beta_f = evol.beta_f
    beta_pts = evol.beta_pts
    ntol = evol.ntol
    deqtol = evol.deqtol
    e_nuc = evol.e_nuc

    # Integrals
    eigs, h1, eri, attrs = evol.loadHDF()
    nso = np.size(h1, axis=0)

    # Transform everything to RHF Basis
    n_rhf = int(nso/2)
    fock = np.diag(eigs[0:n_rhf])
    eri_rhf = np.zeros((n_rhf, n_rhf, n_rhf, n_rhf))
    for p in range(n_rhf):
        for q in range(n_rhf):
            for r in range(n_rhf):
                for s in range(n_rhf):
                    eri_rhf[p, q, r, s] = eri[p, n_rhf+q, r, n_rhf+s]

    n_occ = int(n_elec/2)
    print('Number of Electrons: ', n_elec)
    print('NOcc : ', n_occ)
    print('NAO  : ', n_rhf)
    n_vir = int(n_rhf - n_occ)

    input_time = time.time()

    print('Number of Spin Orbitals:',nso)
    print('--------------------------------------------------------------\n')

    #################################################################
    #               INITIALIZE VARIABLES / PARAMETERS               #
    #       Initializing the arrays for computing various values    #
    #################################################################

    global len_t1
    global len_t2

    len_t1 = int(n_occ * n_vir)
    len_t2 = int(comb(n_occ,2)*comb(n_vir,2))

    # The T and Z vectors
    t1_vec = np.zeros(len_t1)
    t2_vec = np.zeros(len_t2)

    # CC and CI amps
    ci_amps = 1e-2*np.ones(len_t1 + len_t2)

    #################################################################
    #                   Do Fixed Point Iteration                    #
    #################################################################

    energy = cisd_energy(fock, eri_rhf, ci_amps, n_occ)
    print('Abs max of amplitudes before = ',np.max(np.abs(ci_amps)))
    print('Energy unoptimized = ', energy)

    ci_amps_sol = fixed_point(
        lambda x: cisd_equations(fock, eri_rhf, x, n_occ),
        ci_amps,
    )

    # Let's check if the solution worked
    res = cisd_equations(fock, eri_rhf, ci_amps_sol, n_occ)
    print('Final Residual: ', np.linalg.norm(res))

    print('Abs max of amplitudes after = ',np.max(np.abs(ci_amps_sol)))

    e_corr = cisd_energy(fock, eri_rhf, ci_amps_sol, n_occ)
    print('Correlation Energy = ',e_corr)

    init_time = time.time()

    #################################################################
    #                       CLOSE EVERYTHING                        #
    #################################################################

    end_time = init_time

    # Time Profile
    print('----------------------------------------------')
    print('End of Program and Time Profile')
    print('----------------------------------------------')
    print('----------------------------------------------')
    print('Process              Time')
    print('----------------------------------------------')
    print('Reading Input        {}'.format(input_time-start_time))
    print('Init Variables       {}'.format(init_time-input_time))
    print('Writing the Output   {}'.format(end_time-init_time))
    print('Total Time           {}'.format(end_time-start_time))
    print('----------------------------------------------')
    print('----------------------------------------------\n')



    exit()

if __name__ == '__main__':
    main()

