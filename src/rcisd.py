import sys
sys.path.append('../fort_src/')

from scipy.special import comb
from numba import jit, njit

import numpy as np
from CISDEqns import *

# Global parameter
t2_symm_tol = 5e-8

#
# Compression / Decompression functions
#

# @jit(nopython=True)
def _decompress_t2(T2_Compressed, n_occ, n_vir):
    """
    The Idea is to use the antisymmetric property of the T2 tensor:
        T2[p,q,r,s] = -T2[q,p,r,s]
                    = -T2[p,q,s,r]
        i.e. the first 2 and the last 2 indices are anti-symmetric

    This allows a compressed storage of T2's but in passing the T2 amps
    to the FORTRAN SubRoutine, we need to reconstruct the full T2
    which is done by this function
    """

    t2_out = np.zeros((n_vir, n_vir, n_occ, n_occ))
    m = 0

    for i in range(n_vir):
        for j in range(i+1,n_vir):
            for k in range(n_occ):
                for l in range(k+1,n_occ):
                    val = T2_Compressed[m]
                    t2_out[i,j,k,l] = val
                    t2_out[j,i,k,l] = -val
                    t2_out[i,j,l,k] = -val
                    t2_out[j,i,l,k] = val
                    m += 1
    return t2_out

def DecompressT2(T2_Compressed, n_occ, n_vir):
    """
    The Idea is to use the antisymmetric property of the T2 tensor:
        T2[p,q,r,s] = -T2[q,p,r,s]
                    = -T2[p,q,s,r]
        i.e. the first 2 and the last 2 indices are anti-symmetric

    This allows a compressed storage of T2's but in passing the T2 amps
    to the FORTRAN SubRoutine, we need to reconstruct the full T2
    which is done by this function
    """

    if np.size(T2_Compressed) != int(comb(n_occ, 2)*comb(n_vir, 2)):
        raise ValueError('Invalid Size of the compressed T2 array',T2_Compressed)

    t2_out = _decompress_t2(T2_Compressed, n_occ, n_vir)

    return t2_out

# @jit(nopython=True)
def _compress_t2(T2, n_occ, n_vir):
    """
    Numba parallelized function for Compressing 4-fold antisymmetric tensors
    """

    m = 0
    _t2_compressed = np.zeros( int( comb(n_occ, 2)*comb(n_vir, 2) ) )

    for i in range(n_vir):
        for j in range(i+1, n_vir):
            for k in range(n_occ):
                for l in range(k+1, n_occ):
                    val = T2[i,j,k,l]
                    _t2_compressed[m] = val
                    m += 1

    return _t2_compressed

def CompressT2(T2, n_occ, n_vir):
    """
    Here, we also CHECK that the input tensor has the right symmetry
    And proceed directly to COMPRESS in the following way:

        (I,J,K,L) --> I<J and K<L

        Arranged in lexicological order, for instance

            for first pair of indices:  (1,2) comes before (1,3) before (2,3) and so on.

            the actual order is:        (0,1,0,1) -> (0,1,0,2) -> ... -> (0,1,0,NSO-1) ->
                                        (0,1,1,2) -> .... and so on...
    """
    # Check symmetries first
    chk_fail = 0
    if np.max(np.abs(T2 + np.einsum('qprs->pqrs',T2))) > t2_symm_tol:
        chk_fail = 1
        print('first one failed')
    elif np.max(np.abs(T2 + np.einsum('pqsr->pqrs',T2))) > t2_symm_tol:
        chk_fail = 1
        print('second one failed')
    else:
        pass

    if chk_fail:
        raise ValueError('Incorrect Symmetry for the input Tensor')

    T2_compressed = _compress_t2(T2, n_occ, n_vir)

    return T2_compressed


#
# CISD Equations and functions
# 

def cisd_equations(fock, eri, t_amps, n_occ):
    
    #
    # Wrapper function to find root for cisd equations
    # 

    # Dimension parameters
    nao = np.size(fock, axis=0)
    n_vir = nao - n_occ

    # extracting t1 and t2
    t1 = np.reshape(t_amps[0:n_occ*n_vir], (n_vir, n_occ))
    t2 = DecompressT2(t_amps[n_occ*n_vir:], n_occ, n_vir)

    # Get the driver contribution in cisd equations
    r0, r1, r2 = rcisd(fock, eri, t1, t2, n_occ)

    r2_asymm = (
        r2 
        - np.einsum('qprs->pqrs', r2)
        - np.einsum('pqsr->pqrs', r2)
        + np.einsum('qpsr->pqrs', r2)
    )/4.0

    # Compress 'em
    r1_comp = np.reshape(r1, n_occ*n_vir)
    r2_comp = CompressT2(r2_asymm, n_occ, n_vir)

    return np.concatenate((r1_comp, r2_comp))

def cisd_energy(fock, eri, t_amps, n_occ):
    
    #
    # Wrapper function to compute RCISD Energy     
    # 

    # Dimension parameters
    nao = np.size(fock, axis=0)
    n_vir = nao - n_occ

    # extracting t1 and t2
    t1 = np.reshape(t_amps[0:n_occ*n_vir], (n_vir, n_occ))
    t2 = DecompressT2(t_amps[n_occ*n_vir:], n_occ, n_vir)

    # Get the driver contribution in cisd equations
    r0, r1, r2 = rcisd(fock, eri, t1, t2, n_occ)

    return r0
