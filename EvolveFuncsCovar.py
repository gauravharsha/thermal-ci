import numpy as np
from itertools import permutations
from scipy.misc import comb
from PyPTLib import *
from CovMuPT2 import *
from OpDerMP2 import *
# from CovBetaPT2 import *

def parity(perm,perm0):
    """
    Here we assume that the permutation is either 3 or 4 characters long
    """
    if len(perm)!=len(perm0):
        raise(ValueError,'Input strings must of the same length',(perm,perm0))
    elif len(perm)==4:
        loc0 = np.linspace(0,3,4)
        loc1 = loc0
        for i in range(4):
            loc1[i] = perm.find(perm0[i])
        pdiff = np.sum( np.abs(loc1 - loc0) )
        if pdiff%2 == 0:
            return 1
        else:
            return -1
    elif len(perm)==3:
        loc0 = np.linspace(0,2,3)
        loc1 = loc0
        for i in range(3):
            loc1[i] = perm.find(perm0[i])
        pdiff = np.sum( np.abs(loc1 - loc0) )
        if pdiff%2 == 0:
            return -1
        else:
            return 1
    else:
        raise(ValueError,'One or both of the input lengths not supported',(perm,perm0))



def T2_Decompress(T2_Compressed,NSO):
    """
    The Idea is to use the antisymmetric property of the T2 tensor:
        T2[p,q,r,s] = -T2[q,p,r,s]
                    = -T2[p,q,s,r]
        i.e. the first 2 and the last 2 indices are anti-symmetric

    This allows a compressed storage of T2's but in passing the T2 amps
    to the FORTRAN SubRoutine, we need to reconstruct the full T2
    which is done by this function
    """
    if np.size(T2_Compressed) != int(comb(NSO,2)**2):
        raise(ValueError,'Invalid Size of the compressed T2 array',t2_compr)

    t2 = np.zeros((NSO,NSO,NSO,NSO))
    m = 0

    for i in range(NSO):

        for j in range(i+1,NSO):

            for k in range(NSO):

                for l in range(k+1,NSO):

                    val = T2_Compressed[m]
                    t2[i,j,k,l] = val
                    t2[j,i,k,l] = -val
                    t2[i,j,l,k] = -val
                    t2[j,i,l,k] = val

                    m += 1
    return t2


def T2_Compress(T2):
    """
    Here, we also CHECK that the input tensor has the right symmetry
    And proceed directly to COMPRESS in the following way:

        (I,J,K,L) --> I<J  and K<L

        Arranged in lexicological order, for instance

            for first pair of indices:  (1,2) comes before (1,3) before (2,3) and so on.

            the actual order is:        (0,1,0,1) -> (0,1,0,2) -> ... -> (0,1,0,NSO-1) ->
                                        (0,1,1,2) -> .... and so on...
    """
    NSO = np.size(T2,axis=0)
    m = 0
    T2_compressed = np.zeros( int( comb(NSO,2)**2 ) )

    for i in range(NSO):
        for j in range(i+1,NSO):
            for k in range(NSO):
                for l in range(k+1,NSO):

                    val = T2[i,j,k,l]

                    chk1 = np.isclose(T2[j,i,k,l],-val,rtol=1e-7)
                    chk2 = np.isclose(T2[i,j,l,k],-val,rtol=1e-7)
                    chk3 = np.isclose(T2[j,i,l,k],val,rtol=1e-7)

                    if chk1 and chk2 and chk3:
                        T2_compressed[m] = val
                        m += 1
                    else:
                        raise(ValueError,'Incorrect Symmetry for the input Tensor',T2)
    
    return T2_compressed


def T3_Decompress(T3_Compressed,NSO):
    """
    The Idea is to use the antisymmetric property of the T2 tensor:
        T2[p,q,r,s] = -T2[q,p,r,s]
                    = -T2[p,q,s,r]
        i.e. the first 2 and the last 2 indices are anti-symmetric

    This allows a compressed storage of T2's but in passing the T2 amps
    to the FORTRAN SubRoutine, we need to reconstruct the full T2
    which is done by this function
    """
    if np.size(T3_Compressed) != int(comb(NSO,3)**2):
        raise(ValueError,'Invalid Size of the compressed T4 array',T4_Compressed)

    t3 = np.zeros((NSO,NSO,NSO,NSO,NSO,NSO))
    m = 0
    up_str = 'ijk'
    dn_str = 'abc'
    perm_up = [''.join(p) for p in permutations(up_str)]
    perm_dn = [''.join(p) for p in permutations(dn_str)]

    for i in range(NSO):
        for j in range(i+1,NSO):
            for k in range(j+1,NSO):
                for a in range(NSO):
                    for b in range(a+1,NSO):
                        for c in range(b+1,NSO):
                            val = T3_Compressed[m]

                            for p_up in perm_up:
                                for p_dn in perm_dn:
                                    exec(
                                        "t3[%s,%s,%s,%s,%s,%s] = %d*%d*%d",
                                        (p_up[0],p_up[1],p_up[2],p_dn[0],p_dn[1],p_dn[2],
                                            parity(p_up,up_str), parity(p_dn,dn_str), val)
                                    )
                            m += 1

    return t3

def T3_Compress(T3):
    """
    Here, we also CHECK that the input tensor has the right symmetry
    And proceed directly to COMPRESS in the following way:

        (I,J,K,L) --> I<J  and K<L

        Arranged in lexicological order, for instance

            for first pair of indices:  (1,2) comes before (1,3) before (2,3) and so on.

            the actual order is:        (0,1,0,1) -> (0,1,0,2) -> ... -> (0,1,0,NSO-1) ->
                                        (0,1,1,2) -> .... and so on...
    """
    NSO = np.size(T3,axis=0)
    m = 0
    T3_compressed = np.zeros( int( comb(NSO,3)**2 ) )

    for i in range(NSO):
        for j in range(i+1,NSO):
            for k in range(j+1,NSO):
                for a in range(NSO):
                    for b in range(a+1,NSO):
                        for c in range(b+1,NSO):

                            val = T3[i,j,k,a,b,c]

                            chk[0] = np.isclose(T3[j,i,k,a,b,c],-val,rtol=1e-7)
                            chk[1] = np.isclose(T3[i,k,j,a,b,c],-val,rtol=1e-7)
                            chk[2] = np.isclose(T3[k,j,i,a,b,c],-val,rtol=1e-7)

                            chk[3] = np.isclose(T3[i,j,k,b,a,c],-val,rtol=1e-7)
                            chk[4] = np.isclose(T3[i,j,k,a,c,b],-val,rtol=1e-7)
                            chk[5] = np.isclose(T3[i,j,k,c,a,b],-val,rtol=1e-7)

                            if all(chk):
                                T3_compressed[m] = val
                                m += 1
                            else:
                                raise(ValueError,'Incorrect Symmetry for the input Tensor',T3)
    
    return T3_compressed

def T4_Decompress(T4_Compressed,NSO):
    """
    The Idea is to use the antisymmetric property of the T2 tensor:
        T2[p,q,r,s] = -T2[q,p,r,s]
                    = -T2[p,q,s,r]
        i.e. the first 2 and the last 2 indices are anti-symmetric

    This allows a compressed storage of T2's but in passing the T2 amps
    to the FORTRAN SubRoutine, we need to reconstruct the full T2
    which is done by this function
    """
    if np.size(T4_Compressed) != int(comb(NSO,4)**2):
        raise(ValueError,'Invalid Size of the compressed T4 array',T4_Compressed)

    t4 = np.zeros((NSO,NSO,NSO,NSO,NSO,NSO,NSO,NSO))
    m = 0
    up_str = 'ijkl'
    dn_str = 'abcd'
    perm_up = [''.join(p) for p in permutations(up_str)]
    perm_dn = [''.join(p) for p in permutations(dn_str)]

    for i in range(NSO):
        for j in range(i+1,NSO):
            for k in range(j+1,NSO):
                for l in range(k+1,NSO):
                    for a in range(NSO):
                        for b in range(a+1,NSO):
                            for c in range(b+1,NSO):
                                for d in range(c+1,NSO):
                                    val = T4_Compressed[m]

                                    for p_up in perm_up:
                                        for p_dn in perm_dn:
                                            exec(
                                                "t4[%s,%s,%s,%s,%s,%s,%s,%s] = %d*%d*%d",
                                                (p_up[0],p_up[1],p_up[2],p_up[3],p_dn[0],p_dn[1],p_dn[2],p_dn[3],
                                                    parity(p_up,up_str), parity(p_dn,dn_str), val)
                                            )
                                    m += 1

    return t4


def T4_Compress(T4):
    """
    Here, we also CHECK that the input tensor has the right symmetry
    And proceed directly to COMPRESS in the following way:

        (I,J,K,L) --> I<J  and K<L

        Arranged in lexicological order, for instance

            for first pair of indices:  (1,2) comes before (1,3) before (2,3) and so on.

            the actual order is:        (0,1,0,1) -> (0,1,0,2) -> ... -> (0,1,0,NSO-1) ->
                                        (0,1,1,2) -> .... and so on...
    """
    NSO = np.size(T4,axis=0)
    m = 0
    T4_compressed = np.zeros( int( comb(NSO,4)**2 ) )

    chk = [False,False,False,False,False,False]

    for i in range(NSO):
        for j in range(i+1,NSO):
            for k in range(j+1,NSO):
                for l in range(k+1,NSO):
                    for a in range(NSO):
                        for b in range(a+1,NSO):
                            for c in range(b+1,NSO):
                                for d in range(c+1,NSO):

                                    val = T4[i,j,k,a,b,c,d]

                                    chk[0] = np.isclose(T4[j,i,k,l,a,b,c,d],-val,rtol=1e-7)
                                    chk[1] = np.isclose(T4[i,j,l,k,a,b,c,d],-val,rtol=1e-7)
                                    chk[2] = np.isclose(T4[i,k,j,l,a,b,c,d],-val,rtol=1e-7)

                                    chk[3] = np.isclose(T4[i,j,k,l,b,a,c,d],-val,rtol=1e-7)
                                    chk[4] = np.isclose(T4[i,j,k,l,a,b,d,c],-val,rtol=1e-7)
                                    chk[5] = np.isclose(T4[i,j,k,l,a,c,b,d],-val,rtol=1e-7)

                                    if all(chk):
                                        T4_compressed[m] = val
                                        m += 1
                                    else:
                                        raise(ValueError,'Incorrect Symmetry for the input Tensor',T2)
    
    return T4_compressed


"""============================================================================="""
"""         Evolve Functions to be used in ODE evolution                        """
"""============================================================================="""


def mu_evolve(Mu, TSamps, Tau, Alpha, OneH):
    """
    Function that returns the RHS of the Differential equation set up
    Inputs:     
                TSamps   ::  S - amplitudes
                            S0      ::  Constant part in Samps (1)
                            S1      ::  S1 amplitude matrix (Nso*Nso)
                            S2      ::  S2 amplitude matrix (Nso*Nso*Nso*Nso)
                Mu      ::  Chem Pot (off-set taken care of) as an independent
                            parameter without any Beta multiplication
                Tau     ::  Inverse Temperature
                Alpha   ::  Is the initial chemical potential that is used to set the number
                            Alpha = n_elec / (nso + n_elec)
                            u = 1 / np.sqrt( 1 + Alpha )
                OneH    ::  To compute the U and V on the Fly
    Returns:
                Res0    ::  R0
                Ret2    ::  R1
                Res2    ::  R2
    """

    # Number of Spin Orbitals
    Nso = np.size(OneH,axis=0)

    # Computing the HFB parameters
    U = 1/np.sqrt( 1 + np.exp(-Tau*OneH)*Alpha )
    V = np.exp( (-Tau*OneH)/2)*U*np.sqrt(Alpha)

    # lengths of different t and s tensors
    lent1 = int( Nso**2 )
    lent2 = int( comb(Nso,2)**2 )
    lent3 = int( comb(Nso,3)**2 )
    lent4 = int( comb(Nso,4)**2 )

    # Extracting amps from input
    T1 = np.reshape(TSamps[0:lent1], (Nso,Nso))
    S1 = np.reshape(TSamps[lent1:2*lent1], (Nso,Nso))
    T2 = T2_Decompress(TSamps[2*lent1:2*lent1+lent2],Nso)
    S2 = T2_Decompress(TSamps[2*lent1 + lent2:2*lent1 + 2*lent2],Nso)
    S3 = T3_Decompress(TSamps[2*(lent1 + lent2):2*(lent1 + lent2) + lent3],Nso)
    S4 = T4_Decompress(TSamps[2*(lent1 + lent2) + lent3:],Nso)
    
    dt1_dmu, ds1_dmu, ds2_dmu, ds3_dmu = covmupt2(
        T1, T2, U, V
    )

    dt1_dmu *= 1/2
    ds1_dmu *= 1/2
    ds2_dmu *= 1/2
    ds3_dmu *= 1/2
    
    #####################################################
    # TODO: This box is not for the MP - CHANGE         #
    # DERIVATIVE: exp(-T) d exp(T)/dMu                  #
    # We use the following notation                     #
    #                                                   #
    #   F00 = <0| e^(-T) d e^(T) / d Mu |0>             #
    #   F10 = <1| e^(-T) d e^(T) / d Mu |0>             #
    #   F20 = <2| e^(-T) d e^(T) / d Mu |0>             #
    #                                                   #
    #####################################################

    # First two arguments are beta, mu -- here we do not want any effect of either
    t1der,s1der = opdermp2(1, 0, OneH, T2, S2, U, V)[2:]

    dt1_dmu -= t1der 
    ds1_dmu -= s1der 

    # Reshape the array as vectors and compress to send them out.
    dt1_dmu = np.reshape(dt1_dmu,(Nso)**2)
    ds1_dmu = np.reshape(ds1_dmu,(Nso)**2)
    ds2_dmu = T2_Compress(ds2_dmu)
    ds3_dmu = T3_Compress(ds3_dmu)
    dt2_dmu = T2_Compress(T2*0)
    ds4_dmu = T4_Compress(S4*0)

    out = np.concatenate( (dt1_dmu, ds1_dmu, dt2_dmu, ds2_dmu, ds3_dmu, ds4_dmu) )
    return out

def beta_evolve(Tau, TSamps, Mu, Alpha, OneH, Eri):
    """
    Function that returns the RHS of the Differential equation set up
    Inputs:     
                Tau     ::  Imaginary Time
                TSamps  ::  T and S - amplitudes stacked together as follows
                            T is for MP1 and S for MP2
                            T1      ::  T1 amplitude matrix :: (Nso*Nso) elements 
                            S1      ::  S1 amplitude matrix :: (Nso*Nso) elements 
                            T2      ::  T2 amplitude matrix :: (Nso-choose-2)^2 elements
                            T2      ::  S2 amplitude matrix :: (Nso-choose-2)^2 elements
                Mu      ::  Chemical potential (as an independent parameter
                            and not in multiplication with Beta
                Alpha   ::  Is the initial chemical potential that is used to set the number
                            Alpha = n_elec / (nso + n_elec)
                            u = 1 / np.sqrt( 1 + Alpha )
                OneH    ::  One Dimenional Array for expression: Evals
                Eri     ::  Two Elec Integrals in MO basis
    Returns:
                Res0    ::  R0 or the Energy
                Ret2    ::  R1
                Res2    ::  R2
    """
    # Number of Spin Orbitals
    Nso = np.size(Eri,axis=0)

    # Computing the HFB parameters
    U = 1/np.sqrt( 1 + np.exp(-Tau*OneH)*Alpha )
    V = np.exp( (-Tau*OneH)/2)*U*np.sqrt(Alpha)

    # lengths of different t and s tensors
    lent1 = int( Nso**2 )
    lent2 = int( comb(Nso,2)**2 )
    lent3 = int( comb(Nso,3)**2 )
    lent4 = int( comb(Nso,4)**2 )

    # Extracting amps from input
    T1 = np.reshape(TSamps[0:lent1], (Nso,Nso))
    S1 = np.reshape(TSamps[lent1:2*lent1], (Nso,Nso))
    T2 = T2_Decompress(TSamps[2*lent1:2*lent1+lent2],Nso)
    S2 = T2_Decompress(TSamps[2*lent1 + lent2:2*lent1 + 2*lent2],Nso)
    S3 = T3_Decompress(TSamps[2*(lent1 + lent2):2*(lent1 + lent2) + lent3],Nso)
    S4 = T4_Decompress(TSamps[2*(lent1 + lent2) + lent3:],Nso)
    
    dt1_dtau, dt2_dtau, ds1_dtau, ds2_dtau, ds3_dtau, ds4_dtau = covbetapt2(
        OneH, Eri, T1, T2, S1, S2, S3, S4, U, V
    )

    #####################################################
    # TODO: This box is not for the MP - CHANGE         #
    # DERIVATIVE: exp(-T) d exp(T)/dB                   #
    # We use the following notation                     #
    #                                                   #
    #   F00 = <0| e^(-T) d e^(T) / d Beta |0>           #
    #   F10 = <1| e^(-T) d e^(T) / d Beta |0>           #
    #   F20 = <2| e^(-T) d e^(T) / d Beta |0>           #
    #                                                   #
    #####################################################
    # First two arguments are beta, mu -- here we do not want any effect of either
    t1der,s1der = opdermp2(0, 0, OneH, T2, S2, U, V)[0:2]

    dt1_dtau -= t1der 
    ds1_dtau -= s1der 
    
    # Reshape the array as vectors and compress to send them out.
    dt1_dtau = np.reshape(dt1_dtau,(Nso)**2)
    ds1_dtau = np.reshape(ds1_dtau,(Nso)**2)
    dt2_dtau = T2_Compress(dt2_dtau)
    ds2_dtau = T2_Compress(ds2_dtau)
    ds3_dtau = T3_Compress(ds3_dtau)
    ds4_dtau = T4_Compress(ds4_dtau)

    out = np.concatenate( (dt1_dtau, ds1_dtau, dt2_dtau, ds2_dtau, ds3_dtau, ds4_dtau) )
    return out

def EnergyHF(OneH, U, V):
    """
    Returns the Mean-Field Thermal Energy Expectation Value at
    a given temperature BETA and chempot MU.

    The info on BETA and MU is contained in the U and V, the HFB parameters
    """

    # Number of Spin Orbitals
    Nso = np.size(OneH,axis=0)

    # Energy from OneH
    en = np.sum(OneH*V*V)

    return en
