import numpy as np
from scipy.misc import comb
from CovBetaPT2 import *
from CovMuPT2 import *
from OpDerMP2 import *

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



def mu_evolve(Mu, Tamps, Tau, Alpha, OneH):
    """
    Function that returns the RHS of the Differential equation set up
    Inputs:     
                Tamps   ::  S - amplitudes
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

    # Breaking input into the needed pieces
    lent1 = int( Nso**2 )
    lent2 = int( comb(Nso,2)**2 )
    T1 = np.reshape(Tamps[0:lent1], (Nso,Nso))
    S1 = np.reshape(Tamps[lent1:2*lent1], (Nso,Nso))
    T2 = T2_Decompress(Tamps[2*lent1:2*lent1+lent2],Nso)
    S2 = T2_Decompress(Tamps[2*lent1 + lent2:],Nso)
    
    dt1_dmu, ds1_dmu, ds2_dmu = covmupt2(
        T1, T2, U, V
    )
    dt2_dmu = 0*ds2_dmu

    dt1_dmu *= 1/2
    dt2_dmu *= 1/2
    ds1_dmu *= 1/2
    ds2_dmu *= 1/2
    
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
    dt2_dmu = T2_Compress(dt2_dmu)
    ds2_dmu = T2_Compress(ds2_dmu)

    out = np.concatenate( (dt1_dmu, ds1_dmu, dt2_dmu, ds2_dmu) )
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

    # Breaking input into the needed pieces
    lent1 = int( Nso**2 )
    lent2 = int( comb(Nso,2)**2 )
    T1 = np.reshape(TSamps[0:lent1], (Nso,Nso))
    S1 = np.reshape(TSamps[lent1:2*lent1], (Nso,Nso))
    T2 = T2_Decompress(TSamps[2*lent1:2*lent1+lent2],Nso)
    S2 = T2_Decompress(TSamps[2*lent1 + lent2:],Nso)
    
    dt1_dtau, dt2_dtau, ds1_dtau, ds2_dtau = covbetapt2(
        OneH, Eri, T1, T2, S1, S2, U, V
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

    out = np.concatenate( (dt1_dtau, ds1_dtau, dt2_dtau, ds2_dtau) )
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
