import numpy as np

def IntTran4(ERI,Evecs):
    """
        Function to transform the Two Body Hamlitonian Matrix
        using the transformation Evecs
    """
    ERI2 = 0*ERI

    # number of spin orbitals
    nso = np.size(Evecs,axis=0)

    # Get the last index elements
    for p in range(nso):
        for q in range(nso):
            for r in range(nso):
                for a in range(nso):
                    ERI2[p,q,r,a] = ERI[p,q,r,:]@Evecs[:,a]

    # Get the 3rd index transformed
    for p in range(nso):
        for q in range(nso):
            for a in range(nso):
                for s in range(nso):
                    ERI[p,q,a,s] = ERI2[p,q,:,s]@Evecs[:,a]

    # Get the 2nd index transformed
    for p in range(nso):
        for a in range(nso):
            for r in range(nso):
                for s in range(nso):
                    ERI2[p,a,r,s] = ERI[p,:,r,s]@Evecs[:,a]

    # Get the 1st index transformed
    for a in range(nso):
        for q in range(nso):
            for r in range(nso):
                for s in range(nso):
                    ERI[a,q,r,s] = ERI2[:,q,r,s]@Evecs[:,a]

    # all done and return
    return ERI

def IntTran2(OneH):
    """
        Diagonalize the One Body Hamiltonian in the MO basis
        and return both the diagonalized Matrix and the Evecs
    """
    evecs = np.zeros(np.shape(OneH))

    # diagonalize h1

    fock, evecs = np.linalg.eigh(OneH)

    return fock, evecs
