import numpy as np
from scipy.special import comb
from tfdcisd.iofuncs import IOps
from tfdcisd.odefuncs import Evolution, CompressT2, DecompressT2
from tfdcisd.odefuncs import ci_alpha_evolve, ci_beta_evolve


#########################################################################
#                                                                       #
#   Testing the Beta and Alpha ODE drivers and related functiosn        #
#                                                                       #
#########################################################################

def test_compress_decompress_t2():

    #
    # Test the function which compresses and de-compresses T2-like
    # amplitudes based on the symmetries
    #

    # size
    nso = 10

    # t2
    x2 = np.random.rand(nso, nso, nso, nso)

    # anti-symmetrize
    t2 = x2 - np.einsum('qprs->pqrs', x2) \
            - np.einsum('pqsr->pqrs', x2) \
        + np.einsum('qpsr->pqrs', x2)

    t2 /= 4.0

    # Compress
    t2_compr = CompressT2(t2)

    # Check the length of the compressed t2
    assert np.shape(t2_compr) == (int(comb(nso, 2)**2), )

    # DeCompress
    t2_decompr = DecompressT2(t2_compr, nso)

    # Check the shape of the decompressed t2
    assert np.shape(t2_decompr) == (nso, nso, nso, nso)

    # de compressed version should be same as original t2
    assert np.max(np.abs(t2_decompr - t2)) <= 5e-8


def test_evolution_class_attributes():

    #
    # Test the attributes of the `Evolution` class
    #

    # First read the file and the input and do the formalities
    iops = IOps(inp_file='TestInput')
    eigs, h1, eri, attrs = iops.loadHDF()

    fug = iops.n_elec/(iops.nso-iops.n_elec)
    # Define the evolution class instance
    evol = Evolution(inp_file='TestInput')

    assert evol.beta_step == 1e-2
    assert evol.alpha_step == 0.1
    assert evol.nso == iops.nso
    assert evol.n_elec == iops.n_elec
    assert evol.fug == fug
    assert evol.h1.all() == h1.all()
    assert evol.eri.all() == eri.all()
    assert evol.deqtol == iops.deqtol
    assert evol.ntol == iops.ntol

    # amplitudes -- before I set anything
    assert len(evol.ci_amps) == int(evol.nso**2 + 1 + comb(evol.nso, 2)**2)

    # amplitudes -- after changing
    ci_amps = np.random.rand(1+int(iops.nso**2 + comb(iops.nso, 2)**2))
    evol.setAmps(ci_amps=ci_amps)

    assert evol.ci_amps.all() == ci_amps.all()


def test_ci_beta_evolve():

    #
    # Test the integration driver functions
    #

    # Define the evolution class instance
    evol = Evolution(inp_file='TestInput')

    nso = evol.nso
    nso = evol.nso
    fug = evol.fug
    eigs = evol.eigs
    h1 = evol.h1
    eri = evol.eri

    # Construct t1 and t2 = use Tom's MasterCode
    nocc = 6
    t1 = np.loadtxt('T1AMP')
    t1_pre = np.zeros((nso, nso))
    t1_pre[:nocc, nocc:] = t1[:, :]
    t1 = np.einsum('ia->ai', t1_pre)

    t2 = np.loadtxt('T2AMP')
    t2_pre = eri*0
    m = 0
    for i in range(nocc):
        for j in range(nocc):
            for a in range(nocc):
                t2_pre[i, j, nocc+a, nocc:] = t2[m, :]
                m += 1
    t2 = np.einsum('ijab->abij', t2_pre)

    # Other parameters needed for residuals
    y = np.zeros(nso)
    for i in range(nocc):
        y[i] = 1.0

    # Compress and concatenate
    t0 = 0.0
    t1 = np.reshape(t1, int(nso**2))
    t2 = CompressT2(t2)
    ci_amps = np.concatenate(([t0], t1, t2))

    beta_in = 10.0
    alpha_in = beta_in*h1[nocc]

    # Check outputs from CI beta evolution
    yout = ci_beta_evolve(beta_in, ci_amps, alpha_in, fug, eigs, h1, eri)

    r1 = np.reshape(yout[1:1+nso**2], (nso, nso))
    r2 = DecompressT2(yout[1+nso**2:], nso)

    assert np.shape(r1) == (nso, nso)
    assert np.shape(r2) == (nso, nso, nso, nso)

    # Check outputs from CI alpha evolution
    yout = ci_alpha_evolve(alpha_in, ci_amps, beta_in, fug, eigs)

    r1 = np.reshape(yout[1:1+nso**2], (nso, nso))
    r2 = DecompressT2(yout[1+nso**2:], nso)

    assert np.shape(r1) == (nso, nso)
    assert np.shape(r2) == (nso, nso, nso, nso)


def test_evolution_class_beta_integration():

    #
    # Test the attributes of the `Evolution` class
    #

    # Define the evolution class instance
    evol = Evolution(inp_file='TestInput')

    # Do Beta evolution
    evol.DoBetaIntegration()
    ci_amps0 = evol.ci_amps

    assert np.max(np.abs(ci_amps0)) > 0


def test_evolution_class_alpha_integration():

    #
    # Test the attributes of the `Evolution` class
    #

    # Define the evolution class instance
    evol = Evolution(inp_file='TestInput')
    ci_amps0 = evol.ci_amps

    # Do Beta evolution - should not do anything!!
    evol.BisectionAndAlphaIntegrate()

    # amplitudes -- after evolving -- check if it evolve at all
    assert np.max(np.abs(evol.ci_amps - ci_amps0)) < 5e-8
