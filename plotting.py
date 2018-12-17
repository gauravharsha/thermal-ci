import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import numpy as np
import h5py
from Main_Covar import *

global E_NUC, E_CCSD
"""-----------------------------------------------------"""
"""                 Importing the H5 Data               """
"""-----------------------------------------------------"""

fname, E_NUC, E_CCSD = (ParseInput(enuc=True,eccsd=True)[i] for i in (0,4,5))
fname = fname[0:-7]+'thermal_ccsd_out.h5'

f1 = h5py.File(fname,'r')

beta = f1['beta'][:]
ecc = f1['ecc'][:]
chem_pot = f1['chem_pot'][:]
t0 = f1['t0'][:]
t1_rms = f1['t1rms'][:]
t2_rms = f1['t2rms'][:]


f1.close()


"""-----------------------------------------------------"""
"""                 Creating the Plots                  """
"""-----------------------------------------------------"""

fout = fname[0:-7]+'thermal_ccsd_'

plt.figure()
plt.plot(beta,ecc+E_NUC)
plt.xlabel(r'$\beta$')
plt.ylabel('Energy ')
plt.xlim(beta[0],beta[-1])
plt.hlines(E_CCSD,beta[0],beta[-1])
plt.savefig(fout+'energy.pdf')
plt.close()

plt.figure()
plt.plot(beta,t0)
plt.xlabel(r'$\beta$')
plt.ylabel('t0')
plt.xlim(beta[0],beta[-1])
plt.savefig(fout+'t0_amp.pdf')
plt.close()

plt.figure()
plt.plot(beta,t1_rms)
plt.xlabel(r'$\beta$')
plt.ylabel('RMS t1')
plt.xlim(beta[0],beta[-1])
plt.savefig(fout+'t1_amp.pdf')
plt.close()

plt.figure()
plt.plot(beta,t2_rms)
plt.xlabel(r'$\beta$')
plt.ylabel('RMS t2')
plt.xlim(beta[0],beta[-1])
plt.savefig(fout+'t2_amp.pdf')
plt.close()

plt.figure()
plt.plot(beta,chem_pot)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$ \bar{\mu} $')
plt.xlim(beta[0],beta[-1])
plt.savefig(fout+'chempot.pdf')
plt.close()

