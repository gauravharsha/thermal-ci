import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pdb

import numpy as np
import h5py
from Main_Covar import *

global E_NUC, E_CCSD
"""-----------------------------------------------------"""
"""                 Importing the H5 Data               """
"""-----------------------------------------------------"""

fname, E_NUC, E_CCSD = (ParseInput(enuc=True,eccsd=True)[i] for i in (0,4,5))
fname = fname[0:-7]+'thermal_mp_out.h5'

f1 = h5py.File(fname,'r')

beta = np.trim_zeros(f1['beta'][:],'b')
blen = len(beta)

e_hf = f1['e_hf'][:blen]
e_mp1 = f1['e_mp1'][:blen]
e_mp2 = f1['e_mp2'][:blen]
chem_pot = f1['chem_pot'][:blen]
t1_rms = f1['t1rms'][:blen]
t2_rms = f1['t2rms'][:blen]
s1_rms = f1['s1rms'][:blen]
s2_rms = f1['s2rms'][:blen]


f1.close()


"""-----------------------------------------------------"""
"""                 Creating the Plots                  """
"""-----------------------------------------------------"""

fout = fname[0:-7]+'_'

plt.figure()
pdb.set_trace()
plt.plot(beta,e_hf,label='HF')
plt.plot(beta,e_mp1,label='MP1')
plt.plot(beta,e_mp2,label='MP2')
plt.legend()
plt.xlabel(r'$\beta$')
plt.ylabel('Energy ')
plt.xlim(beta[0],beta[-1])
# plt.hlines(E_CCSD,beta[0],beta[-1])
plt.savefig(fout+'energy.pdf')
plt.close()

plt.figure()
plt.plot(beta,t1_rms,label='RMS t1')
plt.plot(beta,s1_rms,label='RMS s1')
plt.legend()
plt.xlabel(r'$\beta$')
plt.ylabel('RMS t1 and s1')
plt.xlim(beta[0],beta[-1])
plt.savefig(fout+'t1_amp.pdf')
plt.close()

plt.figure()
plt.plot(beta,t2_rms,label='RMS t2')
plt.plot(beta,s2_rms,label='RMS s2')
plt.legend()
plt.xlabel(r'$\beta$')
plt.ylabel('RMS t2 and s2')
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

