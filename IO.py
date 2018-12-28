import numpy as np
import h5py

def ParseInput(enuc=False,eccsd=False):
    fin = open('Input')
    line = fin.readline()

    # Reach the line with the file name
    while line[0:3] != 'HDF':
        line = fin.readline()
    pos = line.find(':') + 1
    fname = line[pos:].strip()

    line = fin.readline()
    pos = line.find(':') + 1
    n_elec = int(line[pos:].strip())

    line = fin.readline()
    pos = line.find(':') + 1
    beta_f = float(line[pos:].strip())

    line = fin.readline()
    pos = line.find(':') + 1
    beta_pts = int(line[pos:].strip())

    line = fin.readline()
    pos = line.find(':') + 1
    global ntol
    ntol = float(line[pos:].strip())

    line = fin.readline()
    pos = line.find(':') + 1
    global deqtol
    deqtol = float(line[pos:].strip())

    line = fin.readline()
    pos = line.find(':') + 1
    global E_NUC
    E_NUC = float(line[pos:].strip())

    line = fin.readline()
    pos = line.find(':') + 1
    global E_CCSD
    E_CCSD = float(line[pos:].strip())

    if enuc:
        if eccsd:
            return fname, n_elec, beta_f, beta_pts, E_NUC, E_CCSD
        else:
            return fname, n_elec, beta_f, beta_pts, E_NUC
    else:
        return fname, n_elec, beta_f, beta_pts


def loadHDF(fname):
    f1 = h5py.File(fname,'r')
    en_dset = f1['h1']
    mo_en = en_dset[:]
    attributes = list(en_dset.attrs)
    attr_list1 = tuple((atr,en_dset.attrs[atr]) for atr in attributes)
    attr_list1 = sorted(attr_list1)

    eri_dset = f1['eri']
    mo_eri = eri_dset[:]
    attributes = list(eri_dset.attrs)
    attr_list2 = tuple((atr,en_dset.attrs[atr]) for atr in attributes)
    attr_list2 = sorted(attr_list2)

    if attr_list1 != attr_list2:
        raise ValueError('The attributes of the Energy and ERI data in the file ',fname,' do not match')
    
    # print('The shape of the MO Energy matrix is ',np.shape(mo_en))
    # print('The shape of the MO ERI matrix is ',np.shape(mo_eri))
    print('\n\n\n----------------------------')
    print('    Molecule Parameters')
    print('----------------------------')
    for attr_tup in attr_list1:
        print(attr_tup[0],' = ',attr_tup[1])
    print('----------------------------')

    f1.close()
    
    return mo_en, mo_eri, attr_list1
