import numpy as np
from h5py import File
from os.path import expanduser, expandvars


#
# INPUT / OUTPUT RELATED FUNCTIONS
#

class IOps:
    """
        Input / Output Operations

        When we initialize the class, it automatically looks for the Input File
        and updates the parameters

    """

    input_file = 'Input'
    output_dsets = []  # Eventually it will be a list of H5PY data sets

    def __init__(self, inp_file=None):

        if inp_file is not None:
            self.input_file = inp_file

        # Initialize the variables
        self.fn = ''
        self.nso = 0
        self.n_elec = 0.0
        self.beta_f = 0.0
        self.beta_pts = 0
        self.ntol = 1e-4
        self.deqtol = 1e-6
        self.e_nuc = 0.0

        # Read the Input File and update the parameters
        self.ParseInput(input_file=self.input_file)

        return

    def ParseInput(self, input_file='Input'):
        """
            Read the Input File
            Here is a sample input file highlighting the structure and order:

                  HDF File containing integrals   :   hub_6x1_u2_data.h5
                  Number of electrons             :   6
                  Final Beta value                :   15
                  Number of Beta points           :   151
                  Precision in Number             :   1e-6
                  Precision in Diff Eqn Evolve    :   1e-8
                  Nuclear Repulsion               :   0.0
        """

        fin = open(input_file)
        line = fin.readline()

        # Input HDF file name
        while line[0:3] != 'HDF':
            line = fin.readline()
        pos = line.find(':') + 1
        self.fn = line[pos:].strip()
        self.fn = expandvars(self.fn)
        self.fn = expanduser(self.fn)

        # Number of electrons
        line = fin.readline()
        pos = line.find(':') + 1
        self.n_elec = int(line[pos:].strip())

        # Final Value of Beta evolution
        line = fin.readline()
        pos = line.find(':') + 1
        self.beta_f = float(line[pos:].strip())

        # Number of Grid points in beta evolution
        line = fin.readline()
        pos = line.find(':') + 1
        self.beta_pts = int(line[pos:].strip())

        # Tolerance on the <N> expectation value
        # (Bisection convergence criterion)
        line = fin.readline()
        pos = line.find(':') + 1
        self.ntol = float(line[pos:].strip())

        # Tolerance in the ODE evolutions
        line = fin.readline()
        pos = line.find(':') + 1
        self.deqtol = float(line[pos:].strip())

        # Nuclear repulsion energy (if there should be any)
        line = fin.readline()
        pos = line.find(':') + 1
        self.e_nuc = float(line[pos:].strip())

        return

    def loadHDF(self):
        """
            Function to read the integrals from the HDF file
        """

        # Read the HDF file
        f1 = File(self.fn, 'r')

        # Get the dataset for h1 and its attributes
        h1_dset = f1['h1']
        h1 = h1_dset[:]

        attributes = list(h1_dset.attrs)
        attr_list1 = tuple((atr, h1_dset.attrs[atr]) for atr in attributes)
        attr_list1 = sorted(attr_list1)

        # Get the dataset for eri and its attributes
        eri_dset = f1['eri']
        eri = eri_dset[:]

        attributes = list(eri_dset.attrs)
        attr_list2 = tuple((atr, h1_dset.attrs[atr]) for atr in attributes)
        attr_list2 = sorted(attr_list2)

        # Get the dataset for HF eigenvalues and its attributes
        eigs_dset = f1['eigs']
        eigs = eigs_dset[:]

        attributes = list(eigs_dset.attrs)
        attr_list3 = tuple((atr, eigs_dset.attrs[atr]) for atr in attributes)
        attr_list3 = sorted(attr_list1)

        # Compare the attributes
        if (attr_list1 != attr_list2) and (attr_list1 != attr_list3):
            raise ValueError(
                'The attributes of the Energy and ERI data in the file ',
                self.fn,
                ' do not match'
            )

        # Print the attribute data
        print('\n')
        print('==============================================================')
        print('\tMolecule Parameters from Input File')
        print('==============================================================')
        for attr_tup in attr_list1:
            print(attr_tup[0], ' = ', attr_tup[1])

        # Set the number of Spin orbitals
        self.nso = np.size(h1, axis=0)

        # Close the HDF file
        f1.close()

        # return the integrals and the attribute list
        return eigs, h1, eri, attr_list1

    def createh5(self, fp, dsets, max_size, attrs):
        """
        Function to create and initialize the Output HDF file
        Need: instead of dumping all the data at the end of the program run,
              this function provides a way to continuously update the output
        """
        zarr = np.zeros(max_size)
        for ds in dsets:
            self.output_dsets.append(
                fp.create_dataset(ds, data=zarr)
            )
            for tup in attrs:
                self.output_dsets[-1].attrs[tup[0]] = tup[1]

        return

    def updateh5(self, vals, i_at):
        """
        Function to update the data sets in the output file
        """
        for i in range(len(self.output_dsets)):
            self.output_dsets[i][i_at] = vals[i]

        return
