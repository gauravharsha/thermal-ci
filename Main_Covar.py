import numpy as np
from scipy.integrate import odeint, ode
from scipy.misc import comb
import matplotlib.pyplot as plt
import h5py, time
import pdb

from CovBetaPT import *
from CovMuPT import *
from EnAndOvlp import *
# from BetaDerPT import *
# from MuDerPT import *

from IntTran import *
from EvolveFuncsCovar import *

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

def createh5(fp, dsets, max_size, attrs):
    zarr = np.zeros(max_size)
    dsetlist = []
    for ds in dsets:
        dsetlist.append(
            fp.create_dataset(ds, data=zarr)
        )
        for tup in attrs:
            dsetlist[-1].attrs[tup[0]] = tup[1]
        
    return dsetlist

def updateh5(dsets, vals, i_at):
    for i in range(len(dsets)):
        dsets[i][i_at] = vals[i]

def main():
    start_time = time.time()
    
    #################################################################
    #                       READ INPUT                              #
    #       Get the data first - including the attributes           #
    #################################################################

    print('\n\n----------------------------------------------\n')
    print('Reading Input')
    print('----------------------------------------------\n')

    fn, n_elec, beta_f, beta_pts, e_nuc = ParseInput(enuc=True)

    h1_in, eri_in, attrs  = loadHDF(fname=fn)
    nso = np.size(h1_in,axis=0)
    
    input_time = time.time()

    print('Number of Spin Orbitals:',nso)
    print('-------------------------------------------------------\n')

    #################################################################
    #                   INTEGRAL TRANSFORMATION                     #
    #       Transforming the One and Two Electron Integral          #
    #       so that the OneH or h1 is diagonal                      #
    #################################################################

    # Note that h1 is a 1D array of the eigenvalues
    h1, evecs = IntTran2(h1_in)
    eri = IntTran4(eri_in,evecs)

    #################################################################
    #               INITIAL HFB and OTHER  PARAMETERS               #
    #       Initializing the important x and y HFB parameters       #
    #       and other arrays for a better idea further              #
    #################################################################

    # HFB parameters in the initial thermal vacuum reference at Beta = 0
    # alpha = exp( beta * mu ) 
    # when beta -> 0 and mu -> inf
    alpha = n_elec/(nso-n_elec)

    x = np.ones(nso)/np.sqrt(1+alpha)
    y = np.ones(nso)*np.sqrt(alpha)/np.sqrt(1+alpha)
    
    # Initializing the amplitude vectors - only the unique matrix
    # elements in the t1 and t2 matrices
    len_t1 = int(nso**2)
    len_t2 = int(comb(nso,2)**2)

    t0 = 1.0
    t1_vec = np.zeros(nso**2)
    t2_vec = np.zeros( int( comb(nso,2)**2 ) )

    # These are the actual t1 and t2 matrices
    t1 = np.reshape(t1_vec,(nso,nso))
    t2 = T2_Decompress(t2_vec,nso)

    # Hartree Fock Energy at BETA = 0
    E_hf, o_hf = enandovlp(
        h1, eri, t0, t1, t2, x, y 
    )

    # Beta Grid
    beta_0 = 0
    beta_step = beta_f / (beta_pts-1)
    beta_grid = np.array([beta_0])

    # Alpha Grid
    mu_0 = np.log(alpha)
    mu_step_0 = +1e-1
    mu_f = mu_step_0

    n_data = 1

    #################################################################
    #               SETTING UP THE INITIAL VALUE PROB               #
    #    Compute the Hartree fock energy and define the y-stacks    #
    #    which will eventually be involved in computation           #
    #################################################################

    # initial condition for Tamps at beta = 0 = mu
    y0 = np.concatenate(([t0], t1_vec, t2_vec))
    ystack = np.zeros((n_data,np.size(y0)))
    ystack[0,:] = y0

    # Data to be output at BETA GRID POINT
    e_hf = np.zeros(n_data)
    e_mp2 = np.zeros(n_data)

    mu_cc = np.zeros(n_data)

    n_exp = np.zeros(n_data)

    t0_rms = np.ones(n_data)
    t1_rms = np.zeros(n_data)
    t2_rms = np.zeros(n_data)

    e_hf[0] = e_nuc + E_hf/o_hf
    e_mp2[0] = e_hf[0]

    n_exp[0] = n_elec

    print('Thermal Hartree Fock Energy at T = Inf is :',e_hf[0])
    print('Thermal MP2 Corrected Energy at T = Inf is :',e_mp2[0])
    print('-------------------------------------------------------\n')


    # Time stamp to record that all the initialization process is completed
    init_time = time.time()

    # Print the first things
    num_hf, o_hf = enandovlp(
        h1*0+1, eri*0, t0, t1, t2, x, y
    )
    num_hf /= o_hf

    # Confirm and Print the Number of particles at Beta = 0
    print('Number of Particles = {}'.format(num_hf))
    
    #################################################################
    #                   FILE OUTPUT HANDLE CREATOR                  #
    #   Creating the h5py file headers and update the files after   #
    #   each iteration - so that if the code blows up we know       #
    #   where it went wrong.                                        #
    #################################################################

    fout = fn[0:-7] + 'thermal_mp_out.h5'
    print('Writing output to {}'.format(fout))

    fp1 = h5py.File(fout,'w')

    dsets = ['beta','e_hf','e_mp2','chem_pot','t0rms','t1rms','t2rms','num']

    # Create all but the ystack data sets

    dset_list = createh5(fp1, dsets, beta_pts, attrs)

    dset_tvals = fp1.create_dataset('tvals',data=np.zeros(( beta_pts, len(y0)) ))
    dset_list.append( dset_tvals )

    for tup in attrs:
        dset_tvals.attrs[tup[0]] = tup[1]

    # XXX: Another de-bugger
    # pdb.set_trace()

    vals = [
        0, e_hf[-1], e_mp2[-1], mu_cc[-1], t0_rms[-1], t1_rms[-1],\
        t2_rms[-1], n_exp[-1], ystack[-1]
    ]

    updateh5(dset_list, vals, 0)

    #################################################################
    #                   TIME AND CHEMPOT EVOLUTION                  #
    #   The actual evolution takes place in this section of code.   #
    #   First evolve in Beta Direction to a grid point and then     #
    #   evolve in the Mu Direction to fix the Number of Particles   #
    #################################################################

    beta_solver = ode(beta_evolve).set_integrator('vode',method='bdf',rtol=deqtol)
    mu_solver = ode(mu_evolve).set_integrator('vode',method='bdf',rtol=deqtol)
    # beta_solver = ode(beta_evolve).set_integrator('dopri5',rtol=deqtol)
    # mu_solver = ode(mu_evolve).set_integrator('dopri5',rtol=deqtol)

    j = 1
    # XXX: Checkpoint number 1
    # pdb.set_trace()


    while beta_grid[j-1]<=beta_f:

        beta_2 = beta_grid[j-1] + beta_step
        b_span = [beta_grid[j-1], beta_2]

        ############################### 
        # 1: Update the HFB coefficients
        ############################### 

        x = 1/np.sqrt(1 + np.exp( -b_span[0]*h1 )*alpha )
        y = np.exp( ( -b_span[0]*h1 )/2 )*x*np.sqrt(alpha)

        # Check the number of particles before evolution
        num, ov = enandovlp(
            h1*0+1, eri*0, t0, t1, t2, x, y
        )

        num /= ov
        print('\t\t\tNumber of particles before evolution = {}'.format(num))

        ############################### 
        # 2: Beta evolution 
        ############################### 

        beta_solver.set_initial_value(y0, b_span[0]).set_f_params(alpha, h1, eri)

        # Integrate the BETA_SOLVE
        yf = beta_solver.integrate(b_span[1])

        # Extract the values
        t0 = yf[0]
        t1 = np.reshape(yf[1:1+len_t1],(nso,nso))
        t2 = T2_Decompress(yf[1+len_t1:],nso)

        # Check number after beta evolution
        x = 1/np.sqrt(1 + np.exp( -b_span[1]*h1 )*alpha )
        y = np.exp( ( -b_span[1]*h1 )/2 )*x*np.sqrt(alpha)
        
        num, ov = enandovlp(
            h1*0+1, eri*0, t0, t1, t2, x, y
        )

        num /= ov

        print('\n\t\t\tNew Beta = {}'.format(beta_2))
        print('\t\t\tNumber of particles after evolution = {}'.format(num))
        print('\t\t\tStart value of Mu = {}'.format(mu_0))
        print('\t\t\tEnd value of Mu = {}'.format(mu_f))
        print('\t\t\t----------------------------------------------\n')

        # XXX: Checkpoint number 2
        # pdb.set_trace()

        # Check if the number difference became negative
        ndiff_sgn = np.sign(num - n_elec)
        ndiff_mag = np.abs(num - n_elec)

        ############################### 
        # 3: Mu or ChemPot evolution 
        ############################### 

        # NOTE: The commenting of the next 2 lines is purely experimental - we may resume it back
        # if j==1:
        #    mu_step = mu_step_0
        if ndiff_mag > ntol:
            mu_step = mu_step_0

            # PRE-RUN to obtain the bracket to perform BISECTION

            mu_1 = mu_0
            mu_2 = mu_0 + mu_step

            yf1 = yf
            yf2 = yf

            count = 0
            sp_count = 0
            
            while np.sign(num - n_elec) == ndiff_sgn:
                count += 1
                if count > 1000:
                    print('Could not find the bracket after 6000 steps')
                    count = 0
                    pdb.set_trace()
                    exit()
                    break

                mu_span = [mu_1, mu_2]
                
                # Define the solver
                yf1 = yf2
                mu_solver.set_initial_value(yf1, mu_span[0]).set_f_params(b_span[1], alpha, h1)

                # Evolve to the mu_f
                yf2 = mu_solver.integrate(mu_span[1])

                # Extract the amplitudes
                t0 = yf2[0]
                t1 = np.reshape(yf2[1:1+len_t1],(nso,nso))
                t2 = T2_Decompress(yf2[1+len_t1:],nso)

                # Evaluate the Number Expectation
                num, ov = enandovlp(
                    h1*0+1, eri*0, t0, t1, t2, x, y
                )

                num /= ov

                # Finer grid if we are closer to n_elec
                val = np.abs(num - n_elec) - ndiff_mag
                if (val>0):
                    if val<1e-1:
                        sp_count += 1
                    else:
                        mu_step *= -1
                    
                    if sp_count >= 10:
                        mu_step *= -1
                        sp_count = 0

                ndiff_mag = np.abs(num - n_elec)

                # if np.abs(num - n_elec) < 0.05:
                #     mu_step /= 2

                # Set up for next iteration
                mu_1 = mu_2
                mu_2 = mu_2 + mu_step

                if np.abs(num - n_elec) <= ntol:
                    break

                print('\t\t\tStart value of Mu = {}'.format(mu_span[0]))
                print('\t\t\tEnd value of Mu = {}'.format(mu_span[1]))
                print('\t\t\tNumber of particles after evolution = {}'.format(num))
                print('\t\t\t----------------------------------------------\n')

                # XXX: Checkpoint number 3
                # pdb.set_trace()

            print('Bracket found (unless the previous line says not found)')

            # Bisection bracket
            mu_bisect = mu_span
            mu_mid = mu_bisect[1]

            ndiff_sgn2 = np.sign(num - n_elec)
            ndiff_sgn1 = -ndiff_sgn2

            while np.abs(num - n_elec)>ntol:
                # Set the initial condition for the mu_solver
                mu_solver.set_initial_value(yf1, mu_bisect[0]).set_f_params(b_span[1], alpha, h1)

                # Evolve the ODE to mid point of the bracket
                mu_mid = np.mean(mu_bisect)
                yf_mid = mu_solver.integrate(mu_mid)

                # Extract the amplitudes
                t0 = yf_mid[0]
                t1 = np.reshape(yf_mid[1:1+len_t1],(nso,nso))
                t2 = T2_Decompress(yf_mid[1+len_t1:],nso)

                # Compute the number and update the bracket
                num, ov = enandovlp(
                    h1*0+1, eri*0, t0, t1, t2, x, y
                )

                num /= ov


                if np.sign( num - n_elec ) == ndiff_sgn1:
                    mu_bisect[0] = mu_mid
                    yf1 = yf_mid
                else:
                    mu_bisect[1] = mu_mid
                    yf2 = yf_mid

            print('Bisection converges to mu = {}'.format(mu_mid))

            # Now that we have found the MU_MID, we can use one-shot evolution to avoid any added errors
            mu_f = mu_mid
            mu_solver.set_initial_value(yf, mu_0).set_f_params(b_span[1], alpha, h1)

            yf = mu_solver.integrate(mu_f)


        t0 = yf[0]
        t1 = np.reshape(yf[1:1+len_t1],(nso,nso))
        t2 = T2_Decompress(yf[1+len_t1:],nso)

        # Setting up for next loop
        mu_0 = mu_f
        y0 = yf
        ystack = np.row_stack((ystack,yf))

        # Computing the quantities of interest
        beta_grid = np.append(beta_grid,beta_2)

        E_hf, o_hf = enandovlp(
            h1, eri, t0, t1*0, t2*0, x, y
        )
        e_hf = np.append(
            e_hf,
            e_nuc + E_hf/o_hf
        )

        e_2, o_2 = enandovlp(
            h1, eri, t0, t1, t2, x, y
        )

        # XXX: Checkpoint
        # print('energy correction at order 1: {}'.format(e_1))
        # print('energy correction at order 2: {}'.format(e_2))
        # print('overlap correction at order 3: {}'.format(o_3))
        # print('overlap correction at order 4: {}'.format(o_4))
        # pdb.set_trace()

        e_mp2 = np.append(e_mp2, e_nuc + e_2/o_2 )

        mu_cc = np.append(mu_cc,mu_f)

        n_exp = np.append(n_exp,num)

        t0_rms = np.append(
            t0_rms,
            t0
        )
        t1_rms = np.append(
            t1_rms,
            np.sqrt(np.mean(t1**2))
        )
        t2_rms = np.append(
            t2_rms,
            np.sqrt(np.mean(t2**2))
        )

        vals = [
            beta_grid[-1], e_hf[-1], e_mp2[-1], mu_cc[-1],t0_rms[-1],t1_rms[-1],\
            t2_rms[-1], n_exp[-1], ystack[-1]
        ]

        updateh5(dset_list, vals, j)

        print('At beta = {}'.format(b_span[-1]))
        print('mu = {}'.format(mu_f))
        print('N_elec = {}'.format(n_exp[-1]))
        print('ECC = {}'.format(e_mp2[-1]))
        print('----------------------------------------------\n')

        j += 1
        if j == beta_pts:
            break


    #################################################################
    #                                                               #
    #       CLOSING THE CODE AND PRINT TIMING                       #
    #                                                               #
    #################################################################

    ode_time = time.time()
    fp1.close()
    end_time = time.time()

    print('----------------------------------------------')
    print('End of Program and Time Profile')
    print('----------------------------------------------')
    print('----------------------------------------------')
    print('Process              Time')
    print('----------------------------------------------')
    print('Reading Input        {}'.format(input_time-start_time))
    print('Init Variables       {}'.format(init_time-input_time))
    print('Solve ODE for grid   {}'.format(ode_time-init_time))
    print('Writing the Output   {}'.format(end_time-ode_time))
    print('Total Time           {}'.format(end_time-start_time))
    print('----------------------------------------------')
    print('----------------------------------------------\n')



if __name__ == '__main__':
    main()


