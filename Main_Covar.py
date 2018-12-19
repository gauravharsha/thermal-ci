import numpy as np
from scipy.integrate import odeint, ode
from scipy.misc import comb
import matplotlib.pyplot as plt
import h5py, time
import pdb

# from CovBetaPT2 import *
from PyPTLib import *
from CovMuPT2 import *
from OpDerMP2 import *
# from MPEnergyCov import *

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

def createh5(fp, dsets, max_size):
    zarr = np.zeros(max_size)
    dsetlist = []
    for ds in dsets:
        dsetlist.append(
            fp.create_dataset(ds, data=zarr)
        )
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
    t1_vec = np.zeros(nso**2)
    s1_vec = np.zeros(nso**2)
    t2_vec = np.zeros( int( comb(nso,2)**2 ) )
    s2_vec = np.zeros( int( comb(nso,2)**2 ) )
    s3_vec = np.zeros( int( comb(nso,3)**2 ) )
    s4_vec = np.zeros( int( comb(nso,4)**2 ) )

    # These are the actual t1 and t2 matrices
    t1 = np.reshape(t1_vec,(nso,nso))
    t2 = T2_Decompress(t2_vec,nso)
    s1 = np.reshape(s1_vec,(nso,nso))
    s2 = T2_Decompress(s2_vec,nso)
    s3 = T3_Decompress(s3_vec,nso)
    s4 = T4_Decompress(s4_vec,nso)

    # Hartree Fock Energy at BETA = 0
    E_hf = e_nuc + EnergyHF(
        h1, x, y 
    )

    e_1, e_2, o_1, o_2, o_12 = mpenergycov(
        h1, eri, t1, t2, s1, s2, s3, s4, x, y
    )

    # Beta Grid
    beta_0 = 0
    beta_step = beta_f / (beta_pts-1)
    beta_grid = np.array([beta_0])

    # Alpha Grid
    mu_0 = np.log(alpha)
    mu_step_0 = +1e-2
    mu_f = mu_step_0

    n_data = 1

    #################################################################
    #               SETTING UP THE INITIAL VALUE PROB               #
    #    Compute the Hartree fock energy and define the y-stacks    #
    #    which will eventually be involved in computation           #
    #################################################################

    # initial condition for Tamps at beta = 0 = mu
    # NOTE: there are two singles and doubles parameters here, one each for
    #       the first and second order PT wavefunctions
    #       The order in which they are stacked is as follows:
    #               [t1,s1,t2,s2]
    y0 = np.concatenate((t1_vec,s1_vec,t2_vec,s2_vec,s3_vec,s4_vec))
    ystack = np.zeros((n_data,np.size(y0)))
    ystack[0,:] = y0

    # Data to be output at BETA GRID POINT
    e_hf = np.zeros(n_data)
    e_mp1 = np.zeros(n_data)
    e_mp2 = np.zeros(n_data)

    ov_mp1 = np.zeros(n_data)
    ov_mp2 = np.zeros(n_data)
    ov_mp12 = np.zeros(n_data)

    mu_cc = np.zeros(n_data)

    n_exp = np.zeros(n_data)

    t1_rms = np.zeros(n_data)
    t2_rms = np.zeros(n_data)

    s1_rms = np.zeros(n_data)
    s2_rms = np.zeros(n_data)

    e_hf[0] = E_hf
    e_mp1[0] = (E_hf + e_1)/(1 + o_1)
    e_mp2[0] = (E_hf + e_1 + e_2)/(1 + o_1 + o_12 + o_2)
    n_exp[0] = n_elec

    print('Thermal Hartree Fock Energy at T = Inf is :',E_hf)
    print('Thermal MP1 Energy at T = Inf is :',e_1)
    print('Thermal MP2 Energy at T = Inf is :',e_2)
    print('-------------------------------------------------------\n')


    # Time stamp to record that all the initialization process is completed
    init_time = time.time()

    # Print the first things
    num = EnergyHF(
        h1*0+1, x, y
    )

    # Confirm and Print the Number of particles at Beta = 0
    print('Number of Particles = {}'.format(num))
    
    #################################################################
    #                   FILE OUTPUT HANDLE CREATOR                  #
    #   Creating the h5py file headers and update the files after   #
    #   each iteration - so that if the code blows up we know       #
    #   where it went wrong.                                        #
    #################################################################

    fout = fn[0:-7] + 'thermal_mp_out.h5'
    print('Writing output to {}'.format(fout))

    fp1 = h5py.File(fout,'w')

    dsets = ['beta','e_hf','e_mp1','e_mp2','chem_pot','t1rms','s1rms',\
        't2rms','s2rms','num']

    # Create all but the ystack data sets
    dset_list = createh5(fp1, dsets, beta_pts)
    dset_list.append(
        fp1.create_dataset('tvals',data = np.zeros( (beta_pts, len(y0)) ))
    )

    vals = [
        0, e_hf[-1], e_mp1[-1], e_mp2[-1], mu_cc[-1],t1_rms[-1],\
        s1_rms[-1], t2_rms[-1], s2_rms[-1], n_exp[-1], ystack[-1]
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
    len_t1 = int(nso**2)
    len_t2 = int(comb(nso,2)**2)
    len_t3 = int(comb(nso,3)**2)
    len_t4 = int(comb(nso,4)**2)

    # XXX: Checkpoint number 1
    # pdb.set_trace()


    while beta_grid[j-1]<beta_f:

        beta_2 = beta_grid[j-1] + beta_step
        b_span = [beta_grid[j-1], beta_2]

        ############################### 
        # 1: Update the HFB coefficients
        ############################### 

        x = 1/np.sqrt(1 + np.exp( -b_span[0]*h1 )*alpha )
        y = np.exp( ( -b_span[0]*h1 )/2 )*x*np.sqrt(alpha)

        # Check the number of particles before evolution
        num_hf = EnergyHF(
            h1*0+1, x, y
        )

        n_1, n_2, o_1, o_2, o_12 = mpenergycov(
            h1*0+1, eri*0, t1, t2, s1, s2, s3, s4, x, y
        )

        num = (num_hf + n_1 + n_2)/(1 + o_1 + o_2 + o_12)
        print('\t\t\tNumber of particles before evolution = {}'.format(num))

        ############################### 
        # 2: Beta evolution 
        ############################### 

        beta_solver.set_initial_value(y0, b_span[0]).set_f_params(mu_cc[-1],alpha,h1, eri)

        # Integrate the BETA_SOLVE
        yf = beta_solver.integrate(b_span[1])

        # Extract the values
        t1 = np.reshape(yf[0:len_t1],(nso,nso))
        s1 = np.reshape(yf[len_t1:2*len_t1],(nso,nso))
        t2 = T2_Decompress(yf[2*len_t1:2*len_t1+len_t2],nso)
        s2 = T2_Decompress(yf[2*len_t1+len_t2:2*(len_t1+len_t2)],nso)
        s3 = T3_Decompress(yf[2*(len_t1+len_t2):2*(len_t1+len_t2)+len_t3],nso)
        s4 = T4_Decompress(yf[2*(len_t1+len_t2)+len_t3:2*(len_t1+len_t2)+len_t3+len_t4],nso)

        # Check number after beta evolution
        x = 1/np.sqrt(1 + np.exp( -b_span[1]*h1 )*alpha )
        y = np.exp( ( -b_span[1]*h1 )/2 )*x*np.sqrt(alpha)
        
        num_hf = EnergyHF(
            h1*0+1, x, y
        )
        n_1, n_2, o_1, o_2, o_12 = mpenergycov(
            h1*0+1, eri*0, t1, t2, s1, s2, s3, s4, x, y
        )
        num = (num_hf + n_1 + n_2)/(1 + o_1 + o_2 + o_12)

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
                if count > 500:
                    print('Could not find the bracket after 1000 steps')
                    count = 0
                    exit()
                    break

                mu_span = [mu_1, mu_2]
                
                # Define the solver
                yf1 = yf2
                mu_solver.set_initial_value(yf1, mu_span[0]).set_f_params(b_span[1], alpha, h1)

                # Evolve to the mu_f
                yf2 = mu_solver.integrate(mu_span[1])

                # Extract the amplitudes
                t1 = np.reshape(yf2[0:len_t1],(nso,nso))
                s1 = np.reshape(yf2[len_t1:2*len_t1],(nso,nso))
                t2 = T2_Decompress(yf2[2*len_t1:2*len_t1+len_t2],nso)
                s2 = T2_Decompress(yf2[2*len_t1+len_t2:],nso)
                s3 = T3_Decompress(yf2[2*(len_t1+len_t2):2*(len_t1+len_t2)+len_t3],nso)
                s4 = T4_Decompress(yf2[2*(len_t1+len_t2)+len_t3:2*(len_t1+len_t2)+len_t3+len_t4],nso)

                # Evaluate the Number Expectation
                num_hf = EnergyHF(
                    h1*0+1, x, y
                )
                n_1, n_2, o_1, o_2, o_12 = mpenergycov(
                    h1*0+1, eri*0, t1, t2, s1, s2, s3, s4, x, y
                )
                num = (num_hf + n_1 + n_2)/(1 + o_1 + o_2 + o_12)

                # Finer grid if we are closer to n_elec
                val = np.abs(num - n_elec) - ndiff_mag
                if (val>0):
                    if val<1e-1:
                        sp_count += 1
                    else:
                        mu_step *= -1
                    
                    if sp_count >= 25:
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
                t1 = np.reshape(yf_mid[0:len_t1],(nso,nso))
                s1 = np.reshape(yf_mid[len_t1:2*len_t1],(nso,nso))
                t2 = T2_Decompress(yf_mid[2*len_t1:2*len_t1+len_t2],nso)
                s2 = T2_Decompress(yf_mid[2*len_t1+len_t2:],nso)
                s3 = T3_Decompress(yf_mid[2*(len_t1+len_t2):2*(len_t1+len_t2)+len_t3],nso)
                s4 = T4_Decompress(yf_mid[2*(len_t1+len_t2)+len_t3:2*(len_t1+len_t2)+len_t3+len_t4],nso)

                # Compute the number and update the bracket
                num_hf = EnergyHF(
                    h1*0+1, x, y
                )
                n_1, n_2, o_1, o_2, o_12 = mpenergycov(
                    h1*0+1, eri*0, t1, t2, s1, s2, s3, s4, x, y
                )
                num = (num_hf + n_1 + n_2)/(1 + o_1 + o_2 + o_12)


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


        t1 = np.reshape(yf[0:len_t1],(nso,nso))
        s1 = np.reshape(yf[len_t1:2*len_t1],(nso,nso))
        t2 = T2_Decompress(yf[2*len_t1:2*len_t1+len_t2],nso)
        s2 = T2_Decompress(yf[2*len_t1+len_t2:],nso)
        s3 = T3_Decompress(yf[2*(len_t1+len_t2):2*(len_t1+len_t2)+len_t3],nso)
        s4 = T4_Decompress(yf[2*(len_t1+len_t2)+len_t3:2*(len_t1+len_t2)+len_t3+len_t4],nso)

        # Setting up for next loop
        mu_0 = mu_f
        y0 = yf
        ystack = np.row_stack((ystack,yf))

        # Computing the quantities of interest
        beta_grid = np.append(beta_grid,beta_2)

        e_hf = np.append(
            e_hf,
            e_nuc + EnergyHF(h1, x, y)
        )

        e_1, e_2, o_1, o_2, o_12 = mpenergycov(
            h1, eri, t1, t2, s1, s2, s3, s4, x, y
        )

        e_mp1 = np.append(e_mp1, e_hf[-1] + e_1) / (1+o_1)
        e_mp2 = np.append(e_mp2, e_mp1[-1] + e_2) / (1+o_1+o_12+o_2)

        ov_mp1 = np.append(ov_mp1, o_1)
        ov_mp2 = np.append(ov_mp2, o_2)
        ov_mp12 = np.append(ov_mp2, o_12)
        print('The value of the <1|2> = {}'.format(o_12))

        mu_cc = np.append(mu_cc,mu_f)

        n_exp = np.append(n_exp,num)
        t1_rms = np.append(
            t1_rms,
            np.sqrt(np.mean(t1**2))
        )
        t2_rms = np.append(
            t2_rms,
            np.sqrt(np.mean(t2**2))
        )
        s1_rms = np.append(
            s1_rms,
            np.sqrt(np.mean(s1**2))
        )
        s2_rms = np.append(
            s2_rms,
            np.sqrt(np.mean(s2**2))
        )

        vals = [
            beta_grid[-1], e_hf[-1], e_mp1[-1], e_mp2[-1], mu_cc[-1],t1_rms[-1],\
            s1_rms[-1], t2_rms[-1], s2_rms[-1], n_exp[-1], ystack[-1]
        ]

        updateh5(dset_list, vals, j)

        print('At beta = {}'.format(b_span[-1]))
        print('mu = {}'.format(mu_f))
        print('N_elec = {}'.format(n_exp[j]))
        print('ECC = {}'.format(e_mp2[j]))
        print('----------------------------------------------\n')
        j += 1

    ode_time = time.time()


    # TODO: If the code works, then remove this commented section
    # dset1 = fp1.create_dataset('beta',data=beta_grid)
    # for tup in attrs:
    #     dset1.attrs[tup[0]] = tup[1]
    # 
    # dset2 = fp1.create_dataset('e_hf',data=e_hf)
    # for tup in attrs:
    #     dset2.attrs[tup[0]] = tup[1]

    # dset2a = fp1.create_dataset('e_mp1',data=e_mp1)
    # for tup in attrs:
    #     dset2a.attrs[tup[0]] = tup[1]

    # dset2b = fp1.create_dataset('e_mp2',data=e_mp2)
    # for tup in attrs:
    #     dset2b.attrs[tup[0]] = tup[1]
    #     
    # dset3 = fp1.create_dataset('chem_pot',data=mu_cc)
    # for tup in attrs:
    #     dset3.attrs[tup[0]] = tup[1]

    # dset4 = fp1.create_dataset('t1rms',data=t1_rms)
    # for tup in attrs:
    #     dset4.attrs[tup[0]] = tup[1]

    # dset5 = fp1.create_dataset('s1rms',data=s1_rms)
    # for tup in attrs:
    #     dset5.attrs[tup[0]] = tup[1]

    # dset6 = fp1.create_dataset('t2rms',data=t2_rms)
    # for tup in attrs:
    #     dset6.attrs[tup[0]] = tup[1]

    # dset6a = fp1.create_dataset('s2rms',data=s2_rms)
    # for tup in attrs:
    #     dset6a.attrs[tup[0]] = tup[1]

    # dset7 = fp1.create_dataset('num',data=n_exp)
    # for tup in attrs:
    #     dset7.attrs[tup[0]] = tup[1]

    # dset8 = fp1.create_dataset('tvals',data=ystack)
    # for tup in attrs:
    #     dset8.attrs[tup[0]] = tup[1]

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


