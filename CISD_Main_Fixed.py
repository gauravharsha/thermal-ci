import numpy as np
from scipy.integrate import odeint, ode
from scipy.misc import comb
import matplotlib.pyplot as plt
import h5py, time

from FixBetaCISD import *
from FixMuCISD import *
from FixEnAndOvlp import *

from IntTran import *
from EvolveFuncsFixed import *


#
# GLOBAL VARIABLES
#

mu_step_0 = +5e-2
len_t1 = 0
len_t2 = 0


#
# INPUT / OUTPUT RELATED FUNCTIONS
#


def ParseInput(enuc=False):
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
    global cis_bool
    cis_bool = bool(int(line[pos:].strip()))

    line = fin.readline()
    pos = line.find(':') + 1
    global cisd_bool
    cisd_bool = bool(int(line[pos:].strip()))

    if enuc:
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



#
# ODE SOLVER FUNCTIONS FOR VARIOUS 
# ORDERS OF THE CISD THEORY                              
#


def ci_beta_integrate(integrator, bspan, y_in, x, y, h1, eri):
    """
    Function to perform the task of integrating the BETA evolution equations for zeroth, first
    and second order.

    Inputs:
        integrator      ::  The integrator, an instance of scipy.integrate.ode
                            --this integrator will decide whether it is the cis or cisd
                                level of theory being integrated
                            --the function is therefore general and transparent to the level
                                of theory.
        bspan           ::  Array of length 2 containing initial and final beta value
        y_in            ::  Initial condition at Singles / or Singles and Doubles
        x, y            ::  Fixe-reference's HFB parameters
        h1              ::  One electron integrals
        eri             ::  Two electron integrals

    Outputs:
        fin_ci          ::  Final solution at the respective input order
    """

    ### Check for any possible errors in the input
    if not isinstance(integrator,ode):
        raise ValueError('The argument',integrator,'  needs to be an instance of ',ode)
    if len(bspan) != 2:
        raise ValueError('Inappropriate length of the input array ',bspan)

    # Set the initial conditions
    integrator.set_initial_value(y_in, bspan[0]).set_f_params(x, y, h1, eri)
    fin_ci = integrator.integrate(bspan[1])

    return fin_ci


def mu_find_and_integrate(integrator, mu_in, y_in, nelec, beta, x, y, h1):
    """
    Function that will first findt the chemical potential bracket and then, integrate to the
    correct mu and return the final TAMPS.

    """

    ### Check for any possible errors in the input
    if not isinstance(integrator,ode):
        raise ValueError('The argument',integrator,'  needs to be an instance of ',ode)

    # Global variables
    global mu_step_0
    global ntol

    nso = len(h1)

    # find the initial expectation value, i.e. <N>
    t0 = y_in[0]
    t1 = np.reshape( y_in[1:1+nso**2], (nso,nso) )
    t2 = T2_Decompress( y_in[1+(nso**2):], nso )

    num, ov = fixenandovlp(
        np.ones(nso), t2*0, t0, t1, t2, x, y
    )
    num /= ov

    print('\t\tNumber of particles after the beta evolution = {}'.format(num))

    ndiff_sgn = np.sign( num - nelec )
    ndiff_mag = np.abs( num - nelec )

    yf = y_in
    mu_f = mu_in
    
    # if the number is already converged, then there is no need to do any of the following
    #   and hence we keep an 'if' statement; if the condition evaluates to FALSE, then
    #   the outputs will be yf, mu_f

    mu_0 = mu_in

    if ndiff_mag > ntol:

        mu_step = mu_step_0

        # Obtain the bracket to perform BISECTION

        mu_1 = mu_0
        mu_2 = mu_0 + mu_step

        yf1 = y_in
        yf2 = y_in

        count = 0
        sp_count = 0
        
        while np.sign(num - nelec) == ndiff_sgn:
            count += 1
            if count > 1000:
                print('Could not find the bracket after 6000 steps')
                count = 0
                exit()
                break

            mu_span = [mu_1, mu_2]
            
            # Define the solver

            yf1 = yf2

            integrator.set_initial_value(yf1, mu_span[0])
            integrator.set_f_params(beta, x, y, h1)

            # Evolve to the mu_f
            yf2 = integrator.integrate(mu_span[1])

            # Extract the amplitudes
            t0 = yf2[0]
            t1 = np.reshape(yf2[1:1+nso**2],(nso,nso))
            t2 = T2_Decompress(yf2[1+nso**2:],nso)

            # Evaluate the Number Expectation
            num, ov = fixenandovlp(
                np.ones(nso), t2*0, t0, t1, t2, x, y
            )

            num /= ov

            # Finer grid if we are closer to nelec
            val = np.abs(num - nelec) - ndiff_mag
            if (val>0):
                if val<1e-1:
                    sp_count += 1
                else:
                    mu_step *= -1
                
                if sp_count >= 10:
                    mu_step *= -1
                    sp_count = 0

            ndiff_mag = np.abs(num - nelec)

            # Set up for next iteration
            mu_1 = mu_2
            mu_2 = mu_2 + mu_step

            if np.abs(num - nelec) <= ntol:
                break

            print('\t\t\tStart value of Mu = {}'.format(mu_span[0]))
            print('\t\t\tEnd value of Mu = {}'.format(mu_span[1]))
            print('\t\t\tNumber of particles after evolution = {}'.format(num))
            print('\t\t\t----------------------------------------------\n')

        print('Bracket found between mu = {} and mu = {}'.format(mu_span[0],mu_span[1]))

        # Bisection bracket
        mu_bisect = mu_span
        mu_mid = mu_bisect[1]

        ndiff_sgn2 = np.sign(num - nelec)
        ndiff_sgn1 = -ndiff_sgn2

        while np.abs(num - nelec)>ntol:
            # Set the initial condition for the mu_solver
            integrator.set_initial_value(yf1, mu_bisect[0]).set_f_params(beta, x, y, h1)

            # Evolve the ODE to mid point of the bracket
            mu_mid = np.mean(mu_bisect)
            yf_mid = integrator.integrate(mu_mid)

            # Extract the amplitudes
            t0 = yf_mid[0]
            t1 = np.reshape(yf_mid[1:1+nso**2],(nso,nso))
            t2 = T2_Decompress(yf_mid[1+nso**2:],nso)

            # Compute the number and update the bracket
            num, ov = fixenandovlp(
                h1*0+1, t2*0, t0, t1, t2, x, y
            )

            num /= ov

            if np.sign( num - nelec ) == ndiff_sgn1:
                mu_bisect[0] = mu_mid
                yf1 = yf_mid
            else:
                mu_bisect[1] = mu_mid
                yf2 = yf_mid

        print('Bisection converges to mu = {}'.format(mu_mid))

        # Now that we have found the MU_MID, we can use one-shot evolution to avoid any added errors
        mu_f = mu_mid
        integrator.set_initial_value(yf, mu_0).set_f_params(beta, x, y, h1)

        yf = integrator.integrate(mu_f)

    return yf, mu_f



#
# MAIN
# 

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
    en, ov = fixenandovlp(
        h1, eri, t0, t1, t2, x, y 
    )
    en /= ov

    # Beta Grid
    beta_0 = 0
    beta_step = beta_f / (beta_pts-1)
    beta_grid = np.array([beta_0])

    # Alpha Grid
    mu_0 = np.log(alpha)
    mu_step_0 = -1e-2
    mu_f = mu_step_0

    n_data = 1

    #################################################################
    #               SETTING UP THE INITIAL VALUE PROB               #
    #    Compute the Hartree fock energy and define the y-stacks    #
    #    which will eventually be involved in computation           #
    #################################################################

    # initial condition for Tamps at beta = 0 = mu
    y0_cis = np.concatenate(([t0], t1_vec, t2_vec))
    y0_cisd = np.concatenate(([t0], t1_vec, t2_vec))

    # Data to be output at BETA GRID POINT
    e_cis = np.zeros(n_data)
    e_cisd = np.zeros(n_data)

    mucis = np.zeros(n_data)
    chem_cis = 0.0 # Dummy for mucis
    mucisd = np.zeros(n_data)
    chem_cisd = 0.0 # Dummy for mucisd

    n_exp = np.zeros(n_data)

    t0_cis = np.ones(n_data)
    t1_cis = np.zeros(n_data)
    t2_cis = np.zeros(n_data)

    t0_cisd = np.ones(n_data)
    t1_cisd = np.zeros(n_data)
    t2_cisd = np.zeros(n_data)

    e_cis[0] = e_nuc + en/ov
    e_cisd[0] = e_cis[0]

    print('Thermal CIS Corrected Energy at T = Inf is :',e_cis[0])
    print('Thermal CISD Corrected Energy at T = Inf is :',e_cisd[0])
    print('-------------------------------------------------------\n')


    # Time stamp to record that all the initialization process is completed
    init_time = time.time()

    # Print the first things
    num, ov = fixenandovlp(
        np.ones(nso), eri*0, t0, t1, t2, x, y
    )
    num /= ov

    # Confirm and Print the Number of particles at Beta = 0
    print('Number of Particles = {}'.format(num))
    
    #################################################################
    #                   FILE OUTPUT HANDLE CREATOR                  #
    #   Creating the h5py file headers and update the files after   #
    #   each iteration - so that if the code blows up we know       #
    #   where it went wrong.                                        #
    #################################################################

    fout = fn[0:-7] + 'thermal_cisd_out.h5'
    print('Writing output to {}'.format(fout))

    fp1 = h5py.File(fout,'w')

    dsets = ['beta','e_cis', 'e_cisd','mu_cis','mu_cisd','t0cis','t1icis','t2cis',\
        't0cisd','t1cisd','t2cisd']

    # Create all but the ystack data sets

    dset_list = createh5(fp1, dsets, beta_pts, attrs)

    vals = [
        0, e_cis[-1], e_cisd[-1], chem_cis, chem_cisd, t0_cis[-1], t1_cis[-1],\
        t2_cis[-1], t0_cisd[-1], t1_cisd[-1], t2_cisd[-1]
    ]

    updateh5(dset_list, vals, 0)

    #################################################################
    #                   TIME AND CHEMPOT EVOLUTION                  #
    #   The actual evolution takes place in this section of code.   #
    #   First evolve in Beta Direction to a grid point and then     #
    #   evolve in the Mu Direction to fix the Number of Particles   #
    #################################################################

    beta_solver_cis = ode(beta_cis).set_integrator('vode',method='bdf',rtol=deqtol)
    beta_solver_cisd = ode(beta_evolve).set_integrator('vode',method='bdf',rtol=deqtol)

    mu_solver_cis = ode(mu_cis).set_integrator('vode',method='bdf',rtol=deqtol)
    mu_solver_cisd = ode(mu_evolve).set_integrator('vode',method='bdf',rtol=deqtol)

    j = 1

    while beta_grid[j-1]<=beta_f:

        beta_2 = beta_grid[j-1] + beta_step
        b_span = [beta_grid[j-1], beta_2]

        ############################### 
        # 1: Beta evolution 
        ############################### 

        print('\t\t\t Beta-Evolution for CI-Singles')
        yf_cis = ci_beta_integrate(
            beta_solver_cis, b_span, y0_cis, x, y, h1, eri
        )
        print('\t\t\t Beta-Evolution for CI-Singles and Doubles')
        yf_cisd = ci_beta_integrate(
            beta_solver_cisd, b_span, y0_cisd, x, y, h1, eri
        )

        y0_cis = yf_cis
        y0_cisd = yf_cisd

        print('\n\t\t\tNew Beta = {}'.format(beta_2))

        ############################### 
        # 2: Mu or ChemPot evolution 
        ############################### 

        print('\t\t\t Mu-Evolution for CI-Singles')
        yf_cis, chem_cis = mu_find_and_integrate(
            mu_solver_cis, chem_cis, y0_cis, n_elec, b_span[1], x, y, h1
        )

        print('\t\t\t Mu-Evolution for CI-Singles and Doubles')
        yf_cisd, chem_cisd = mu_find_and_integrate(
            mu_solver_cisd, chem_cisd, y0_cisd, n_elec, b_span[1], x, y, h1
        )

        ###############################################
        # 3: Update H5 and set up for the next loop
        ###############################################

        # Setting up for next loop
        y0_cis = yf_cis
        y0_cisd = yf_cisd

        # Computing the quantities of interest
        beta_grid = np.append(beta_grid,b_span[1])

        # for cis
        t0 = y0_cis[0]
        t1 = np.reshape( y0_cis[1:1+nso**2], (nso, nso) )
        t2 = T2_Decompress( y0_cis[1+nso**2:], nso)

        en, ov = fixenandovlp(
            h1, eri, t0, t1, t2, x, y
        )
        en /= ov

        e_cis = np.append( e_cis, en )
        t0_cis = np.append( t0_cis, t0 )
        t1_cis = np.append( t1_cis, np.sqrt( np.mean( t1**2 ) ) )
        t2_cis = np.append( t2_cis, np.sqrt( np.mean( t2**2 ) ) )

        # for cisd
        t0 = y0_cisd[0]
        t1 = np.reshape( y0_cisd[1:1+nso**2], (nso, nso) )
        t2 = T2_Decompress( y0_cisd[1+nso**2:], nso)

        en, ov = fixenandovlp(
            h1, eri, t0, t1, t2, x, y
        )
        en /= ov

        e_cisd = np.append( e_cisd, en )
        t0_cisd = np.append( t0_cisd, t0 )
        t1_cisd = np.append( t1_cisd, np.sqrt( np.mean( t1**2 ) ) )
        t2_cisd = np.append( t2_cisd, np.sqrt( np.mean( t2**2 ) ) )

        # chem pots
        mucis = np.append(mucis, chem_cis)
        mucisd = np.append(mucisd, chem_cisd)

        vals = [
            beta_grid[-1], e_cis[-1], e_cisd[-1], mucis[-1], mucisd[-1],\
            t0_cis[-1], t1_cis[-1], t2_cis[-1], t0_cisd[-1], t1_cisd[-1], t2_cisd[-1],
        ]

        updateh5(dset_list, vals, j)

        print('At beta = {}'.format(b_span[-1]))
        print('CIS mu = {}'.format(chem_cis))
        print('CISD mu = {}'.format(chem_cisd))
        print('CIS En = {}'.format(e_cis[-1]))
        print('CISD En = {}'.format(e_cisd[-1]))
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


