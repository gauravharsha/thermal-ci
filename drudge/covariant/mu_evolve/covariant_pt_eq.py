"""
    Date: Dec 2, 2018
    Modified: Dec 12, 2018
    Python Script to Carry Out the Algebraic Calculations for Ab-Initial Thermal Perturbation Theory
    This is covariant version - i.e. the reference keeps evolving
"""

from pyspark import SparkContext
# from dummy_spark import SparkContext
from drudge import *
from ThermofieldDrudge import *
# from sympy import Symbol, symbols, IndexedBase, sqrt, init_printing, KroneckerDelta
from sympy import *
from gristmill import *
import time

delK = KroneckerDelta

#Start Time
start_time = time.time()

"""-------------------------------------------------------------------------"""
"""                        Initializing the Environments                    """
"""-------------------------------------------------------------------------"""

#ctx = SparkContext('local[*]','thermal-coupled-cluster')
ctx = SparkContext()

# Thermofield Particle Hole Drudge
dr = ThermofieldDrudge(ctx)
dr2 = dr._sp_dr

nam1 = dr.names
nam2 = dr2.names
# Gen operators
d_ = nam1.d_
d_dag = nam1.d_dag
# Spin operators - for thermal state kill
c_ = nam2.c_
c_dag = nam2.c_dag

a, b, c, d, i, j, k, l, p, q, r, s = nam2.a, nam2.b, nam2.c, nam2.d, nam2.i, nam2.j, nam2.k, nam2.l, nam2.p, nam2.q, nam2.r, nam2.s

# Full Hamiltonian
full_ham = dr.simplify(dr.ham)

ham1 = dr.einst(
    d_dag[p]*d_[p]
)

full_ham = dr.simplify(ham1)

# Orig Hamiltonian Construct Time
ham_time = time.time()

"""-------------------------------------------------------------------------"""
"""   Encoding the Translation to the new representation of the             """
"""   Thermofield Spin Drudge and forming the Thermal Hamiltonian           """
"""-------------------------------------------------------------------------"""


x = IndexedBase('x')
y = IndexedBase('y')

d_def = dr.define(d_[a],
    (x[a]*c_[a,UP] + y[a]*c_dag[a,DOWN])
)
ddag_def = dr.define(d_dag[a],
    (x[a]*c_dag[a,UP] + y[a]*c_[a,DOWN])
)

tfd_defs = [
    d_def,
    ddag_def
]

ham_th = Tensor(
    dr2,
    full_ham.subst_all(tfd_defs).terms
)

ham_th = dr2.simplify(dr2.simplify(ham_th))
print('\nThermal Hamiltonian has {} terms'.format(ham_th.n_terms))

###############################################################################
# Splitting the Hamiltonian into 2, 4 and higher body terms                   #
# for ease of computation                                                     #
###############################################################################

ham_th_1body = ham_th.filter(lambda x: len(x.vecs)<=2)
print('\nThermal Hamiltonian has {} 2 or less body terms'.format(
    ham_th_1body.n_terms
    )
)
ham_th_2body = ham_th.filter(lambda x: len(x.vecs)>2)
print('\nThermal Hamiltonian has {} 3 or more body terms'.format(
    ham_th_2body.n_terms
    )
)

en_rhf = dr2.simplify(dr2.eval_phys_vev( ham_th ))

# Thermal Hamiltonian Time
hamth_time = time.time()

"""-------------------------------------------------------------------------"""
"""   Perturbation Theory - Building 1st order wavefunction correction      """
"""-------------------------------------------------------------------------"""

# First order perturbation wavefunction description
t1 = IndexedBase('t1')
t2 = IndexedBase('t2')

# Second order perturbation wavefunction description
s1 = IndexedBase('s1')
s2 = IndexedBase('s2')
s3 = IndexedBase('s3')
s4 = IndexedBase('s4')

dr2.set_symm(t2,
    Perm([1,0,2,3],NEG),
    Perm([0,1,3,2],NEG),
)

dr2.set_symm(s2,
    Perm([1,0,2,3],NEG),
    Perm([0,1,3,2],NEG),
)

dr2.set_symm(s3,
    Perm([1,0,2,3,4,5],NEG),
    Perm([0,2,1,3,4,5],NEG),
    Perm([2,0,1,3,4,5],IDENT),
    Perm([0,1,2,4,3,5],NEG),
    Perm([0,1,2,3,5,4],NEG),
    Perm([0,1,2,5,4,3],IDENT),
)

dr2.set_symm(s4,
    Perm([1,0,2,3,4,5,6,7],NEG),
    Perm([0,2,1,3,4,5,6,7],NEG),
    Perm([0,1,3,2,4,5,6,7],NEG),
    Perm([0,1,2,3,5,4,6,7],NEG),
    Perm([0,1,2,3,4,6,5,7],NEG),
    Perm([0,1,2,3,4,5,7,6],NEG),
)

T1 = dr2.einst(
    t1[a, i] * ( (c_dag[a,UP] * c_dag[i,DOWN]) )
)

T2 = dr2.einst(
    t2[a, b, i, j] * (c_dag[a,UP] * c_dag[b,UP] * c_dag[j,DOWN] \
        * c_dag[i,DOWN]) / 4 
)

S1 = dr2.einst(
    s1[a, i] * ( (c_dag[a,UP] * c_dag[i,DOWN]) )
)

S2 = dr2.einst(
    s2[a, b, i, j] * (c_dag[a,UP] * c_dag[b,UP] * c_dag[j,DOWN] \
        * c_dag[i,DOWN]) / 4 
)

S3 = dr2.einst(
    s3[a, b, c, i, j, k] * (c_dag[a,UP] * c_dag[b,UP] * c_dag[c,UP] \
        * c_dag[k,DOWN] * c_dag[j,DOWN] * c_dag[i,DOWN]) / 36
)

Tvec = dr2.simplify(T1 + T2)
Svec = dr2.simplify(S1 + S2 + S3)

#########################################################################
##      Perturbation Theory equation looks like                        ##
##      TODO: THE FOLLOWING TWO EQNS ARE WRONG HERE                    ##
##      ( d/dMu - 0.5 * Beta * N ) | 1 > = 0                           ##
##      ( d/dMu - 0.5 * Beta * N ) | 2 > = 0                           ##
##                                                                     ##
##      So here, we will generate the RHS components                   ##
##      The ``0.5 * Beta'' factor will be taken care of in odesolve    ##
#########################################################################

# First order theory
mp1_rhs_op = dr2.simplify( (ham_th) )

# Second order theory
mp2_rhs_op = dr2.simplify( (ham_th * Tvec) )

print('\n\n--------------------------------------------------------------------------------')
print('RHS operator terms evaluated')
print('--------------------------------------------------------------------------------')

# Hbar Computing Time
hbar_time = time.time()

# Projection operators
proj_t1 = (
    c_[a,DOWN] * c_[p,UP]
)

proj_t2 = (
    c_[a,DOWN] * c_[b,DOWN] * c_[q,UP] * c_[p,UP]
)

proj_s1 = (
    c_[a,DOWN] * c_[p,UP]
)

proj_s2 = (
    c_[a,DOWN] * c_[b,DOWN] * c_[q,UP] * c_[p,UP]
)

proj_s3 = (
    c_[a,DOWN] * c_[b,DOWN] * c_[c,DOWN] * c_dag[r,UP] * c_dag[q,UP] * c_dag[p,UP]
)

rt1_th = dr2.memoize(lambda: dr2.simplify(proj_t1 * (mp1_rhs_op)),'RTpa.pickle')
rt2_th = dr2.memoize(lambda: dr2.simplify(proj_t2 * (mp1_rhs_op)),'RTpqab.pickle')

rs1_th = dr2.memoize(lambda: dr2.simplify(proj_s1 * (mp2_rhs_op)),'RSpa.pickle')
rs2_th = dr2.memoize(lambda: dr2.simplify(proj_s2 * (mp2_rhs_op)),'RSpqab.pickle')
rs3_th = dr2.memoize(lambda: dr2.simplify(proj_s3 * (mp2_rhs_op)),'RSpqrabc.pickle')

# Derive the working equations by projection
rt_1body = dr2.simplify(dr2.eval_phys_vev( dr2.simplify( rt1_th )))
rt_2body = dr2.simplify(dr2.eval_phys_vev( dr2.simplify( rt2_th )))

rs_1body = dr2.simplify(dr2.eval_phys_vev( dr2.simplify( rs1_th )))
rs_2body = dr2.simplify(dr2.eval_phys_vev( dr2.simplify( rs2_th )))
rs_3body = dr2.simplify(dr2.eval_phys_vev( dr2.simplify( rs3_th )))

print('Equations obtained')

# Confirming that the contractions are unity
t1dag_t = dr2.simplify( proj_t1 * ( Tvec ) )
t2dag_t = dr2.simplify( proj_t2 * ( Tvec ) )

s1dag_s = dr2.simplify( proj_s1 * ( Svec ) )
s2dag_s = dr2.simplify( proj_s2 * ( Svec ) )
s3dag_s = dr2.simplify( proj_s3 * ( Svec ) )

t1_dag_t_exp = dr2.simplify( dr2.eval_phys_vev( t1dag_t ) )
t2_dag_t_exp = dr2.simplify( dr2.eval_phys_vev( t2dag_t ) )

s1_dag_s_exp = dr2.simplify( dr2.eval_phys_vev( s1dag_s ) )
s2_dag_s_exp = dr2.simplify( dr2.eval_phys_vev( s2dag_s ) )
s3_dag_s_exp = dr2.simplify( dr2.eval_phys_vev( s1dag_s ) )

# T2 equation time
t2eqn_time = time.time()

"""-------------------------------------------------------------------------"""
"""   Optimizing the Working Eqn and setting up the Fortran Printing        """
"""-------------------------------------------------------------------------"""

###############################################################################
# Some Pre processing to convert everything into a                            # 
# single drudge type data for gristmill                                       #
###############################################################################


rt1_eqn_in_dr = Tensor(
    dr,
    rt_1body.terms
)
rt2_eqn_in_dr = Tensor(
    dr,
    rt_2body.terms
)

rs1_eqn_in_dr = Tensor(
    dr,
    rs_1body.terms
)
rs2_eqn_in_dr = Tensor(
    dr,
    rs_2body.terms
)
rs3_eqn_in_dr = Tensor(
    dr,
    rs_3body.terms
)

# Now the actual optimization process
rt1 = IndexedBase('rt1')
rt2 = IndexedBase('rt2')

rs1 = IndexedBase('rs1')
rs2 = IndexedBase('rs2')

work_eqn = [
    dr.define(rt1[p,a], rt1_eqn_in_dr),
    dr.define(rs1[p,a], rs1_eqn_in_dr),
    dr.define(rs2[p,q,a,b], rs2_eqn_in_dr),
]

# Original Un-optimized cost
orig_cost = get_flop_cost(work_eqn, leading=True)
init_printing()
print(orig_cost)

# Optimizing
eval_seq = optimize(
    work_eqn,
    interm_fmt = 'tau{}'
)
# eval_seq = work_eqn

print(len(eval_seq))
opt_cost = get_flop_cost(eval_seq, leading=True)
print(opt_cost)

# Code Generation
fort_print = FortranPrinter(default_type='Real (Kind=8)', explicit_bounds=True)
code = fort_print.doprint(eval_seq, separate_decls=False)

with open('CovMuPT2.f90','w') as fp:
    print(code, file=fp)



"""-------------------------------------------------------------------------"""
"""                        Testing the Stuff Generated                      """
"""-------------------------------------------------------------------------"""
with dr.report('number.html','Thermal Perturbation Theory') as rep:
    rep.add('Full Hamiltonian',full_ham)
    rep.add('Thermal Ham',ham_th)
    rep.add('Tvec',Tvec)
    rep.add('Svec',Svec)
    rep.add('The Spin Rep for Tilde and non Tilde',c_[a,UP] + c_[b,DOWN])
    rep.add('MP1 1 body eqn',rt_1body)
    rep.add('MP1 2 body eqn',rt_2body)
    rep.add('MP2 1 body eqn',rs_1body)
    rep.add('MP2 2 body eqn',rs_2body)
    rep.add('MP2 3 body eqn',rs_3body)
    rep.add('Tdag1 * T',t1_dag_t_exp)
    rep.add('Tdag2 * T',t2_dag_t_exp)
    rep.add('Sdag1 * S',s1_dag_s_exp)
    rep.add('Sdag2 * S',s2_dag_s_exp)
    rep.add('Sdag3 * S',s3_dag_s_exp)
    rep.add('Symm Check',dr2.simplify(dr2.einst(
        (t2[a,b] + t2[b,a])*c_dag[a,UP]*c_dag[b,UP]
        )))
    i = 0
    # for eq in eval_seq:
    #     rep.add('Intermediate Eqn {}'.format(i),eq)
    #     i += 1

# End_time
end_time = time.time()

print('-----------------------------------------------------')
print('-----------------------------------------------------')
print('Time Profiling of the Script Run')
print('-----------------------------------------------------')
print('-----------------------------------------------------')
print('Process              Time')
print('-----------------------------------------------------')
print('Ham Init             {}'.format(ham_time-start_time))
print('Thermal Ham          {}'.format(hamth_time-ham_time))
print('Hbar Comput          {}'.format(hbar_time-hamth_time))
print('T2eq Comput          {}'.format(t2eqn_time-hbar_time))
print('Opt and Print        {}'.format(end_time-t2eqn_time))
print('Total Time           {}'.format(end_time-start_time))
print('-----------------------------------------------------')


