"""
    Date: Dec 4, 2018
    Modified: Dec 25, 2018
    Python Script to Carry Out the Algebraic Calculations for Ab-Initial Thermal Perturbation Theory
    This is covariant version - i.e. the reference keeps evolving

    Here, we compute the derivatives of the operator parts of the cluster operator
"""

from pyspark import SparkContext
# from dummy_spark import SparkContext
from drudge import *
from ThermofieldDrudge import *
from sympy import Symbol, symbols, IndexedBase, sqrt, init_printing, KroneckerDelta
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

# # Full Hamiltonian
# full_ham = dr.simplify(dr.ham)
# 
e0 = IndexedBase('e0')
mu = Symbol(r'\mu')
Beta = Symbol(r'\beta')

"""-------------------------------------------------------------------------"""
"""   Encoding the Translation to the new representation of the             """
"""   Thermofield Spin Drudge and forming the Thermal Hamiltonian           """
"""-------------------------------------------------------------------------"""


x = IndexedBase('x')
y = IndexedBase('y')


"""-------------------------------------------------------------------------"""
"""   Perturbation Theory - LHS == Derivatives of the T Operator            """
"""-------------------------------------------------------------------------"""

# First order perturbation wavefunction description
t1 = IndexedBase('t1')
t2 = IndexedBase('t2')

dr2.set_symm(t2,
    Perm([1,0,2,3],NEG),
    Perm([0,1,3,2],NEG),
)

T1 = dr2.einst(
    t1[a, b] * ( (c_dag[a,UP] * c_dag[b,DOWN]) )
)

T2 = dr2.einst(
    t2[a, b, c, d] * (c_dag[a,UP] * c_dag[b,UP] * c_dag[d,DOWN] * c_dag[c,DOWN]) / 4 
)

Tvec = dr2.simplify(T1 + T2)

# Derivatives of the T1 and S1 operator part
dT1_dBeta = dr2.einst( t1[a, b] *
    ( 
        ( e0[a] - mu ) * x[a]*y[a] * c_[a,DOWN] * c_dag[b,DOWN] / 2 -\
        ( e0[b] - mu ) * x[b]*y[b] * c_dag[a,UP] * c_[b,UP] / 2
    )
)
dT1_dBeta = dr2.simplify( dT1_dBeta )

dT1_dMu = dr2.einst( t1[a, b] *
    ( 
        ( -Beta ) * x[a]*y[a] * c_[a,DOWN] * c_dag[b,DOWN] / 2 +\
        ( Beta ) * x[b]*y[b] * c_dag[a,UP] * c_[b,UP] / 2
    )
)
dT1_dMu = dr2.simplify( dT1_dMu )

# Derivatives of the T2 and S2 operator part
dT2_dBeta = dr2.einst( t2[a, b, c, d] *
    ( 
        ( e0[a] - mu ) * x[a]*y[a] * c_[a,DOWN] * c_dag[b,UP] * c_dag[d,DOWN] * c_dag[c,DOWN] / 2 + 
        ( e0[b] - mu ) * x[b]*y[b] * c_dag[a,UP] * c_[b,DOWN] * c_dag[d,DOWN] * c_dag[c,DOWN] / 2 - 
        ( e0[d] - mu ) * x[d]*y[d] * c_dag[a,UP] * c_dag[b,UP] * c_[d,UP] * c_dag[c,DOWN] / 2 - 
        ( e0[c] - mu ) * x[c]*y[c] * c_dag[a,UP] * c_dag[b,UP] * c_dag[d,DOWN] * c_[c,UP] / 2 
    )/4
)
dT2_dBeta = dr2.simplify( dT2_dBeta )

dT2_dMu = dr2.einst( t2[a, b, c, d] *
    ( 
        ( -Beta ) * x[a]*y[a] * c_[a,DOWN] * c_dag[b,UP] * c_dag[d,DOWN] * c_dag[c,DOWN] / 2 +
        ( -Beta ) * x[b]*y[b] * c_dag[a,UP] * c_[b,DOWN] * c_dag[d,DOWN] * c_dag[c,DOWN] / 2 +
        ( Beta ) * x[d]*y[d] * c_dag[a,UP] * c_dag[b,UP] * c_[d,UP] * c_dag[c,DOWN] / 2 +
        ( Beta ) * x[c]*y[c] * c_dag[a,UP] * c_dag[b,UP] * c_dag[d,DOWN] * c_[c,UP] / 2
    )/4
)
dT2_dMu = dr2.simplify( dT2_dMu )

# Combining the T' derivatives

dT_dBeta = dr2.simplify( dT1_dBeta + dT2_dBeta )
dT_dMu = dr2.simplify( dT1_dMu + dT2_dMu )

# Various Tbar Computing Time
Tder_bar_time = time.time()

"""-------------------------------------------------------------------------"""
"""   Obtaining the Necessary equations out of the Operators                """
"""-------------------------------------------------------------------------"""

# Projection operators
proj_t1 = (
    c_[q,DOWN] * c_[p,UP]
)

proj_t2 = (
    c_[r,DOWN] * c_[s,DOWN] * c_[q,UP] * c_[p,UP]
)

# Projected with proj_t1

dT_dBeta_1body = dr2.simplify( proj_t1 * dT_dBeta )
dT_dMu_1body = dr2.simplify( proj_t1 * dT_dMu )

# Projected with proj_t2

dT_dBeta_2body = dr2.simplify( proj_t2 * dT_dBeta )
dT_dMu_2body = dr2.simplify( proj_t2 * dT_dMu )

# Calculating the vev

dT_dBeta_0body_eqn = dr2.simplify( dr2.eval_phys_vev( dT_dBeta ) )
dT_dBeta_1body_eqn = dr2.simplify( dr2.eval_phys_vev( dT_dBeta_1body ) )
dT_dBeta_2body_eqn = dr2.simplify( dr2.eval_phys_vev( dT_dBeta_2body ) )

dT_dMu_0body_eqn = dr2.simplify( dr2.eval_phys_vev( dT_dMu ) )
dT_dMu_1body_eqn = dr2.simplify( dr2.eval_phys_vev( dT_dMu_1body ) )
dT_dMu_2body_eqn = dr2.simplify( dr2.eval_phys_vev( dT_dMu_2body ) )

print('Equations obtained')

t1dag_t = dr2.simplify( proj_t1 * ( Tvec ) )
t2dag_t = dr2.simplify( proj_t2 * ( Tvec ) )

t1_dag_t_exp = dr2.simplify( dr2.eval_phys_vev( t1dag_t ) )
t2_dag_t_exp = dr2.simplify( dr2.eval_phys_vev( t2dag_t ) )


# T2 equation time
t2eqn_time = time.time()

# Some Pre processing to convert everything into a single drudge type data for gristmill

dT_dBeta_0body_eqn_in_dr = Tensor(
    dr,
    dT_dBeta_0body_eqn.terms
)
dT_dBeta_1body_eqn_in_dr = Tensor(
    dr,
    dT_dBeta_1body_eqn.terms
)
dT_dBeta_2body_eqn_in_dr = Tensor(
    dr,
    dT_dBeta_2body_eqn.terms
)

dT_dMu_0body_eqn_in_dr = Tensor(
    dr,
    dT_dMu_0body_eqn.terms
)
dT_dMu_1body_eqn_in_dr = Tensor(
    dr,
    dT_dMu_1body_eqn.terms
)
dT_dMu_2body_eqn_in_dr = Tensor(
    dr,
    dT_dMu_2body_eqn.terms
)

# Now the actual optimization process
tbeta0 = Symbol('tb0')
tbeta1 = IndexedBase('tb1')
tbeta2 = IndexedBase('tb2')
tmu0 = Symbol('tm0')
tmu1 = IndexedBase('tm1')
tmu2 = IndexedBase('tm2')

work_beta = [
    dr.define(tbeta0, dT_dBeta_0body_eqn_in_dr),
    dr.define(tbeta1[p,q], dT_dBeta_1body_eqn_in_dr),
]

work_mu = [
    dr.define(tmu0, dT_dMu_0body_eqn_in_dr),
    dr.define(tmu1[p,q], dT_dMu_1body_eqn_in_dr),
]

#####################################################################
#       BETA equations - code optimization                          #
#####################################################################

# Original Un-optimized cost
orig_cost = get_flop_cost(work_beta, leading=True)
init_printing()
print('Original cost of equations')
print(orig_cost)

# optimizing
try:
    eval_beta = optimize(
        work_beta,
        interm_fmt = 'tau{}'
    )
except ValueError:
    print('\n\nOne or more Beta equations were zeros on the RHS\n\n')
    eval_beta = work_beta


opt_cost = get_flop_cost(eval_beta, leading=True)
print('Optimized cost of equations')
print(opt_cost)
print('Length of eval seq')
print(len(eval_beta))

# Code Generation
fort_print = FortranPrinter(default_type='Real (Kind=8)', explicit_bounds=True)
code = fort_print.doprint(eval_beta, separate_decls=False)

# ein_print = EinsumPrinter()
# code = ein_print.doprint(eval_beta, separate_decls=False)

with open('BetaDerPT.f90','w') as fp:
    print(code, file=fp)


#####################################################################
#       MU equations - code optimization                            #
#####################################################################

# Original Un-optimized cost
orig_cost = get_flop_cost(work_mu, leading=True)
init_printing()
print('Original cost of equations')
print(orig_cost)

# optimizing
try:
    eval_mu = optimize(
        work_mu,
        interm_fmt = 'tau{}'
    )
except ValueError:
    print('\n\nOne or more Mu equations were zeros on the RHS\n\n')
    eval_mu = work_mu


opt_cost = get_flop_cost(eval_mu, leading=True)
print('Optimized cost of equations')
print(opt_cost)
print('Length of eval seq')
print(len(eval_mu))

# Code Generation
fort_print = FortranPrinter(default_type='Real (Kind=8)', explicit_bounds=True)
code = fort_print.doprint(eval_mu, separate_decls=False)

# code = ein_print.doprint(eval_mu, separate_decls=False)

with open('MuDerPT.f90','w') as fp:
    print(code, file=fp)


"""-------------------------------------------------------------------------"""
"""                        Testing the Stuff Generated                      """
"""-------------------------------------------------------------------------"""
with dr.report('derivatives.html','Covariant Derivatives') as rep:
    rep.add('Tvec',Tvec)
    # rep.add('dT1_dBeta',dT1_dBeta)
    # rep.add('dT2_dBeta',dT2_dBeta)
    # rep.add('dS1_dBeta',dS1_dBeta)
    # rep.add('dS2_dBeta',dS2_dBeta)
    # rep.add('dS3_dBeta',dS3_dBeta)
    # rep.add('dT1_dMu',dT1_dMu)
    # rep.add('dT2_dMu',dT2_dMu)
    # rep.add('dS1_dMu',dS1_dMu)
    # rep.add('dS2_dMu',dS2_dMu)
    # rep.add('dS3_dMu',dS3_dMu)
    rep.add('dT_dBeta',dT_dBeta)
    rep.add('dT_dMu',dT_dMu)
    rep.add('dT_dBeta_1body',dT_dBeta_1body_eqn_in_dr)
    rep.add('dT_dBeta_2body',dT_dBeta_2body_eqn_in_dr)
    rep.add('dT_dMu_1body',dT_dMu_1body_eqn_in_dr)
    rep.add('dT_dMu_2body',dT_dMu_2body_eqn_in_dr)
    # i = 0
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
print('Tder_bar comput      {}'.format(Tder_bar_time-start_time))
print('T eqn Comput         {}'.format(t2eqn_time-Tder_bar_time))
print('Opt and Print        {}'.format(end_time-t2eqn_time))
print('Total Time           {}'.format(end_time-start_time))
print('-----------------------------------------------------')


