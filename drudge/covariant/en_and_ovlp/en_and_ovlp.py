"""
    Date: Dec 4, 2018
    Modified: Dec 20, 2018
    Python Script to Carry Out the Algebraic Calculations for Ab-Initial Thermal Perturbation Theory
    This is covariant version - i.e. the reference keeps evolving
    Location: Documents/projects/Molecules/thermal_covariant/drudge/ccd/

    Here, we are only exporting the energy and the overlap
"""

# from pyspark import SparkContext
from dummy_spark import SparkContext
from drudge import *
from ThermofieldDrudge import *
from sympy import *
#Symbol, symbols, IndexedBase, sqrt, init_printing, KroneckerDelta
from gristmill import *
import time
import pdb

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
full_ham2 = dr.simplify(dr.ham)

e0 = IndexedBase('e0')
ham1 = dr.einst(
    e0[p]*d_dag[p]*d_[p]
)

ham2 = full_ham2.filter(lambda x: len(x.vecs)>2)

full_ham = dr.simplify(ham1 + ham2)

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

ham_th1 = Tensor(
    dr2,
    ham1.subst_all(tfd_defs).terms
)

ham_th2 = Tensor(
    dr2,
    ham2.subst_all(tfd_defs).terms
)

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
"""   Constructing the 1st order wavefunction correction                    """
"""-------------------------------------------------------------------------"""

# First order perturbation wavefunction description
t0 = Symbol('t0')
t1 = IndexedBase('t1')
t2 = IndexedBase('t2')

# Second order perturbation wavefunction description
s0 = Symbol('s0')
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
    Perm([0,1,2,4,3,5],NEG),
    Perm([0,1,2,3,5,4],NEG),
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
    t1[a, b] * ( (c_dag[a,UP] * c_dag[b,DOWN]) )
)

T2 = dr2.einst(
    t2[a, b, c, d] * (c_dag[a,UP] * c_dag[b,UP] * c_dag[d,DOWN] \
        * c_dag[c,DOWN]) / 4 
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

S4 = dr2.einst(
    s4[a, b, c, d, i, j, k, l] * (c_dag[a,UP] * c_dag[b,UP] * c_dag[c,UP] * c_dag[d,UP] \
        * c_dag[l,DOWN] * c_dag[k,DOWN] * c_dag[j,DOWN] * c_dag[i,DOWN]) / 576
)

T1dag = dr2.einst(
    t1[a,b] * (c_[b,DOWN] * c_[a,UP] )
)

T2dag = dr2.einst(
    t2[a,b,c,d] * ( c_[c,DOWN] * c_[d,DOWN] * c_[b,UP] * c_[a,UP] ) / 4
)

S1dag = dr2.einst(
    s1[a, i] * ( (c_[i,DOWN] * c_[a,UP]) )
)

S2dag = dr2.einst(
    s2[a, b, i, j] * (c_[i,DOWN] * c_[j,DOWN] \
        * c_[b,UP] * c_[a,UP]) / 4 
)

S3dag = dr2.einst(
    s3[a, b, c, i, j, k] * ( c_[i,DOWN] * c_[j,DOWN] * c_[k,DOWN]\
        * c_[c,UP] * c_[b,UP] * c_[a,UP]) / 36
)

S4dag = dr2.einst(
    s4[a, b, c, d, i, j, k, l] * (c_[i,DOWN] * c_[j,DOWN] * c_[k,DOWN] * c_[l,DOWN] \
        * c_[d,UP] * c_[c,UP] * c_[b,UP] * c_[a,UP]) / 576
)

Tvec = dr2.simplify(t0 + T1 + T2)
Tdagvec = dr2.simplify(t0 + T1dag + T2dag)

Svec = dr2.simplify(s0 + S1 + S2 + S3 + S4)
Sdagvec = dr2.simplify(s0 + S1dag + S2dag + S3dag + S4dag)


"""-------------------------------------------------------------------------"""
"""     Computing the energy and the overlap operators at given T           """
"""-------------------------------------------------------------------------"""

# Overlap Operators
mp2_ovlp_op = dr2.memoize(lambda: dr2.simplify( Tdagvec * Tvec ), 'mp2_ovlp_order.pickle')
mp3_ovlp_op = dr2.memoize(lambda: dr2.simplify( Tdagvec * Svec + Sdagvec * Tvec), 'mp3_ovlp_order.pickle')
mp4_ovlp_op = dr2.memoize(lambda: dr2.simplify( Sdagvec * Svec ), 'mp4_ovlp_order.pickle')

# Energy Operators (NOTE that this is just the numerator of the energy)
mp1_enum_op = dr2.memoize(
    lambda: dr2.simplify(
        Tdagvec*ham_th1 + ham_th1*Tvec + ham_th2
    ),
    'mp1_enum_op.pickle'
)
mp2_enum_op = dr2.memoize(
    lambda: dr2.simplify(
        Tdagvec*ham_th1*Tvec + Sdagvec*ham_th1 + ham_th1*Svec + ham_th2*Tvec + Tdagvec*ham_th2
    ),
    'mp2_enum_op.pickle'
)

# Hbar Computing Time
hbar_time = time.time()

# Now, Overlap expectation values
mp2_ovlp = dr2.memoize( lambda: dr2.simplify( dr2.simplify( dr2.eval_phys_vev( dr2.simplify( mp2_ovlp_op ))) ), 'mp2_ovlp_eqn.pickle')
mp3_ovlp = dr2.memoize( lambda: dr2.simplify( dr2.simplify( dr2.eval_phys_vev( dr2.simplify( mp3_ovlp_op ))) ), 'mp3_ovlp_eqn.pickle')
mp4_ovlp = dr2.memoize( lambda: dr2.simplify( dr2.simplify( dr2.eval_phys_vev( dr2.simplify( mp4_ovlp_op ))) ), 'mp4_ovlp_eqn.pickle')

# and Energy-numerator expectation values
mp1_enum = dr2.memoize( lambda: dr2.simplify( dr2.simplify( dr2.eval_phys_vev( dr2.simplify( mp1_enum_op ))) ), 'mp1_enum_eqn.pickle')
mp2_enum = dr2.memoize( lambda: dr2.simplify( dr2.simplify( dr2.eval_phys_vev( dr2.simplify( mp2_enum_op ))) ), 'mp2_enum_eqn.pickle')

print('Equations obtained')

# T2 equation time
t2eqn_time = time.time()

"""-------------------------------------------------------------------------"""
"""   Optimizing the Working Eqn and setting up the Fortran Printing        """
"""-------------------------------------------------------------------------"""

###############################################################################
# Some Pre processing to convert everything into a                            # 
# single drudge type data for gristmill                                       #
###############################################################################


mp2_ovlp_in_dr = Tensor(
    dr,
    mp2_ovlp.terms
)
mp3_ovlp_in_dr = Tensor(
    dr,
    mp3_ovlp.terms
)
mp4_ovlp_in_dr = Tensor(
    dr,
    mp4_ovlp.terms
)
mp1_enum_in_dr = Tensor(
    dr,
    mp1_enum.terms
)
mp2_enum_in_dr = Tensor(
    dr,
    mp2_enum.terms
)

# Now the actual optimization process
work_eqn = [
    dr.define(Symbol('e1'), mp1_enum_in_dr),
    dr.define(Symbol('e2'), mp2_enum_in_dr),
    dr.define(Symbol('ov2'), mp2_ovlp_in_dr),
    dr.define(Symbol('ov3'), mp3_ovlp_in_dr),
    dr.define(Symbol('ov4'), mp4_ovlp_in_dr),
]

# Original Un-optimized cost
orig_cost = get_flop_cost(work_eqn, leading=True)
init_printing()
print(orig_cost)

# Optimizing
try:
    eval_seq = optimize(
        work_eqn,
        interm_fmt = 'tau{}'
    )
except ValueError:
    print("Encountered Value Error, so proceeding with the 'work_eqn' instead")
    eval_seq = work_eqn

print(len(eval_seq))
opt_cost = get_flop_cost(eval_seq, leading=True)
print(opt_cost)

# Code Generation
# fort_print = C99CodePrinter(default_type='Real (Kind=8)', explicit_bounds=True)
# code = fort_print.doprint(eval_seq, separate_decls=False)
# c_print = C99CodePrinter()
# pdb.set_trace()

ein_print = EinsumPrinter()
code = ein_print.doprint(eval_seq, separate_decls=False)

with open('ovlp.py','w') as fp:
    print(code, file=fp)



"""-------------------------------------------------------------------------"""
"""                        Testing the Stuff Generated                      """
"""-------------------------------------------------------------------------"""
with dr.report('en_and_ovlp_cov_pt.html','Thermal Perturbation Theory') as rep:
    rep.add('Full Hamiltonian',full_ham)
    rep.add('Thermal Ham',ham_th)
    rep.add('Tvec',Tvec)
    rep.add('Svec',Svec)
    rep.add('The Spin Rep for Tilde and non Tilde',c_[a,UP] + c_[b,DOWN])
    rep.add('Ovlp at order 2',mp2_ovlp_in_dr)
    rep.add('Ovlp at order 3',mp3_ovlp_in_dr)
    rep.add('Ovlp at order 4',mp4_ovlp_in_dr)
    rep.add('MP1 Energy',mp1_enum_in_dr)
    rep.add('MP2 Energy',mp2_enum_in_dr)
    i = 0
    for eq in eval_seq:
        rep.add('Intermediate Eqn {}'.format(i),eq)
        i += 1

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


