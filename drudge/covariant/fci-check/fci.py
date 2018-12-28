"""
    Date: Dec 2, 2018
    Modified: Dec 25, 2018
    Python Script to Carry Out the Algebraic Calculations for Ab-Initial Thermal Perturbation Theory
    This is covariant version - i.e. the reference keeps evolving

    Here, we expand the wavefunction as a CI-like series expansion and 
    then form the working equations (following Griffiths)
"""

# from pyspark import SparkContext
from dummy_spark import SparkContext
from drudge import *
from ThermofieldDrudge import *
# from sympy import Symbol, symbols, IndexedBase, sqrt, init_printing, KroneckerDelta
from sympy import *
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
u = Symbol('u')

"""========================================================================="""
"""         Defining the thermal operators                                  """
"""========================================================================="""

adag1 = dr2.simplify( ( c_dag[1,UP] - c_[1,DOWN] ) / sqrt(2) )
adag2 = dr2.simplify( ( c_dag[2,UP] - c_[2,DOWN] ) / sqrt(2) )
bdag1 = dr2.simplify( ( c_dag[1,DOWN] + c_[1,UP] ) / sqrt(2) )
bdag2 = dr2.simplify( ( c_dag[2,DOWN] + c_[2,UP] ) / sqrt(2) )

# Hamiltonian - in the On-Site Basis
full_ham = dr2.simplify( u * c_dag[1,UP] * c_[1,UP] * c_dag[2,UP] * c_[2,UP] )

"""========================================================================="""
"""   FCI Thermal State : Building in terms of the on-site basis            """
"""========================================================================="""

th_zero = dr2.simplify( ( 1 + c_dag[1,UP] * c_dag[1,DOWN] ) * ( 1 + c_dag[2,UP] * c_dag[2,DOWN] ) )

th_one = dr2.simplify( dr2.simplify( adag1 * bdag1 * th_zero ) )

th_two = dr2.simplify( dr2.simplify( adag2 * bdag2 * th_zero ) )

th_onetwo = dr2.simplify( dr2.simplify( adag1 * adag2 * bdag1 * bdag2 * th_zero ) )

"""========================================================================="""
"""                        Testing the Stuff Generated                      """
"""========================================================================="""
with dr2.report('fci.html','Thermal Perturbation Theory') as rep:
    rep.add('Full Hamiltonian',full_ham)
    rep.add('adag1',adag1)
    rep.add('bdag1',bdag1)
    rep.add('|-,->',th_zero)
    rep.add('adag1*bdag1',dr2.simplify(adag1*bdag1))
    rep.add('|1,1>',th_one)
    rep.add('adag2*bdag2',dr2.simplify(adag2*bdag2))
    rep.add('|2,2>',th_two)
    rep.add('|12,12>',th_onetwo)
