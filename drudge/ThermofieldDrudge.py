from sympy import IndexedBase, Symbol
from drudge.term import Range
from drudge import GenMBDrudge, SpinOneHalfGenDrudge


class ThermofieldDrudge(GenMBDrudge):
    """Drudge for the fermionic problems for thermal wave function theories in
    the spin-orbital picture.

    The 'd' labelled creation and annihilation operators represent ZERO
    TEMPERATURE fermionic operators.

    The 'c' labelled creation and annihilation operators represent
    FINITE-TEMPERATURE creation annihilation operators.

    Spin UP can be considered as the physical system while spin 'DOWN' as the
    auxiliary one.
    """

    DEFAULT_ORB_DUMMS = tuple(Symbol(i) for i in 'abcdijklpqrs') + tuple(
        Symbol('p{}'.format(i)) for i in range(50)
    )

    def __init__(
            self, *args, op_label='d', one_body=IndexedBase('h'),
            orb=((Range('L', 0, Symbol('na')), DEFAULT_ORB_DUMMS),), spin=(),
            dbbar=True, **kwargs
    ):
        """Initialize the particle-hole drudge."""

        super().__init__(
            *args, op_label=op_label, one_body=one_body, spin=spin,
            dbbar=dbbar, orb=orb, **kwargs
        )

        self._sp_dr = SpinOneHalfGenDrudge(
            *args, orb=orb, dbbar=dbbar, **kwargs
        )

        c_ = self._sp_dr.names.c_
        c_dag = self._sp_dr.names.c_dag

        self.c_ = c_
        self.c_dag = c_dag

        self.set_name(**{
            c_.label[0]+'_': c_,
            c_.label[0]+'_dag': c_dag
        })
