import abc
import collections
import functools
import itertools
import types
import typing

from gristmill import *
from sympy import (
    Expr, Mul, Pow, Integer, Rational, Add, Indexed, IndexedBase, Symbol
)

from sympy.printing.ccode import C99CodePrinter
from sympy.printing.printer import Printer


def print_c_indexed(base, indices):
    """Print indexed objects according to the C syntax.

    The indexed will be printed as multi-dimensional array.
    """
    return base + ''.join(
        '[{}]'.format(i.index) for i in indices
    )


class CPrinter(NaiveCodePrinter):
    """Naive C code printer.

    In this class, just some parameters for the C programming language is fixed
    relative to the base :py:class:`NaiveCodePrinter`.
    """

    def __init__(self, print_indexed_cb=print_c_indexed, **kwargs):
        """Initialize a C code printer.

        The printer class, the name of the template, the line continuation
        symbol, and the statement ending will be set automatically.
        """

        super().__init__(
            C99CodePrinter(), print_indexed_cb=print_indexed_cb,
            line_cont='\\', stmt_end=';', **kwargs
        )

    def form_loop_open(self, ctx):
        """Form the loop opening for C.
        """
        return 'for({index}={lower}; {index}<{upper}, {index}++)'.format(
            index=ctx.index, lower=ctx.lower, upper=ctx.upper
        ) + ' {'

    def form_loop_close(self, _):
        """Form the loop closing for C.
        """
        return '}'

    #
    # Other abstract methods.
    #

    def print_decl(self, event: TensorDecl):
        """Print declaration of an intermediate tensor.

        Here we simply write a declaration for an automatic array.
        """
        ctx = event.comput.ctx

        return '{} {}{}'.format('double', ctx.base, ''.join(
            '[{}]'.format(i.size) for i in ctx.indices
        ))

    def print_begin_body(self, event: BeginBody):
        """Do nothing.
        """
        return None

    def print_out_of_use(self, event: OutOfUse):
        """Do nothing.
        """
        return None

    def print_end_body(self, event: EndBody):
        """Do nothing.
        """
        return None

