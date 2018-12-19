import abc
import collections
import functools
import itertools
import types
import typing

from sympy import (
    Expr, Mul, Pow, Integer, Rational, Add, Indexed, IndexedBase, Symbol
)

from gristmill import *

from sympy.printing.ccode import C99CodePrinter
from sympy.printing.printer import Printer


def print_c_indexed(base, indices):
    """Print indexed objects according to the C syntax.

    The indexed will be printed as multi-dimensional array.
    """
    return base + ''.join(
        '[{}]'.format(i.index) for i in indices
    )

class TensorDecl(typing.NamedTuple):
    """Events for declaration of intermediate tensors.
    """

    comput: TensorComp

    def __repr__(self):
        """Form a string, mostly for debugging.
        """
        return '_TensorDecl({!s})'.format(self.comput)


class BeginBody:
    """Events for the beginning of the main computational body.
    """
    __slots__ = []

    def __repr__(self):
        """Form a string representation.
        """
        return 'BeginBody()'


class BeforeComp(typing.NamedTuple):
    """Events that come before the first computation of any tensor.

    Normally, it should be rendered as memory allocation or initialization to
    zero.
    """

    comput: TensorComp

    def __repr__(self):
        """Form a string, mostly for debugging."""
        return '_BeforeCompute({!s})'.format(self.comput)


class CompTerm(typing.NamedTuple):
    """Events for the computation of a term in a tensor.

    This events attempt to add a term to the LHS of the computation.

    Attributes
    ----------

    comput
       The actual full computation.

    term_idx
        Index to the term being computed.

    term_ctx
        The rendering context for the current term.

    """

    comput: TensorComp
    term_idx: int
    term_ctx: types.SimpleNamespace

    def __repr__(self):
        """Format a string, mostly for debugging."""
        return "_ComputeTerm({!s}, {!s})".format(
            self.comput, self.term_ctx.orig_term
        )


class OutOfUse(typing.NamedTuple):
    """Events after intermediate tensors are no longer in use.

    This event can be used for the freeing of the associated computer memory.
    """

    comput: TensorComp

    def __repr__(self):
        """Form a string, mostly for debugging."""
        return '_NoLongerInUse({!s})'.format(self.comput)


class EndBody:
    """Events for the end of the main computational body.
    """
    __slots__ = []

    def __repr__(self):
        """Form a string representation.
        """
        return 'EndBody()'


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

