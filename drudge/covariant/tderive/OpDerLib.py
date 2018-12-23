from numpy import zeros, einsum

def BetaDerMP(e0, mu, t2, s2, s4, x, y):
    na = len(e0)

    tb1 = zeros((na, na))

    tb1 += einsum(
        "a,a,a,apaq->pq", e0, x, y, t2
    ) / 2

    tb1 -= mu * einsum(
        "a,a,apaq->pq", x, y, t2
    ) / 2

    sb1 = zeros((na, na))

    sb1 += einsum(
        "a,a,a,apaq->pq", e0, x, y, s2
    ) / 2

    sb1 -= mu * einsum(
        "a,a,apaq->pq", x, y, s2
    ) / 2

    sb3 = zeros((na, na, na, na, na, na))

    sb3 += mu * einsum(
        "d,d,dpqrabcd->pqrabc", x, y, s4
    )

    sb3 -= einsum(
        "d,d,d,dpqrabcd->pqrabc", e0, x, y, s4
    )

    return tb1, sb1, sb3

def MuDerMP(e0, beta, t2, s2, s4, x, y):
    na = len(x)

    tm1 = zeros((na, na))

    tm1 -= beta * einsum(
        "a,a,apaq->pq", x, y, t2
    ) / 2

    sm1 = zeros((na, na))

    sm1 -= beta * einsum(
        "a,a,apaq->pq", x, y, s2
    ) / 2

    sm3 = zeros((na, na, na, na, na, na))

    sm3 += beta * einsum(
        "d,d,dpqrabcd->pqrabc", x, y, s4
    )

    return tm1, sm1, sm3
