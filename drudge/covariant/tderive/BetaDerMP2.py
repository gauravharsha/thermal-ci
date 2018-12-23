    tb0 = 0

    tb0 += einsum(
        "a,a,a,aa->", e0, x, y, t1
    ) / 2

    tb0 -= mu * einsum(
        "a,a,aa->", x, y, t1
    ) / 2

    tb1 = zeros((na, na))

    tb1 += einsum(
        "a,a,a,apaq->pq", e0, x, y, t2
    ) / 2

    tb1 -= mu * einsum(
        "a,a,apaq->pq", x, y, t2
    ) / 2

    sb0 = 0

    sb0 += einsum(
        "a,a,a,aa->", e0, x, y, s1
    ) / 2

    sb0 -= mu * einsum(
        "a,a,aa->", x, y, s1
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

