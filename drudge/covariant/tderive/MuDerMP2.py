    tm0 = 0

    tm0 -= beta * einsum(
        "a,a,aa->", x, y, t1
    ) / 2

    tm1 = zeros((na, na))

    tm1 -= beta * einsum(
        "a,a,apaq->pq", x, y, t2
    ) / 2

    sm0 = 0

    sm0 -= beta * einsum(
        "a,a,aa->", x, y, s1
    ) / 2

    sm1 = zeros((na, na))

    sm1 -= beta * einsum(
        "a,a,apaq->pq", x, y, s2
    ) / 2

    sm2 = zeros((na, na, na, na))

    sm2 -= beta * einsum(
        "a,a,apqars->pqrs", x, y, s3
    ) / 2

    sm3 = zeros((na, na, na, na, na, na))

    sm3 += beta * einsum(
        "d,d,dpqrabcd->pqrabc", x, y, s4
    ) / 2

