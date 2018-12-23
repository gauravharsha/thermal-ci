    tau0 = zeros((na))

    tau0 += einsum(
        "b,baba->a", y**2, u
    )

    tau1 = zeros((na))

    tau1 += einsum(
        "a->a", tau0
    )

    tau13 = zeros((na))

    tau13 += t0 * einsum(
        "a->a", tau0
    )

    del tau0

    tau1 += 4*t0 * einsum(
        "a->a", e0
    )

    e1 = 0

    e1 += einsum(
        "a,a->", y**2, tau1
    ) / 2

    del tau1

    tau2 = zeros((na))

    tau2 += einsum(
        "ab->a", t1**2
    )

    tau11 = zeros((na))

    tau11 += 2 * einsum(
        "a,a->a", x**2, tau2
    )

    del tau2

    tau3 = zeros((na))

    tau3 += einsum(
        "ba->a", t1**2
    )

    tau11 -= 2 * einsum(
        "a,a->a", y**2, tau3
    )

    del tau3

    tau4 = zeros((na))

    tau4 += einsum(
        "abcd->a", t2**2
    )

    tau11 += einsum(
        "a,a->a", x**2, tau4
    )

    del tau4

    tau5 = zeros((na))

    tau5 += einsum(
        "bcad->a", t2**2
    )

    tau11 -= einsum(
        "a,a->a", y**2, tau5
    )

    del tau5

    tau6 = zeros((na))

    tau6 += einsum(
        "cb,acab->a", t1, t2
    )

    tau7 = zeros((na))

    tau7 += einsum(
        "a,a->a", tau6, y
    )

    del tau6

    tau7 += einsum(
        "a,aa->a", y, s1
    )

    tau7 += t0 * einsum(
        "a,aa->a", y, t1
    )

    tau11 += 4 * einsum(
        "a,a->a", tau7, x
    )

    del tau7

    tau8 = 0

    tau8 += einsum(
        "ab->", t1**2
    )

    tau10 = 0

    tau10 += 4*tau8 

    ov2 = 0

    ov2 += tau8

    del tau8

    tau9 = 0

    tau9 += einsum(
        "abcd->", t2**2
    )

    tau10 += tau9

    tau11 += tau10 * einsum(
        "a->a", y**2
    ) / 2

    del tau10

    ov2 += tau9 / 4

    del tau9

    tau11 += (4*s0 + 2*t0**2) * einsum(
        "a->a", y**2
    )

    e2 = 0

    e2 += einsum(
        "a,a->", e0, tau11
    ) / 2

    del tau11

    tau12 = zeros((na, na, na))

    tau12 += einsum(
        "abac->abc", u
    )

    tau12 += einsum(
        "acab->abc", u
    )

    tau13 += einsum(
        "b,c,bc,abc->a", x, y, t1, tau12
    )

    del tau12

    e2 += einsum(
        "a,a->", y**2, tau13
    )

    del tau13

    tau14 = zeros((na, na, na, na))

    tau14 += einsum(
        "abcd,badc->abcd", t2, u
    )

    tau14 += einsum(
        "abdc,cdba->abcd", t2, u
    )

    tau15 = zeros((na))

    tau15 += einsum(
        "b,c,d,bcda->a", x, x, y, tau14
    )

    del tau14

    e2 += einsum(
        "a,a->", tau15, y
    ) / 4

    del tau15

    e1 += 2 * einsum(
        "a,a,a,aa->", e0, x, y, t1
    )

    ov2 += t0**2

    ov3 = 0

    ov3 += 2*s0*t0

    ov3 += 2 * einsum(
        "ba,ba->", s1, t1
    )

    ov3 += einsum(
        "cbda,cbda->", s2, t2
    ) / 2

    ov4 = 0

    ov4 += s0**2

    ov4 += einsum(
        "abcd->", s2**2
    ) / 4

    ov4 += einsum(
        "abcdij->", s3**2
    ) / 36

    ov4 += einsum(
        "abcdijkl->", s4**2
    ) / 576

    ov4 += einsum(
        "ab->", s1**2
    )

