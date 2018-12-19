from numpy import einsum

def MPEnergyCov(e0, eri, t1, t2, s1, s2, s3, s4, x, y):
    e1 = 0

    e1 += einsum(
        "a,b,abab->", y**2, y**2, eri
    ) / 2

    e1 += 2 * einsum(
        "a,a,a,aa->", e0, x, y, t1
    )

    e2 = 0

    e2 += einsum(
        "a,b,c,d,abcd,abcd->", x, x, y, y, t2, eri
    )

    e2 += einsum(
        "a,b,c,d,abcd,cdab->", x, x, y, y, t2, eri
    ) / 4

    e2 += 2 * einsum(
        "a,a,a,aa->", e0, x, y, s1
    )

    e2 += einsum(
        "c,a,b,ab,acbc->", y**2, x, y, t1, eri
    )

    e2 += einsum(
        "c,a,b,ab,bcac->", y**2, x, y, t1, eri
    )

    ov1 = 0

    ov1 += einsum(
        "abcd->", t2**2
    )

    ov1 += einsum(
        "ab->", t1**2
    )

    ov2 = 0

    ov2 += einsum(
        "abcdijkl->", s4**2
    ) / 576

    ov2 += einsum(
        "abcd->", s2**2
    ) / 4

    ov2 += einsum(
        "abcdij->", s3**2
    ) / 36

    ov2 += einsum(
        "ab->", s1**2
    )

    return e1, e2, ov1, ov2

