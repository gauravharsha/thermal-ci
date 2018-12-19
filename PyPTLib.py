from numpy import einsum, zeros

def mpenergycov(e0, eri, t1, t2, s1, s2, s3, s4, x, y):
    """
    Function to return the 1st and 2nd order correction to energy
    and corresponding 2nd and 4th order correction to the overlaps
    
    Input arguments
        e0      ::      array of dimension (nso,1) containing matrix 
                        elements of the 1-body part of Hamiltonian 
                        (which we assume is in the diagonalized form)
        eri     ::      array of dimension (nso,nso,nso,nso) containing
                        matrix elements for 2-body Hamiltonian
        t1, t2  ::      arrays of dimension (nso,nso) and (nso,nso,nso,nso)
                        respectively and containing info about 1st order
                        correction to the thermal state
        s1--s4  ::      arrays of dimension (nso,nso), (nso,nso,nso,nso)
                        and so on, and containing info about 2nd order
                        correction to the thermal state
        x,y     ::      HFB parameters that define the thermal state or
                        the temperature

    Output parameters
        e1      ::      1st order correction to the < psi | H | psi >
        e2      ::      2nd order correction to the < psi | H | psi >
        ov1     ::      1st order correction to the < psi | psi >
        ov2     ::      2nd order correction to the < psi | psi >

    """
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
    ) / 4

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

    ov12 = 0

    ov12 += einsum(
        "abcd,abcd->", t2, s2
    ) / 4

    ov12 += einsum(
        "ab,ab->", t1, s1
    )

    return e1, e2, ov1, ov2, ov12


def covbetapt2(e0, eri, t1, t2, s1, s2, s3, s4, x, y):

    na = len(x)

    tau0 = zeros((na, na))

    tau0 += einsum(
        "c,cacb->ab", y**2, eri
    )

    tau3 = zeros((na, na))

    tau3 += einsum(
        "c,ca,bc->ab", x, t1, tau0
    )

    rs1 = zeros((na, na))

    rs1 += -0.5 * einsum(
        "a,ba->ab", x, tau3
    )

    del tau3

    tau5 = zeros((na, na, na))

    tau5 += 1.0 * einsum(
        "ab,bc->abc", t1, tau0
    )

    tau6 = zeros((na, na, na, na))

    tau6 += einsum(
        "d,i,di,iabc->abcd", x, x, tau0, t2
    )

    rs2 = zeros((na, na, na, na))

    rs2 += 0.5 * einsum(
        "acdb->abcd", tau6
    )

    rs2 += -0.5 * einsum(
        "bcda->abcd", tau6
    )

    del tau6

    tau7 = zeros((na, na, na, na))

    tau7 += einsum(
        "d,i,id,abci->abcd", y, y, tau0, t2
    )

    rs2 += 0.5 * einsum(
        "abcd->abcd", tau7
    )

    rs2 += -0.5 * einsum(
        "abdc->abcd", tau7
    )

    del tau7

    tau12 = zeros((na, na, na, na))

    tau12 += 1.0 * einsum(
        "d,ab,cd->abcd", y, t1, tau0
    )

    tau12 += -1.0 * einsum(
        "b,ad,cb->abcd", y, t1, tau0
    )

    rs2 += -0.5 * einsum(
        "b,acbd->abcd", x, tau12
    )

    del tau12

    tau13 = zeros((na, na, na, na))

    tau13 += -1.0 * einsum(
        "d,ab,cd->abcd", y, t1, tau0
    )

    tau13 += 1.0 * einsum(
        "b,ad,cb->abcd", y, t1, tau0
    )

    rs2 += -0.5 * einsum(
        "a,bcad->abcd", x, tau13
    )

    del tau13

    tau20 = zeros((na, na, na, na, na, na))

    tau20 += 1.0 * einsum(
        "dj,acbi->abcdij", tau0, t2
    )

    tau21 = zeros((na, na, na, na, na, na))

    tau21 += 1.0 * einsum(
        "i,ij,abcd->abcdij", x, tau0, t2
    )

    tau21 += -1.0 * einsum(
        "b,bj,aicd->abcdij", x, tau0, t2
    )

    tau21 += 1.0 * einsum(
        "a,aj,bicd->abcdij", x, tau0, t2
    )

    rs3 = zeros((na, na, na, na, na, na))

    rs3 += 0.5 * einsum(
        "i,abdjci->abcdij", y, tau21
    )

    rs3 += -0.5 * einsum(
        "j,abdicj->abcdij", y, tau21
    )

    del tau21

    tau22 = zeros((na, na, na, na, na, na))

    tau22 += 1.0 * einsum(
        "i,ij,abcd->abcdij", x, tau0, t2
    )

    tau22 += -1.0 * einsum(
        "b,bj,aicd->abcdij", x, tau0, t2
    )

    rs3 += -0.5 * einsum(
        "d,abijcd->abcdij", y, tau22
    )

    del tau22

    rt1 = zeros((na, na))

    rt1 += -0.5 * einsum(
        "a,b,ab->ab", x, y, tau0
    )

    rs1 += 0.5 * einsum(
        "c,d,dc,cabd->ab", x, y, tau0, t2
    )

    del tau0

    tau1 = zeros((na))

    tau1 += einsum(
        "b,abab->a", y**2, eri
    )

    tau2 = 0

    tau2 += einsum(
        "a,a->", y**2, tau1
    )

    del tau1

    rs1 += -0.25*tau2 * einsum(
        "ab->ab", t1
    )

    rs2 += -0.25*tau2 * einsum(
        "abcd->abcd", t2
    )

    del tau2

    tau4 = zeros((na, na, na, na))

    tau4 += 1.0 * einsum(
        "d,ab,cbda->abcd", y, t1, eri
    )

    tau4 += -0.5 * einsum(
        "i,iadb,cbai->abcd", x, t2, eri
    )

    rs1 += -0.5 * einsum(
        "a,c,d,cdab->ab", x, x, y, tau4
    )

    del tau4

    tau5 += 0.5 * einsum(
        "d,i,daib,bicd->abc", x, y, t2, eri
    )

    rs1 += 0.5 * einsum(
        "b,c,acb->ab", y, y, tau5
    )

    del tau5

    tau8 = zeros((na, na, na, na, na))

    tau8 += -0.5 * einsum(
        "j,jabc,idaj->abcdi", x, t2, eri
    )

    tau8 += -1.0 * einsum(
        "c,ab,idca->abcdi", y, t1, eri
    )

    tau8 += 1.0 * einsum(
        "b,ac,idba->abcdi", y, t1, eri
    )

    rs2 += 0.5 * einsum(
        "a,b,i,icdab->abcd", x, x, x, tau8
    )

    del tau8

    tau9 = zeros((na, na, na, na, na))

    tau9 += -0.5 * einsum(
        "j,abjc,cjid->abcdi", y, t2, eri
    )

    tau9 += 1.0 * einsum(
        "b,ac,bcid->abcdi", x, t1, eri
    )

    tau9 += -1.0 * einsum(
        "a,bc,acid->abcdi", x, t1, eri
    )

    rs2 += -0.5 * einsum(
        "c,d,i,abidc->abcd", y, y, y, tau9
    )

    del tau9

    tau10 = zeros((na, na, na, na, na, na))

    tau10 += -1.0 * einsum(
        "j,abcd,idja->abcdij", y, t2, eri
    )

    tau10 += 1.0 * einsum(
        "c,abjd,idca->abcdij", y, t2, eri
    )

    rs2 += -0.5 * einsum(
        "a,i,j,ibdjac->abcd", x, x, y, tau10
    )

    del tau10

    tau11 = zeros((na, na, na, na, na, na))

    tau11 += 1.0 * einsum(
        "j,abcd,idja->abcdij", y, t2, eri
    )

    tau11 += -1.0 * einsum(
        "c,abjd,idca->abcdij", y, t2, eri
    )

    rs2 += -0.5 * einsum(
        "b,i,j,iadjbc->abcd", x, x, y, tau11
    )

    del tau11

    tau14 = zeros((na, na, na, na, na, na, na))

    tau14 += 1.0 * einsum(
        "k,abcd,jika->abcdijk", y, t2, eri
    )

    tau14 += -1.0 * einsum(
        "d,abck,jida->abcdijk", y, t2, eri
    )

    tau14 += 1.0 * einsum(
        "c,abdk,jica->abcdijk", y, t2, eri
    )

    rs3 += 0.5 * einsum(
        "a,b,k,kcdiabj->abcdij", x, x, x, tau14
    )

    del tau14

    tau15 = zeros((na, na, na, na, na, na, na))

    tau15 += -1.0 * einsum(
        "k,abcd,jika->abcdijk", y, t2, eri
    )

    tau15 += 1.0 * einsum(
        "d,abck,jida->abcdijk", y, t2, eri
    )

    tau15 += 1.0 * einsum(
        "c,abkd,jica->abcdijk", y, t2, eri
    )

    rs3 += -0.5 * einsum(
        "a,c,k,kbdjaci->abcdij", x, x, x, tau15
    )

    rs3 += -0.5 * einsum(
        "b,c,k,kadjcbi->abcdij", x, x, x, tau15
    )

    del tau15

    tau16 = zeros((na, na, na, na, na, na, na))

    tau16 += -1.0 * einsum(
        "i,abcd,dikj->abcdijk", x, t2, eri
    )

    tau16 += 1.0 * einsum(
        "b,aicd,dbkj->abcdijk", x, t2, eri
    )

    tau16 += -1.0 * einsum(
        "a,bicd,dakj->abcdijk", x, t2, eri
    )

    tau17 = zeros((na, na, na, na, na, na))

    tau17 += 0.5 * einsum(
        "a,b,k,cdikjba->abcdij", y, y, y, tau16
    )

    del tau16

    rs3 += 1.0 * einsum(
        "djabic->abcdij", tau17
    )

    rs3 += -1.0 * einsum(
        "diabjc->abcdij", tau17
    )

    rs3 += -1.0 * einsum(
        "ijabdc->abcdij", tau17
    )

    del tau17

    tau18 = zeros((na, na, na, na, na, na))

    tau18 += 1.0 * einsum(
        "i,ab,dcij->abcdij", y, t1, eri
    )

    tau18 += -1.0 * einsum(
        "b,ai,dcbj->abcdij", y, t1, eri
    )

    rs3 += -0.5 * einsum(
        "a,b,j,cdbaij->abcdij", x, x, y, tau18
    )

    rs3 += -0.5 * einsum(
        "a,c,i,bdcaji->abcdij", x, x, y, tau18
    )

    rs3 += -0.5 * einsum(
        "b,c,j,adcbij->abcdij", x, x, y, tau18
    )

    del tau18

    tau19 = zeros((na, na, na, na, na, na))

    tau19 += 1.0 * einsum(
        "c,ab,dcji->abcdij", x, t1, eri
    )

    tau19 += -1.0 * einsum(
        "a,cb,daji->abcdij", x, t1, eri
    )

    rs3 += 0.5 * einsum(
        "b,d,i,cjabid->abcdij", x, y, y, tau19
    )

    del tau19

    tau20 += -1.0 * einsum(
        "c,i,ab,dcij->abcdij", x, y, t1, eri
    )

    rs3 += -0.5 * einsum(
        "a,d,bicajd->abcdij", x, y, tau20
    )

    del tau20

    tau23 = zeros((na, na, na, na, na, na, na, na))

    tau23 += 1.0 * einsum(
        "k,abcd,ijkl->abcdijkl", y, t2, eri
    )

    tau23 += -1.0 * einsum(
        "d,abck,ijdl->abcdijkl", y, t2, eri
    )

    tau23 += 1.0 * einsum(
        "c,abdk,ijcl->abcdijkl", y, t2, eri
    )

    rs4 = zeros((na, na, na, na, na, na, na, na))

    rs4 += 0.5 * einsum(
        "a,b,i,cdjkabli->abcdijkl", x, x, y, tau23
    )

    del tau23

    tau24 = zeros((na, na, na, na, na, na, na, na))

    tau24 += 1.0 * einsum(
        "i,abcd,jilk->abcdijkl", x, t2, eri
    )

    tau24 += -1.0 * einsum(
        "b,aicd,jblk->abcdijkl", x, t2, eri
    )

    tau24 += 1.0 * einsum(
        "a,bicd,jalk->abcdijkl", x, t2, eri
    )

    rs4 += 0.5 * einsum(
        "a,j,l,bcikdalj->abcdijkl", x, y, y, tau24
    )

    rs4 += -0.5 * einsum(
        "a,j,k,bcildakj->abcdijkl", x, y, y, tau24
    )

    rs4 += -0.5 * einsum(
        "a,k,l,bcijdalk->abcdijkl", x, y, y, tau24
    )

    del tau24

    tau25 = zeros((na, na, na, na, na, na, na, na))

    tau25 += 1.0 * einsum(
        "k,abcd,ijkl->abcdijkl", y, t2, eri
    )

    tau25 += 1.0 * einsum(
        "d,abkc,ijdl->abcdijkl", y, t2, eri
    )

    tau25 += -1.0 * einsum(
        "c,abkd,ijcl->abcdijkl", y, t2, eri
    )

    rs4 += -0.5 * einsum(
        "a,c,i,bdklacji->abcdijkl", x, x, y, tau25
    )

    rs4 += 0.5 * einsum(
        "b,d,i,ackldbji->abcdijkl", x, x, y, tau25
    )

    del tau25

    tau26 = zeros((na, na, na, na, na, na, na, na))

    tau26 += -1.0 * einsum(
        "k,abcd,ijkl->abcdijkl", y, t2, eri
    )

    tau26 += 1.0 * einsum(
        "d,abck,ijdl->abcdijkl", y, t2, eri
    )

    tau26 += 1.0 * einsum(
        "c,abkd,ijcl->abcdijkl", y, t2, eri
    )

    rs4 += 0.5 * einsum(
        "a,d,i,bcjladki->abcdijkl", x, x, y, tau26
    )

    rs4 += -0.5 * einsum(
        "b,c,i,adjlcbki->abcdijkl", x, x, y, tau26
    )

    rs4 += -0.5 * einsum(
        "c,d,i,abjldcki->abcdijkl", x, x, y, tau26
    )

    del tau26

    tau27 = zeros((na, na, na, na, na, na, na, na))

    tau27 += 1.0 * einsum(
        "i,abcd,jilk->abcdijkl", x, t2, eri
    )

    tau27 += -1.0 * einsum(
        "b,aicd,jblk->abcdijkl", x, t2, eri
    )

    rs4 += -0.5 * einsum(
        "b,j,l,acikdblj->abcdijkl", x, y, y, tau27
    )

    rs4 += 0.5 * einsum(
        "c,j,k,adilbckj->abcdijkl", x, y, y, tau27
    )

    rs4 += 0.5 * einsum(
        "c,k,l,adijbclk->abcdijkl", x, y, y, tau27
    )

    del tau27

    tau28 = zeros((na, na, na, na, na, na, na, na))

    tau28 += 1.0 * einsum(
        "k,abcd,ijkl->abcdijkl", y, t2, eri
    )

    tau28 += -1.0 * einsum(
        "d,abck,ijdl->abcdijkl", y, t2, eri
    )

    rs4 += -0.5 * einsum(
        "b,d,k,acildbjk->abcdijkl", x, x, y, tau28
    )

    del tau28

    rt1 += 0.5 * einsum(
        "b,b,ab->ab", y**2, e0, t1
    )

    rt1 += -0.5 * einsum(
        "a,a,ab->ab", x**2, e0, t1
    )

    rt1 += 0.5 * einsum(
        "c,c,c,cabc->ab", e0, x, y, t2
    )

    rt2 = zeros((na, na, na, na))

    rt2 += 0.5 * einsum(
        "c,c,abcd->abcd", y**2, e0, t2
    )

    rt2 += 0.5 * einsum(
        "d,d,abcd->abcd", y**2, e0, t2
    )

    rt2 += -0.5 * einsum(
        "a,a,abcd->abcd", x**2, e0, t2
    )

    rt2 += -0.5 * einsum(
        "b,b,abcd->abcd", x**2, e0, t2
    )

    rt2 += -0.5 * einsum(
        "a,b,c,d,badc->abcd", x, x, y, y, eri
    )

    rs1 += 0.5 * einsum(
        "b,b,ab->ab", y**2, e0, s1
    )

    rs1 += -0.5 * einsum(
        "a,a,ab->ab", x**2, e0, s1
    )

    rs1 += 0.5 * einsum(
        "c,c,c,cabc->ab", e0, x, y, s2
    )

    rs2 += 0.5 * einsum(
        "c,c,abcd->abcd", y**2, e0, s2
    )

    rs2 += 0.5 * einsum(
        "d,d,abcd->abcd", y**2, e0, s2
    )

    rs2 += -0.5 * einsum(
        "a,a,abcd->abcd", x**2, e0, s2
    )

    rs2 += -0.5 * einsum(
        "b,b,abcd->abcd", x**2, e0, s2
    )

    rs2 += -0.5 * einsum(
        "i,i,i,iabcdi->abcd", e0, x, y, s3
    )

    rs3 += 0.5 * einsum(
        "d,d,abcdij->abcdij", y**2, e0, s3
    )

    rs3 += 0.5 * einsum(
        "i,i,abcdij->abcdij", y**2, e0, s3
    )

    rs3 += 0.5 * einsum(
        "j,j,abcdij->abcdij", y**2, e0, s3
    )

    rs3 += -0.5 * einsum(
        "a,a,abcdij->abcdij", x**2, e0, s3
    )

    rs3 += -0.5 * einsum(
        "b,b,abcdij->abcdij", x**2, e0, s3
    )

    rs3 += -0.5 * einsum(
        "c,c,abcdij->abcdij", x**2, e0, s3
    )

    rs3 += 0.5 * einsum(
        "k,k,k,kabcdijk->abcdij", e0, x, y, s4
    )

    rs4 += 0.5 * einsum(
        "i,i,abcdijkl->abcdijkl", y**2, e0, s4
    )

    rs4 += 0.5 * einsum(
        "j,j,abcdijkl->abcdijkl", y**2, e0, s4
    )

    rs4 += 0.5 * einsum(
        "k,k,abcdijkl->abcdijkl", y**2, e0, s4
    )

    rs4 += 0.5 * einsum(
        "l,l,abcdijkl->abcdijkl", y**2, e0, s4
    )

    rs4 += -0.5 * einsum(
        "a,a,abcdijkl->abcdijkl", x**2, e0, s4
    )

    rs4 += -0.5 * einsum(
        "b,b,abcdijkl->abcdijkl", x**2, e0, s4
    )

    rs4 += -0.5 * einsum(
        "c,c,abcdijkl->abcdijkl", x**2, e0, s4
    )

    rs4 += -0.5 * einsum(
        "d,d,abcdijkl->abcdijkl", x**2, e0, s4
    )

    rs4 += -0.5 * einsum(
        "c,d,j,l,abik,cdlj->abcdijkl", x, x, y, y, t2, eri
    )

    return rt1, rt2, rs1, rs2, rs3, rs4
