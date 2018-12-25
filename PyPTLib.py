from numpy import einsum, zeros

def covbetapt2(e0, eri, t0, t1, t2, s0, s1, s2, s3, s4, x, y):
    """
    Function to return the residuals at various ranks corresponding 
    to 1st and 2nd order correction in the wavefunction
    
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
        rt1     ::      residual that drives t1 evolution 
        rt2     ::      residual that drives t2 evolution 
        rs1     ::      residual that drives s1 evolution 
        rs2     ::      residual that drives s2 evolution 
        rs3     ::      residual that drives s3 evolution 
        rs4     ::      residual that drives s4 evolution 
    """
    na = len(x)

    tau0 = zeros((na))

    tau0 += einsum(
        "b,baba->a", y**2, eri
    )

    tau1 = 0

    tau1 += einsum(
        "a,a->", y**2, tau0
    )

    rt0 = 0

    rt0 -= tau1 / 4

    rs1 = zeros((na, na))

    rs1 -= tau1 * einsum(
        "ab->ab", t1
    ) / 4

    rs2 = zeros((na, na, na, na))

    rs2 -= tau1 * einsum(
        "abcd->abcd", t2
    ) / 4

    del tau1

    tau5 = zeros((na))

    tau5 += t0 * einsum(
        "a->a", tau0
    ) / 2

    del tau0

    tau2 = zeros((na, na))

    tau2 += einsum(
        "c,cacb->ab", y**2, eri
    )

    tau6 = zeros((na, na))

    tau6 += einsum(
        "c,ca,bc->ab", x, t1, tau2
    )

    rs1 -= einsum(
        "a,ba->ab", x, tau6
    ) / 2

    del tau6

    tau8 = zeros((na, na, na))

    tau8 += 2 * einsum(
        "ab,bc->abc", t1, tau2
    )

    tau9 = zeros((na, na, na, na))

    tau9 += einsum(
        "d,i,di,iabc->abcd", x, x, tau2, t2
    )

    rs2 += einsum(
        "acdb->abcd", tau9
    ) / 2

    rs2 -= einsum(
        "bcda->abcd", tau9
    ) / 2

    del tau9

    tau10 = zeros((na, na, na, na))

    tau10 += einsum(
        "d,i,id,abci->abcd", y, y, tau2, t2
    )

    rs2 += einsum(
        "abcd->abcd", tau10
    ) / 2

    rs2 -= einsum(
        "abdc->abcd", tau10
    ) / 2

    del tau10

    tau15 = zeros((na, na, na, na))

    tau15 += einsum(
        "d,ab,cd->abcd", y, t1, tau2
    )

    tau15 -= einsum(
        "b,ad,cb->abcd", y, t1, tau2
    )

    rs2 -= einsum(
        "b,acbd->abcd", x, tau15
    ) / 2

    rs2 -= einsum(
        "a,bdac->abcd", x, tau15
    ) / 2

    del tau15

    tau25 = zeros((na, na, na, na, na, na))

    tau25 += einsum(
        "dj,acbi->abcdij", tau2, t2
    )

    tau26 = zeros((na, na, na, na, na, na))

    tau26 += einsum(
        "i,ij,abcd->abcdij", x, tau2, t2
    )

    tau26 += einsum(
        "b,bj,iacd->abcdij", x, tau2, t2
    )

    tau26 -= einsum(
        "a,aj,ibcd->abcdij", x, tau2, t2
    )

    rs3 = zeros((na, na, na, na, na, na))

    rs3 += einsum(
        "i,bcdjai->abcdij", y, tau26
    ) / 2

    del tau26

    tau27 = zeros((na, na, na, na, na, na))

    tau27 -= einsum(
        "i,ij,abcd->abcdij", x, tau2, t2
    )

    tau27 += einsum(
        "b,bj,aicd->abcdij", x, tau2, t2
    )

    tau27 += einsum(
        "a,aj,ibcd->abcdij", x, tau2, t2
    )

    rs3 -= einsum(
        "j,acdibj->abcdij", y, tau27
    ) / 2

    del tau27

    tau28 = zeros((na, na, na, na, na, na))

    tau28 += einsum(
        "i,ij,abcd->abcdij", x, tau2, t2
    )

    tau28 -= einsum(
        "b,bj,aicd->abcdij", x, tau2, t2
    )

    rs3 -= einsum(
        "d,abijcd->abcdij", y, tau28
    ) / 2

    del tau28

    rt1 = zeros((na, na))

    rt1 -= einsum(
        "a,b,ab->ab", x, y, tau2
    ) / 2

    rs1 += einsum(
        "c,d,dc,cabd->ab", x, y, tau2, t2
    ) / 2

    del tau2

    tau3 = zeros((na, na))

    tau3 += einsum(
        "c,d,abcd,dcba->ab", y, y, t2, eri
    )

    tau4 = zeros((na))

    tau4 += einsum(
        "b,ab->a", x, tau3
    )

    del tau3

    tau4 += 4 * einsum(
        "a,a,aa->a", e0, y, s1
    )

    rs0 = 0

    rs0 -= einsum(
        "a,a->", tau4, x
    ) / 8

    del tau4

    tau5 += einsum(
        "b,c,bc,acab->a", x, y, t1, eri
    )

    rs0 -= einsum(
        "a,a->", y**2, tau5
    ) / 2

    del tau5

    tau7 = zeros((na, na, na, na))

    tau7 += 2 * einsum(
        "d,ab,cbda->abcd", y, t1, eri
    )

    tau7 -= einsum(
        "i,iadb,cbai->abcd", x, t2, eri
    )

    rs1 -= einsum(
        "a,c,d,cdab->ab", x, x, y, tau7
    ) / 4

    del tau7

    tau8 += einsum(
        "d,i,daib,bicd->abc", x, y, t2, eri
    )

    rs1 += einsum(
        "b,c,acb->ab", y, y, tau8
    ) / 4

    del tau8

    tau11 = zeros((na, na, na, na, na))

    tau11 -= einsum(
        "j,jabc,idaj->abcdi", x, t2, eri
    )

    tau11 -= 2 * einsum(
        "c,ab,idca->abcdi", y, t1, eri
    )

    tau11 += 2 * einsum(
        "b,ac,idba->abcdi", y, t1, eri
    )

    rs2 += einsum(
        "a,b,i,icdab->abcd", x, x, x, tau11
    ) / 4

    del tau11

    tau12 = zeros((na, na, na, na, na))

    tau12 -= einsum(
        "j,abjc,cjid->abcdi", y, t2, eri
    )

    tau12 += 2 * einsum(
        "b,ac,bcid->abcdi", x, t1, eri
    )

    tau12 -= 2 * einsum(
        "a,bc,acid->abcdi", x, t1, eri
    )

    rs2 -= einsum(
        "c,d,i,abidc->abcd", y, y, y, tau12
    ) / 4

    del tau12

    tau13 = zeros((na, na, na, na, na, na))

    tau13 += einsum(
        "j,abcd,idja->abcdij", y, t2, eri
    )

    tau13 -= einsum(
        "c,abjd,idca->abcdij", y, t2, eri
    )

    rs2 -= einsum(
        "a,i,j,ibcjad->abcd", x, x, y, tau13
    ) / 2

    del tau13

    tau14 = zeros((na, na, na, na, na, na))

    tau14 -= einsum(
        "j,abcd,idja->abcdij", y, t2, eri
    )

    tau14 += einsum(
        "c,abjd,idca->abcdij", y, t2, eri
    )

    rs2 -= einsum(
        "b,i,j,iacjbd->abcd", x, x, y, tau14
    ) / 2

    del tau14

    tau16 = zeros((na, na, na, na, na, na, na))

    tau16 += einsum(
        "k,abcd,jika->abcdijk", y, t2, eri
    )

    tau16 += einsum(
        "d,abkc,jida->abcdijk", y, t2, eri
    )

    tau16 -= einsum(
        "c,abkd,jica->abcdijk", y, t2, eri
    )

    rs3 += einsum(
        "a,b,k,kcijabd->abcdij", x, x, x, tau16
    ) / 2

    del tau16

    tau17 = zeros((na, na, na, na, na, na, na))

    tau17 += einsum(
        "k,abcd,jika->abcdijk", y, t2, eri
    )

    tau17 -= einsum(
        "d,abck,jida->abcdijk", y, t2, eri
    )

    tau17 += einsum(
        "c,abdk,jica->abcdijk", y, t2, eri
    )

    rs3 -= einsum(
        "a,c,k,kbdiacj->abcdij", x, x, x, tau17
    ) / 2

    del tau17

    tau18 = zeros((na, na, na, na, na, na, na))

    tau18 -= einsum(
        "k,abcd,jika->abcdijk", y, t2, eri
    )

    tau18 += einsum(
        "d,abck,jida->abcdijk", y, t2, eri
    )

    tau18 += einsum(
        "c,abkd,jica->abcdijk", y, t2, eri
    )

    rs3 -= einsum(
        "b,c,k,kadjcbi->abcdij", x, x, x, tau18
    ) / 2

    del tau18

    tau19 = zeros((na, na, na, na, na, na, na))

    tau19 -= einsum(
        "i,abcd,dikj->abcdijk", x, t2, eri
    )

    tau19 += einsum(
        "b,aicd,dbkj->abcdijk", x, t2, eri
    )

    tau19 -= einsum(
        "a,bicd,dakj->abcdijk", x, t2, eri
    )

    tau20 = zeros((na, na, na, na, na, na))

    tau20 += einsum(
        "a,b,k,cdikjba->abcdij", y, y, y, tau19
    )

    del tau19

    rs3 += einsum(
        "djabic->abcdij", tau20
    ) / 2

    rs3 -= einsum(
        "diabjc->abcdij", tau20
    ) / 2

    del tau20

    tau21 = zeros((na, na, na, na, na, na, na))

    tau21 += einsum(
        "i,abcd,dikj->abcdijk", x, t2, eri
    )

    tau21 -= einsum(
        "b,aicd,dbkj->abcdijk", x, t2, eri
    )

    tau21 -= einsum(
        "a,ibcd,dakj->abcdijk", x, t2, eri
    )

    rs3 -= einsum(
        "i,j,k,acdkbji->abcdij", y, y, y, tau21
    ) / 2

    del tau21

    tau22 = zeros((na, na, na, na, na, na))

    tau22 -= einsum(
        "i,ab,dcij->abcdij", y, t1, eri
    )

    tau22 += einsum(
        "b,ai,dcbj->abcdij", y, t1, eri
    )

    rs3 -= einsum(
        "a,b,j,cibadj->abcdij", x, x, y, tau22
    ) / 2

    rs3 -= einsum(
        "b,c,j,aicbdj->abcdij", x, x, y, tau22
    ) / 2

    del tau22

    tau23 = zeros((na, na, na, na, na, na))

    tau23 += einsum(
        "i,ab,dcij->abcdij", y, t1, eri
    )

    tau23 -= einsum(
        "b,ai,dcbj->abcdij", y, t1, eri
    )

    rs3 -= einsum(
        "a,c,i,bdcaji->abcdij", x, x, y, tau23
    ) / 2

    del tau23

    tau24 = zeros((na, na, na, na, na, na))

    tau24 += einsum(
        "c,ab,dcji->abcdij", x, t1, eri
    )

    tau24 -= einsum(
        "a,cb,daji->abcdij", x, t1, eri
    )

    rs3 += einsum(
        "b,d,i,cjabid->abcdij", x, y, y, tau24
    ) / 2

    del tau24

    tau25 -= einsum(
        "c,i,ab,dcij->abcdij", x, y, t1, eri
    )

    rs3 -= einsum(
        "a,d,bicajd->abcdij", x, y, tau25
    ) / 2

    del tau25

    tau29 = zeros((na, na, na, na, na, na, na, na))

    tau29 += einsum(
        "k,abcd,ijkl->abcdijkl", y, t2, eri
    )

    tau29 += einsum(
        "d,abkc,ijdl->abcdijkl", y, t2, eri
    )

    tau29 -= einsum(
        "c,abkd,ijcl->abcdijkl", y, t2, eri
    )

    rs4 = zeros((na, na, na, na, na, na, na, na))

    rs4 += einsum(
        "a,b,i,cdklabji->abcdijkl", x, x, y, tau29
    ) / 2

    rs4 -= einsum(
        "b,c,i,adklcbji->abcdijkl", x, x, y, tau29
    ) / 2

    rs4 += einsum(
        "b,d,i,ackldbji->abcdijkl", x, x, y, tau29
    ) / 2

    rs4 -= einsum(
        "c,d,i,abkldcji->abcdijkl", x, x, y, tau29
    ) / 2

    del tau29

    tau30 = zeros((na, na, na, na, na, na, na, na))

    tau30 += einsum(
        "i,abcd,jilk->abcdijkl", x, t2, eri
    )

    tau30 -= einsum(
        "b,aicd,jblk->abcdijkl", x, t2, eri
    )

    tau30 += einsum(
        "a,bicd,jalk->abcdijkl", x, t2, eri
    )

    rs4 += einsum(
        "a,j,l,bcikdalj->abcdijkl", x, y, y, tau30
    ) / 2

    rs4 -= einsum(
        "a,j,k,bcildakj->abcdijkl", x, y, y, tau30
    ) / 2

    del tau30

    tau31 = zeros((na, na, na, na, na, na, na, na))

    tau31 += einsum(
        "k,abcd,ijkl->abcdijkl", y, t2, eri
    )

    tau31 -= einsum(
        "d,abck,ijdl->abcdijkl", y, t2, eri
    )

    tau31 += einsum(
        "c,abdk,ijcl->abcdijkl", y, t2, eri
    )

    rs4 -= einsum(
        "a,c,i,bdjkacli->abcdijkl", x, x, y, tau31
    ) / 2

    del tau31

    tau32 = zeros((na, na, na, na, na, na, na, na))

    tau32 -= einsum(
        "i,abcd,jilk->abcdijkl", x, t2, eri
    )

    tau32 += einsum(
        "b,aicd,jblk->abcdijkl", x, t2, eri
    )

    tau32 += einsum(
        "a,ibcd,jalk->abcdijkl", x, t2, eri
    )

    rs4 -= einsum(
        "a,k,l,bdijcalk->abcdijkl", x, y, y, tau32
    ) / 2

    del tau32

    tau33 = zeros((na, na, na, na, na, na, na, na))

    tau33 -= einsum(
        "k,abcd,ijkl->abcdijkl", y, t2, eri
    )

    tau33 += einsum(
        "d,abck,ijdl->abcdijkl", y, t2, eri
    )

    tau33 += einsum(
        "c,abkd,ijcl->abcdijkl", y, t2, eri
    )

    rs4 += einsum(
        "a,d,i,bcjladki->abcdijkl", x, x, y, tau33
    ) / 2

    del tau33

    tau34 = zeros((na, na, na, na, na, na, na, na))

    tau34 -= einsum(
        "i,abcd,jilk->abcdijkl", x, t2, eri
    )

    tau34 += einsum(
        "b,aicd,jblk->abcdijkl", x, t2, eri
    )

    rs4 -= einsum(
        "b,j,l,adikcblj->abcdijkl", x, y, y, tau34
    ) / 2

    rs4 += einsum(
        "c,k,l,abijdclk->abcdijkl", x, y, y, tau34
    ) / 2

    del tau34

    tau35 = zeros((na, na, na, na, na, na, na, na))

    tau35 += einsum(
        "i,abcd,jilk->abcdijkl", x, t2, eri
    )

    tau35 -= einsum(
        "b,aicd,jblk->abcdijkl", x, t2, eri
    )

    rs4 += einsum(
        "c,j,k,adilbckj->abcdijkl", x, y, y, tau35
    ) / 2

    del tau35

    tau36 = zeros((na, na, na, na, na, na, na, na))

    tau36 -= einsum(
        "k,abcd,ijkl->abcdijkl", y, t2, eri
    )

    tau36 += einsum(
        "d,abck,ijdl->abcdijkl", y, t2, eri
    )

    rs4 -= einsum(
        "b,d,k,acijdblk->abcdijkl", x, x, y, tau36
    ) / 2

    del tau36

    rt0 -= einsum(
        "a,a,a,aa->", e0, x, y, t1
    ) / 2

    rt1 += einsum(
        "b,b,ab->ab", y**2, e0, t1
    ) / 2

    rt1 -= einsum(
        "a,a,ab->ab", x**2, e0, t1
    ) / 2

    rt1 += einsum(
        "c,c,c,cabc->ab", e0, x, y, t2
    ) / 2

    rt2 = zeros((na, na, na, na))

    rt2 += einsum(
        "c,c,abcd->abcd", y**2, e0, t2
    ) / 2

    rt2 += einsum(
        "d,d,abcd->abcd", y**2, e0, t2
    ) / 2

    rt2 -= einsum(
        "a,a,abcd->abcd", x**2, e0, t2
    ) / 2

    rt2 -= einsum(
        "b,b,abcd->abcd", x**2, e0, t2
    ) / 2

    rt2 -= einsum(
        "a,b,c,d,badc->abcd", x, x, y, y, eri
    ) / 2

    rs1 += einsum(
        "b,b,ab->ab", y**2, e0, s1
    ) / 2

    rs1 -= einsum(
        "a,a,ab->ab", x**2, e0, s1
    ) / 2

    rs1 += einsum(
        "c,c,c,cabc->ab", e0, x, y, s2
    ) / 2

    rs2 += einsum(
        "c,c,abcd->abcd", y**2, e0, s2
    ) / 2

    rs2 += einsum(
        "d,d,abcd->abcd", y**2, e0, s2
    ) / 2

    rs2 -= einsum(
        "a,a,abcd->abcd", x**2, e0, s2
    ) / 2

    rs2 -= einsum(
        "b,b,abcd->abcd", x**2, e0, s2
    ) / 2

    rs2 -= einsum(
        "i,i,i,iabcdi->abcd", e0, x, y, s3
    ) / 2

    rs3 += einsum(
        "d,d,abcdij->abcdij", y**2, e0, s3
    ) / 2

    rs3 += einsum(
        "i,i,abcdij->abcdij", y**2, e0, s3
    ) / 2

    rs3 += einsum(
        "j,j,abcdij->abcdij", y**2, e0, s3
    ) / 2

    rs3 -= einsum(
        "a,a,abcdij->abcdij", x**2, e0, s3
    ) / 2

    rs3 -= einsum(
        "b,b,abcdij->abcdij", x**2, e0, s3
    ) / 2

    rs3 -= einsum(
        "c,c,abcdij->abcdij", x**2, e0, s3
    ) / 2

    rs3 += einsum(
        "k,k,k,kabcdijk->abcdij", e0, x, y, s4
    ) / 2

    rs4 += einsum(
        "i,i,abcdijkl->abcdijkl", y**2, e0, s4
    ) / 2

    rs4 += einsum(
        "j,j,abcdijkl->abcdijkl", y**2, e0, s4
    ) / 2

    rs4 += einsum(
        "k,k,abcdijkl->abcdijkl", y**2, e0, s4
    ) / 2

    rs4 += einsum(
        "l,l,abcdijkl->abcdijkl", y**2, e0, s4
    ) / 2

    rs4 -= einsum(
        "a,a,abcdijkl->abcdijkl", x**2, e0, s4
    ) / 2

    rs4 -= einsum(
        "b,b,abcdijkl->abcdijkl", x**2, e0, s4
    ) / 2

    rs4 -= einsum(
        "c,c,abcdijkl->abcdijkl", x**2, e0, s4
    ) / 2

    rs4 -= einsum(
        "d,d,abcdijkl->abcdijkl", x**2, e0, s4
    ) / 2

    rs4 -= einsum(
        "c,d,j,l,abik,cdlj->abcdijkl", x, x, y, y, t2, eri
    ) / 2

    return rt0, rt1, rt2, rs0, rs1, rs2, rs3, rs4


def mpenergycov(e0, eri, t0, t1, t2, s0, s1, s2, s3, s4, x, y):
    """
    Function to return the 1st and 2nd order correction to energy
    and corresponding 2nd, 3rd and 4th order correction to the overlaps
    
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
        ov2     ::      2nd order correction to the < psi | psi >
        ov3     ::      3rd order correction to the < psi | psi >
        ov4     ::      4th order correction to the < psi | psi > <<--- NOTE: this is not exact
    """
    na = len(x)

    tau0 = zeros((na))

    tau0 += einsum(
        "b,baba->a", y**2, eri
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

    ov1 = 2*t0

    ov2 = 2*s0

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
        "abac->abc", eri
    )

    tau12 += einsum(
        "acab->abc", eri
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
        "abcd,badc->abcd", t2, eri
    )

    tau14 += einsum(
        "abdc,cdba->abcd", t2, eri
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

    return e1, e2, ov1, ov2, ov3, ov4

