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
    tau0 = zeros((na, na))

    tau0 += einsum(
        "c,cacb->ab", y**2,eri
    )

    tau3 = zeros((na, na))

    tau3 += einsum(
        "c,ca,bc->ab", x, t1, tau0
    )

    rs1 = zeros((na, na))

    rs1 -= einsum(
        "a,ba->ab", x, tau3
    ) / 2

    del tau3

    tau5 = zeros((na, na, na))

    tau5 += 2 * einsum(
        "ab,bc->abc", t1, tau0
    )

    tau6 = zeros((na, na, na, na))

    tau6 += einsum(
        "d,i,di,iabc->abcd", x, x, tau0, t2
    )

    rs2 = zeros((na, na, na, na))

    rs2 += einsum(
        "acdb->abcd", tau6
    ) / 2

    rs2 -= einsum(
        "bcda->abcd", tau6
    ) / 2

    del tau6

    tau7 = zeros((na, na, na, na))

    tau7 += einsum(
        "d,i,id,abci->abcd", y, y, tau0, t2
    )

    rs2 += einsum(
        "abcd->abcd", tau7
    ) / 2

    rs2 -= einsum(
        "abdc->abcd", tau7
    ) / 2

    del tau7

    tau12 = zeros((na, na, na, na))

    tau12 += einsum(
        "d,ab,cd->abcd", y, t1, tau0
    )

    tau12 -= einsum(
        "b,ad,cb->abcd", y, t1, tau0
    )

    rs2 -= einsum(
        "b,acbd->abcd", x, tau12
    ) / 2

    rs2 -= einsum(
        "a,bdac->abcd", x, tau12
    ) / 2

    del tau12

    tau20 = zeros((na, na, na, na, na, na))

    tau20 += einsum(
        "dj,acbi->abcdij", tau0, t2
    )

    tau21 = zeros((na, na, na, na, na, na))

    tau21 -= einsum(
        "i,ij,abcd->abcdij", x, tau0, t2
    )

    tau21 += einsum(
        "b,bj,aicd->abcdij", x, tau0, t2
    )

    tau21 += einsum(
        "a,aj,ibcd->abcdij", x, tau0, t2
    )

    rs3 = zeros((na, na, na, na, na, na))

    rs3 += einsum(
        "i,acdjbi->abcdij", y, tau21
    ) / 2

    rs3 -= einsum(
        "j,acdibj->abcdij", y, tau21
    ) / 2

    del tau21

    tau22 = zeros((na, na, na, na, na, na))

    tau22 -= einsum(
        "i,ij,abcd->abcdij", x, tau0, t2
    )

    tau22 += einsum(
        "b,bj,aicd->abcdij", x, tau0, t2
    )

    rs3 -= einsum(
        "d,acijbd->abcdij", y, tau22
    ) / 2

    del tau22

    rt1 = zeros((na, na))

    rt1 -= einsum(
        "a,b,ab->ab", x, y, tau0
    ) / 2

    rs1 += einsum(
        "c,d,dc,cabd->ab", x, y, tau0, t2
    ) / 2

    del tau0

    tau1 = zeros((na))

    tau1 += einsum(
        "b,abab->a", y**2,eri
    )

    tau2 = 0

    tau2 += einsum(
        "a,a->", y**2, tau1
    )

    del tau1

    rs1 -= tau2 * einsum(
        "ab->ab", t1
    ) / 4

    rs2 -= tau2 * einsum(
        "abcd->abcd", t2
    ) / 4

    del tau2

    tau4 = zeros((na, na, na, na))

    tau4 += 2 * einsum(
        "d,ab,cbda->abcd", y, t1,eri
    )

    tau4 -= einsum(
        "i,iadb,cbai->abcd", x, t2,eri
    )

    rs1 -= einsum(
        "a,c,d,cdab->ab", x, x, y, tau4
    ) / 4

    del tau4

    tau5 += einsum(
        "d,i,daib,bicd->abc", x, y, t2,eri
    )

    rs1 += einsum(
        "b,c,acb->ab", y, y, tau5
    ) / 4

    del tau5

    tau8 = zeros((na, na, na, na, na))

    tau8 -= einsum(
        "j,jabc,idaj->abcdi", x, t2,eri
    )

    tau8 -= 2 * einsum(
        "c,ab,idca->abcdi", y, t1,eri
    )

    tau8 += 2 * einsum(
        "b,ac,idba->abcdi", y, t1,eri
    )

    rs2 += einsum(
        "a,b,i,icdab->abcd", x, x, x, tau8
    ) / 4

    del tau8

    tau9 = zeros((na, na, na, na, na))

    tau9 -= einsum(
        "j,abjc,cjid->abcdi", y, t2,eri
    )

    tau9 += 2 * einsum(
        "b,ac,bcid->abcdi", x, t1,eri
    )

    tau9 -= 2 * einsum(
        "a,bc,acid->abcdi", x, t1,eri
    )

    rs2 -= einsum(
        "c,d,i,abidc->abcd", y, y, y, tau9
    ) / 4

    del tau9

    tau10 = zeros((na, na, na, na, na, na))

    tau10 += einsum(
        "j,abcd,idja->abcdij", y, t2,eri
    )

    tau10 -= einsum(
        "c,abjd,idca->abcdij", y, t2,eri
    )

    rs2 -= einsum(
        "a,i,j,ibcjad->abcd", x, x, y, tau10
    ) / 2

    del tau10

    tau11 = zeros((na, na, na, na, na, na))

    tau11 -= einsum(
        "j,abcd,idja->abcdij", y, t2,eri
    )

    tau11 += einsum(
        "c,abjd,idca->abcdij", y, t2,eri
    )

    rs2 -= einsum(
        "b,i,j,iacjbd->abcd", x, x, y, tau11
    ) / 2

    del tau11

    tau13 = zeros((na, na, na, na, na, na, na))

    tau13 -= einsum(
        "k,abcd,jika->abcdijk", y, t2,eri
    )

    tau13 += einsum(
        "d,abck,jida->abcdijk", y, t2,eri
    )

    tau13 += einsum(
        "c,abkd,jica->abcdijk", y, t2,eri
    )

    tau14 = zeros((na, na, na, na, na, na))

    tau14 += einsum(
        "a,b,k,kcdiabj->abcdij", x, x, x, tau13
    )

    rs3 += einsum(
        "abcdji->abcdij", tau14
    ) / 2

    rs3 -= einsum(
        "acbdji->abcdij", tau14
    ) / 2

    del tau14

    rs3 -= einsum(
        "b,c,k,kadjcbi->abcdij", x, x, x, tau13
    ) / 2

    del tau13

    tau15 = zeros((na, na, na, na, na, na, na))

    tau15 -= einsum(
        "i,abcd,dikj->abcdijk", x, t2,eri
    )

    tau15 -= einsum(
        "b,iacd,dbkj->abcdijk", x, t2,eri
    )

    tau15 += einsum(
        "a,ibcd,dakj->abcdijk", x, t2,eri
    )

    rs3 += einsum(
        "d,j,k,bcikajd->abcdij", y, y, y, tau15
    ) / 2

    del tau15

    tau16 = zeros((na, na, na, na, na, na, na))

    tau16 += einsum(
        "i,abcd,dikj->abcdijk", x, t2,eri
    )

    tau16 -= einsum(
        "b,aicd,dbkj->abcdijk", x, t2,eri
    )

    tau16 -= einsum(
        "a,ibcd,dakj->abcdijk", x, t2,eri
    )

    tau17 = zeros((na, na, na, na, na, na))

    tau17 += einsum(
        "a,b,k,cdikjba->abcdij", y, y, y, tau16
    )

    del tau16

    rs3 -= einsum(
        "diacjb->abcdij", tau17
    ) / 2

    rs3 -= einsum(
        "ijacdb->abcdij", tau17
    ) / 2

    del tau17

    tau18 = zeros((na, na, na, na, na, na))

    tau18 -= einsum(
        "i,ab,dcij->abcdij", y, t1,eri
    )

    tau18 += einsum(
        "b,ai,dcbj->abcdij", y, t1,eri
    )

    rs3 -= einsum(
        "a,b,j,cibadj->abcdij", x, x, y, tau18
    ) / 2

    rs3 -= einsum(
        "a,c,i,bjcadi->abcdij", x, x, y, tau18
    ) / 2

    rs3 -= einsum(
        "b,c,j,aicbdj->abcdij", x, x, y, tau18
    ) / 2

    del tau18

    tau19 = zeros((na, na, na, na, na, na))

    tau19 += einsum(
        "c,ab,dcji->abcdij", x, t1,eri
    )

    tau19 -= einsum(
        "a,cb,daji->abcdij", x, t1,eri
    )

    rs3 += einsum(
        "b,d,i,cjabid->abcdij", x, y, y, tau19
    ) / 2

    del tau19

    tau20 -= einsum(
        "c,i,ab,dcij->abcdij", x, y, t1,eri
    )

    rs3 -= einsum(
        "a,d,bicajd->abcdij", x, y, tau20
    ) / 2

    del tau20

    tau23 = zeros((na, na, na, na, na, na, na, na))

    tau23 += einsum(
        "k,abcd,ijkl->abcdijkl", y, t2,eri
    )

    tau23 -= einsum(
        "d,abck,ijdl->abcdijkl", y, t2,eri
    )

    tau23 += einsum(
        "c,abdk,ijcl->abcdijkl", y, t2,eri
    )

    rs4 = zeros((na, na, na, na, na, na, na, na))

    rs4 += einsum(
        "a,b,i,cdjkabli->abcdijkl", x, x, y, tau23
    ) / 2

    rs4 -= einsum(
        "b,c,i,adjkcbli->abcdijkl", x, x, y, tau23
    ) / 2

    rs4 += einsum(
        "b,d,i,acjkdbli->abcdijkl", x, x, y, tau23
    ) / 2

    del tau23

    tau24 = zeros((na, na, na, na, na, na, na, na))

    tau24 += einsum(
        "i,abcd,jilk->abcdijkl", x, t2,eri
    )

    tau24 += einsum(
        "b,iacd,jblk->abcdijkl", x, t2,eri
    )

    tau24 -= einsum(
        "a,ibcd,jalk->abcdijkl", x, t2,eri
    )

    rs4 += einsum(
        "a,j,l,cdikbalj->abcdijkl", x, y, y, tau24
    ) / 2

    rs4 -= einsum(
        "a,j,k,cdilbakj->abcdijkl", x, y, y, tau24
    ) / 2

    rs4 -= einsum(
        "a,k,l,cdijbalk->abcdijkl", x, y, y, tau24
    ) / 2

    del tau24

    tau25 = zeros((na, na, na, na, na, na, na, na))

    tau25 -= einsum(
        "k,abcd,ijkl->abcdijkl", y, t2,eri
    )

    tau25 += einsum(
        "d,abck,ijdl->abcdijkl", y, t2,eri
    )

    tau25 += einsum(
        "c,abkd,ijcl->abcdijkl", y, t2,eri
    )

    rs4 -= einsum(
        "a,c,i,bdjlacki->abcdijkl", x, x, y, tau25
    ) / 2

    rs4 += einsum(
        "a,d,i,bcjladki->abcdijkl", x, x, y, tau25
    ) / 2

    del tau25

    tau26 = zeros((na, na, na, na, na, na, na, na))

    tau26 += einsum(
        "k,abcd,ijkl->abcdijkl", y, t2,eri
    )

    tau26 += einsum(
        "d,abkc,ijdl->abcdijkl", y, t2,eri
    )

    tau26 -= einsum(
        "c,abkd,ijcl->abcdijkl", y, t2,eri
    )

    rs4 -= einsum(
        "c,d,i,abkldcji->abcdijkl", x, x, y, tau26
    ) / 2

    del tau26

    tau27 = zeros((na, na, na, na, na, na, na, na))

    tau27 -= einsum(
        "i,abcd,jilk->abcdijkl", x, t2,eri
    )

    tau27 += einsum(
        "b,aicd,jblk->abcdijkl", x, t2,eri
    )

    rs4 -= einsum(
        "b,j,l,adikcblj->abcdijkl", x, y, y, tau27
    ) / 2

    rs4 += einsum(
        "c,j,k,abildckj->abcdijkl", x, y, y, tau27
    ) / 2

    rs4 += einsum(
        "c,k,l,abijdclk->abcdijkl", x, y, y, tau27
    ) / 2

    del tau27

    tau28 = zeros((na, na, na, na, na, na, na, na))

    tau28 -= einsum(
        "k,abcd,ijkl->abcdijkl", y, t2,eri
    )

    tau28 += einsum(
        "d,abck,ijdl->abcdijkl", y, t2,eri
    )

    rs4 -= einsum(
        "b,d,k,acijdblk->abcdijkl", x, x, y, tau28
    ) / 2

    del tau28

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
        "a,b,c,d,badc->abcd", x, x, y, y,eri
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
        "c,d,j,l,abik,cdlj->abcdijkl", x, x, y, y, t2,eri
    ) / 2

    return rt1, rt2, rs1, rs2, rs3, rs4

