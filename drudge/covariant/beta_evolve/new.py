from numpy import zeros, einsum

def
    tau0 = zeros((na, na))

    tau0 += einsum(
        "c,cacb->ab", y**2, u
    )

    tau3 = zeros((na, na))

    tau3 += einsum(
        "c,ca,bc->ab", x, t1, tau0
    )

    rs1 = zeros((na, na))

    rs1 += Float('-0.5', precision=53) * einsum(
        "a,ba->ab", x, tau3
    )

    del tau3

    tau5 = zeros((na, na, na))

    tau5 += Float('1.0', precision=53) * einsum(
        "ab,bc->abc", t1, tau0
    )

    tau6 = zeros((na, na, na, na))

    tau6 += einsum(
        "d,i,di,iabc->abcd", x, x, tau0, t2
    )

    rs2 = zeros((na, na, na, na))

    rs2 += Float('0.5', precision=53) * einsum(
        "acdb->abcd", tau6
    )

    rs2 += Float('-0.5', precision=53) * einsum(
        "bcda->abcd", tau6
    )

    del tau6

    tau7 = zeros((na, na, na, na))

    tau7 += einsum(
        "d,i,id,abci->abcd", y, y, tau0, t2
    )

    rs2 += Float('0.5', precision=53) * einsum(
        "abcd->abcd", tau7
    )

    rs2 += Float('-0.5', precision=53) * einsum(
        "abdc->abcd", tau7
    )

    del tau7

    tau12 = zeros((na, na, na, na))

    tau12 += Float('1.0', precision=53) * einsum(
        "d,ab,cd->abcd", y, t1, tau0
    )

    tau12 += Float('-1.0', precision=53) * einsum(
        "b,ad,cb->abcd", y, t1, tau0
    )

    rs2 += Float('-0.5', precision=53) * einsum(
        "b,acbd->abcd", x, tau12
    )

    del tau12

    tau13 = zeros((na, na, na, na))

    tau13 += Float('-1.0', precision=53) * einsum(
        "d,ab,cd->abcd", y, t1, tau0
    )

    tau13 += Float('1.0', precision=53) * einsum(
        "b,ad,cb->abcd", y, t1, tau0
    )

    rs2 += Float('-0.5', precision=53) * einsum(
        "a,bcad->abcd", x, tau13
    )

    del tau13

    tau20 = zeros((na, na, na, na, na, na))

    tau20 += Float('1.0', precision=53) * einsum(
        "dj,acbi->abcdij", tau0, t2
    )

    tau21 = zeros((na, na, na, na, na, na))

    tau21 += Float('1.0', precision=53) * einsum(
        "i,ij,abcd->abcdij", x, tau0, t2
    )

    tau21 += Float('-1.0', precision=53) * einsum(
        "b,bj,aicd->abcdij", x, tau0, t2
    )

    tau21 += Float('1.0', precision=53) * einsum(
        "a,aj,bicd->abcdij", x, tau0, t2
    )

    rs3 = zeros((na, na, na, na, na, na))

    rs3 += Float('0.5', precision=53) * einsum(
        "i,abdjci->abcdij", y, tau21
    )

    rs3 += Float('-0.5', precision=53) * einsum(
        "j,abdicj->abcdij", y, tau21
    )

    del tau21

    tau22 = zeros((na, na, na, na, na, na))

    tau22 += Float('1.0', precision=53) * einsum(
        "i,ij,abcd->abcdij", x, tau0, t2
    )

    tau22 += Float('-1.0', precision=53) * einsum(
        "b,bj,aicd->abcdij", x, tau0, t2
    )

    rs3 += Float('-0.5', precision=53) * einsum(
        "d,abijcd->abcdij", y, tau22
    )

    del tau22

    rt1 = zeros((na, na))

    rt1 += Float('-0.5', precision=53) * einsum(
        "a,b,ab->ab", x, y, tau0
    )

    rs1 += Float('0.5', precision=53) * einsum(
        "c,d,dc,cabd->ab", x, y, tau0, t2
    )

    del tau0

    tau1 = zeros((na))

    tau1 += einsum(
        "b,abab->a", y**2, u
    )

    tau2 = 0

    tau2 += einsum(
        "a,a->", y**2, tau1
    )

    del tau1

    rs1 += -Float('0.25', precision=53)*tau2 * einsum(
        "ab->ab", t1
    )

    rs2 += -Float('0.25', precision=53)*tau2 * einsum(
        "abcd->abcd", t2
    )

    del tau2

    tau4 = zeros((na, na, na, na))

    tau4 += Float('1.0', precision=53) * einsum(
        "d,ab,cbda->abcd", y, t1, u
    )

    tau4 += Float('-0.5', precision=53) * einsum(
        "i,iadb,cbai->abcd", x, t2, u
    )

    rs1 += Float('-0.5', precision=53) * einsum(
        "a,c,d,cdab->ab", x, x, y, tau4
    )

    del tau4

    tau5 += Float('0.5', precision=53) * einsum(
        "d,i,daib,bicd->abc", x, y, t2, u
    )

    rs1 += Float('0.5', precision=53) * einsum(
        "b,c,acb->ab", y, y, tau5
    )

    del tau5

    tau8 = zeros((na, na, na, na, na))

    tau8 += Float('-0.5', precision=53) * einsum(
        "j,jabc,idaj->abcdi", x, t2, u
    )

    tau8 += Float('-1.0', precision=53) * einsum(
        "c,ab,idca->abcdi", y, t1, u
    )

    tau8 += Float('1.0', precision=53) * einsum(
        "b,ac,idba->abcdi", y, t1, u
    )

    rs2 += Float('0.5', precision=53) * einsum(
        "a,b,i,icdab->abcd", x, x, x, tau8
    )

    del tau8

    tau9 = zeros((na, na, na, na, na))

    tau9 += Float('-0.5', precision=53) * einsum(
        "j,abjc,cjid->abcdi", y, t2, u
    )

    tau9 += Float('1.0', precision=53) * einsum(
        "b,ac,bcid->abcdi", x, t1, u
    )

    tau9 += Float('-1.0', precision=53) * einsum(
        "a,bc,acid->abcdi", x, t1, u
    )

    rs2 += Float('-0.5', precision=53) * einsum(
        "c,d,i,abidc->abcd", y, y, y, tau9
    )

    del tau9

    tau10 = zeros((na, na, na, na, na, na))

    tau10 += Float('-1.0', precision=53) * einsum(
        "j,abcd,idja->abcdij", y, t2, u
    )

    tau10 += Float('1.0', precision=53) * einsum(
        "c,abjd,idca->abcdij", y, t2, u
    )

    rs2 += Float('-0.5', precision=53) * einsum(
        "a,i,j,ibdjac->abcd", x, x, y, tau10
    )

    del tau10

    tau11 = zeros((na, na, na, na, na, na))

    tau11 += Float('1.0', precision=53) * einsum(
        "j,abcd,idja->abcdij", y, t2, u
    )

    tau11 += Float('-1.0', precision=53) * einsum(
        "c,abjd,idca->abcdij", y, t2, u
    )

    rs2 += Float('-0.5', precision=53) * einsum(
        "b,i,j,iadjbc->abcd", x, x, y, tau11
    )

    del tau11

    tau14 = zeros((na, na, na, na, na, na, na))

    tau14 += Float('1.0', precision=53) * einsum(
        "k,abcd,jika->abcdijk", y, t2, u
    )

    tau14 += Float('-1.0', precision=53) * einsum(
        "d,abck,jida->abcdijk", y, t2, u
    )

    tau14 += Float('1.0', precision=53) * einsum(
        "c,abdk,jica->abcdijk", y, t2, u
    )

    rs3 += Float('0.5', precision=53) * einsum(
        "a,b,k,kcdiabj->abcdij", x, x, x, tau14
    )

    del tau14

    tau15 = zeros((na, na, na, na, na, na, na))

    tau15 += Float('-1.0', precision=53) * einsum(
        "k,abcd,jika->abcdijk", y, t2, u
    )

    tau15 += Float('1.0', precision=53) * einsum(
        "d,abck,jida->abcdijk", y, t2, u
    )

    tau15 += Float('1.0', precision=53) * einsum(
        "c,abkd,jica->abcdijk", y, t2, u
    )

    rs3 += Float('-0.5', precision=53) * einsum(
        "a,c,k,kbdjaci->abcdij", x, x, x, tau15
    )

    rs3 += Float('-0.5', precision=53) * einsum(
        "b,c,k,kadjcbi->abcdij", x, x, x, tau15
    )

    del tau15

    tau16 = zeros((na, na, na, na, na, na, na))

    tau16 += Float('-1.0', precision=53) * einsum(
        "i,abcd,dikj->abcdijk", x, t2, u
    )

    tau16 += Float('1.0', precision=53) * einsum(
        "b,aicd,dbkj->abcdijk", x, t2, u
    )

    tau16 += Float('-1.0', precision=53) * einsum(
        "a,bicd,dakj->abcdijk", x, t2, u
    )

    tau17 = zeros((na, na, na, na, na, na))

    tau17 += Float('0.5', precision=53) * einsum(
        "a,b,k,cdikjba->abcdij", y, y, y, tau16
    )

    del tau16

    rs3 += Float('1.0', precision=53) * einsum(
        "djabic->abcdij", tau17
    )

    rs3 += Float('-1.0', precision=53) * einsum(
        "diabjc->abcdij", tau17
    )

    rs3 += Float('-1.0', precision=53) * einsum(
        "ijabdc->abcdij", tau17
    )

    del tau17

    tau18 = zeros((na, na, na, na, na, na))

    tau18 += Float('1.0', precision=53) * einsum(
        "i,ab,dcij->abcdij", y, t1, u
    )

    tau18 += Float('-1.0', precision=53) * einsum(
        "b,ai,dcbj->abcdij", y, t1, u
    )

    rs3 += Float('-0.5', precision=53) * einsum(
        "a,b,j,cdbaij->abcdij", x, x, y, tau18
    )

    rs3 += Float('-0.5', precision=53) * einsum(
        "a,c,i,bdcaji->abcdij", x, x, y, tau18
    )

    rs3 += Float('-0.5', precision=53) * einsum(
        "b,c,j,adcbij->abcdij", x, x, y, tau18
    )

    del tau18

    tau19 = zeros((na, na, na, na, na, na))

    tau19 += Float('1.0', precision=53) * einsum(
        "c,ab,dcji->abcdij", x, t1, u
    )

    tau19 += Float('-1.0', precision=53) * einsum(
        "a,cb,daji->abcdij", x, t1, u
    )

    rs3 += Float('0.5', precision=53) * einsum(
        "b,d,i,cjabid->abcdij", x, y, y, tau19
    )

    del tau19

    tau20 += Float('-1.0', precision=53) * einsum(
        "c,i,ab,dcij->abcdij", x, y, t1, u
    )

    rs3 += Float('-0.5', precision=53) * einsum(
        "a,d,bicajd->abcdij", x, y, tau20
    )

    del tau20

    tau23 = zeros((na, na, na, na, na, na, na, na))

    tau23 += Float('1.0', precision=53) * einsum(
        "k,abcd,ijkl->abcdijkl", y, t2, u
    )

    tau23 += Float('-1.0', precision=53) * einsum(
        "d,abck,ijdl->abcdijkl", y, t2, u
    )

    tau23 += Float('1.0', precision=53) * einsum(
        "c,abdk,ijcl->abcdijkl", y, t2, u
    )

    rs4 = zeros((na, na, na, na, na, na, na, na))

    rs4 += Float('0.5', precision=53) * einsum(
        "a,b,i,cdjkabli->abcdijkl", x, x, y, tau23
    )

    del tau23

    tau24 = zeros((na, na, na, na, na, na, na, na))

    tau24 += Float('1.0', precision=53) * einsum(
        "i,abcd,jilk->abcdijkl", x, t2, u
    )

    tau24 += Float('-1.0', precision=53) * einsum(
        "b,aicd,jblk->abcdijkl", x, t2, u
    )

    tau24 += Float('1.0', precision=53) * einsum(
        "a,bicd,jalk->abcdijkl", x, t2, u
    )

    rs4 += Float('0.5', precision=53) * einsum(
        "a,j,l,bcikdalj->abcdijkl", x, y, y, tau24
    )

    rs4 += Float('-0.5', precision=53) * einsum(
        "a,j,k,bcildakj->abcdijkl", x, y, y, tau24
    )

    rs4 += Float('-0.5', precision=53) * einsum(
        "a,k,l,bcijdalk->abcdijkl", x, y, y, tau24
    )

    del tau24

    tau25 = zeros((na, na, na, na, na, na, na, na))

    tau25 += Float('1.0', precision=53) * einsum(
        "k,abcd,ijkl->abcdijkl", y, t2, u
    )

    tau25 += Float('1.0', precision=53) * einsum(
        "d,abkc,ijdl->abcdijkl", y, t2, u
    )

    tau25 += Float('-1.0', precision=53) * einsum(
        "c,abkd,ijcl->abcdijkl", y, t2, u
    )

    rs4 += Float('-0.5', precision=53) * einsum(
        "a,c,i,bdklacji->abcdijkl", x, x, y, tau25
    )

    rs4 += Float('0.5', precision=53) * einsum(
        "b,d,i,ackldbji->abcdijkl", x, x, y, tau25
    )

    del tau25

    tau26 = zeros((na, na, na, na, na, na, na, na))

    tau26 += Float('-1.0', precision=53) * einsum(
        "k,abcd,ijkl->abcdijkl", y, t2, u
    )

    tau26 += Float('1.0', precision=53) * einsum(
        "d,abck,ijdl->abcdijkl", y, t2, u
    )

    tau26 += Float('1.0', precision=53) * einsum(
        "c,abkd,ijcl->abcdijkl", y, t2, u
    )

    rs4 += Float('0.5', precision=53) * einsum(
        "a,d,i,bcjladki->abcdijkl", x, x, y, tau26
    )

    rs4 += Float('-0.5', precision=53) * einsum(
        "b,c,i,adjlcbki->abcdijkl", x, x, y, tau26
    )

    rs4 += Float('-0.5', precision=53) * einsum(
        "c,d,i,abjldcki->abcdijkl", x, x, y, tau26
    )

    del tau26

    tau27 = zeros((na, na, na, na, na, na, na, na))

    tau27 += Float('1.0', precision=53) * einsum(
        "i,abcd,jilk->abcdijkl", x, t2, u
    )

    tau27 += Float('-1.0', precision=53) * einsum(
        "b,aicd,jblk->abcdijkl", x, t2, u
    )

    rs4 += Float('-0.5', precision=53) * einsum(
        "b,j,l,acikdblj->abcdijkl", x, y, y, tau27
    )

    rs4 += Float('0.5', precision=53) * einsum(
        "c,j,k,adilbckj->abcdijkl", x, y, y, tau27
    )

    rs4 += Float('0.5', precision=53) * einsum(
        "c,k,l,adijbclk->abcdijkl", x, y, y, tau27
    )

    del tau27

    tau28 = zeros((na, na, na, na, na, na, na, na))

    tau28 += Float('1.0', precision=53) * einsum(
        "k,abcd,ijkl->abcdijkl", y, t2, u
    )

    tau28 += Float('-1.0', precision=53) * einsum(
        "d,abck,ijdl->abcdijkl", y, t2, u
    )

    rs4 += Float('-0.5', precision=53) * einsum(
        "b,d,k,acildbjk->abcdijkl", x, x, y, tau28
    )

    del tau28

    rt1 += Float('0.5', precision=53) * einsum(
        "b,b,ab->ab", y**2, e0, t1
    )

    rt1 += Float('-0.5', precision=53) * einsum(
        "a,a,ab->ab", x**2, e0, t1
    )

    rt1 += Float('0.5', precision=53) * einsum(
        "c,c,c,cabc->ab", e0, x, y, t2
    )

    rt2 = zeros((na, na, na, na))

    rt2 += Float('0.5', precision=53) * einsum(
        "c,c,abcd->abcd", y**2, e0, t2
    )

    rt2 += Float('0.5', precision=53) * einsum(
        "d,d,abcd->abcd", y**2, e0, t2
    )

    rt2 += Float('-0.5', precision=53) * einsum(
        "a,a,abcd->abcd", x**2, e0, t2
    )

    rt2 += Float('-0.5', precision=53) * einsum(
        "b,b,abcd->abcd", x**2, e0, t2
    )

    rt2 += Float('-0.5', precision=53) * einsum(
        "a,b,c,d,badc->abcd", x, x, y, y, u
    )

    rs1 += Float('0.5', precision=53) * einsum(
        "b,b,ab->ab", y**2, e0, s1
    )

    rs1 += Float('-0.5', precision=53) * einsum(
        "a,a,ab->ab", x**2, e0, s1
    )

    rs1 += Float('0.5', precision=53) * einsum(
        "c,c,c,cabc->ab", e0, x, y, s2
    )

    rs2 += Float('0.5', precision=53) * einsum(
        "c,c,abcd->abcd", y**2, e0, s2
    )

    rs2 += Float('0.5', precision=53) * einsum(
        "d,d,abcd->abcd", y**2, e0, s2
    )

    rs2 += Float('-0.5', precision=53) * einsum(
        "a,a,abcd->abcd", x**2, e0, s2
    )

    rs2 += Float('-0.5', precision=53) * einsum(
        "b,b,abcd->abcd", x**2, e0, s2
    )

    rs2 += Float('-0.5', precision=53) * einsum(
        "i,i,i,iabcdi->abcd", e0, x, y, s3
    )

    rs3 += Float('0.5', precision=53) * einsum(
        "d,d,abcdij->abcdij", y**2, e0, s3
    )

    rs3 += Float('0.5', precision=53) * einsum(
        "i,i,abcdij->abcdij", y**2, e0, s3
    )

    rs3 += Float('0.5', precision=53) * einsum(
        "j,j,abcdij->abcdij", y**2, e0, s3
    )

    rs3 += Float('-0.5', precision=53) * einsum(
        "a,a,abcdij->abcdij", x**2, e0, s3
    )

    rs3 += Float('-0.5', precision=53) * einsum(
        "b,b,abcdij->abcdij", x**2, e0, s3
    )

    rs3 += Float('-0.5', precision=53) * einsum(
        "c,c,abcdij->abcdij", x**2, e0, s3
    )

    rs3 += Float('0.5', precision=53) * einsum(
        "k,k,k,kabcdijk->abcdij", e0, x, y, s4
    )

    rs4 += Float('0.5', precision=53) * einsum(
        "i,i,abcdijkl->abcdijkl", y**2, e0, s4
    )

    rs4 += Float('0.5', precision=53) * einsum(
        "j,j,abcdijkl->abcdijkl", y**2, e0, s4
    )

    rs4 += Float('0.5', precision=53) * einsum(
        "k,k,abcdijkl->abcdijkl", y**2, e0, s4
    )

    rs4 += Float('0.5', precision=53) * einsum(
        "l,l,abcdijkl->abcdijkl", y**2, e0, s4
    )

    rs4 += Float('-0.5', precision=53) * einsum(
        "a,a,abcdijkl->abcdijkl", x**2, e0, s4
    )

    rs4 += Float('-0.5', precision=53) * einsum(
        "b,b,abcdijkl->abcdijkl", x**2, e0, s4
    )

    rs4 += Float('-0.5', precision=53) * einsum(
        "c,c,abcdijkl->abcdijkl", x**2, e0, s4
    )

    rs4 += Float('-0.5', precision=53) * einsum(
        "d,d,abcdijkl->abcdijkl", x**2, e0, s4
    )

    rs4 += Float('-0.5', precision=53) * einsum(
        "c,d,j,l,abik,cdlj->abcdijkl", x, x, y, y, t2, u
    )

