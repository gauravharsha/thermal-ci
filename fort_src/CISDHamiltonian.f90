    Subroutine CISDHamiltonian(Fock, ERI, T1, T2, Energy, NAO, NOcc, NVir)
        Implicit None
        Integer, parameter  :: pr = Selected_Real_Kind(15,307)

        Integer, Intent(In) :: NAO, NOcc, Nvir
        Real (Kind=pr), Intent(In)  ::  Fock(NAO, NAO), ERI(NAO, NAO, NAO, NAO)
        Real (Kind=pr), Intent(In)  ::  T1(NOcc+1:NAO, NOcc)
        Real (Kind=pr), Intent(In)  ::  T2(NOcc+1:NAO, NOcc+1:NAO, NOcc, NOcc)
        Real (Kind=pr), Intent(Out) ::  Energy
        Real (Kind=pr), Allocatable ::  Ham
        Integer         ::  a, b, c, d, i, j, k, l, NDet

        !Get number of determinants and allocate the Hamiltonian
        NDet = 1 + NOcc*NVir + (NOcc*NVir*(NOcc-1)*(Nvir-1))/2
        Allocate(Ham(1:NDet,1:NDet))

        Call MakeHamSS(Ham, T1, T2, NAO, NOcc, NVir, NDet)

    End Subroutine CIResiduals

    Subroutine MakeHamSS(Ham, T1, T2, NAO, NOcc, NVir, NDet)
        Implicit None
        Integer, parameter  :: pr = Selected_Real_Kind(15,307)

        Integer, Intent(In) :: NAO, NOcc, Nvir, NDet
        Real (Kind=pr), Intent(In)  ::  Fock(NAO, NAO), ERI(NAO, NAO, NAO, NAO)
        Real (Kind=pr), Intent(In)  ::  T1(NOcc+1:NAO, NOcc)
        Real (Kind=pr), Intent(In)  ::  T2(NOcc+1:NAO, NOcc+1:NAO, NOcc, NOcc)
        Real (Kind=pr), Intent(Out) ::  Ham(NDet, NDet)
        Integer     ::  a, b, i, j, x, y

        ! The expression is simple
        !     < phi_j^b | H | phi_i^a > = delK(a,b)F(i,j) - delK(i,j)F(a,b) +
        !             2*u(bija) - u(biaj)

        Do i = 1, NOcc




    Subroutine A_dot_B(A,B,C,MA,NA,NB)
        ! Subroutine performs the following Matrix Multiplication
        !
        !      C = C + A . B
        !
        !  where A is of dimension (MA,NA)
        !    and B is of dimension (MB,NB)
        ! Obviously, NB = NA
        Implicit None
        Integer, parameter  :: pr = Selected_Real_Kind(15,307)
        Integer, Intent(In)           ::  MA,NA,NB
        Real (Kind=pr), Intent(In)    ::  A(MA,NA), B(NA,NB)
        Real (Kind=pr), Intent(InOut) ::  C(MA,NB)

        C = C + MatMul(A,B)
        ! Call DGEMM('N','N',MA,NB,NA,1.0_pr,A,MA,B,NA,1.0_pr,C,MA)

        Return

    End Subroutine A_dot_B
