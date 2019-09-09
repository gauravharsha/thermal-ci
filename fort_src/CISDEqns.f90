    Subroutine RCISD(Fock, ERI, T1, T2, En, R1, R2, NAO, NOcc)
        Implicit None
        Integer, parameter  :: pr = Selected_Real_Kind(15,307)

        Integer, Intent(In) :: NAO, NOcc
        Real (Kind=pr), Intent(In)  ::  Fock(NAO, NAO), ERI(NAO, NAO, NAO, NAO)
        Real (Kind=pr), Intent(In)  ::  T1(NOcc+1:NAO, NOcc)
        Real (Kind=pr), Intent(In)  ::  T2(NOcc+1:NAO, NOcc+1:NAO, NOcc, NOcc)
        Real (Kind=pr), Intent(Out) ::  En
        Real (Kind=pr), Intent(Out) ::  R1(NOcc, NOcc+1:NAO)
        Real (Kind=pr), Intent(Out) ::  R2(NOcc, NOcc, NOcc+1:NAO, NOcc+1:NAO)
        Real (Kind=pr)              ::  Ov1(NOcc, NOcc+1:NAO)
        Real (Kind=pr)              ::  Ov2(NOcc, NOcc, NOcc+1:NAO, NOcc+1:NAO)
        Integer         ::  a, b, c, d, i, j, k, l

        !print *, MaxVal(t2)
        !print *, MaxVal(abs(fock))
        !print *, MaxVal(abs(eri))
        !$omp parallel default(shared)

        !$omp single
        
        en = 0.0_pr
        
        !$omp end single
        

        !$omp do schedule(static) reduction(+:en)
        
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        en = en + ( &
                            3.0_pr * t2(a, b, i, j) * eri(i, j, a, b)&
                        )
                    end do
                end do
            end do
        end do
        
        !$omp end do

        !$omp do schedule(static) reduction(+:en)
        
        do i=1, nocc
            do a=nocc + 1, nao
                en = en + ( &
                    2.0_pr * fock(i, a) * t1(a, i)&
                )
            end do
        end do
        
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do a=nocc + 1, nao
                R1(i, a) = 0.0_pr
            end do
        end do
        !$omp end do
        

        !$omp do schedule(static)
        do i=1, nocc
            do a=nocc + 1, nao
                do j=1, nocc
                    do k=1, nocc
                        do b=nocc + 1, nao
                            R1(i, a) = R1(i, a) - ( &
                                6 * t2(a, b, j, k) * eri(j, k, i, b)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do a=nocc + 1, nao
        
                R1(i, a) = R1(i, a) + ( &
                    2 * fock(a, i)&
                )
        
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do a=nocc + 1, nao
                do j=1, nocc
                    R1(i, a) = R1(i, a) - ( &
                        2 * fock(j, i) * t1(a, j)&
                    )
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do a=nocc + 1, nao
                do j=1, nocc
                    do b=nocc + 1, nao
                        R1(i, a) = R1(i, a) - ( &
                            2 * t1(b, j) * eri(a, j, b, i)&
                        )
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do a=nocc + 1, nao
                do j=1, nocc
                    do b=nocc + 1, nao
                        R1(i, a) = R1(i, a) + ( &
                            4 * t1(b, j) * eri(a, j, i, b)&
                        )
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do a=nocc + 1, nao
                do j=1, nocc
                    do b=nocc + 1, nao
                        R1(i, a) = R1(i, a) + ( &
                            6 * fock(j, b) * t2(a, b, i, j)&
                        )
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do a=nocc + 1, nao
                do j=1, nocc
                    do b=nocc + 1, nao
                        do c=nocc + 1, nao
                            R1(i, a) = R1(i, a) + ( &
                                6 * t2(b, c, i, j) * eri(a, j, b, c)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do a=nocc + 1, nao
                do b=nocc + 1, nao
                    R1(i, a) = R1(i, a) + ( &
                        2 * fock(a, b) * t1(b, i)&
                    )
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        R2(i, j, a, b) = 0.0_pr
                    end do
                end do
            end do
        end do
        !$omp end do
        

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
        
                        R2(i, j, a, b) = R2(i, j, a, b) - ( &
                            2 * eri(a, b, j, i)&
                        )
        
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
        
                        R2(i, j, a, b) = R2(i, j, a, b) + ( &
                            4 * eri(a, b, i, j)&
                        )
        
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
        
                        R2(i, j, a, b) = R2(i, j, a, b) - ( &
                            2 * fock(a, j) * t1(b, i)&
                        )
        
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
        
                        R2(i, j, a, b) = R2(i, j, a, b) - ( &
                            2 * fock(b, i) * t1(a, j)&
                        )
        
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
        
                        R2(i, j, a, b) = R2(i, j, a, b) + ( &
                            4 * fock(a, i) * t1(b, j)&
                        )
        
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
        
                        R2(i, j, a, b) = R2(i, j, a, b) + ( &
                            4 * fock(b, j) * t1(a, i)&
                        )
        
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            R2(i, j, a, b) = R2(i, j, a, b) - ( &
                                6 * fock(k, i) * t2(a, b, k, j)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            R2(i, j, a, b) = R2(i, j, a, b) - ( &
                                6 * fock(k, j) * t2(a, b, i, k)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            R2(i, j, a, b) = R2(i, j, a, b) - ( &
                                4 * t1(a, k) * eri(b, k, j, i)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            R2(i, j, a, b) = R2(i, j, a, b) - ( &
                                4 * t1(b, k) * eri(a, k, i, j)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            R2(i, j, a, b) = R2(i, j, a, b) + ( &
                                2 * t1(a, k) * eri(b, k, i, j)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            R2(i, j, a, b) = R2(i, j, a, b) + ( &
                                2 * t1(b, k) * eri(a, k, j, i)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            do c=nocc + 1, nao
                                R2(i, j, a, b) = R2(i, j, a, b) - ( &
                                    12 * t2(b, c, k, j) * eri(a, k, i, c)&
                                )
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            do c=nocc + 1, nao
                                R2(i, j, a, b) = R2(i, j, a, b) - ( &
                                    6 * t2(a, c, i, k) * eri(b, k, c, j)&
                                )
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            do c=nocc + 1, nao
                                R2(i, j, a, b) = R2(i, j, a, b) - ( &
                                    6 * t2(a, c, k, j) * eri(b, k, c, i)&
                                )
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            do c=nocc + 1, nao
                                R2(i, j, a, b) = R2(i, j, a, b) - ( &
                                    6 * t2(b, c, i, k) * eri(a, k, j, c)&
                                )
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            do c=nocc + 1, nao
                                R2(i, j, a, b) = R2(i, j, a, b) + ( &
                                    6 * t2(a, c, k, j) * eri(b, k, i, c)&
                                )
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            do c=nocc + 1, nao
                                R2(i, j, a, b) = R2(i, j, a, b) + ( &
                                    6 * t2(b, c, i, k) * eri(a, k, c, j)&
                                )
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            do c=nocc + 1, nao
                                R2(i, j, a, b) = R2(i, j, a, b) + ( &
                                    6 * t2(b, c, k, j) * eri(a, k, c, i)&
                                )
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            do c=nocc + 1, nao
                                R2(i, j, a, b) = R2(i, j, a, b) + ( &
                                    12 * t2(a, c, i, k) * eri(b, k, j, c)&
                                )
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do c=nocc + 1, nao
                            do d=nocc + 1, nao
                                R2(i, j, a, b) = R2(i, j, a, b) + ( &
                                    6 * t2(c, d, i, j) * eri(a, b, c, d)&
                                )
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do c=nocc + 1, nao
                            R2(i, j, a, b) = R2(i, j, a, b) - ( &
                                6 * fock(a, c) * t2(b, c, i, j)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do c=nocc + 1, nao
                            R2(i, j, a, b) = R2(i, j, a, b) - ( &
                                2 * t1(c, i) * eri(a, b, j, c)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do c=nocc + 1, nao
                            R2(i, j, a, b) = R2(i, j, a, b) - ( &
                                2 * t1(c, j) * eri(a, b, c, i)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do c=nocc + 1, nao
                            R2(i, j, a, b) = R2(i, j, a, b) + ( &
                                4 * t1(c, i) * eri(a, b, c, j)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do c=nocc + 1, nao
                            R2(i, j, a, b) = R2(i, j, a, b) + ( &
                                4 * t1(c, j) * eri(a, b, i, c)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do c=nocc + 1, nao
                            R2(i, j, a, b) = R2(i, j, a, b) + ( &
                                6 * fock(b, c) * t2(a, c, i, j)&
                            )
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        do k=1, nocc
                            do l=1, nocc
                                R2(i, j, a, b) = R2(i, j, a, b) + ( &
                                    6 * t2(a, b, k, l) * eri(k, l, i, j)&
                                )
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do a=nocc + 1, nao
                Ov1(i, a) = 0.0_pr
            end do
        end do
        !$omp end do
        

        !$omp do schedule(static)
        do i=1, nocc
            do a=nocc + 1, nao
        
                Ov1(i, a) = Ov1(i, a) + ( &
                    2 * t1(a, i)&
                    )
        
            end do
        end do
        !$omp end do

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
                        Ov2(i, j, a, b) = 0.0_pr
                    end do
                end do
            end do
        end do
        !$omp end do
        

        !$omp do schedule(static)
        do i=1, nocc
            do j=1, nocc
                do a=nocc + 1, nao
                    do b=nocc + 1, nao
        
                        Ov2(i, j, a, b) = Ov2(i, j, a, b) + ( &
                            6 * t2(a, b, i, j)&
                            )
        
                    end do
                end do
            end do
        end do
        !$omp end do

        !$omp end parallel

        R1 = R1/( 2.0_pr*En )
        R2 = R2/( 6.0_pr*En )
        !R1 = R1 - En*Ov1
        !R2 = R2 - En*Ov2

    End Subroutine RCISD

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
