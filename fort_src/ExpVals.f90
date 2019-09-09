      Subroutine EvalEnergy(OneH, ERI, T1, T2, NA, X, Y, Energy)
          Implicit None

          Integer, parameter  :: pr = Selected_Real_Kind(15,307)
          Integer, Intent(In) :: NA
          Real (Kind=pr), Intent(In) :: X(NA), Y(NA)
          Real (Kind=pr), Intent(In) :: T1(NA,NA), T2(NA,NA,NA,NA)
          Real (Kind=pr), Intent(In) :: OneH(NA,NA),  ERI(NA,NA,NA,NA)
          Real (Kind=pr), Intent(Out)   ::  Energy
          Real (Kind=pr)    ::  R0, Ovlp
          Real (Kind=pr)    ::  R1(NA,NA)
          Real (Kind=pr)    ::  R2(NA,NA,NA,NA)
          !Real (Kind=pr)    ::  Amat(NA,NA,NA,NA)
          Integer :: a, b, c, d, p, q, r, s

          !Real (Kind=pr) :: tau0
          !Real (Kind=pr) :: tau1
          !Real (Kind=pr), dimension(:), allocatable :: tau2
          !Real (Kind=pr) :: tau3
          !Real (Kind=pr) :: tau4
          !Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau5
          !Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau7
          !Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau8
          !Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau9
          !Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau10
          !Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau11

          ! Pre Processing
          ! Defining the Hamiltonian Matrix Elements

          Real (Kind=pr) ::  h0, h20(na,na), h02(na,na), h11(na,na)
          Real (Kind=pr) ::  h40(na,na,na,na), h04(na,na,na,na)
          Real (Kind=pr) ::  h31(na,na,na,na), h13(na,na,na,na)
          Real (Kind=pr) ::  h221(na,na,na,na)
          Real (Kind=pr) ::  h222(na,na,na,na)
          Real (Kind=pr) ::  scr1(na), scr2(na,na), delK(na,na)

          h0 = 0.0_pr
          delK = 0.0_pr

          do a=1, na
              h0 = h0 + y(a)*y(a)*oneh(a,a)
              delK(a,a) = 1.0_pr
              do b=1, na
                  scr2(a,b) = 0.0_pr
                  do c=1,na
                      scr1(c) = eri(a,c,b,c)
                  end do
                  scr2(a,b) = Sum(y*y*scr1) + oneh(a,b)
                  h0 = h0 + ( y(a)**2 * y(b)**2 * eri(a,b,a,b) )/2
              end do
          end do

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  h20(a,b) = x(a)*x(b)*scr2(a,b)
                  h02(a,b) = -y(a)*y(b)*scr2(a,b)
                  h11(a,b) = x(a)*y(b)*scr2(a,b)
                  do c=1, na
                      do d=1, na
                          h40(a,b,c,d) = eri(a,b,c,d)*x(a)*x(b)*x(c)*x(d)/4.0
                          h04(a,b,c,d) = eri(c,d,a,b)*y(a)*y(b)*y(c)*y(d)/4.0
                          h31(a,b,c,d) = -eri(a,b,c,d)*x(a)*x(b)*y(c)*x(d)/2.0
                          h13(a,b,c,d) = -eri(a,d,b,c)*x(a)*y(b)*y(c)*y(d)/2.0
                          h221(a,b,c,d) = eri(a,b,c,d)*x(a)*x(b)*y(c)*y(d)/4.0
                          h222(a,b,c,d) = eri(a,d,b,c)*x(a)*x(c)*y(b)*y(d)
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp single
          
          R0 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:R0)
          
          do a=1, na
              do b=1, na
                  R0 = R0 + ( &
                      h11(a, b) * t1(a, b)&
                  )
              end do
          end do
          
          !$omp end do

          !$omp single
          
          
          R0 = R0 + ( &
              h0&
          )
          
          
          !$omp end single

          !$omp do schedule(static) reduction(+:R0)
          
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          R0 = R0 + ( &
                              h221(a, b, c, d) * t2(a, b, c, d)&
                          )
                      end do
                  end do
              end do
          end do
          
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  R1(p, q) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          R1(p, q) = R1(p, q) + ( &
                              h11(a, b) * t2(a, p, b, q)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          R1(p, q) = R1(p, q) + ( &
                              t1(a, b) * h222(p, q, a, b)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
          
                  R1(p, q) = R1(p, q) + ( &
                      h0 * t1(p, q)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
          
                  R1(p, q) = R1(p, q) + ( &
                      h11(p, q)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              R1(p, q) = R1(p, q) + ( &
                                  h13(a, b, c, q) * t2(a, p, b, c)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              R1(p, q) = R1(p, q) - ( &
                                  h31(a, b, c, p) * t2(a, b, c, q)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      R1(p, q) = R1(p, q) + ( &
                          h02(a, q) * t1(p, a)&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      R1(p, q) = R1(p, q) + ( &
                          h20(a, p) * t1(a, q)&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          R2(p, q, r, s) = 0.0
                      end do
                  end do
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              do b=1, na
                                  R2(p, q, r, s) = R2(p, q, r, s) + ( &
                                      h222(p, r, a, b) * t2(a, q, b, s)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              do b=1, na
                                  R2(p, q, r, s) = R2(p, q, r, s) + ( &
                                      h222(q, s, a, b) * t2(a, p, b, r)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              do b=1, na
                                  R2(p, q, r, s) = R2(p, q, r, s) - ( &
                                      h222(p, s, a, b) * t2(a, q, b, r)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              do b=1, na
                                  R2(p, q, r, s) = R2(p, q, r, s) - ( &
                                      h222(q, r, a, b) * t2(a, p, b, s)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              do b=1, na
                                  R2(p, q, r, s) = R2(p, q, r, s) + ( &
                                      2 * h04(a, b, r, s) * t2(p, q, a, b)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              do b=1, na
                                  R2(p, q, r, s) = R2(p, q, r, s) + ( &
                                      2 * h40(a, b, p, q) * t2(a, b, r, s)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
          
                          R2(p, q, r, s) = R2(p, q, r, s) + ( &
                              4 * h221(p, q, r, s)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
          
                          R2(p, q, r, s) = R2(p, q, r, s) + ( &
                              h0 * t2(p, q, r, s)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
          
                          R2(p, q, r, s) = R2(p, q, r, s) + ( &
                              h11(p, r) * t1(q, s)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
          
                          R2(p, q, r, s) = R2(p, q, r, s) + ( &
                              h11(q, s) * t1(p, r)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
          
                          R2(p, q, r, s) = R2(p, q, r, s) - ( &
                              h11(p, s) * t1(q, r)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
          
                          R2(p, q, r, s) = R2(p, q, r, s) - ( &
                              h11(q, r) * t1(p, s)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              R2(p, q, r, s) = R2(p, q, r, s) + ( &
                                  h02(a, r) * t2(p, q, a, s)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              R2(p, q, r, s) = R2(p, q, r, s) + ( &
                                  h20(a, p) * t2(a, q, r, s)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              R2(p, q, r, s) = R2(p, q, r, s) - ( &
                                  h02(a, s) * t2(p, q, a, r)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              R2(p, q, r, s) = R2(p, q, r, s) - ( &
                                  h20(a, q) * t2(a, p, r, s)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              R2(p, q, r, s) = R2(p, q, r, s) - ( &
                                  2 * t1(a, s) * h31(p, q, r, a)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              R2(p, q, r, s) = R2(p, q, r, s) - ( &
                                  2 * t1(p, a) * h13(q, r, s, a)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              R2(p, q, r, s) = R2(p, q, r, s) + ( &
                                  2 * t1(a, r) * h31(p, q, s, a)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do s=1, na
                          do a=1, na
                              R2(p, q, r, s) = R2(p, q, r, s) + ( &
                                  2 * t1(q, a) * h13(p, r, s, a)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp end parallel



          Ovlp = 1.0_pr + Sum(T1*T1) + Sum(T2*T2)/4.0_pr
          Energy = (R0 + Sum(T1*R1) + Sum(T2*R2)/4.0_pr)/Ovlp

      End Subroutine EvalEnergy


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

           !C = C + MatMul(A,B)
           Call DGEMM('N','N',MA,NB,NA,1.0_pr,A,MA,B,NA,1.0_pr,C,MA)

           Return

      End Subroutine A_dot_B

      Subroutine A_dot_Btran(A,B,C,MA,NA,MB)
          ! Subroutine performs the following Matrix Multiplication
          !
          !      C = C + A . B^transpose
          ! 
          !  where A is of dimension (MA,NA)
          !    and B is of dimension (MB,NB)
          ! Obviously, NB = NA
           Implicit None
           Integer, parameter  :: pr = Selected_Real_Kind(15,307)
           Integer, Intent(In)           ::  MA,NA,MB
           Real (Kind=pr), Intent(In)    ::  A(MA,NA), B(MB,NA)
           Real (Kind=pr), Intent(InOut) ::  C(MA,MB)

           !C = C + MatMul(A,Transpose(B))
           Call DGEMM('N','T',MA,MB,NA,1.0_pr,A,MA,B,MB,1.0_pr,C,MA)

           Return

      End Subroutine A_dot_Btran 

      Subroutine Atran_dot_B(A,B,C,NA,MA,NB)
          ! Subroutine performs the following Matrix Multiplication
          !
          !      C = C + A^transpose . B
          ! 
          !  where A is of dimension (NA,MA)
          !    and B is of dimension (MB,NB)
          ! Obviously, MB = NA
           Implicit None
           Integer, parameter  :: pr = Selected_Real_Kind(15,307)
           Integer, Intent(In)           ::  MA,NA,NB
           Real (Kind=pr), Intent(In)    ::  A(NA,MA), B(NA,NB)
           Real (Kind=pr), Intent(InOut) ::  C(MA,NB)

           !C = C + MatMul(Transpose(A),B)
           Call DGEMM('T','N',MA,NB,NA,1.0_pr,A,NA,B,NA,1.0_pr,C,MA)

           Return

      End Subroutine Atran_dot_B

      Subroutine Atran_dot_Btran(A,B,C,NA,MA,MB)
          ! Subroutine performs the following Matrix Multiplication
          !
          !      C = C + A^transpose . B
          ! 
          !  where A is of dimension (NA,MA)
          !    and B is of dimension (MB,NB)
          ! Obviously, NB = NA
           Implicit None
           Integer, parameter  :: pr = Selected_Real_Kind(15,307)
           Integer, Intent(In)           ::  MA,NA,MB
           Real (Kind=pr), Intent(In)    ::  A(NA,MA), B(MB,NA)
           Real (Kind=pr), Intent(InOut) ::  C(MA,MB)

           !C = C + MatMul(Transpose(A),Transpose(B))
           Call DGEMM('T','T',MA,MB,NA,1.0_pr,A,NA,B,MB,1.0_pr,C,MA)

           Return

      End Subroutine Atran_dot_Btran
