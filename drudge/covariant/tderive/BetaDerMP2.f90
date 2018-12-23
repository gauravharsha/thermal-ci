      Subroutine BetaDerMP2(MU,E0,T2,S2,S4,NA,x,y,TB1,SB1,SB3)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: MU
          Real (Kind=8), Intent(In) :: x(NA), y(NA), E0(NA)
          Real (Kind=8), Intent(In) :: T2(NA,NA,NA,NA), S2(NA,NA,NA,NA)
          Real (Kind=8), Intent(In) :: S4(NA,NA,NA,NA,NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: TB1(NA,NA), SB1(NA,NA)
          Real (Kind=8), Intent(Out) :: SB3(NA,NA,NA,NA,NA,NA)
          Integer :: p,q,r,s, a, b, c, d

          Real (Kind=8), dimension(:), allocatable :: tau0

          !$omp parallel default(shared)

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  tb1(p, q) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      tb1(p, q) = tb1(p, q) + ( &
                          e0(a) * x(a) * y(a) * t2(a, p, a, q) / 2&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      tb1(p, q) = tb1(p, q) - ( &
                          mu * x(a) * y(a) * t2(a, p, a, q) / 2&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  sb1(p, q) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      sb1(p, q) = sb1(p, q) + ( &
                          e0(a) * x(a) * y(a) * s2(a, p, a, q) / 2&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      sb1(p, q) = sb1(p, q) - ( &
                          mu * x(a) * y(a) * s2(a, p, a, q) / 2&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do r=1, na
                      do a=1, na
                          do b=1, na
                              do c=1, na
                                  sb3(p, q, r, a, b, c) = 0.0
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
                      do a=1, na
                          do b=1, na
                              do c=1, na
                                  do d=1, na
                                      sb3(p, q, r, a, b, c) = sb3(p, q, r, a, b, c) + ( &
                                          mu * x(d) * y(d) * s4(d, p, q, r, a, b, c, d)&
                                      )
                                  end do
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
                      do a=1, na
                          do b=1, na
                              do c=1, na
                                  do d=1, na
                                      sb3(p, q, r, a, b, c) = sb3(p, q, r, a, b, c) - ( &
                                          e0(d) * x(d) * y(d) * s4(d, p, q, r, a, b, c, d)&
                                      )
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp end parallel

      End Subroutine BetaDerMP2
