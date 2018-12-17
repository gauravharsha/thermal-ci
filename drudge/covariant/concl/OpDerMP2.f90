      Subroutine OpDerMP2(BETA,MU,E0,T2,S2,NA,x,y,TB1,SB1,TM1,SM1)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: MU, BETA
          Real (Kind=8), Intent(In) :: x(NA), y(NA), E0(NA)
          Real (Kind=8), Intent(In) :: T2(NA,NA,NA,NA), S2(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: TB1(NA,NA), SB1(NA,NA)
          Real (Kind=8), Intent(Out) :: TM1(NA,NA), SM1(NA,NA)
          Integer :: a, b, c

          Real (Kind=8), dimension(:), allocatable :: tau0

          !$omp parallel default(shared)

          allocate(tau0(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau0(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
          
              tau0(a) = tau0(a) + ( &
                  mu * y(a)&
              )
          
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
          
              tau0(a) = tau0(a) - ( &
                  e0(a) * y(a)&
              )
          
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tb1(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tb1(a, b) = tb1(a, b) - ( &
                          tau0(c) * x(c) * t2(c, a, c, b) / 2&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  sb1(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      sb1(a, b) = sb1(a, b) - ( &
                          tau0(c) * x(c) * s2(c, a, c, b) / 2&
                      )
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau0)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tm1(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tm1(a, b) = tm1(a, b) - ( &
                          beta * x(c) * y(c) * t2(c, a, c, b) / 2&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  sm1(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      sm1(a, b) = sm1(a, b) - ( &
                          beta * x(c) * y(c) * s2(c, a, c, b) / 2&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp end parallel

      End Subroutine OpDerMP2
