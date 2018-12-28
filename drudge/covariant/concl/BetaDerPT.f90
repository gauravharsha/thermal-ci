      Subroutine BetaDerPT(E0,Mu,T1,T2,NA,X,Y,TB0,TB1)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: E0(NA), Mu
          Real (Kind=8), Intent(In) :: T1(NA,NA), T2(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: TB0
          Real (Kind=8), Intent(Out) :: TB1(NA,NA)
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

          !$omp single
          
          tb0 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:tb0)
          
          do a=1, na
              tb0 = tb0 - ( &
                  tau0(a) * x(a) * t1(a, a) / 2&
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

          deallocate(tau0)

          !$omp end parallel

      End Subroutine BetaDerPT
