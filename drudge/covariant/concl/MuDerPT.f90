      Subroutine MuDerPT(Beta,T1,T2,NA,X,Y,TM0,TM1)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: Beta, X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: T1(NA,NA), T2(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: TM0
          Real (Kind=8), Intent(Out) :: TM1(NA,NA)
          Integer :: a, b, c

          !$omp parallel default(shared)

          !$omp single
          
          tm0 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:tm0)
          
          do a=1, na
              tm0 = tm0 - ( &
                  beta * x(a) * y(a) * t1(a, a) / 2&
              )
          end do
          
          !$omp end do

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

          !$omp end parallel

      End Subroutine MuDerPT
