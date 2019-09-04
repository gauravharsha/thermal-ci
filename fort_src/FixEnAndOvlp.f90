      Subroutine FixEnAndOvlp(E0,ERI,T0,T1,T2,NA,X,Y,E,Ov)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: T0
          Real (Kind=8), Intent(In) :: T1(NA,NA), T2(NA,NA,NA,NA)
          Real (Kind=8), Intent(In) :: E0(NA), ERI(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: E, Ov
          Integer :: a, b, c, d

          Real (Kind=8), dimension(:), allocatable :: tau0
          Real (Kind=8), dimension(:), allocatable :: tau1
          Real (Kind=8), dimension(:, :), allocatable :: tau2
          Real (Kind=8), dimension(:), allocatable :: tau3

          !$omp parallel default(shared)

          allocate(tau0(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau0(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau0(a) = tau0(a) + ( &
                      y(b)**2 * eri(b, a, b, a)&
                  )
              end do
          end do
          !$omp end do

          !$omp single
          
          e = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:e)
          
          do a=1, na
              e = e + ( &
                  t0 * y(a)**2 * tau0(a) / 2&
              )
          end do
          
          !$omp end do

          deallocate(tau0)

          allocate(tau1(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau1(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau1(a) = tau1(a) + ( &
                          x(b) * y(c) * t1(b, c) * eri(a, c, a, b)&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static) reduction(+:e)
          
          do a=1, na
              e = e + ( &
                  y(a)**2 * tau1(a)&
              )
          end do
          
          !$omp end do

          deallocate(tau1)

          allocate(tau2(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau2(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau2(a, b) = tau2(a, b) - ( &
                              x(c) * y(d) * t2(a, c, d, b) * eri(b, d, c, a)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau3(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau3(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau3(a) = tau3(a) - ( &
                      x(b) * tau2(b, a)&
                  )
              end do
          end do
          !$omp end do

          deallocate(tau2)

          !$omp do schedule(static)
          do a=1, na
          
              tau3(a) = tau3(a) + ( &
                  4 * e0(a) * x(a) * t1(a, a)&
              )
          
          end do
          !$omp end do

          !$omp do schedule(static) reduction(+:e)
          
          do a=1, na
              e = e + ( &
                  tau3(a) * y(a) / 4&
              )
          end do
          
          !$omp end do

          deallocate(tau3)

          !$omp do schedule(static) reduction(+:e)
          
          do a=1, na
              e = e + ( &
                  t0 * y(a)**2 * e0(a)&
              )
          end do
          
          !$omp end do

          !$omp single
          
          ov = 0.0
          
          !$omp end single
          

          !$omp single
          
          
          ov = ov + t0
          
          
          !$omp end single

          !$omp end parallel

      End Subroutine FixEnAndOvlp
