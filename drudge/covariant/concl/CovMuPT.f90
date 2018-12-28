      Subroutine CovMuPT(T0,T1,T2,NA,X,Y,RT0,RT1,RT2)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: T0, T1(NA,NA), T2(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: RT0
          Real (Kind=8), Intent(Out) :: RT1(NA,NA)
          Real (Kind=8), Intent(Out) :: RT2(NA,NA,NA,NA)
          Integer :: a, b, c, d

          Real (Kind=8) :: tau0
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau1

          !$omp parallel default(shared)

          !$omp single
          
          tau0 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:tau0)
          
          do a=1, na
              tau0 = tau0 + ( &
                  y(a)**2&
              )
          end do
          
          !$omp end do

          !$omp single
          
          rt0 = 0.0
          
          !$omp end single
          

          !$omp single
          
          
          rt0 = rt0 + ( &
              t0*tau0 / 2&
          )
          
          
          !$omp end single

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  rt1(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rt1(a, b) = rt1(a, b) + ( &
                      tau0 * t1(a, b) / 2&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          rt2(a, b, c, d) = 0.0
                      end do
                  end do
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) + ( &
                              tau0 * t2(a, b, c, d) / 2&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau1(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau1(a, b, c, d) = 0.0
                      end do
                  end do
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do c=1, na
                  do d=1, na
          
                      tau1(a, a, c, d) = tau1(a, a, c, d) + t1(c, d)
          
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
          
                      tau1(a, b, c, a) = tau1(a, b, c, a) - t1(c, b)
          
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) - ( &
                              x(a) * y(a) * tau1(a, d, b, c) / 2&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) - ( &
                              x(b) * y(b) * tau1(b, c, a, d) / 2&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau1)

          !$omp do schedule(static) reduction(+:rt0)
          
          do a=1, na
              rt0 = rt0 + ( &
                  x(a) * y(a) * t1(a, a) / 2&
              )
          end do
          
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      rt1(a, b) = rt1(a, b) - ( &
                          x(c) * y(c) * t2(c, a, b, c) / 2&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rt1(a, b) = rt1(a, b) + ( &
                      x(a)**2 * t1(a, b) / 2&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rt1(a, b) = rt1(a, b) - ( &
                      y(b)**2 * t1(a, b) / 2&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
          
              rt1(a, a) = rt1(a, a) + ( t0 * x(a) * y(a) / 2 )
          
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) + ( &
                              x(a)**2 * t2(a, b, c, d) / 2&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) + ( &
                              x(b)**2 * t2(a, b, c, d) / 2&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) - ( &
                              y(c)**2 * t2(a, b, c, d) / 2&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) - ( &
                              y(d)**2 * t2(a, b, c, d) / 2&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp end parallel

      End Subroutine CovMuPT
