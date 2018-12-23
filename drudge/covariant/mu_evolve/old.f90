      Subroutine CovMuPT2(T1,T2,NA,X,Y,RT1,RS1,RS2)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: T1(NA,NA)
          Real (Kind=8), Intent(In) :: T2(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: RT1(NA,NA), RS1(NA,NA)
          Real (Kind=8), Intent(Out) :: RS2(NA,NA,NA,NA)
          Integer :: a, b, c, d

          Real (Kind=8) :: tau0
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau1
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau2

          !=================================================================!
          !     NOTE: This code is for Mu Evolution such that the Mu        !
          !     dependence is totally ignored in the reference state        !
          !=================================================================!


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

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  rs1(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rs1(a, b) = rs1(a, b) + ( &
                      tau0 * t1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          rs2(a, b, c, d) = 0.0
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
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) + ( &
                              tau0 * t2(a, b, c, d)&
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
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) - ( &
                              x(a) * y(a) * tau1(a, d, b, c)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau1)

          allocate(tau2(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau2(a, b, c, d) = 0.0
                      end do
                  end do
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do c=1, na
                  do d=1, na
          
                      tau2(a, a, c, d) = tau2(a, a, c, d) - t1(c, d)
          
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
          
                      tau2(a, b, c, a) = tau2(a, b, c, a) + t1(c, b)
          
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) - ( &
                              x(b) * y(b) * tau2(b, d, a, c)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau2)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  rt1(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
          
              rt1(a, a) = rt1(a, a) + ( x(a) * y(a) )
          
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rs1(a, b) = rs1(a, b) + ( &
                      x(a)**2 * t1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rs1(a, b) = rs1(a, b) - ( &
                      y(b)**2 * t1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      rs1(a, b) = rs1(a, b) - ( &
                          x(c) * y(c) * t2(c, a, b, c)&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) + ( &
                              x(a)**2 * t2(a, b, c, d)&
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
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) + ( &
                              x(b)**2 * t2(a, b, c, d)&
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
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) - ( &
                              y(c)**2 * t2(a, b, c, d)&
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
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) - ( &
                              y(d)**2 * t2(a, b, c, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp end parallel

      End Subroutine CovMuPT2
