      Subroutine CovMuPT2(T0,T1,T2,NA,X,Y,RT0,RT1,RS0,RS1,RS2,RS3)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: T0
          Real (Kind=8), Intent(In) :: T1(NA,NA)
          Real (Kind=8), Intent(In) :: T2(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: RT0, RS0
          Real (Kind=8), Intent(Out) :: RT1(NA,NA), RS1(NA,NA)
          Real (Kind=8), Intent(Out) :: RS2(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: RS3(NA,NA,NA,NA,NA,NA)
          Integer :: a, b, c, d, i, j

          Real (Kind=8) :: tau0
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau1
          Real (Kind=8), dimension(:, :, :, :, :, :), allocatable :: tau2
          Real (Kind=8), dimension(:, :, :, :, :, :), allocatable :: tau3
          Real (Kind=8), dimension(:, :, :, :, :, :), allocatable :: tau4

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
              tau0 / 2&
          )
          
          
          !$omp end single

          !$omp single
          
          rs0 = 0.0
          
          !$omp end single
          

          !$omp single
          
          
          rs0 = rs0 + ( &
              t0*tau0 / 2&
          )
          
          
          !$omp end single

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
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) - ( &
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
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) - ( &
                              x(b) * y(b) * tau1(b, c, a, d) / 2&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau1)

          allocate(tau2(1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  tau2(a, b, c, d, i, j) = 0.0
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do c=1, na
                  do d=1, na
                      do i=1, na
                          do j=1, na
          
                              tau2(a, a, c, d, i, j) = tau2(a, a, c, d, i, j) + t2(c, d, i, j)
          
                          end do
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
                          do j=1, na
          
                              tau2(a, b, c, d, a, j) = tau2(a, b, c, d, a, j) + t2(c, d, j, b)
          
                          end do
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
                          do i=1, na
          
                              tau2(a, b, c, d, i, a) = tau2(a, b, c, d, i, a) - t2(c, d, i, b)
          
                          end do
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
                          do i=1, na
                              do j=1, na
                                  rs3(a, b, c, d, i, j) = 0.0
                              end do
                          end do
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
                          do i=1, na
                              do j=1, na
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      x(a) * y(a) * tau2(a, j, b, c, d, i) / 2&
                                  )
          
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau2)

          allocate(tau3(1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  tau3(a, b, c, d, i, j) = 0.0
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do c=1, na
                  do d=1, na
                      do i=1, na
                          do j=1, na
          
                              tau3(a, a, c, d, i, j) = tau3(a, a, c, d, i, j) - t2(c, d, i, j)
          
                          end do
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
                          do j=1, na
          
                              tau3(a, b, c, d, a, j) = tau3(a, b, c, d, a, j) + t2(c, d, b, j)
          
                          end do
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
                          do i=1, na
          
                              tau3(a, b, c, d, i, a) = tau3(a, b, c, d, i, a) + t2(c, d, i, b)
          
                          end do
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
                          do i=1, na
                              do j=1, na
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) - ( &
                                      x(b) * y(b) * tau3(b, i, a, c, d, j) / 2&
                                  )
          
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau3)

          allocate(tau4(1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  tau4(a, b, c, d, i, j) = 0.0
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do c=1, na
                  do d=1, na
                      do i=1, na
                          do j=1, na
          
                              tau4(a, a, c, d, i, j) = tau4(a, a, c, d, i, j) + t2(c, d, i, j)
          
                          end do
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
                          do j=1, na
          
                              tau4(a, b, c, d, a, j) = tau4(a, b, c, d, a, j) - t2(c, d, b, j)
          
                          end do
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
                          do i=1, na
          
                              tau4(a, b, c, d, i, a) = tau4(a, b, c, d, i, a) + t2(c, d, b, i)
          
                          end do
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
                          do i=1, na
                              do j=1, na
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      x(c) * y(c) * tau4(c, d, a, b, i, j) / 2&
                                  )
          
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau4)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  rt1(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
          
              rt1(a, a) = rt1(a, a) + ( x(a) * y(a) / 2 )
          
          end do
          !$omp end do

          !$omp do schedule(static) reduction(+:rs0)
          
          do a=1, na
              rs0 = rs0 + ( &
                  x(a) * y(a) * t1(a, a) / 2&
              )
          end do
          
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rs1(a, b) = rs1(a, b) + ( &
                      x(a)**2 * t1(a, b) / 2&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rs1(a, b) = rs1(a, b) - ( &
                      y(b)**2 * t1(a, b) / 2&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      rs1(a, b) = rs1(a, b) - ( &
                          x(c) * y(c) * t2(c, a, b, c) / 2&
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
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) + ( &
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
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) - ( &
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
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) - ( &
                              y(d)**2 * t2(a, b, c, d) / 2&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp end parallel

      End Subroutine CovMuPT2
