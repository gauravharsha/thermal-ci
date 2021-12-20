      Subroutine FixBetaCISD(E0,ERI,T0,T1,T2,NA,X,Y,RT0,RT1,RT2)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: T0
          Real (Kind=8), Intent(In) :: T1(NA,NA), T2(NA,NA,NA,NA)
          Real (Kind=8), Intent(In) :: E0(NA), ERI(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: RT0
          Real (Kind=8), Intent(Out) :: RT1(NA,NA)
          Real (Kind=8), Intent(Out) :: RT2(NA,NA,NA,NA)
          Integer :: a, b, c, d, i, j

          Real (Kind=8), dimension(:), allocatable :: tau0
          Real (Kind=8), dimension(:), allocatable :: tau1
          Real (Kind=8), dimension(:, :), allocatable :: tau2
          Real (Kind=8), dimension(:), allocatable :: tau3
          Real (Kind=8), dimension(:, :), allocatable :: tau4
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau5
          Real (Kind=8), dimension(:, :, :), allocatable :: tau6
          Real (Kind=8), dimension(:, :), allocatable :: tau7
          Real (Kind=8), dimension(:), allocatable :: tau8
          Real (Kind=8) :: tau9
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau10
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau11
          Real (Kind=8), dimension(:, :, :, :, :), allocatable :: tau12
          Real (Kind=8), dimension(:, :, :, :, :), allocatable :: tau13
          Real (Kind=8), dimension(:, :, :, :, :, :), allocatable :: tau14
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau15
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau16
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau17
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau18
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau19

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

          allocate(tau1(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau1(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
          
              tau1(a) = tau1(a) + ( &
                  t0 * tau0(a) / 2&
              )
          
          end do
          !$omp end do

          allocate(tau8(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau8(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
          
              tau8(a) = tau8(a) + ( &
                  tau0(a)&
              )
          
          end do
          !$omp end do

          deallocate(tau0)

          !$omp do schedule(static)
          do a=1, na
          
              tau1(a) = tau1(a) + ( &
                  t0 * e0(a)&
              )
          
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

          !$omp single
          
          rt0 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:rt0)
          
          do a=1, na
              rt0 = rt0 + ( &
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

          !$omp do schedule(static) reduction(+:rt0)
          
          do a=1, na
              rt0 = rt0 + ( &
                  tau3(a) * y(a) / 4&
              )
          end do
          
          !$omp end do

          deallocate(tau3)

          allocate(tau4(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau4(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau4(a, b) = tau4(a, b) + ( &
                          y(c)**2 * eri(c, a, c, b)&
                      )
                  end do
              end do
          end do
          !$omp end do

          allocate(tau6(1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau6(a, b, c) = 0.0
                  end do
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
          
                      tau6(a, b, c) = tau6(a, b, c) + ( &
                          2 * t1(a, b) * tau4(b, c)&
                      )
          
                  end do
              end do
          end do
          !$omp end do

          allocate(tau7(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau7(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  tau7(a, b) = tau7(a, b) + ( &
                      t0 * y(b) * tau4(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau7(a, b) = tau7(a, b) + ( &
                          x(c) * t1(c, b) * tau4(a, c)&
                      )
                  end do
              end do
          end do
          !$omp end do

          allocate(tau10(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau10(a, b, c, d) = 0.0
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
                              tau10(a, b, c, d) = tau10(a, b, c, d) + ( &
                                  x(d) * x(i) * tau4(d, i) * t2(i, a, b, c)&
                              )
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
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) - ( &
                              tau10(a, c, d, b)&
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
                              tau10(b, c, d, a)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau10)

          allocate(tau11(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau11(a, b, c, d) = 0.0
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
                              tau11(a, b, c, d) = tau11(a, b, c, d) + ( &
                                  y(d) * y(i) * tau4(i, d) * t2(a, b, c, i)&
                              )
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
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) - ( &
                              tau11(a, b, c, d)&
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
                              tau11(a, b, d, c)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau11)

          allocate(tau16(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau16(a, b, c, d) = 0.0
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
          
                          tau16(a, b, c, d) = tau16(a, b, c, d) + ( &
                              t1(a, c) * tau4(b, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau18(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau18(a, b, c, d) = 0.0
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
          
                          tau18(a, b, c, d) = tau18(a, b, c, d) - ( &
                              y(d) * t1(c, b) * tau4(a, d)&
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
          
                          tau18(a, b, c, d) = tau18(a, b, c, d) + ( &
                              y(b) * t1(c, d) * tau4(a, b)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

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
                  do c=1, na
                      do d=1, na
                          rt1(a, b) = rt1(a, b) - ( &
                              x(c) * y(d) * tau4(d, c) * t2(c, a, b, d)&
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
                              x(a) * y(c) * t1(b, d) * tau4(a, c)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau4)

          allocate(tau5(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau5(a, b, c, d) = 0.0
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
          
                          tau5(a, b, c, d) = tau5(a, b, c, d) + ( &
                              2 * y(d) * t1(a, b) * eri(c, b, d, a)&
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
                          do i=1, na
                              tau5(a, b, c, d) = tau5(a, b, c, d) - ( &
                                  x(i) * t2(i, a, d, b) * eri(c, b, a, i)&
                              )
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
                          rt1(a, b) = rt1(a, b) + ( &
                              x(a) * x(c) * y(d) * tau5(c, d, a, b) / 2&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau5)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              tau6(a, b, c) = tau6(a, b, c) + ( &
                                  x(d) * y(i) * t2(d, a, i, b) * eri(b, i, c, d)&
                              )
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
                      rt1(a, b) = rt1(a, b) - ( &
                          y(b) * y(c) * tau6(a, c, b) / 2&
                      )
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau6)

          !$omp do schedule(static)
          do a=1, na
          
              tau7(a, a) = tau7(a, a) + ( t0 * e0(a) * y(a) )
          
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rt1(a, b) = rt1(a, b) + ( &
                      x(a) * tau7(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          deallocate(tau7)

          !$omp do schedule(static)
          do a=1, na
          
              tau8(a) = tau8(a) + ( &
                  2 * e0(a)&
              )
          
          end do
          !$omp end do

          !$omp single
          
          tau9 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:tau9)
          
          do a=1, na
              tau9 = tau9 + ( &
                  y(a)**2 * tau8(a)&
              )
          end do
          
          !$omp end do

          deallocate(tau8)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rt1(a, b) = rt1(a, b) + ( &
                      tau9 * t1(a, b) / 2&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) + ( &
                              tau9 * t2(a, b, c, d) / 2&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau12(1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              tau12(a, b, c, d, i) = 0.0
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
                                  tau12(a, b, c, d, i) = tau12(a, b, c, d, i) - ( &
                                      x(j) * t2(j, a, b, c) * eri(i, d, a, j)&
                                  )
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
          
                              tau12(a, b, c, d, i) = tau12(a, b, c, d, i) - ( &
                                  2 * y(c) * t1(a, b) * eri(i, d, c, a)&
                              )
          
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
          
                              tau12(a, b, c, d, i) = tau12(a, b, c, d, i) + ( &
                                  2 * y(b) * t1(a, c) * eri(i, d, b, a)&
                              )
          
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
                              rt2(a, b, c, d) = rt2(a, b, c, d) - ( &
                                  x(a) * x(b) * x(i) * tau12(i, c, d, a, b) / 2&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau12)

          allocate(tau13(1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              tau13(a, b, c, d, i) = 0.0
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
                                  tau13(a, b, c, d, i) = tau13(a, b, c, d, i) - ( &
                                      y(j) * t2(a, b, j, c) * eri(c, j, i, d)&
                                  )
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
          
                              tau13(a, b, c, d, i) = tau13(a, b, c, d, i) + ( &
                                  2 * x(b) * t1(a, c) * eri(b, c, i, d)&
                              )
          
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
          
                              tau13(a, b, c, d, i) = tau13(a, b, c, d, i) - ( &
                                  2 * x(a) * t1(b, c) * eri(a, c, i, d)&
                              )
          
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
                              rt2(a, b, c, d) = rt2(a, b, c, d) + ( &
                                  y(c) * y(d) * y(i) * tau13(a, b, i, d, c) / 2&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau13)

          allocate(tau14(1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  tau14(a, b, c, d, i, j) = 0.0
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
          
                                  tau14(a, b, c, d, i, j) = tau14(a, b, c, d, i, j) + ( &
                                      y(j) * t2(a, b, c, d) * eri(i, d, j, a)&
                                  )
          
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
          
                                  tau14(a, b, c, d, i, j) = tau14(a, b, c, d, i, j) - ( &
                                      y(c) * t2(a, b, j, d) * eri(i, d, c, a)&
                                  )
          
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau15(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau15(a, b, c, d) = 0.0
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
                                  tau15(a, b, c, d) = tau15(a, b, c, d) + ( &
                                      x(c) * x(i) * y(j) * tau14(i, a, b, j, c, d)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau14)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) - ( &
                              tau15(b, d, a, c)&
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
                              tau15(a, c, b, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau15)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          tau16(a, b, c, d) = tau16(a, b, c, d) + ( &
                              t0 * x(a) * y(c) * eri(a, b, c, d)&
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
                              x(a) * y(d) * tau16(b, a, c, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau16)

          allocate(tau17(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau17(a, b, c, d) = 0.0
                      end do
                  end do
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do c=1, na
                  do d=1, na
          
                      tau17(a, a, c, d) = tau17(a, a, c, d) + t1(c, d)
          
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
          
                      tau17(a, b, c, a) = tau17(a, b, c, a) - t1(c, b)
          
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          tau18(a, b, c, d) = tau18(a, b, c, d) - ( &
                              e0(a) * y(a) * tau17(a, d, c, b)&
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
                              x(b) * tau18(b, c, a, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau18)

          allocate(tau19(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau19(a, b, c, d) = 0.0
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
          
                          tau19(a, b, c, d) = tau19(a, b, c, d) - ( &
                              x(a) * y(a) * tau17(a, d, b, c)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau17)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          tau19(a, b, c, d) = tau19(a, b, c, d) + ( &
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
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) + ( &
                              e0(a) * tau19(a, b, c, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau19)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rt1(a, b) = rt1(a, b) + ( &
                      x(a)**2 * e0(a) * t1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rt1(a, b) = rt1(a, b) - ( &
                      y(b)**2 * e0(b) * t1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      rt1(a, b) = rt1(a, b) - ( &
                          e0(c) * x(c) * y(c) * t2(c, a, b, c)&
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
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) + ( &
                              x(b)**2 * e0(b) * t2(a, b, c, d)&
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
                              y(c)**2 * e0(c) * t2(a, b, c, d)&
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
                              y(d)**2 * e0(d) * t2(a, b, c, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp end parallel

      End Subroutine FixBetaCISD


      Subroutine FixMuCISD(T0,T1,T2,NA,X,Y,RT0,RT1,RT2)
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
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau2

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
              t0*tau0&
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
          
                      tau1(a, a, c, d) = tau1(a, a, c, d) - t1(c, d)
          
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
          
                      tau1(a, b, c, a) = tau1(a, b, c, a) + t1(c, b)
          
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
                              x(a) * y(a) * tau1(a, c, b, d)&
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
          
                      tau2(a, a, c, d) = tau2(a, a, c, d) + t1(c, d)
          
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
          
                      tau2(a, b, c, a) = tau2(a, b, c, a) - t1(c, b)
          
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
                              x(b) * y(b) * tau2(b, c, a, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau2)

          !$omp do schedule(static) reduction(+:rt0)
          
          do a=1, na
              rt0 = rt0 + ( &
                  x(a) * y(a) * t1(a, a)&
              )
          end do
          
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      rt1(a, b) = rt1(a, b) - ( &
                          x(c) * y(c) * t2(c, a, b, c)&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rt1(a, b) = rt1(a, b) + ( &
                      x(a)**2 * t1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rt1(a, b) = rt1(a, b) - ( &
                      y(b)**2 * t1(a, b)&
                  )
          
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
                  do c=1, na
                      do d=1, na
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) + ( &
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
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) + ( &
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
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) - ( &
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
          
                          rt2(a, b, c, d) = rt2(a, b, c, d) - ( &
                              y(d)**2 * t2(a, b, c, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp end parallel

      End Subroutine FixMuCISD
