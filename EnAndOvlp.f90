      Subroutine EnAndOvlp(E0,ERI,T0,T1,T2,NA,X,Y,E,Ov)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: T0
          Real (Kind=8), Intent(In) :: T1(NA,NA), T2(NA,NA,NA,NA)
          Real (Kind=8), Intent(In) :: E0(NA), ERI(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: E, Ov
          Integer :: a, b, c, d, i, j

          Real (Kind=8), dimension(:, :), allocatable :: tau0
          Real (Kind=8), dimension(:, :), allocatable :: tau1
          Real (Kind=8), dimension(:, :), allocatable :: tau2
          Real (Kind=8), dimension(:, :), allocatable :: tau3
          Real (Kind=8), dimension(:, :), allocatable :: tau4
          Real (Kind=8) :: tau5
          Real (Kind=8) :: tau6
          Real (Kind=8) :: tau7
          Real (Kind=8), dimension(:), allocatable :: tau8
          Real (Kind=8), dimension(:), allocatable :: tau9
          Real (Kind=8), dimension(:), allocatable :: tau10
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau11
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau12
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau13
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau14
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau15
          Real (Kind=8), dimension(:), allocatable :: tau16
          Real (Kind=8), dimension(:), allocatable :: tau17
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau18
          Real (Kind=8), dimension(:, :), allocatable :: tau19
          Real (Kind=8), dimension(:), allocatable :: tau20
          Real (Kind=8), dimension(:), allocatable :: tau21

          !$omp parallel default(shared)

          allocate(tau0(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau0(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau0(a, b) = tau0(a, b) + ( &
                              t1(c, d) * t2(c, a, d, b)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau1(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau1(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  tau1(a, b) = tau1(a, b) + ( &
                      tau0(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          deallocate(tau0)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  tau1(a, b) = tau1(a, b) + ( &
                      t0 * t1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

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
          
                  tau4(a, b) = tau4(a, b) + ( &
                      x(b) * y(a) * tau1(b, a)&
                  )
          
              end do
          end do
          !$omp end do

          allocate(tau10(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau10(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau10(a) = tau10(a) + ( &
                          x(c) * y(b) * tau1(c, b) * eri(a, b, a, c)&
                      )
                  end do
              end do
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
                      tau2(a, b) = tau2(a, b) + ( &
                          2 * t1(c, a) * t1(c, b)&
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
                          do i=1, na
                              tau2(a, b) = tau2(a, b) + ( &
                                  t2(i, d, c, a) * t2(i, d, c, b)&
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
          
                  tau4(a, b) = tau4(a, b) - ( &
                      y(a) * y(b) * tau2(b, a) / 2&
                  )
          
              end do
          end do
          !$omp end do

          deallocate(tau2)

          allocate(tau3(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau3(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau3(a, b) = tau3(a, b) + ( &
                          2 * t1(a, c) * t1(b, c)&
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
                          do i=1, na
                              tau3(a, b) = tau3(a, b) + ( &
                                  t2(i, a, d, c) * t2(i, b, d, c)&
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
          
                  tau4(a, b) = tau4(a, b) + ( &
                      x(a) * x(b) * tau3(b, a) / 2&
                  )
          
              end do
          end do
          !$omp end do

          deallocate(tau3)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau10(a) = tau10(a) + ( &
                          tau4(c, b) * eri(a, b, a, c)&
                      )
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau4)

          !$omp single
          
          tau5 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:tau5)
          
          do a=1, na
              do b=1, na
                  tau5 = tau5 + ( &
                      t1(a, b)**2&
                  )
              end do
          end do
          
          !$omp end do

          !$omp single
          
          tau7 = 0.0
          
          !$omp end single
          

          !$omp single
          
          
          tau7 = tau7 + ( &
              4*tau5&
          )
          
          
          !$omp end single

          !$omp single
          
          ov = 0.0
          
          !$omp end single
          

          !$omp single
          
          
          ov = ov + ( &
              tau5&
          )
          
          
          !$omp end single

          !$omp single
          
          tau6 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:tau6)
          
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau6 = tau6 + ( &
                              t2(a, b, c, d)**2&
                          )
                      end do
                  end do
              end do
          end do
          
          !$omp end do

          !$omp single
          
          
          tau7 = tau7 + ( &
              tau6&
          )
          
          
          !$omp end single

          allocate(tau8(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau8(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
          
              tau8(a) = tau8(a) + ( &
                  tau7 * y(a)**2 / 4&
              )
          
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
          
              tau10(a) = tau10(a) + ( &
                  tau7 * e0(a) / 4&
              )
          
          end do
          !$omp end do

          !$omp single
          
          
          ov = ov + ( &
              tau6 / 4&
          )
          
          
          !$omp end single

          !$omp do schedule(static)
          do a=1, na
          
              tau8(a) = tau8(a) + ( &
                  t0**2 * y(a)**2&
              )
          
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau10(a) = tau10(a) + ( &
                      tau8(b) * eri(b, a, b, a) / 2&
                  )
              end do
          end do
          !$omp end do

          deallocate(tau8)

          allocate(tau9(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau9(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau9(a) = tau9(a) + ( &
                      2 * t1(b, a)**2&
                  )
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau9(a) = tau9(a) + ( &
                              t2(b, c, a, d)**2&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
          
              tau10(a) = tau10(a) - ( &
                  e0(a) * tau9(a) / 2&
              )
          
          end do
          !$omp end do

          deallocate(tau9)

          !$omp do schedule(static)
          do a=1, na
          
              tau10(a) = tau10(a) + ( &
                  t0**2 * e0(a)&
              )
          
          end do
          !$omp end do

          !$omp single
          
          e = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:e)
          
          do a=1, na
              e = e + ( &
                  y(a)**2 * tau10(a)&
              )
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
          
                          tau11(a, b, c, d) = tau11(a, b, c, d) + ( &
                              eri(b, a, d, c)&
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
          
                          tau11(a, b, c, d) = tau11(a, b, c, d) + ( &
                              eri(d, c, b, a)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau17(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau17(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              tau17(a) = tau17(a) - ( &
                                  x(c) * x(d) * y(i) * t1(a, b) * t2(d, c, b, i) * tau11(c, d, a, i) / 4&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau20(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau20(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              tau20(a) = tau20(a) - ( &
                                  4 * x(c) * y(d) * y(i) * t1(b, a) * t2(b, c, i, d) * tau11(d, i, a, c)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau11)

          allocate(tau12(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau12(a, b, c, d) = 0.0
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
                                  tau12(a, b, c, d) = tau12(a, b, c, d) + ( &
                                      t2(a, c, i, j) * t2(b, d, i, j)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau14(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau14(a, b, c, d) = 0.0
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
          
                          tau14(a, b, c, d) = tau14(a, b, c, d) - ( &
                              x(a) * x(c) * x(d) * tau12(b, a, c, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau12)

          allocate(tau13(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau13(a, b, c, d) = 0.0
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
          
                          tau13(a, b, c, d) = tau13(a, b, c, d) + ( &
                              t1(a, c) * t1(b, d)&
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
                              do j=1, na
                                  tau13(a, b, c, d) = tau13(a, b, c, d) + ( &
                                      t2(j, a, i, c) * t2(j, b, i, d)&
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
          
                          tau14(a, b, c, d) = tau14(a, b, c, d) + ( &
                              8 * x(a) * y(c) * y(d) * tau13(b, a, d, c)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau13)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau17(a) = tau17(a) - ( &
                              tau14(d, a, b, c) * eri(b, a, c, d) / 16&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau14)

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
          
                          tau15(a, b, c, d) = tau15(a, b, c, d) + ( &
                              t2(a, b, c, d) * eri(b, a, d, c)&
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
          
                          tau15(a, b, c, d) = tau15(a, b, c, d) + ( &
                              t2(b, a, c, d) * eri(d, c, a, b)&
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
                          tau17(a) = tau17(a) + ( &
                              t0 * x(b) * y(c) * y(d) * tau15(b, a, d, c) / 8&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau15)

          allocate(tau16(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau16(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
          
              tau16(a) = tau16(a) + ( &
                  t0 * t1(a, a)&
              )
          
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau16(a) = tau16(a) + ( &
                          t1(b, c) * t2(a, b, a, c)&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
          
              tau17(a) = tau17(a) + ( &
                  e0(a) * tau16(a) * y(a)&
              )
          
          end do
          !$omp end do

          deallocate(tau16)

          !$omp do schedule(static) reduction(+:e)
          
          do a=1, na
              e = e + ( &
                  2 * tau17(a) * x(a)&
              )
          end do
          
          !$omp end do

          deallocate(tau17)

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
                          do i=1, na
                              do j=1, na
                                  tau18(a, b, c, d) = tau18(a, b, c, d) + ( &
                                      t2(j, i, a, c) * t2(j, i, b, d)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau19(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau19(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau19(a, b) = tau19(a, b) + ( &
                              y(c) * y(d) * tau18(b, a, d, c) * eri(c, a, d, b)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau18)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau20(a) = tau20(a) - ( &
                      y(b) * tau19(b, a)&
                  )
              end do
          end do
          !$omp end do

          deallocate(tau19)

          !$omp do schedule(static) reduction(+:e)
          
          do a=1, na
              e = e - ( &
                  tau20(a) * y(a) / 8&
              )
          end do
          
          !$omp end do

          deallocate(tau20)

          allocate(tau21(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau21(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau21(a) = tau21(a) + ( &
                      2 * t1(a, b)**2&
                  )
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau21(a) = tau21(a) + ( &
                              t2(a, b, c, d)**2&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static) reduction(+:e)
          
          do a=1, na
              e = e + ( &
                  x(a)**2 * e0(a) * tau21(a) / 2&
              )
          end do
          
          !$omp end do

          deallocate(tau21)

          !$omp single
          
          
          ov = ov + ( &
              t0**2&
          )
          
          
          !$omp end single

          !$omp end parallel

      End Subroutine EnAndOvlp
