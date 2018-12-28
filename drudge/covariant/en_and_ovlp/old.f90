      Subroutine EnAndOvlp(E0,ERI,T0,T1,T2,NA,X,Y,E1,E2,Ov1,Ov2)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: T0
          Real (Kind=8), Intent(In) :: T1(NA,NA), T2(NA,NA,NA,NA)
          Real (Kind=8), Intent(In) :: E0(NA), ERI(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: E1, E2, Ov1, Ov2
          Integer :: a, b, c, d, i, j

          Real (Kind=8), dimension(:), allocatable :: tau0
          Real (Kind=8), dimension(:), allocatable :: tau1
          Real (Kind=8) :: tau2
          Real (Kind=8), dimension(:, :), allocatable :: tau3
          Real (Kind=8), dimension(:, :), allocatable :: tau4
          Real (Kind=8), dimension(:, :), allocatable :: tau5
          Real (Kind=8), dimension(:, :), allocatable :: tau6
          Real (Kind=8), dimension(:), allocatable :: tau7
          Real (Kind=8), dimension(:), allocatable :: tau8
          Real (Kind=8), dimension(:, :), allocatable :: tau9
          Real (Kind=8), dimension(:), allocatable :: tau10
          Real (Kind=8), dimension(:), allocatable :: tau11
          Real (Kind=8), dimension(:), allocatable :: tau12
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau13
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau14
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau15
          Real (Kind=8), dimension(:, :), allocatable :: tau16
          Real (Kind=8), dimension(:), allocatable :: tau17
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau18
          Real (Kind=8), dimension(:, :), allocatable :: tau19
          Real (Kind=8), dimension(:), allocatable :: tau20
          Real (Kind=8) :: tau21
          Real (Kind=8), dimension(:), allocatable :: tau22
          Real (Kind=8) :: tau23
          Real (Kind=8), dimension(:), allocatable :: tau24
          Real (Kind=8), dimension(:, :), allocatable :: tau25
          Real (Kind=8), dimension(:, :), allocatable :: tau26
          Real (Kind=8), dimension(:, :), allocatable :: tau27
          Real (Kind=8), dimension(:, :), allocatable :: tau28
          Real (Kind=8), dimension(:, :), allocatable :: tau29
          Real (Kind=8), dimension(:), allocatable :: tau30

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
                      t1(a, b)**2&
                  )
              end do
          end do
          !$omp end do

          !$omp single
          
          e1 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:e1)
          
          do a=1, na
              e1 = e1 + ( &
                  x(a)**2 * e0(a) * tau0(a)&
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
                  tau1(a) = tau1(a) + ( &
                      t1(b, a)**2&
                  )
              end do
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
          
              tau8(a) = tau8(a) - ( &
                  e0(a) * tau1(a)&
              )
          
          end do
          !$omp end do

          deallocate(tau1)

          !$omp single
          
          tau2 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:tau2)
          
          do a=1, na
              do b=1, na
                  tau2 = tau2 + ( &
                      t1(a, b)**2&
                  )
              end do
          end do
          
          !$omp end do

          allocate(tau7(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau7(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
          
              tau7(a) = tau7(a) + ( &
                  tau2 * y(a)**2&
              )
          
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
          
              tau8(a) = tau8(a) + ( &
                  tau2 * e0(a)&
              )
          
          end do
          !$omp end do

          !$omp single
          
          ov1 = 0.0
          
          !$omp end single
          

          !$omp single
          
          
          ov1 = ov1 + ( &
              tau2&
          )
          
          
          !$omp end single

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
                          t1(c, a) * t1(c, b)&
                      )
                  end do
              end do
          end do
          !$omp end do

          allocate(tau6(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau6(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  tau6(a, b) = tau6(a, b) - ( &
                      y(a) * y(b) * tau3(a, b)&
                  )
          
              end do
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
                          t1(a, c) * t1(b, c)&
                      )
                  end do
              end do
          end do
          !$omp end do

          allocate(tau5(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau5(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  tau5(a, b) = tau5(a, b) + ( &
                      x(b) * tau4(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          deallocate(tau4)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  tau5(a, b) = tau5(a, b) + ( &
                      t0 * y(b) * t1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  tau6(a, b) = tau6(a, b) + ( &
                      x(b) * tau5(b, a)&
                  )
          
              end do
          end do
          !$omp end do

          deallocate(tau5)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau8(a) = tau8(a) + ( &
                          tau6(b, c) * eri(a, c, a, b)&
                      )
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau6)

          !$omp do schedule(static)
          do a=1, na
          
              tau7(a) = tau7(a) + ( &
                  t0**2 * y(a)**2&
              )
          
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau8(a) = tau8(a) + ( &
                      tau7(b) * eri(b, a, b, a) / 2&
                  )
              end do
          end do
          !$omp end do

          deallocate(tau7)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau8(a) = tau8(a) + ( &
                          t0 * x(b) * y(c) * t1(b, c) * eri(a, c, a, b)&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static) reduction(+:e1)
          
          do a=1, na
              e1 = e1 + ( &
                  y(a)**2 * tau8(a)&
              )
          end do
          
          !$omp end do

          deallocate(tau8)

          allocate(tau9(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau9(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau9(a, b) = tau9(a, b) - ( &
                              x(c) * y(d) * t1(c, d) * eri(a, c, d, b)&
                          )
                      end do
                  end do
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
                  tau10(a) = tau10(a) + ( &
                      x(b) * t1(b, a) * tau9(a, b)&
                  )
              end do
          end do
          !$omp end do

          deallocate(tau9)

          !$omp do schedule(static)
          do a=1, na
          
              tau10(a) = tau10(a) + ( &
                  2*t0 * e0(a) * x(a) * t1(a, a)&
              )
          
          end do
          !$omp end do

          !$omp do schedule(static) reduction(+:e1)
          
          do a=1, na
              e1 = e1 + ( &
                  tau10(a) * y(a)&
              )
          end do
          
          !$omp end do

          deallocate(tau10)

          allocate(tau11(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau11(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau11(a) = tau11(a) + ( &
                              t2(a, b, c, d)**2&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp single
          
          e2 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:e2)
          
          do a=1, na
              e2 = e2 + ( &
                  x(a)**2 * e0(a) * tau11(a) / 2&
              )
          end do
          
          !$omp end do

          deallocate(tau11)

          allocate(tau12(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau12(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau12(a) = tau12(a) + ( &
                          t1(c, b) * t2(a, c, a, b)&
                      )
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
          
              tau17(a) = tau17(a) + ( &
                  e0(a) * tau12(a) * x(a)&
              )
          
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
          
                          tau13(a, b, c, d) = tau13(a, b, c, d) + ( &
                              eri(d, c, b, a)&
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
                              tau17(a) = tau17(a) + ( &
                                  x(c) * y(d) * y(i) * t1(b, a) * t2(b, c, i, d) * tau13(a, c, d, i) / 4&
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
                          tau17(a) = tau17(a) + ( &
                              t0 * x(b) * x(c) * y(d) * t2(b, c, a, d) * tau13(a, d, b, c) / 8&
                          )
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
                                  4 * x(c) * x(d) * y(i) * t1(a, b) * t2(d, c, b, i) * tau13(c, d, a, i)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau13)

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
                          do i=1, na
                              do j=1, na
                                  tau14(a, b, c, d) = tau14(a, b, c, d) + ( &
                                      t2(j, i, a, c) * t2(j, i, b, d)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau16(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau16(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau16(a, b) = tau16(a, b) - ( &
                              y(c) * y(d) * tau14(a, c, b, d) * eri(d, c, b, a)&
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
                          do i=1, na
                              do j=1, na
                                  tau15(a, b, c, d) = tau15(a, b, c, d) + ( &
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
                          tau16(a, b) = tau16(a, b) + ( &
                              8 * x(c) * x(d) * tau15(d, c, a, b) * eri(a, c, b, d)&
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
                  tau17(a) = tau17(a) - ( &
                      y(b) * tau16(a, b) / 16&
                  )
              end do
          end do
          !$omp end do

          deallocate(tau16)

          !$omp do schedule(static) reduction(+:e2)
          
          do a=1, na
              e2 = e2 + ( &
                  2 * tau17(a) * y(a)&
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
                                      t2(a, c, i, j) * t2(b, d, i, j)&
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
                              x(c) * x(d) * tau18(c, a, d, b) * eri(b, a, d, c)&
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
                  tau20(a) = tau20(a) + ( &
                      x(b) * tau19(b, a)&
                  )
              end do
          end do
          !$omp end do

          deallocate(tau19)

          !$omp do schedule(static) reduction(+:e2)
          
          do a=1, na
              e2 = e2 + ( &
                  tau20(a) * x(a) / 8&
              )
          end do
          
          !$omp end do

          deallocate(tau20)

          !$omp single
          
          tau21 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:tau21)
          
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau21 = tau21 + ( &
                              t2(a, b, c, d)**2&
                          )
                      end do
                  end do
              end do
          end do
          
          !$omp end do

          !$omp single
          
          ov2 = 0.0
          
          !$omp end single
          

          !$omp single
          
          
          ov2 = ov2 + ( &
              tau21 / 4&
          )
          
          
          !$omp end single

          allocate(tau22(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau22(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
          
              tau22(a) = tau22(a) + ( &
                  2 * e0(a)&
              )
          
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau22(a) = tau22(a) + ( &
                      y(b)**2 * eri(a, b, a, b)&
                  )
              end do
          end do
          !$omp end do

          !$omp single
          
          tau23 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:tau23)
          
          do a=1, na
              tau23 = tau23 + ( &
                  y(a)**2 * tau22(a)&
              )
          end do
          
          !$omp end do

          deallocate(tau22)

          !$omp single
          
          
          e2 = e2 + ( &
              tau21*tau23 / 8&
          )
          
          
          !$omp end single

          allocate(tau24(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau24(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau24(a) = tau24(a) + ( &
                              t2(b, c, a, d)**2&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau30(1:na))
          !$omp do schedule(static)
          do a=1, na
              tau30(a) = 0.0
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
          
              tau30(a) = tau30(a) - ( &
                  e0(a) * tau24(a)&
              )
          
          end do
          !$omp end do

          deallocate(tau24)

          allocate(tau25(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau25(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau25(a, b) = tau25(a, b) + ( &
                              t1(d, c) * t2(d, a, c, b)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau28(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau28(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  tau28(a, b) = tau28(a, b) + ( &
                      2 * y(b) * tau25(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau30(a) = tau30(a) + ( &
                          2 * x(b) * y(c) * tau25(b, c) * eri(a, c, a, b)&
                      )
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau25)

          allocate(tau26(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau26(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              tau26(a, b) = tau26(a, b) + ( &
                                  t2(i, c, d, a) * t2(i, c, d, b)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau29(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau29(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  tau29(a, b) = tau29(a, b) - ( &
                      y(a) * y(b) * tau26(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          deallocate(tau26)

          allocate(tau27(1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau27(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              tau27(a, b) = tau27(a, b) + ( &
                                  t2(i, a, c, d) * t2(i, b, c, d)&
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
          
                  tau28(a, b) = tau28(a, b) + ( &
                      x(b) * tau27(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          deallocate(tau27)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  tau29(a, b) = tau29(a, b) + ( &
                      x(a) * tau28(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          deallocate(tau28)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau30(a) = tau30(a) + ( &
                          tau29(b, c) * eri(a, b, a, c)&
                      )
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau29)

          !$omp do schedule(static) reduction(+:e2)
          
          do a=1, na
              e2 = e2 + ( &
                  y(a)**2 * tau30(a) / 2&
              )
          end do
          
          !$omp end do

          deallocate(tau30)

          !$omp end parallel

      End Subroutine EnAndOvlp
