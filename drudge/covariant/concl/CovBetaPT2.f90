      Subroutine CovBetaPT2(E0,ERI,T1,T2,S1,S2,NA,X,Y,RT1,RT2,RS1,RS2)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: E0(NA), ERI(NA,NA,NA,NA)
          Real (Kind=8), Intent(In) :: X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: T1(NA,NA), S1(NA,NA)
          Real (Kind=8), Intent(In) :: T2(NA,NA,NA,NA), S2(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: RT1(NA,NA), RS1(NA,NA)
          Real (Kind=8), Intent(Out) :: RT2(NA,NA,NA,NA), RS2(NA,NA,NA,NA)
          Integer :: i, j, a, b, c, d

          Real (Kind=8), dimension(:, :), allocatable :: tau0
          Real (Kind=8), dimension(:), allocatable :: tau1
          Real (Kind=8) :: tau2
          Real (Kind=8), dimension(:, :), allocatable :: tau3
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau4
          Real (Kind=8), dimension(:, :, :), allocatable :: tau5
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau6
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau7
          Real (Kind=8), dimension(:, :, :, :, :), allocatable :: tau8
          Real (Kind=8), dimension(:, :, :, :, :), allocatable :: tau9
          Real (Kind=8), dimension(:, :, :, :, :, :), allocatable :: tau10
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau11
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau12
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau13

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
                      tau0(a, b) = tau0(a, b) + ( &
                          y(c)**2 * eri(c, a, c, b)&
                      )
                  end do
              end do
          end do
          !$omp end do

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
                          x(c) * t1(c, a) * tau0(b, c)&
                      )
                  end do
              end do
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
                      -0.5 * x(a) * tau3(b, a)&
                  )
          
              end do
          end do
          !$omp end do

          deallocate(tau3)

          allocate(tau5(1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      tau5(a, b, c) = 0.0
                  end do
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
          
                      tau5(a, b, c) = tau5(a, b, c) + ( &
                          t1(a, b) * tau0(b, c)&
                      )
          
                  end do
              end do
          end do
          !$omp end do

          allocate(tau6(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau6(a, b, c, d) = 0.0
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
                              tau6(a, b, c, d) = tau6(a, b, c, d) + ( &
                                  x(d) * x(i) * tau0(d, i) * t2(i, a, b, c)&
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
                              0.5 * tau6(a, c, d, b)&
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
                              -0.5 * tau6(b, c, d, a)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau6)

          allocate(tau7(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau7(a, b, c, d) = 0.0
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
                              tau7(a, b, c, d) = tau7(a, b, c, d) + ( &
                                  y(d) * y(i) * tau0(i, d) * t2(a, b, c, i)&
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
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) + ( &
                              0.5 * tau7(a, b, c, d)&
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
                              -0.5 * tau7(a, b, d, c)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau7)

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
          
                          tau12(a, b, c, d) = tau12(a, b, c, d) + ( &
                              y(d) * t1(a, b) * tau0(c, d)&
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
          
                          tau12(a, b, c, d) = tau12(a, b, c, d) - ( &
                              y(b) * t1(a, d) * tau0(c, b)&
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
                              -0.5 * x(b) * tau12(a, c, b, d)&
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
          
                          tau13(a, b, c, d) = tau13(a, b, c, d) - ( &
                              y(d) * t1(a, b) * tau0(c, d)&
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
                              y(b) * t1(a, d) * tau0(c, b)&
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
                              -0.5 * x(a) * tau13(b, c, a, d)&
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
                  rt1(a, b) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rt1(a, b) = rt1(a, b) + ( &
                      -0.5 * x(a) * y(b) * tau0(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          rs1(a, b) = rs1(a, b) + ( &
                              0.5 * x(c) * y(d) * tau0(d, c) * t2(c, a, b, d)&
                          )
                      end do
                  end do
              end do
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
                      y(b)**2 * eri(a, b, a, b)&
                  )
              end do
          end do
          !$omp end do

          !$omp single
          
          tau2 = 0.0
          
          !$omp end single
          

          !$omp do schedule(static) reduction(+:tau2)
          
          do a=1, na
              tau2 = tau2 + ( &
                  y(a)**2 * tau1(a)&
              )
          end do
          
          !$omp end do

          deallocate(tau1)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rs1(a, b) = rs1(a, b) + ( &
                      -0.25*tau2 * t1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) + ( &
                              -0.25*tau2 * t2(a, b, c, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau4(1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau4(a, b, c, d) = 0.0
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
          
                          tau4(a, b, c, d) = tau4(a, b, c, d) + ( &
                              y(d) * t1(a, b) * eri(c, b, d, a)&
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
                              tau4(a, b, c, d) = tau4(a, b, c, d) + ( &
                                  -0.5 * x(i) * t2(i, a, d, b) * eri(c, b, a, i)&
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
                          rs1(a, b) = rs1(a, b) + ( &
                              -0.5 * x(a) * x(c) * y(d) * tau4(c, d, a, b)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau4)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              tau5(a, b, c) = tau5(a, b, c) + ( &
                                  0.5 * x(d) * y(i) * t2(d, a, i, b) * eri(b, i, c, d)&
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
                      rs1(a, b) = rs1(a, b) + ( &
                          0.5 * y(b) * y(c) * tau5(a, c, b)&
                      )
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau5)

          allocate(tau8(1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              tau8(a, b, c, d, i) = 0.0
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
                                  tau8(a, b, c, d, i) = tau8(a, b, c, d, i) + ( &
                                      -0.5 * x(j) * t2(j, a, b, c) * eri(i, d, a, j)&
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
          
                              tau8(a, b, c, d, i) = tau8(a, b, c, d, i) - ( &
                                  y(c) * t1(a, b) * eri(i, d, c, a)&
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
          
                              tau8(a, b, c, d, i) = tau8(a, b, c, d, i) + ( &
                                  y(b) * t1(a, c) * eri(i, d, b, a)&
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
                              rs2(a, b, c, d) = rs2(a, b, c, d) + ( &
                                  0.5 * x(a) * x(b) * x(i) * tau8(i, c, d, a, b)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau8)

          allocate(tau9(1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              tau9(a, b, c, d, i) = 0.0
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
                                  tau9(a, b, c, d, i) = tau9(a, b, c, d, i) + ( &
                                      -0.5 * y(j) * t2(a, b, j, c) * eri(c, j, i, d)&
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
          
                              tau9(a, b, c, d, i) = tau9(a, b, c, d, i) + ( &
                                  x(b) * t1(a, c) * eri(b, c, i, d)&
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
          
                              tau9(a, b, c, d, i) = tau9(a, b, c, d, i) - ( &
                                  x(a) * t1(b, c) * eri(a, c, i, d)&
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
                              rs2(a, b, c, d) = rs2(a, b, c, d) + ( &
                                  -0.5 * y(c) * y(d) * y(i) * tau9(a, b, i, d, c)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau9)

          allocate(tau10(1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  tau10(a, b, c, d, i, j) = 0.0
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
          
                                  tau10(a, b, c, d, i, j) = tau10(a, b, c, d, i, j) + ( &
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
          
                                  tau10(a, b, c, d, i, j) = tau10(a, b, c, d, i, j) - ( &
                                      y(c) * t2(a, b, j, d) * eri(i, d, c, a)&
                                  )
          
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

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
                              do j=1, na
                                  tau11(a, b, c, d) = tau11(a, b, c, d) + ( &
                                      0.5 * x(c) * x(i) * y(j) * tau10(i, a, b, j, c, d)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau10)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
          
                          rs2(a, b, c, d) = rs2(a, b, c, d) - ( &
                              tau11(b, c, a, d)&
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
                              tau11(a, d, b, c)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau11)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rt1(a, b) = rt1(a, b) + ( &
                      0.5 * y(b)**2 * e0(b) * t1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rt1(a, b) = rt1(a, b) + ( &
                      -0.5 * x(a)**2 * e0(a) * t1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      rt1(a, b) = rt1(a, b) + ( &
                          0.5 * e0(c) * x(c) * y(c) * t2(c, a, b, c)&
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
                              0.5 * y(c)**2 * e0(c) * t2(a, b, c, d)&
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
                              0.5 * y(d)**2 * e0(d) * t2(a, b, c, d)&
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
                              -0.5 * x(a)**2 * e0(a) * t2(a, b, c, d)&
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
                              -0.5 * x(b)**2 * e0(b) * t2(a, b, c, d)&
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
                              -0.5 * x(a) * x(b) * y(c) * y(d) * eri(b, a, d, c)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rs1(a, b) = rs1(a, b) + ( &
                      0.5 * y(b)**2 * e0(b) * s1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rs1(a, b) = rs1(a, b) + ( &
                      -0.5 * x(a)**2 * e0(a) * s1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      rs1(a, b) = rs1(a, b) + ( &
                          0.5 * e0(c) * x(c) * y(c) * s2(c, a, b, c)&
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
                              0.5 * y(c)**2 * e0(c) * s2(a, b, c, d)&
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
                              0.5 * y(d)**2 * e0(d) * s2(a, b, c, d)&
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
                              -0.5 * x(a)**2 * e0(a) * s2(a, b, c, d)&
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
                              -0.5 * x(b)**2 * e0(b) * s2(a, b, c, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp end parallel

      End Subroutine CovBetaPT2
