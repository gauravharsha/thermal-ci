      Subroutine CovBetaPT2(E0,ERI,T1,T2,S1,S2,S3,S4,NA,X,Y,RT1,RT2,RS1,RS2,RS3,RS4)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: E0(NA), ERI(NA,NA,NA,NA)
          Real (Kind=8), Intent(In) :: X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: T1(NA,NA), S1(NA,NA)
          Real (Kind=8), Intent(In) :: T2(NA,NA,NA,NA), S2(NA,NA,NA,NA)
          Real (Kind=8), Intent(In) :: S3(NA,NA,NA,NA,NA,NA)
          Real (Kind=8), Intent(In) :: S4(NA,NA,NA,NA,NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: RT1(NA,NA), RS1(NA,NA)
          Real (Kind=8), Intent(Out) :: RT2(NA,NA,NA,NA), RS2(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: RS3(NA,NA,NA,NA,NA,NA)
          Real (Kind=8), Intent(In) :: RS4(NA,NA,NA,NA,NA,NA,NA,NA)
          Integer :: i, j, k, l, a, b, c, d

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
          Real (Kind=8), dimension(:, :, :, :, :, :), allocatable :: tau11
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau12
          Real (Kind=8), dimension(:, :, :, :), allocatable :: tau13
          Real (Kind=8), dimension(:, :, :, :, :, :, :), allocatable :: tau14
          Real (Kind=8), dimension(:, :, :, :, :, :, :), allocatable :: tau15
          Real (Kind=8), dimension(:, :, :, :, :, :, :), allocatable :: tau16
          Real (Kind=8), dimension(:, :, :, :, :, :, :), allocatable :: tau17
          Real (Kind=8), dimension(:, :, :, :, :, :, :), allocatable :: tau18
          Real (Kind=8), dimension(:, :, :, :, :, :), allocatable :: tau19
          Real (Kind=8), dimension(:, :, :, :, :, :), allocatable :: tau20
          Real (Kind=8), dimension(:, :, :, :, :, :), allocatable :: tau21
          Real (Kind=8), dimension(:, :, :, :, :, :), allocatable :: tau22
          Real (Kind=8), dimension(:, :, :, :, :, :), allocatable :: tau23
          Real (Kind=8), dimension(:, :, :, :, :, :, :, :), allocatable :: tau24
          Real (Kind=8), dimension(:, :, :, :, :, :, :, :), allocatable :: tau25
          Real (Kind=8), dimension(:, :, :, :, :, :, :, :), allocatable :: tau26
          Real (Kind=8), dimension(:, :, :, :, :, :, :, :), allocatable :: tau27
          Real (Kind=8), dimension(:, :, :, :, :, :, :, :), allocatable :: tau28
          Real (Kind=8), dimension(:, :, :, :, :, :, :, :), allocatable :: tau29
          Real (Kind=8), dimension(:, :, :, :, :, :, :, :), allocatable :: tau30
          Real (Kind=8), dimension(:, :, :, :, :, :, :, :), allocatable :: tau31

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
                          1.0 * t1(a, b) * tau0(b, c)&
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
                              -1.0 * y(d) * t1(a, b) * tau0(c, d)&
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
          
                          tau12(a, b, c, d) = tau12(a, b, c, d) + ( &
                              1.0 * y(b) * t1(a, d) * tau0(c, b)&
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
                              -0.5 * x(b) * tau12(a, d, b, c)&
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
                              1.0 * y(d) * t1(a, b) * tau0(c, d)&
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
                              -1.0 * y(b) * t1(a, d) * tau0(c, b)&
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
                              -0.5 * x(a) * tau13(b, d, a, c)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau13)

          allocate(tau21(1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  tau21(a, b, c, d, i, j) = 0.0
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
          
                                  tau21(a, b, c, d, i, j) = tau21(a, b, c, d, i, j) + ( &
                                      1.0 * tau0(d, j) * t2(a, c, b, i)&
                                  )
          
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          allocate(tau22(1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  tau22(a, b, c, d, i, j) = 0.0
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
          
                                  tau22(a, b, c, d, i, j) = tau22(a, b, c, d, i, j) + ( &
                                      1.0 * x(i) * tau0(i, j) * t2(a, b, c, d)&
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
          
                                  tau22(a, b, c, d, i, j) = tau22(a, b, c, d, i, j) + ( &
                                      1.0 * x(b) * tau0(b, j) * t2(i, a, c, d)&
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
          
                                  tau22(a, b, c, d, i, j) = tau22(a, b, c, d, i, j) + ( &
                                      -1.0 * x(a) * tau0(a, j) * t2(i, b, c, d)&
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
                                      0.5 * y(i) * tau22(b, c, d, j, a, i)&
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
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      -0.5 * y(j) * tau22(b, c, d, i, a, j)&
                                  )
          
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau22)

          allocate(tau23(1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  tau23(a, b, c, d, i, j) = 0.0
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
          
                                  tau23(a, b, c, d, i, j) = tau23(a, b, c, d, i, j) + ( &
                                      1.0 * x(i) * tau0(i, j) * t2(a, b, c, d)&
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
          
                                  tau23(a, b, c, d, i, j) = tau23(a, b, c, d, i, j) + ( &
                                      -1.0 * x(b) * tau0(b, j) * t2(a, i, c, d)&
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
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      -0.5 * y(d) * tau23(a, b, i, j, c, d)&
                                  )
          
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau23)

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
                      -0.25d0*tau2 * t1(a, b)&
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
                              -0.25d0*tau2 * t2(a, b, c, d)&
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
                              1.0 * y(d) * t1(a, b) * eri(c, b, d, a)&
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
          
                              tau8(a, b, c, d, i) = tau8(a, b, c, d, i) + ( &
                                  -1.0 * y(c) * t1(a, b) * eri(i, d, c, a)&
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
                                  1.0 * y(b) * t1(a, c) * eri(i, d, b, a)&
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
                                  1.0 * x(b) * t1(a, c) * eri(b, c, i, d)&
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
          
                              tau9(a, b, c, d, i) = tau9(a, b, c, d, i) + ( &
                                  -1.0 * x(a) * t1(b, c) * eri(a, c, i, d)&
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
                                      1.0 * y(j) * t2(a, b, c, d) * eri(i, d, j, a)&
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
          
                                  tau10(a, b, c, d, i, j) = tau10(a, b, c, d, i, j) + ( &
                                      -1.0 * y(c) * t2(a, b, j, d) * eri(i, d, c, a)&
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
                                  rs2(a, b, c, d) = rs2(a, b, c, d) + ( &
                                      -0.5 * x(a) * x(i) * y(j) * tau10(i, b, c, j, a, d)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau10)

          allocate(tau11(1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  tau11(a, b, c, d, i, j) = 0.0
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
          
                                  tau11(a, b, c, d, i, j) = tau11(a, b, c, d, i, j) + ( &
                                      -1.0 * y(j) * t2(a, b, c, d) * eri(i, d, j, a)&
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
          
                                  tau11(a, b, c, d, i, j) = tau11(a, b, c, d, i, j) + ( &
                                      1.0 * y(c) * t2(a, b, j, d) * eri(i, d, c, a)&
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
                                  rs2(a, b, c, d) = rs2(a, b, c, d) + ( &
                                      -0.5 * x(b) * x(i) * y(j) * tau11(i, a, c, j, b, d)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau11)

          allocate(tau14(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      tau14(a, b, c, d, i, j, k) = 0.0
                                  end do
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
                                  do k=1, na
          
                                      tau14(a, b, c, d, i, j, k) = tau14(a, b, c, d, i, j, k) + ( &
                                          -1.0 * y(k) * t2(a, b, c, d) * eri(j, i, k, a)&
                                      )
          
                                  end do
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
                                  do k=1, na
          
                                      tau14(a, b, c, d, i, j, k) = tau14(a, b, c, d, i, j, k) + ( &
                                          1.0 * y(d) * t2(a, b, c, k) * eri(j, i, d, a)&
                                      )
          
                                  end do
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
                                  do k=1, na
          
                                      tau14(a, b, c, d, i, j, k) = tau14(a, b, c, d, i, j, k) + ( &
                                          1.0 * y(c) * t2(a, b, k, d) * eri(j, i, c, a)&
                                      )
          
                                  end do
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
                                  do k=1, na
                                      rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                          0.5 * x(a) * x(b) * x(k) * tau14(k, c, d, j, a, b, i)&
                                      )
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau14)

          allocate(tau15(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      tau15(a, b, c, d, i, j, k) = 0.0
                                  end do
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
                                  do k=1, na
          
                                      tau15(a, b, c, d, i, j, k) = tau15(a, b, c, d, i, j, k) + ( &
                                          1.0 * y(k) * t2(a, b, c, d) * eri(j, i, k, a)&
                                      )
          
                                  end do
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
                                  do k=1, na
          
                                      tau15(a, b, c, d, i, j, k) = tau15(a, b, c, d, i, j, k) + ( &
                                          -1.0 * y(d) * t2(a, b, c, k) * eri(j, i, d, a)&
                                      )
          
                                  end do
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
                                  do k=1, na
          
                                      tau15(a, b, c, d, i, j, k) = tau15(a, b, c, d, i, j, k) + ( &
                                          1.0 * y(c) * t2(a, b, d, k) * eri(j, i, c, a)&
                                      )
          
                                  end do
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
                                  do k=1, na
                                      rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                          -0.5 * x(a) * x(c) * x(k) * tau15(k, b, d, i, a, c, j)&
                                      )
                                  end do
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
                                  do k=1, na
                                      rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                          -0.5 * x(b) * x(c) * x(k) * tau15(k, a, d, i, c, b, j)&
                                      )
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau15)

          allocate(tau16(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      tau16(a, b, c, d, i, j, k) = 0.0
                                  end do
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
                                  do k=1, na
          
                                      tau16(a, b, c, d, i, j, k) = tau16(a, b, c, d, i, j, k) + ( &
                                          -1.0 * x(i) * t2(a, b, c, d) * eri(d, i, k, j)&
                                      )
          
                                  end do
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
                                  do k=1, na
          
                                      tau16(a, b, c, d, i, j, k) = tau16(a, b, c, d, i, j, k) + ( &
                                          -1.0 * x(b) * t2(i, a, c, d) * eri(d, b, k, j)&
                                      )
          
                                  end do
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
                                  do k=1, na
          
                                      tau16(a, b, c, d, i, j, k) = tau16(a, b, c, d, i, j, k) + ( &
                                          1.0 * x(a) * t2(i, b, c, d) * eri(d, a, k, j)&
                                      )
          
                                  end do
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
                                  do k=1, na
                                      rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                          0.5 * y(d) * y(j) * y(k) * tau16(b, c, i, k, a, j, d)&
                                      )
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau16)

          allocate(tau17(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      tau17(a, b, c, d, i, j, k) = 0.0
                                  end do
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
                                  do k=1, na
          
                                      tau17(a, b, c, d, i, j, k) = tau17(a, b, c, d, i, j, k) + ( &
                                          1.0 * x(i) * t2(a, b, c, d) * eri(d, i, k, j)&
                                      )
          
                                  end do
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
                                  do k=1, na
          
                                      tau17(a, b, c, d, i, j, k) = tau17(a, b, c, d, i, j, k) + ( &
                                          -1.0 * x(b) * t2(a, i, c, d) * eri(d, b, k, j)&
                                      )
          
                                  end do
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
                                  do k=1, na
          
                                      tau17(a, b, c, d, i, j, k) = tau17(a, b, c, d, i, j, k) + ( &
                                          -1.0 * x(a) * t2(i, b, c, d) * eri(d, a, k, j)&
                                      )
          
                                  end do
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
                                  do k=1, na
                                      rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                          -0.5 * y(d) * y(i) * y(k) * tau17(a, c, j, k, b, i, d)&
                                      )
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau17)

          allocate(tau18(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      tau18(a, b, c, d, i, j, k) = 0.0
                                  end do
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
                                  do k=1, na
          
                                      tau18(a, b, c, d, i, j, k) = tau18(a, b, c, d, i, j, k) + ( &
                                          -1.0 * x(i) * t2(a, b, c, d) * eri(d, i, k, j)&
                                      )
          
                                  end do
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
                                  do k=1, na
          
                                      tau18(a, b, c, d, i, j, k) = tau18(a, b, c, d, i, j, k) + ( &
                                          1.0 * x(b) * t2(a, i, c, d) * eri(d, b, k, j)&
                                      )
          
                                  end do
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
                                  do k=1, na
          
                                      tau18(a, b, c, d, i, j, k) = tau18(a, b, c, d, i, j, k) + ( &
                                          -1.0 * x(a) * t2(b, i, c, d) * eri(d, a, k, j)&
                                      )
          
                                  end do
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
                                  do k=1, na
                                      rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                          -0.5 * y(i) * y(j) * y(k) * tau18(a, b, d, k, c, j, i)&
                                      )
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau18)

          allocate(tau19(1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  tau19(a, b, c, d, i, j) = 0.0
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
          
                                  tau19(a, b, c, d, i, j) = tau19(a, b, c, d, i, j) + ( &
                                      -1.0 * y(i) * t1(a, b) * eri(d, c, i, j)&
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
          
                                  tau19(a, b, c, d, i, j) = tau19(a, b, c, d, i, j) + ( &
                                      1.0 * y(b) * t1(a, i) * eri(d, c, b, j)&
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
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      -0.5 * x(a) * x(b) * y(j) * tau19(c, i, b, a, d, j)&
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
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      -0.5 * x(a) * x(c) * y(i) * tau19(b, j, c, a, d, i)&
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
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      -0.5 * x(b) * x(c) * y(j) * tau19(a, i, c, b, d, j)&
                                  )
          
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau19)

          allocate(tau20(1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  tau20(a, b, c, d, i, j) = 0.0
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
          
                                  tau20(a, b, c, d, i, j) = tau20(a, b, c, d, i, j) + ( &
                                      1.0 * x(c) * t1(a, b) * eri(d, c, j, i)&
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
          
                                  tau20(a, b, c, d, i, j) = tau20(a, b, c, d, i, j) + ( &
                                      -1.0 * x(a) * t1(c, b) * eri(d, a, j, i)&
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
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      0.5 * x(b) * y(d) * y(i) * tau20(c, j, a, b, i, d)&
                                  )
          
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau20)

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
          
                                  tau21(a, b, c, d, i, j) = tau21(a, b, c, d, i, j) + ( &
                                      -1.0 * x(c) * y(i) * t1(a, b) * eri(d, c, i, j)&
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
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      -0.5 * x(a) * y(d) * tau21(b, i, c, a, j, d)&
                                  )
          
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau21)

          allocate(tau24(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      do l=1, na
                                          tau24(a, b, c, d, i, j, k, l) = 0.0
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau24(a, b, c, d, i, j, k, l) = tau24(a, b, c, d, i, j, k, l) + ( &
                                              -1.0 * y(k) * t2(a, b, c, d) * eri(i, j, k, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau24(a, b, c, d, i, j, k, l) = tau24(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * y(d) * t2(a, b, c, k) * eri(i, j, d, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau24(a, b, c, d, i, j, k, l) = tau24(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * y(c) * t2(a, b, k, d) * eri(i, j, c, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
                                          rs4(a, b, c, d, i, j, k, l) = 0.0
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              0.5 * x(a) * x(b) * y(i) * tau24(c, d, j, l, a, b, k, i)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              0.5 * x(a) * x(d) * y(i) * tau24(b, c, j, l, a, d, k, i)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              -0.5 * x(b) * x(c) * y(i) * tau24(a, d, j, l, c, b, k, i)&
                                          )
          
                                      end do
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau24)

          allocate(tau25(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      do l=1, na
                                          tau25(a, b, c, d, i, j, k, l) = 0.0
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau25(a, b, c, d, i, j, k, l) = tau25(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * x(i) * t2(a, b, c, d) * eri(j, i, l, k)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau25(a, b, c, d, i, j, k, l) = tau25(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * x(b) * t2(i, a, c, d) * eri(j, b, l, k)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau25(a, b, c, d, i, j, k, l) = tau25(a, b, c, d, i, j, k, l) + ( &
                                              -1.0 * x(a) * t2(i, b, c, d) * eri(j, a, l, k)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              0.5 * x(a) * y(j) * y(l) * tau25(c, d, i, k, b, a, l, j)&
                                          )
          
                                      end do
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau25)

          allocate(tau26(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      do l=1, na
                                          tau26(a, b, c, d, i, j, k, l) = 0.0
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau26(a, b, c, d, i, j, k, l) = tau26(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * y(k) * t2(a, b, c, d) * eri(i, j, k, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau26(a, b, c, d, i, j, k, l) = tau26(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * y(d) * t2(a, b, k, c) * eri(i, j, d, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau26(a, b, c, d, i, j, k, l) = tau26(a, b, c, d, i, j, k, l) + ( &
                                              -1.0 * y(c) * t2(a, b, k, d) * eri(i, j, c, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              -0.5 * x(a) * x(c) * y(i) * tau26(b, d, k, l, a, c, j, i)&
                                          )
          
                                      end do
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau26)

          allocate(tau27(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      do l=1, na
                                          tau27(a, b, c, d, i, j, k, l) = 0.0
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau27(a, b, c, d, i, j, k, l) = tau27(a, b, c, d, i, j, k, l) + ( &
                                              -1.0 * x(i) * t2(a, b, c, d) * eri(j, i, l, k)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau27(a, b, c, d, i, j, k, l) = tau27(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * x(b) * t2(a, i, c, d) * eri(j, b, l, k)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau27(a, b, c, d, i, j, k, l) = tau27(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * x(a) * t2(i, b, c, d) * eri(j, a, l, k)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              -0.5 * x(a) * y(j) * y(k) * tau27(b, d, i, l, c, a, k, j)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              -0.5 * x(a) * y(k) * y(l) * tau27(b, d, i, j, c, a, l, k)&
                                          )
          
                                      end do
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau27)

          allocate(tau28(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      do l=1, na
                                          tau28(a, b, c, d, i, j, k, l) = 0.0
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau28(a, b, c, d, i, j, k, l) = tau28(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * y(k) * t2(a, b, c, d) * eri(i, j, k, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau28(a, b, c, d, i, j, k, l) = tau28(a, b, c, d, i, j, k, l) + ( &
                                              -1.0 * y(d) * t2(a, b, c, k) * eri(i, j, d, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau28(a, b, c, d, i, j, k, l) = tau28(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * y(c) * t2(a, b, d, k) * eri(i, j, c, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              0.5 * x(b) * x(d) * y(i) * tau28(a, c, j, k, d, b, l, i)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              -0.5 * x(c) * x(d) * y(i) * tau28(a, b, j, k, d, c, l, i)&
                                          )
          
                                      end do
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau28)

          allocate(tau29(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      do l=1, na
                                          tau29(a, b, c, d, i, j, k, l) = 0.0
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau29(a, b, c, d, i, j, k, l) = tau29(a, b, c, d, i, j, k, l) + ( &
                                              -1.0 * x(i) * t2(a, b, c, d) * eri(j, i, l, k)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau29(a, b, c, d, i, j, k, l) = tau29(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * x(b) * t2(a, i, c, d) * eri(j, b, l, k)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              -0.5 * x(b) * y(j) * y(l) * tau29(a, d, i, k, c, b, l, j)&
                                          )
          
                                      end do
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau29)

          allocate(tau30(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      do l=1, na
                                          tau30(a, b, c, d, i, j, k, l) = 0.0
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau30(a, b, c, d, i, j, k, l) = tau30(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * x(i) * t2(a, b, c, d) * eri(j, i, l, k)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau30(a, b, c, d, i, j, k, l) = tau30(a, b, c, d, i, j, k, l) + ( &
                                              -1.0 * x(b) * t2(a, i, c, d) * eri(j, b, l, k)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              0.5 * x(c) * y(j) * y(k) * tau30(a, d, i, l, b, c, k, j)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              0.5 * x(c) * y(k) * y(l) * tau30(a, d, i, j, b, c, l, k)&
                                          )
          
                                      end do
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau30)

          allocate(tau31(1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na, 1:na))
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              do j=1, na
                                  do k=1, na
                                      do l=1, na
                                          tau31(a, b, c, d, i, j, k, l) = 0.0
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau31(a, b, c, d, i, j, k, l) = tau31(a, b, c, d, i, j, k, l) + ( &
                                              -1.0 * y(k) * t2(a, b, c, d) * eri(i, j, k, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          tau31(a, b, c, d, i, j, k, l) = tau31(a, b, c, d, i, j, k, l) + ( &
                                              1.0 * y(d) * t2(a, b, c, k) * eri(i, j, d, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              -0.5 * x(b) * x(d) * y(k) * tau31(a, c, i, j, d, b, l, k)&
                                          )
          
                                      end do
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          deallocate(tau31)

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

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              rs2(a, b, c, d) = rs2(a, b, c, d) + ( &
                                  -0.5 * e0(i) * x(i) * y(i) * s3(i, a, b, c, d, i)&
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
                              do j=1, na
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      0.5 * y(d)**2 * e0(d) * s3(a, b, c, d, i, j)&
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
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      0.5 * y(i)**2 * e0(i) * s3(a, b, c, d, i, j)&
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
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      0.5 * y(j)**2 * e0(j) * s3(a, b, c, d, i, j)&
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
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      -0.5 * x(a)**2 * e0(a) * s3(a, b, c, d, i, j)&
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
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      -0.5 * x(b)**2 * e0(b) * s3(a, b, c, d, i, j)&
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
          
                                  rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                      -0.5 * x(c)**2 * e0(c) * s3(a, b, c, d, i, j)&
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
                                  do k=1, na
                                      rs3(a, b, c, d, i, j) = rs3(a, b, c, d, i, j) + ( &
                                          0.5 * e0(k) * x(k) * y(k) * s4(k, a, b, c, d, i, j, k)&
                                      )
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              0.5 * y(i)**2 * e0(i) * s4(a, b, c, d, i, j, k, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              0.5 * y(j)**2 * e0(j) * s4(a, b, c, d, i, j, k, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              0.5 * y(k)**2 * e0(k) * s4(a, b, c, d, i, j, k, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              0.5 * y(l)**2 * e0(l) * s4(a, b, c, d, i, j, k, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              -0.5 * x(a)**2 * e0(a) * s4(a, b, c, d, i, j, k, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              -0.5 * x(b)**2 * e0(b) * s4(a, b, c, d, i, j, k, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              -0.5 * x(c)**2 * e0(c) * s4(a, b, c, d, i, j, k, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              -0.5 * x(d)**2 * e0(d) * s4(a, b, c, d, i, j, k, l)&
                                          )
          
                                      end do
                                  end do
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
                                  do k=1, na
                                      do l=1, na
          
                                          rs4(a, b, c, d, i, j, k, l) = rs4(a, b, c, d, i, j, k, l) + ( &
                                              -0.5 * x(c) * x(d) * y(j) * y(l) * t2(a, b, i, k) * eri(c, d, l, j)&
                                          )
          
                                      end do
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp end parallel

      End Subroutine CovBetaPT2
