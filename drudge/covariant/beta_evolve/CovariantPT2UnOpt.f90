      Subroutine CovariantPT2UnOpt(T1,T2,S1,S2,NA,X,Y,RT1,RT2,RS1,RS2)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: T1(NA,NA), S1(NA,NA)
          Real (Kind=8), Intent(In) :: T2(NA,NA,NA,NA), S2(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: RT1(NA,NA), RS1(NA,NA)
          Real (Kind=8), Intent(Out) :: RT2(NA,NA,NA,NA), RS2(NA,NA,NA,NA)
          Integer :: p, q, i, a, b, c, d

          !$omp parallel default(shared)

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  rt1(p, a) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              rt1(p, a) = rt1(p, a) + ( &
                                  y(c)**2 * x(a) * y(b) * t2(a, p, b, q) * u(b, c, a, c)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              rt1(p, a) = rt1(p, a) + ( &
                                  x(a) * x(b) * x(p) * y(c) * t2(a, b, c, q) * u(c, p, a, b) / 2&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              rt1(p, a) = rt1(p, a) - ( &
                                  x(a) * y(b) * y(c) * y(q) * t2(a, p, b, c) * u(b, c, a, q) / 2&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
          
                  rt1(p, a) = rt1(p, a) + ( &
                      y(q)**2 * e0(q) * t1(p, q)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
          
                  rt1(p, a) = rt1(p, a) - ( &
                      x(p)**2 * e0(p) * t1(p, q)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do a=1, na
                      do b=1, na
                          rt1(p, a) = rt1(p, a) + ( &
                              y(a)**2 * y(b)**2 * t1(p, q) * u(a, b, a, b) / 2&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do a=1, na
                      do b=1, na
                          rt1(p, a) = rt1(p, a) + ( &
                              y(b)**2 * y(a) * y(q) * t1(p, a) * u(a, b, b, q)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do a=1, na
                      do b=1, na
                          rt1(p, a) = rt1(p, a) - ( &
                              y(b)**2 * x(a) * x(p) * t1(a, q) * u(b, p, a, b)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do a=1, na
                      do b=1, na
                          rt1(p, a) = rt1(p, a) + ( &
                              x(a) * x(p) * y(b) * y(q) * t1(a, b) * u(b, p, a, q)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do a=1, na
                      rt1(p, a) = rt1(p, a) - ( &
                          e0(a) * x(a) * y(a) * t2(a, p, a, q)&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          rt2(p, q, a, b) = 0.0
                      end do
                  end do
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
          
                          rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                              y(r)**2 * e0(r) * t2(p, q, r, s)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
          
                          rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                              y(s)**2 * e0(s) * t2(p, q, r, s)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
          
                          rt2(p, q, a, b) = rt2(p, q, a, b) - ( &
                              x(p)**2 * e0(p) * t2(p, q, r, s)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
          
                          rt2(p, q, a, b) = rt2(p, q, a, b) - ( &
                              x(q)**2 * e0(q) * t2(p, q, r, s)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              do b=1, na
                                  rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                                      y(a)**2 * y(b)**2 * t2(p, q, r, s) * u(a, b, a, b) / 2&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              do b=1, na
                                  rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                                      y(b)**2 * x(a) * x(q) * t2(a, p, r, s) * u(b, q, a, b)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              do b=1, na
                                  rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                                      y(b)**2 * y(a) * y(r) * t2(p, q, a, s) * u(a, b, b, r)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              do b=1, na
                                  rt2(p, q, a, b) = rt2(p, q, a, b) - ( &
                                      y(b)**2 * x(a) * x(p) * t2(a, q, r, s) * u(b, p, a, b)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              do b=1, na
                                  rt2(p, q, a, b) = rt2(p, q, a, b) - ( &
                                      y(b)**2 * y(a) * y(s) * t2(p, q, a, r) * u(a, b, b, s)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              do b=1, na
                                  rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                                      x(a) * x(p) * y(b) * y(r) * t2(a, q, b, s) * u(b, p, a, r)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              do b=1, na
                                  rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                                      x(a) * x(q) * y(b) * y(s) * t2(a, p, b, r) * u(b, q, a, s)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              do b=1, na
                                  rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                                      x(a) * x(b) * x(p) * x(q) * t2(a, b, r, s) * u(p, q, a, b) / 2&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              do b=1, na
                                  rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                                      y(a) * y(b) * y(r) * y(s) * t2(p, q, a, b) * u(a, b, r, s) / 2&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              do b=1, na
                                  rt2(p, q, a, b) = rt2(p, q, a, b) - ( &
                                      x(a) * x(p) * y(b) * y(s) * t2(a, q, b, r) * u(b, p, a, s)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              do b=1, na
                                  rt2(p, q, a, b) = rt2(p, q, a, b) - ( &
                                      x(a) * x(q) * y(b) * y(r) * t2(a, p, b, s) * u(b, q, a, r)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                                  y(a)**2 * x(p) * y(r) * t1(q, s) * u(a, p, a, r)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                                  y(a)**2 * x(q) * y(s) * t1(p, r) * u(a, q, a, s)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              rt2(p, q, a, b) = rt2(p, q, a, b) - ( &
                                  y(a)**2 * x(p) * y(s) * t1(q, r) * u(a, p, a, s)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              rt2(p, q, a, b) = rt2(p, q, a, b) - ( &
                                  y(a)**2 * x(q) * y(r) * t1(p, s) * u(a, q, a, r)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                                  x(a) * x(p) * x(q) * y(s) * t1(a, r) * u(p, q, a, s)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              rt2(p, q, a, b) = rt2(p, q, a, b) + ( &
                                  x(p) * y(a) * y(r) * y(s) * t1(q, a) * u(a, p, r, s)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              rt2(p, q, a, b) = rt2(p, q, a, b) - ( &
                                  x(a) * x(p) * x(q) * y(r) * t1(a, s) * u(p, q, a, r)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do a=1, na
                              rt2(p, q, a, b) = rt2(p, q, a, b) - ( &
                                  x(q) * y(a) * y(r) * y(s) * t1(p, a) * u(a, q, r, s)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  rs1(p, a) = 0.0
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
          
                  rs1(p, a) = rs1(p, a) + ( &
                      y(a)**2 * e0(a) * s1(p, a)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
          
                  rs1(p, a) = rs1(p, a) - ( &
                      x(p)**2 * e0(p) * s1(p, a)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do b=1, na
                      do c=1, na
                          rs1(p, a) = rs1(p, a) + ( &
                              y(b)**2 * y(c)**2 * t1(p, a) * u(b, c, b, c) / 2&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do b=1, na
                      do c=1, na
                          rs1(p, a) = rs1(p, a) - ( &
                              y(c)**2 * x(b) * x(p) * t1(b, a) * u(c, p, b, c)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do b=1, na
                      do c=1, na
                          rs1(p, a) = rs1(p, a) - ( &
                              y(c)**2 * y(a) * y(b) * t1(p, b) * u(b, c, a, c)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do b=1, na
                      do c=1, na
                          rs1(p, a) = rs1(p, a) - ( &
                              x(b) * x(p) * y(a) * y(c) * t1(b, c) * u(c, p, a, b)&
                          )
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do b=1, na
                      do c=1, na
                          do d=1, na
                              rs1(p, a) = rs1(p, a) - ( &
                                  y(d)**2 * x(b) * y(c) * t2(b, p, a, c) * u(c, d, b, d)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do b=1, na
                      do c=1, na
                          do d=1, na
                              rs1(p, a) = rs1(p, a) + ( &
                                  x(b) * y(a) * y(c) * y(d) * t2(b, p, c, d) * u(c, d, a, b) / 2&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do b=1, na
                      do c=1, na
                          do d=1, na
                              rs1(p, a) = rs1(p, a) - ( &
                                  x(b) * x(c) * x(p) * y(d) * t2(b, c, a, d) * u(d, p, b, c) / 2&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do a=1, na
                  do b=1, na
                      rs1(p, a) = rs1(p, a) + ( &
                          e0(b) * x(b) * y(b) * s2(b, p, a, b)&
                      )
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          rs2(p, q, a, b) = 0.0
                      end do
                  end do
              end do
          end do
          !$omp end do
          

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              do d=1, na
                                  rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                                      y(c)**2 * y(d)**2 * t2(p, q, a, b) * u(c, d, c, d) / 2&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              do d=1, na
                                  rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                                      y(d)**2 * x(c) * x(q) * t2(c, p, a, b) * u(d, q, c, d)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              do d=1, na
                                  rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                                      y(d)**2 * y(a) * y(c) * t2(p, q, b, c) * u(c, d, a, d)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              do d=1, na
                                  rs2(p, q, a, b) = rs2(p, q, a, b) - ( &
                                      y(d)**2 * x(c) * x(p) * t2(c, q, a, b) * u(d, p, c, d)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              do d=1, na
                                  rs2(p, q, a, b) = rs2(p, q, a, b) - ( &
                                      y(d)**2 * y(b) * y(c) * t2(p, q, a, c) * u(c, d, b, d)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              do d=1, na
                                  rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                                      x(c) * x(p) * y(a) * y(d) * t2(c, q, b, d) * u(d, p, a, c)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              do d=1, na
                                  rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                                      x(c) * x(q) * y(b) * y(d) * t2(c, p, a, d) * u(d, q, b, c)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              do d=1, na
                                  rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                                      x(c) * x(d) * x(p) * x(q) * t2(c, d, a, b) * u(p, q, c, d) / 2&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              do d=1, na
                                  rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                                      y(a) * y(b) * y(c) * y(d) * t2(p, q, c, d) * u(c, d, a, b) / 2&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              do d=1, na
                                  rs2(p, q, a, b) = rs2(p, q, a, b) - ( &
                                      x(c) * x(p) * y(b) * y(d) * t2(c, q, a, d) * u(d, p, b, c)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              do d=1, na
                                  rs2(p, q, a, b) = rs2(p, q, a, b) - ( &
                                      x(c) * x(q) * y(a) * y(d) * t2(c, p, b, d) * u(d, q, a, c)&
                                  )
                              end do
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
          
                          rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                              y(a)**2 * e0(a) * s2(p, q, a, b)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
          
                          rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                              y(b)**2 * e0(b) * s2(p, q, a, b)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
          
                          rs2(p, q, a, b) = rs2(p, q, a, b) - ( &
                              x(p)**2 * e0(p) * s2(p, q, a, b)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
          
                          rs2(p, q, a, b) = rs2(p, q, a, b) - ( &
                              x(q)**2 * e0(q) * s2(p, q, a, b)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              rs2(p, q, a, b) = rs2(p, q, a, b) - ( &
                                  e0(c) * x(c) * y(c) * s3(c, p, q, a, b, c)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                                  y(c)**2 * x(p) * y(b) * t1(q, a) * u(c, p, b, c)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                                  y(c)**2 * x(q) * y(a) * t1(p, b) * u(c, q, a, c)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              rs2(p, q, a, b) = rs2(p, q, a, b) - ( &
                                  y(c)**2 * x(p) * y(a) * t1(q, b) * u(c, p, a, c)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              rs2(p, q, a, b) = rs2(p, q, a, b) - ( &
                                  y(c)**2 * x(q) * y(b) * t1(p, a) * u(c, q, b, c)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                                  x(c) * x(p) * x(q) * y(a) * t1(c, b) * u(p, q, a, c)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              rs2(p, q, a, b) = rs2(p, q, a, b) + ( &
                                  x(p) * y(a) * y(b) * y(c) * t1(q, c) * u(c, p, a, b)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              rs2(p, q, a, b) = rs2(p, q, a, b) - ( &
                                  x(c) * x(p) * x(q) * y(b) * t1(c, a) * u(p, q, b, c)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  do a=1, na
                      do b=1, na
                          do c=1, na
                              rs2(p, q, a, b) = rs2(p, q, a, b) - ( &
                                  x(q) * y(a) * y(b) * y(c) * t1(p, c) * u(c, q, a, b)&
                              )
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp end parallel

      End Subroutine CovariantPT2UnOpt
