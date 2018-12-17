      Subroutine CovMuPT2(T1,T2,S1,S2,NA,X,Y,RT1,RT2,RS1,RS2)
          Implicit None
          Integer, Intent(In) :: NA
          Real (Kind=8), Intent(In) :: X(NA), Y(NA)
          Real (Kind=8), Intent(In) :: T1(NA,NA), S1(NA,NA)
          Real (Kind=8), Intent(In) :: T2(NA,NA,NA,NA), S2(NA,NA,NA,NA)
          Real (Kind=8), Intent(Out) :: RT1(NA,NA), RS1(NA,NA)
          Real (Kind=8), Intent(Out) :: RT2(NA,NA,NA,NA), RS2(NA,NA,NA,NA)
          Integer :: a, b, c, d

          !$omp parallel default(shared)

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
                      x(a)**2 * s1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
          
                  rs1(a, b) = rs1(a, b) - ( &
                      y(b)**2 * s1(a, b)&
                  )
          
              end do
          end do
          !$omp end do

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      rs1(a, b) = rs1(a, b) - ( &
                          x(c) * y(c) * s2(c, a, b, c)&
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
                              x(a)**2 * s2(a, b, c, d)&
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
                              x(b)**2 * s2(a, b, c, d)&
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
                              y(c)**2 * s2(a, b, c, d)&
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
                              y(d)**2 * s2(a, b, c, d)&
                          )
          
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp end parallel

      End Subroutine CovMuPT2
