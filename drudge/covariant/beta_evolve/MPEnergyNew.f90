    Real (Kind=8), dimension(:), allocatable :: tau0

    Real (Kind=8), dimension(:, :, :, :), allocatable :: tau1

    Real (Kind=8), dimension(:), allocatable :: tau2

    Real (Kind=8), dimension(:, :, :), allocatable :: tau3

    Real (Kind=8), dimension(:), allocatable :: tau4

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
                y(b)**2 * u(b, a, b, a)&
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
            y(a)**2 * tau0(a) / 2&
        )
    end do
    
    !$omp end do

    deallocate(tau0)

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
        do b=1, na
            do c=1, na
                do d=1, na
    
                    tau1(a, b, c, d) = tau1(a, b, c, d) + ( &
                        4 * t2(a, b, c, d) * u(b, a, d, c)&
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
    
                    tau1(a, b, c, d) = tau1(a, b, c, d) + ( &
                        t2(a, b, d, c) * u(c, d, b, a)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    allocate(tau2(1:na))
    !$omp do schedule(static)
    do a=1, na
        tau2(a) = 0.0
    end do
    !$omp end do
    

    !$omp do schedule(static)
    do a=1, na
        do b=1, na
            do c=1, na
                do d=1, na
                    tau2(a) = tau2(a) + ( &
                        x(b) * x(c) * y(d) * tau1(b, c, a, d)&
                    )
                end do
            end do
        end do
    end do
    !$omp end do

    deallocate(tau1)

    !$omp do schedule(static)
    do a=1, na
    
        tau2(a) = tau2(a) + ( &
            8 * e0(a) * x(a) * s1(a, a)&
        )
    
    end do
    !$omp end do

    !$omp single
    
    e2 = 0.0
    
    !$omp end single
    

    !$omp do schedule(static) reduction(+:e2)
    
    do a=1, na
        e2 = e2 + ( &
            tau2(a) * y(a) / 4&
        )
    end do
    
    !$omp end do

    deallocate(tau2)

    allocate(tau3(1:na, 1:na, 1:na))
    !$omp do schedule(static)
    do a=1, na
        do b=1, na
            do c=1, na
                tau3(a, b, c) = 0.0
            end do
        end do
    end do
    !$omp end do
    

    !$omp do schedule(static)
    do a=1, na
        do b=1, na
            do c=1, na
    
                tau3(a, b, c) = tau3(a, b, c) + ( &
                    u(a, b, a, c)&
                )
    
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do a=1, na
        do b=1, na
            do c=1, na
    
                tau3(a, b, c) = tau3(a, b, c) + ( &
                    u(a, c, a, b)&
                )
    
            end do
        end do
    end do
    !$omp end do

    allocate(tau4(1:na))
    !$omp do schedule(static)
    do a=1, na
        tau4(a) = 0.0
    end do
    !$omp end do
    

    !$omp do schedule(static)
    do a=1, na
        do b=1, na
            do c=1, na
                tau4(a) = tau4(a) + ( &
                    x(b) * y(c) * t1(b, c) * tau3(a, b, c)&
                )
            end do
        end do
    end do
    !$omp end do

    deallocate(tau3)

    !$omp do schedule(static) reduction(+:e2)
    
    do a=1, na
        e2 = e2 + ( &
            y(a)**2 * tau4(a)&
        )
    end do
    
    !$omp end do

    deallocate(tau4)

    !$omp do schedule(static) reduction(+:e1)
    
    do a=1, na
        e1 = e1 + ( &
            2 * e0(a) * x(a) * y(a) * t1(a, a)&
        )
    end do
    
    !$omp end do

    !$omp single
    
    ov1 = 0.0
    
    !$omp end single
    

    !$omp do schedule(static) reduction(+:ov1)
    
    do a=1, na
        do b=1, na
            do c=1, na
                do d=1, na
                    ov1 = ov1 + ( &
                        t2(a, b, c, d)**2&
                    )
                end do
            end do
        end do
    end do
    
    !$omp end do

    !$omp do schedule(static) reduction(+:ov1)
    
    do a=1, na
        do b=1, na
            ov1 = ov1 + ( &
                t1(a, b)**2&
            )
        end do
    end do
    
    !$omp end do

    !$omp single
    
    ov2 = 0.0
    
    !$omp end single
    

    !$omp do schedule(static) reduction(+:ov2)
    
    do a=1, na
        do b=1, na
            do c=1, na
                do d=1, na
                    do i=1, na
                        do j=1, na
                            ov2 = ov2 + ( &
                                s3(a, b, c, d, i, j)**2 / 36&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    !$omp end do

    !$omp do schedule(static) reduction(+:ov2)
    
    do a=1, na
        do b=1, na
            ov2 = ov2 + ( &
                s1(a, b)**2&
            )
        end do
    end do
    
    !$omp end do

    !$omp do schedule(static) reduction(+:ov2)
    
    do a=1, na
        do b=1, na
            do c=1, na
                do d=1, na
                    ov2 = ov2 + ( &
                        s2(a, b, c, d)**2 / 4&
                    )
                end do
            end do
        end do
    end do
    
    !$omp end do

    !$omp do schedule(static) reduction(+:ov2)
    
    do a=1, na
        do b=1, na
            do c=1, na
                do d=1, na
                    do i=1, na
                        do j=1, na
                            do k=1, na
                                do l=1, na
                                    ov2 = ov2 + ( &
                                        s4(a, b, c, d, i, j, k, l)**2 / 576&
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

