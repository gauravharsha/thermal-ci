    !$omp parallel default(shared)

    !$omp do schedule(static)
    do p=1, na
        do q=1, na
            tm1(p, q) = 0.0
        end do
    end do
    !$omp end do
    

    !$omp do schedule(static)
    do p=1, na
        do q=1, na
            do a=1, na
                tm1(p, q) = tm1(p, q) + ( &
                    e0(a) * x(a) * y(a) * t2(a, p, a, q) / 2&
                )
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do p=1, na
        do q=1, na
            do a=1, na
                tm1(p, q) = tm1(p, q) - ( &
                    \mu * x(a) * y(a) * t2(a, p, a, q) / 2&
                )
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do p=1, na
        do q=1, na
            sm1(p, q) = 0.0
        end do
    end do
    !$omp end do
    

    !$omp do schedule(static)
    do p=1, na
        do q=1, na
            do a=1, na
                sm1(p, q) = sm1(p, q) + ( &
                    e0(a) * x(a) * y(a) * s2(a, p, a, q) / 2&
                )
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do p=1, na
        do q=1, na
            do a=1, na
                sm1(p, q) = sm1(p, q) - ( &
                    \mu * x(a) * y(a) * s2(a, p, a, q) / 2&
                )
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do p=1, na
        do q=1, na
            do r=1, na
                do a=1, na
                    do b=1, na
                        do c=1, na
                            sm3(p, q, r, a, b, c) = 0.0
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
            do r=1, na
                do a=1, na
                    do b=1, na
                        do c=1, na
                            do d=1, na
                                sm3(p, q, r, a, b, c) = sm3(p, q, r, a, b, c) + ( &
                                    \mu * x(d) * y(d) * s4(d, p, q, r, a, b, c, d)&
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
    do p=1, na
        do q=1, na
            do r=1, na
                do a=1, na
                    do b=1, na
                        do c=1, na
                            do d=1, na
                                sm3(p, q, r, a, b, c) = sm3(p, q, r, a, b, c) - ( &
                                    e0(d) * x(d) * y(d) * s4(d, p, q, r, a, b, c, d)&
                                )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp end parallel

