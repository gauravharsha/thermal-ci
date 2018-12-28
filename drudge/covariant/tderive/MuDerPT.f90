    !$omp parallel default(shared)

    !$omp single
    
    tm0 = 0.0
    
    !$omp end single
    

    !$omp do schedule(static) reduction(+:tm0)
    
    do a=1, na
        tm0 = tm0 - ( &
            \beta * x(a) * y(a) * t1(a, a) / 2&
        )
    end do
    
    !$omp end do

    !$omp do schedule(static)
    do a=1, na
        do b=1, na
            tm1(a, b) = 0.0
        end do
    end do
    !$omp end do
    

    !$omp do schedule(static)
    do a=1, na
        do b=1, na
            do c=1, na
                tm1(a, b) = tm1(a, b) - ( &
                    \beta * x(c) * y(c) * t2(c, a, c, b) / 2&
                )
            end do
        end do
    end do
    !$omp end do

    !$omp end parallel

