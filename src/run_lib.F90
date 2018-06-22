module run_lib
    contains
    subroutine run_mm(N)
    use matrix_lib, only : mm

    real(kind= 8), allocatable :: first(:, :), second(:, :), result(:, :)
    real :: start, finish
    integer(kind=4) :: i, j
    integer(kind=4), intent(in) :: N

    allocate(first(n, n))
    allocate(second(n, n))
    allocate(result(n, n))

    do i=1,N
        do j=1,N
            first(i, j) = i*j + 10
            second(i, j) = i**j + 1
        enddo
    enddo

    call cpu_time(start)
    call mm(first, second, result, n, n, n)
    call cpu_time(finish)
    print '(i8, f10.4)',N, finish-start

    end subroutine

    subroutine run_gauss(N)
    use matrix_lib, only : gauss_elimination

    real(kind= 8), allocatable :: A(:, :), X(:)
    real :: start, finish
    integer(kind=4) :: i
    integer(kind=4), intent(in) :: N
    real(kind = 8) :: P1, P2, h2

    allocate(A(n, n))
    allocate(X(n))

    h2 = 1.0 / real(N * N)
    P1 = 1.0 / h2
    P2 = -2.0 / h2

    A = 0
    do i = 1,N
        A(i, i) = P2
        if (i .NE. 1) then 
            A(i, i - 1) = P1 
        endif
        if (i .NE. N) then 
            A(i, i + 1) = P1 
        endif
    enddo

    X(:) = 0.0
    X(N) = 1

    call cpu_time(start)
    call gauss_elimination(A, X, N)
    call cpu_time(finish)
    print '(i8, f10.4)',N, finish-start
    X = X * real(N * (N + 1) * (-1))
    write (*, *) X

    end subroutine

end module