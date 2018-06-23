subroutine run_mm(N)
    use matrix_lib_coarr, only : mm

    real(kind= 8), allocatable, codimension[:, :], dimension(:, :) :: first, second, result
    integer :: start, finish, count_rate, count_max
    integer(kind=4), intent(in) :: N

    allocate(first(n, n)[num_images(), *])
    allocate(second(n, n)[num_images(), *])
    allocate(result(n, n)[num_images(), *])

    first = 1
    second = 2

    call system_clock(start, count_rate, count_max)
    call mm(first, second, result, n, n, n)
    call system_clock(finish, count_rate, count_max)

    if (this_image() .EQ. 1) then
        print '(i8, f10.4)',N, real(finish-start)/real(count_rate)
    endif

end subroutine

subroutine run_gauss(N)
    use matrix_lib_coarr, only : gauss_elimination

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

end subroutine

program main_coarr

    integer(kind=4) :: n, mode
    integer(kind=4) :: l, s
    character(len=10) :: buffer

    call get_command_argument(1, buffer, l, s)
    read(buffer,*) n

    call get_command_argument(2, buffer, l, s)
    read(buffer,*) mode

    select case (mode)
        case (0)
            ! write (*, *) "mm"
            call run_mm(n)
        case (1)
            write (*, *) "gauss"
            call run_gauss(n)
        case default
            write (*, *) "incorrect mode"
    end select

end program main_coarr