program main

    use run_lib, only : run_mm, run_gauss

    integer(kind=4) :: i, n, m, k, mode
    integer(kind=4) :: l, s
    character(len=10) :: buffer

    call get_command_argument(1, buffer, l, s)
    read(buffer,*) n

    call get_command_argument(2, buffer, l, s)
    read(buffer,*) m

    call get_command_argument(3, buffer, l, s)
    read(buffer,*) k

    call get_command_argument(4, buffer, l, s)
    read(buffer,*) mode

    do i=n, m, k
    select case (mode)
        case (0)
            write (*, *) "mm"
            call run_mm(i)
        case (1)
            write (*, *) "gauss"
            call run_gauss(i)
        case default
            write (*, *) "incorrect mode"
    end select
    enddo

end program main