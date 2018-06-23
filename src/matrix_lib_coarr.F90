module matrix_lib_coarr
contains

subroutine mm(first, second, multiply, n, m, o)
    implicit none
    integer(kind= 4), intent(in) :: n, m, o
    real(kind= 8),intent(in), codimension[num_images(), *], dimension(n,m) :: first
    real(kind= 8),intent(in), codimension[num_images(), *] :: second(m ,o)
    real(kind= 8),intent(out), codimension[num_images(), *] :: multiply(n,o)
    integer, codimension[:], allocatable :: obeg, oend
    integer(kind= 4) :: i, j, k, jj, kk, ichunk, p

    allocate(obeg[*])
    allocate(oend[*])

    ichunk = 180

    p = ceiling(real(o)/real(num_images()))
    obeg = (p*(this_image()-1)+1)
    oend = (p*this_image())
    if (this_image() .EQ. num_images()) then
        oend = o
    endif
    multiply(:, obeg:oend) = 0

    do jj = obeg, oend, ichunk
        do kk = 1, m, ichunk
            do j = jj, min(jj + ichunk - 1, oend)
                do k = kk, min(kk + ichunk - 1, m)
                    do i = 1,n
                        multiply(i, j) = multiply(i, j) + first(i, k) * second(k, j)
                    end do
                end do
            end do
        end do
    end do

    syncall()

    if (this_image() .EQ. 1) then
        write (*, *) ichunk
        do i=2, num_images()
            multiply(:, obeg[i]:oend[i]) = multiply(:, obeg[i]:oend[i])[i, 1]
        enddo
    endif
end subroutine

subroutine gauss_elimination(A, X, N)
    integer(kind=4), intent(in) :: N
    real(kind= 8), intent(inout) :: A(N, N), X(N)
    integer(kind= 4) :: i, j
    real(kind= 8) :: C

    do I = 1,N
        do J = 1,N
            if (I .NE. J) then
                C = A(I, J) / A(I, I)
                if (C .NE. 0) then
                    A(:, J) = A(:, J) - C * A(:, I)
                    X(J) = X(J) - C*X(I)
                endif
                if (A(I, I) .NE. 1) then
                    X(I) = X(I) / A(I, I)
                    A(:, I) = A(:, I) / A(I, I)
                endif
            endif
        enddo
    enddo
    X = X * real(N * (N + 1) * (-1))
end subroutine

end module matrix_lib_coarr