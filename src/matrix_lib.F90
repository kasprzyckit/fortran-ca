module matrix_lib        
contains


!> Performs matrix multiplication.
!! This subroutine utilises loop nest optimization. If supported,
!! put into the 'ichunk' variable value fitting to your CPU cache size.
!! 
!! @param first first matrix
!! @param second seocnd matrix
!! @param multiply result matrix
!! @param n row size od the first matrix, assert this is true: n = shape(first, 1) = shape(multiply, 1)
!! @param m column size od the first matrix, assert this is true: n = shape(first, 2) = shape(second, 1)
!! @param o column size od the second matrix, assert this is true: n = shape(second, 2) = shape(multiply, 2)
subroutine mm(first, second, multiply, n, m, o)
    implicit none
    integer(kind= 4), intent(in) :: n, m, o
    real(kind= 8),intent(in):: first(n,m)
    real(kind= 8),intent(in):: second(m ,o)
    real(kind= 8),intent(out):: multiply(n,o)
    integer(kind= 4) :: i, j, k, jj, kk, ichunk

    !f2py intent(in) :: first, second
    !f2py intent(out) :: multiply, status

    ichunk = 180 !< Adjust to your own cache size

    multiply = 0

    do jj = 1, o, ichunk
        do kk = 1, m, ichunk
            do j = jj, min(jj + ichunk - 1, o)
                do k = kk, min(kk + ichunk - 1, m)
                    do i = 1,n
                        multiply(i, j) = multiply(i, j) + first(i, k) * second(k, j)
                    end do
                end do
            end do
        end do
    end do

end subroutine

!> Performs gaussian elimination
!! 
!! @param A matrix of coefficients
!! @param X result column
!! @param N size of the matrix
subroutine gauss_elimination(A, X, N)
    integer(kind=4), intent(in) :: N
    real(kind= 8), intent(inout) :: A(N, N), X(N)
    integer(kind= 4) :: i, j
    real(kind= 8) :: C

    !f2py intent(in, out) :: A, X

    do I = 1,N
        do J = 1,N
            if (I .NE. J) then
                if (A(I, J) .NE. 0) then
                    C = A(I, J) / A(I, I)
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

end module matrix_lib