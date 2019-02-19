module linalgebra
use const
use random
implicit none
contains
    function lin_ones(rows, cols) result(M)
    implicit none
        integer, intent(in) :: rows, cols
        real(8), dimension(rows,cols) :: M
        !-- allocate(M(rows,cols))
        M = 1.0
    end function lin_ones

    function lin_zeros(rows, cols) result(M)
    implicit none
        integer, intent(in) :: rows, cols
        real(8), dimension(rows,cols) :: M
        !allocate(M(rows,cols))
        M = 0.0
    end function lin_zeros
    
    function lin_eye(msize, ival) result(M)
    implicit none
        integer, intent(in) :: msize
        real(8), intent(in), optional :: ival
        real(8), dimension(msize,msize) :: M
        real(8) :: val
        integer :: i
        if(present(ival)) then
            val = ival
        else
            val = 1.0
        endif
        M=0.0
        do i=1,msize
            M(i,i) = val
        enddo
        !-- print *, "test on la : eye", M
    end function lin_eye
    
    function lin_dot(X, Y) result(M)
    implicit none
        real(8), dimension(:,:), intent(in) :: X,Y
        real(8), dimension(size(X,dim=1), size(Y,dim=2)) :: M
        if( size(X,dim=2) .ne. size(Y,dim=1) )  then
            stop "matmul error in linagebra"
        else
            M = matmul(X,Y)
        endif
    end function lin_dot
    
    function lin_outer(X, Y) result(M)
    implicit none
        real(8), dimension(:), intent(in) :: X,Y
        real(8), dimension(size(X),size(Y)) :: M
        integer :: i,j
        if( size(X).ne. size(Y)) then
            stop "matmul error in linagebra"
        else
            do i=1,size(X)
                do j=1,size(Y)
                    M(i,j)=X(i)*Y(j)
                enddo
            enddo
        endif
    end function lin_outer
    
    function lin_trace(X) result(Z)
    implicit none
        real(8), dimension(:,:), intent(in) :: X
        real(8) :: Z
        integer :: i
        if( size(X,dim=1) .ne. size(X,dim=2) )  then
            stop "matmul error in linagebra"
        else
            Z = 0.0
            do i=1,size(X,dim=1)
                Z = Z + X(i,i)
            enddo
        endif
    end function lin_trace
    
    function lin_eps(X, eps) result(M)
    implicit none
        real(8), dimension(:,:) :: X
        real(8), intent(in) :: eps
        real(8), dimension(size(X,dim=1),size(X,dim=2)) :: M
        integer :: i,j
        M=0.0
        do i=1,size(X,dim=1)
            do j=1,size(X,dim=2)
                if(abs(X(i,j)) > eps) M(i,j) = X(i,j)
            enddo
        enddo
        !-- print *, "test on la : eye", M
    end function lin_eps
    
    function lin_diag(X) result(M)
    implicit none
        real(8), dimension(:,:), intent(in) :: X
        real(8), dimension(size(X,dim=1),size(X,dim=2)) :: M
        integer :: i
        if( size(X,dim=1) .ne. size(X,dim=2) )  then
            stop "matmul error in linagebra"
        else
            M = 0.0
            do i=1,size(X,dim=1)
                M(i,i) = X(i,i)
            enddo
        endif
    end function lin_diag
    
    function lin_darray(X) result(M)
    implicit none
        real(8), dimension(:), intent(in) :: X
        real(8), dimension(size(X,dim=1),size(X,dim=1)) :: M
        integer :: i
        M = 0.0
        do i=1,size(X,dim=1)
            M(i,i) = X(i)
        enddo
    end function lin_darray
    
    function lin_offdiag(X) result(M)
    implicit none
        real(8), dimension(:,:), intent(in) :: X
        real(8), dimension(size(X,dim=1),size(X,dim=2)) :: M
        integer :: i,j
        if( size(X,dim=1) .ne. size(X,dim=2) )  then
            stop "matmul error in linagebra"
        else
            M = 0.0
            do i=1,size(X,dim=1)
                do j=1,size(X,dim=2)
                    if(i.eq. j) cycle
                    M(i,j) = X(i,j)
                enddo
            enddo
        endif
    end function lin_offdiag
    
    function lin_expdiag(X) result(M)
    implicit none
        real(8), dimension(:,:), intent(in) :: X
        real(8), dimension(size(X,dim=1),size(X,dim=2)) :: M
        integer :: i
        if( size(X,dim=1) .ne. size(X,dim=2) )  then
            stop "matmul error in linagebra"
        else
            M = 0.0
            do i=1,size(X,dim=1)
                M(i,i) = dexp(X(i,i))
            enddo
        endif
    end function lin_expdiag
    
    function lin_cumul(P) result(M)
    implicit none
        real(8), dimension(:), intent(in) ::  P
        real(8), dimension(size(P)) :: M
        integer ::  k
        M = (/(sum(P(1:k)),k=1,size(P))/)
    end function lin_cumul
    
    function lin_op_power(X, n) result(M)
    implicit none
        real(8), dimension(:,:), intent(in) :: X
        integer, intent(in) :: n
        integer :: k
        real(8), dimension(size(X,dim=1),size(X,dim=2)) :: M
        M = lin_eye(size(X,dim=1))
        do k=1,n
            M = matmul(X,M)
        enddo
    end function lin_op_power
    
    function lin_powerdiag(X, n) result(M)
    implicit none
        real(8), dimension(:,:), intent(in) :: X
        integer, intent(in) :: n
        integer :: k
        real(8), dimension(size(X,dim=1),size(X,dim=2)) :: M
        M = 0.0
        do k=1,size(X,dim=1)
            M(k,k) = X(k,k) ** n
        enddo
    end function lin_powerdiag
    
    function lin_op_exp(X, prec) result(M)
        real(8), dimension(:,:), intent(in) :: X
        real(8), dimension(size(X,dim=1),size(X,dim=2)) :: M, P, Y
        integer, intent(in), optional :: prec
        integer :: k, n
        if(present(prec) .and. prec > 50) then
            n = prec
        else
            n = 100
        endif
        M = lin_eye(size(X,dim=1))
        P = lin_eye(size(X,dim=1))
        Y = X/real(n,kind=8)
        do k=1,100
            P = matmul(Y,P)/real(k,kind=8)
            M = M + P
        enddo
        M = lin_op_power(M, 100)
    end function lin_op_exp
    
    function lin_op_log(X) result(M)
        real(8), dimension(:,:), intent(in) :: X
        real(8), dimension(size(X,dim=1),size(X,dim=2)) :: M, P
        integer :: k
        M = lin_eye(size(X,dim=1))
        P = - lin_eye(size(X,dim=1))
        do k=1,100
            P = - matmul(X,P)
            M = M + P/real(k,kind=8)
        enddo
    end function lin_op_log
    
    recursive function lin_op_cos(X) result(M)
        real(8), dimension(:,:), intent(in) :: X
        real(8), dimension(size(X,dim=1),size(X,dim=2)) :: M, P, Y
        integer :: k
        if(maxval(abs(X)) > 0.05) then
            Y = lin_op_cos(X/real(3.0_8))
            M = 4.0_8 * lin_op_power(Y,3) - 3.0_8 * Y
        else
            M = lin_eye(size(X,dim=1))
            P = -lin_eye(size(X,dim=1))
            Y = matmul(X,X)
            do k=1,100
                P = -matmul(Y,P)/real(k*(k+1),kind=8)
                M = M + P
            enddo
        endif
    end function lin_op_cos
    
    recursive function lin_op_sin(X) result(M)
        real(8), dimension(:,:), intent(in) :: X
        real(8), dimension(size(X,dim=1),size(X,dim=2)) :: M, P, Y
        integer :: k
        if(maxval(abs(X)) > 0.05) then
            Y = lin_op_sin(X/real(3.0_8))
            M = -4.0_8 * lin_op_power(Y,3) + 3.0_8 * Y
        else
            M = X
            P = X
            Y = matmul(X,X)
            do k=2,100
                P = -matmul(Y,P)/real(k*(k+1),kind=8)
                M = M + P
            enddo
        endif
    end function lin_op_sin
    
    function lin_op_inv(A) result(Ainv)
        real(8), dimension(:,:), intent(in) :: A
        real(8), dimension(size(A,1),size(A,2)) :: Ainv

        real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info
        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)
        if (info /= 0) then
            stop 'Matrix is numerically singular!'
        end if
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)

        if (info /= 0) then
           stop 'Matrix inversion failed!'
        end if
    end function lin_op_inv
    
    function lin_random_norm(rows, cols) result(M)
    implicit none
        integer, intent(in) :: rows, cols
        real(8), dimension(rows,cols) :: M
        integer :: i,j
        real(8) :: irand_norm
        do i=1,rows
            do j=1,cols
                call random_normal(irand_norm)
                M(i,j) = irand_norm
            enddo
        enddo
    end function lin_random_norm
    
    subroutine lin_ev(Es, Vs, A)
        real(8), dimension(:,:), intent(in) :: A
        real(8), dimension(size(A,1),size(A,2)), intent(out) :: Vs
        real(8), dimension(size(A,1)), intent(out) :: Es
        integer :: LDA, n, lwork, info
        integer, parameter :: lwmax = 1000
        real(8) :: work( lwmax )
        ! External procedures defined in LAPACK
        external DSYEV
        ! Store A
        Vs = A
        n = size(A,1)
        LDA = n
        
        ! Query the optimal workspace.
        lwork = -1
        CALL DSYEV( 'Vectors', 'Upper', n, Vs, LDA, Es, work, lwork, info )
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
        
        ! Solve eigenproblem.
        CALL DSYEV( 'Vectors', 'Upper', n, vs, LDA, Es, work, lwork, info )

        ! Check for convergence.
        if( info .GT. 0 ) then
            print *, 'The algorithm failed to compute eigenvalues.'
            stop
        endif
    end subroutine lin_ev

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- subroutines
    
end module linalgebra










