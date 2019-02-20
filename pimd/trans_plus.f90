!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- PIMD Simulation procedure
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- transformation code
!-- primitive, staging & normal-mode
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module trans_plus
use const
implicit none
private
    public :: m2mb, x2ks, ks2x, fx2fks, set_trans, del_trans
    
    !-- orthogonal matrix
    real(dp), dimension(:,:), allocatable :: OM
    !-- transfrom methods
    !-- val = 0 : primtive
    !-- val = 1 : staging
    !-- val = 2 : centroid normal mode
    integer :: trans_mod = 1 !-- by default
    
contains

subroutine set_trans(imod, P)
    integer, intent(in) :: imod, P
    trans_mod = imod
    if( trans_mod .eq. 2 ) then
        if( mod(P,2) .ne. 0) stop 'P should even for NM-trans'
        allocate( OM(P,P) )
    endif
end subroutine set_trans


subroutine del_trans()
    if(trans_mod .eq. 2) deallocate(OM)
end subroutine del_trans


subroutine m2mb(mb, m, P)
implicit none
    integer, intent(in) :: P
    real(dp), dimension(P), intent(in) :: m
    real(dp), dimension(P), intent(out) :: mb
    integer :: i,j
    
    select case (trans_mod)
        case (0)
            mb = m
        case (1)
            mb(1) = m(1)
            if(P.eq.1) return
            do i=2, P
                mb(i) = m(i) * real(i)/real(i-1)
            enddo
        case (2)
            mb(1) = real(P,dp) * m(1)
            OM(1,:) = 1.0_dp/ SQRT( real(P,dp) )
            do i=2,P/2
                mb(2*i-2) = 2.0_dp*( 1 - DCOS(twopi * real(i-1)/real(P,dp) ) )*real(P,dp) * m(2*i-2)
                mb(2*i-1) = 2.0_dp*( 1 - DCOS(twopi * real(i-1)/real(P,dp) ) )*real(P,dp) * m(2*i-1)
                do j=1,P
                    OM(2*i-2,j) = DCOS( twopi* real(i-1)*real(j-1)/real(P,dp) ) * SQRT(2.0_dp/real(P,dp))
                    OM(2*i-1,j) =-DSIN( twopi* real(i-1)*real(j-1)/real(P,dp) ) * SQRT(2.0_dp/real(P,dp))
                end do
            end do
            mb(P) = 4.0_dp*real(P,dp) * m(P)
            OM(P,1::2) = 1.0_dp / SQRT( real(P,dp) )
            OM(P,2::2) =-1.0_dp / SQRT( real(P,dp) )
    endselect
    
end subroutine m2mb

subroutine x2ks(ks, x, P)
implicit none
    integer, intent(in) :: P
    real(dp), dimension(P), intent(in) :: x
    real(dp), dimension(P), intent(out) :: ks
    integer :: i,j
    
    select case (trans_mod)
        case (0)
            ks = x
        case (1)
            ks(1) = x(1)
            if(P.eq.1) return
            ks(P) = x(P) - x(1)
            do i=2, P-1
                ks(i) = x(i) - ( x(i+1) * real(i - 1) + x(1) ) / real(i)
            enddo
        case (2)
            do i=1,P
                ks(i) = 0.0_dp
                do j=1,P
                    ks(i) = ks(i) + x(j) * OM(i,j)
                end do
                ks(i) = ks(i) / SQRT(real(P,dp))
            end do
    endselect
    
end subroutine x2ks

subroutine ks2x(x, ks, P)
implicit none
    integer, intent(in) :: P
    real(dp), dimension(P), intent(in) :: ks
    real(dp), dimension(P), intent(out) :: x
    integer :: i,j
    select case (trans_mod)
        case (0)
            x = ks
        case (1)      
            x(1) = ks(1)
            if(P .eq.1) return
            x(P) = ks(P) + ks(1)
            do i=P-1,2,-1
                x(i) = ks(i) + ( real( i - 1 ) * x(i+1) + ks(1) ) / real(i)
            enddo
        case (2)
            do i=1,P
                x(i) = 0.0_dp
                do j=1,P
                    x(i) = x(i) + ks(j) * OM(j,i)
                end do
                x(i) = x(i) * SQRT(real(P,dp))
            end do
    endselect
end subroutine ks2x

subroutine fx2fks(fks, fx, P)
implicit none
    integer, intent(in) :: P
    real(dp), dimension(P), intent(in) :: fx
    real(dp), dimension(P), intent(out) :: fks
    integer :: i,j
    select case (trans_mod)
        case (0)
            fks = fx
        case (1)
            fks(1) = 0.0_dp
            do i=1,P
                fks(1) = fks(1) + fx(i)
            enddo
            
            do i=2,P,1
                fks(i) = fx(i) + fks(i-1) * real(i-2)/real(i-1)
            enddo
            fks = fks / real(P,dp)
        case (2)
            do i=1,P
                fks(i) = 0.0_dp
                do j=1,P
                    fks(i) = fks(i) + fx(j) * OM(i,j)
                end do
                fks(i) = fks(i) * SQRT( real(P,dp) )
            end do
            fks = fks / SQRT(real(P,dp))
    endselect
end subroutine fx2fks

end module trans_plus


