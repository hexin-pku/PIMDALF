!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- PIMD Simulation procedure
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- for staging transformation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module staging
implicit none
private
    public :: m2mb, x2ks, ks2x, fx2fks
contains

subroutine m2mb(mb, m, P)
implicit none
    integer, intent(in) :: P
    real(8), dimension(P), intent(in) :: m
    real(8), dimension(P), intent(out) :: mb
    integer :: i
    mb(1) = m(1)
    if(P.eq.1) return
    do i=2, P
        mb(i) = m(i) * real(i)/real(i-1)
    enddo
end subroutine m2mb

subroutine x2ks(ks, x, P)
implicit none
    integer, intent(in) :: P
    real(8), dimension(P), intent(in) :: x
    real(8), dimension(P), intent(out) :: ks
    integer :: i
    ks(1) = x(1)
    if(P.eq.1) return
    ks(P) = x(P) - x(1)
    do i=2, P-1
        ks(i) = x(i) - ( x(i+1) * real(i - 1) + x(1) ) / real(i)
    enddo
end subroutine x2ks

subroutine ks2x(x, ks, P)
implicit none
    integer, intent(in) :: P
    real(8), dimension(P), intent(in) :: ks
    real(8), dimension(P), intent(out) :: x
    integer :: i
    x(1) = ks(1)
    if(P .eq.1) return
    x(P) = ks(P) + ks(1)
    do i=P-1,2,-1
        x(i) = ks(i) + ( real( i - 1 ) * x(i+1) + ks(1) ) / real(i)
    enddo
end subroutine ks2x

subroutine fx2fks(fks, fx, P)
implicit none
    integer, intent(in) :: P
    real(8), dimension(P), intent(in) :: fx
    real(8), dimension(P), intent(out) :: fks
    integer :: i

    fks(1) = 0.0
    do i=1,P
        fks(1) = fks(1) + fx(i)
    enddo
    
    do i=2,P,1
        fks(i) = fx(i) + fks(i-1) * real(i-2)/real(i-1)
    enddo
    fks = fks /real(P)
end subroutine fx2fks

end module staging


