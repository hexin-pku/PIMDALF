!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- PIMD Simulation procedure
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Basic model --- 1-D QPW model
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module model
use const
implicit none
    integer, parameter :: ndof = 1
    real(8), dimension(ndof), parameter :: M = 1
contains

subroutine model_V0(v,x)
    real(8), intent(inout) ::v, x(:)
    v = 0.25_8 * sum(x**4)
end subroutine model_V0

subroutine model_V1(v,dv,x,iflag)
    real(dp), intent(inout) :: v, dv(:), x(:)
    integer, optional :: iflag
    v = 0.25_8 * sum(x**4)
    if( present(iflag) .and. iflag < 1) return
    dv = x**3
end subroutine model_V1

subroutine model_V2(v,dv,ddv,x,iflag)
    real(dp), intent(inout) :: v, dv(:), ddv(:,:), x(:)
    integer, optional :: iflag
    integer :: i
    v = 0.25_8 * sum(x*x)
    if( present(iflag) .and. iflag < 1) return
    dv = x**3
    if( present(iflag) .and. iflag < 2) return
    ddv = 0.0_8
    do i=1,size(x)
        ddv(i,i) = 3*x(i)*x(i)
    enddo
end subroutine model_V2

end module model
