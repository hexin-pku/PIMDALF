!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- PIMD Simulation procedure
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Harmonic Oscillitor Model
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module model
use const
implicit none
    integer, parameter :: ndof = 1
    real(dp), dimension(ndof), parameter :: M = 1
contains

subroutine model_V0(v,x)
    real(dp), intent(inout) :: v,x(:)
    v = 0.5_8 * sum(x*x)
end subroutine model_V0

subroutine model_V1(v,dv,x,iflag)
    real(dp), intent(inout) :: v, dv(:), x(:)
    integer, optional :: iflag
    v = 0.5_8 * sum(x*x)
    if( present(iflag) .and. iflag < 1) return
    dv = x
end subroutine model_V1

subroutine model_V2(v,dv,ddv,x,iflag)
    real(dp), intent(inout) :: v, dv(:), ddv(:,:), x(:)
    integer, optional :: iflag
    v = 0.5_8 * sum(x*x)
    if( present(iflag) .and. iflag < 1) return
    dv = x
    if( present(iflag) .and. iflag < 2) return
    ddv = 1.0_8
end subroutine model_V2

end module model



