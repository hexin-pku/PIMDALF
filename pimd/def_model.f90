!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- PIMD Simulation procedure
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module def_model
implicit none
    integer, parameter :: ndof = 1
    real(8), dimension(ndof), parameter :: M = 1
contains

elemental subroutine model_V(v,x)
    real(8), intent(inout) :: v,x
    v = 0.5_8 * x*x
end subroutine model_V

elemental subroutine model_dV(dv,x)
    real(8), intent(inout) :: dv,x
    dv = x
end subroutine model_dV

end module def_model
