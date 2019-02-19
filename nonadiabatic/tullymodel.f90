!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Tully models
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module tullymodel
implicit none
private
    real(8) :: V0, V1, E0, a, b
    integer, parameter :: ndof = 1
    integer, parameter :: vdim = 2
    real(8), dimension(ndof), parameter :: M = 2000
public :: M, ndof, vdim
public :: tully_read, tully_v1, tully_v2, tully_v3, tully_dv1, tully_dv2, tully_dv3

contains
subroutine tully_read(parmsfile)
    character(*), intent(in) :: parmsfile
    logical :: my_exist
    inquire(file=parmsfile, exist=my_exist ) 
    if(.not. my_exist) stop "parms lossing"
    open(unit=20, file=parmsfile, status="old")
    read(20,*)
    read(20,*) V0
    read(20,*)
    read(20,*) V1
    read(20,*)
    read(20,*) E0
    read(20,*)
    read(20,*) a
    read(20,*)
    read(20,*) b
    close(unit=20)
end subroutine tully_read

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Tully model I : avoid crossing
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine tully_v1(V, R)
    real(8), dimension(2,2), intent(out) :: V
    real(8), dimension(1), intent(in) :: R
    
    if(R(1)>0) then
        V(1,1) = V0 * ( 1.0 - dexp(-a*R(1)) )
    else
        V(1,1) = -V0 * ( 1.0 - dexp(a*R(1)) )
    endif
    V(2,2) = -V(1,1)
    V(1,2) = V1 * dexp(- b * R(1)*R(1))
    V(2,1) = V(1,2)
    return
end subroutine tully_v1


subroutine tully_dv1(dV, R)
    real(8), dimension(1,2,2), intent(out) :: dV
    real(8), dimension(1), intent(in) :: R
    
    if(R(1)>0) then
        dV(1,1,1) = V0 * a * dexp(-a*R(1))
    else
        dV(1,1,1) = V0 * a * dexp(a*R(1))
    endif
    dV(1,2,2) = -dV(1,1,1)
    dV(1,1,2) = -2.0*b*R(1)* V1 * dexp(- b * R(1)*R(1))
    dV(1,2,1) = dV(1,1,2)
    return
end subroutine tully_dv1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Tully model II : double crossing
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine tully_v2(V, R)
    real(8), dimension(2,2), intent(out) :: V
    real(8), dimension(1), intent(in) :: R
    
    V(1,1) = 0
    V(2,2) = -V0*dexp(-a*R(1)*R(1)) + E0
    V(1,2) = V1 * dexp(- b * R(1)*R(1))
    V(2,1) = V(1,2)
    return
end subroutine tully_v2

subroutine tully_dv2(dV, R)
    real(8), dimension(1,2,2), intent(out) :: dV
    real(8), dimension(1), intent(in) :: R
    
    dV(1,1,1) = 0
    dV(1,2,2) = 2*a*R(1) * V0 * dexp(-a*R(1)*R(1))
    dV(1,1,2) = -2.0*b*R(1)* V1 * dexp(- b * R(1)*R(1))
    dV(1,2,1) = dV(1,1,2)
    return
end subroutine tully_dv2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Tully model III : avoid crossing
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine tully_v3(V, R)
    real(8), dimension(2,2), intent(out) :: V
    real(8), dimension(1), intent(in) :: R
    
    V(1,1) = -V0
    V(2,2) = V0
    if(R(1) < 0) then
        V(1,2) = V1 * dexp(b * R(1))
    else
        V(1,2) = V1 * (2.0 - dexp(-b * R(1)) )
    endif
    V(2,1) = V(1,2)
    return
end subroutine tully_v3

subroutine tully_dv3(dV, R)
    real(8), dimension(1,2,2), intent(out) :: dV
    real(8), dimension(1), intent(in) :: R
    
    dV(1,1,1) = 0
    dV(1,2,2) = 0
    if(R(1) < 0) then
        dV(1,1,2) = V1 * b * dexp(b * R(1))
    else
        dV(1,1,2) = V1 * b * dexp(-b * R(1))
    endif
    dV(1,2,1) = dV(1,1,2)
    return
end subroutine tully_dv3

end module tullymodel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Test
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




