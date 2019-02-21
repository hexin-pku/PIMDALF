!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- mapping variables
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module map_core
implicit none
    integer, public :: map4_type = 0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Mapping variables type    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type mapvar
    integer :: sz   !-- size of diabatic basis
    integer :: cdim !-- c number dimension
    real(8), dimension(:,:), allocatable :: c !-- coefficients
    real(8), dimension(:), allocatable :: p   !-- populations
end type mapvar

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Allocate mapvar object
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine alloc_mapvar(imv, isz, icdim)
    type(mapvar), intent(inout) :: imv
    integer, intent(in) :: isz, icdim
    if(allocated(imv%c)) then
        stop "re-allocate warning"
    else
        imv%sz = isz
        imv%cdim = icdim
        allocate(imv%c(isz,icdim))
        allocate(imv%p(isz))
    endif
end subroutine alloc_mapvar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Deallocate mapvar object
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine del_mapvar(imv)
    type(mapvar), intent(inout) :: imv
    if(allocated(imv%c)) then
        deallocate(imv%c)
        deallocate(imv%p)
    endif
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- sampling mapvar object from random gaussian
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine sample_mapvar(imv)
use random
    type(mapvar), intent(inout) :: imv
    integer :: i
    if(imv%cdim .ne. 2) stop 'imv type error'
    imv%c(:,1) = gaussian_rns(imv%sz)
    imv%c(:,2) = gaussian_rns(imv%sz)
end subroutine sample_mapvar


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- initialize mapvar object for specific x-p
!-- Position-Momentum variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine init_mapvar_xp(imv, occ, x, p, iflag)
    type(mapvar), intent(inout) :: imv
    integer, intent(in) :: occ
    real(8), intent(in) :: x,p
    integer, intent(in), optional :: iflag
    
    if(present(iflag)) then
        if(iflag.eq. 0 ) imv%c = 0.0
    endif
    
    if(imv%cdim .eq. 2) then
        imv%c(occ,1) = x
        imv%c(occ,2) = p
    else
        stop "wrong cdim"
    endif
end subroutine init_mapvar_xp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- initialize mapvar object for specific a-a
!-- Action-Angle variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine init_mapvar_aa(imv, occ, act, ang, iflag)
    type(mapvar), intent(inout) :: imv
    integer, intent(in) :: occ
    real(8), intent(in) :: act, ang
    integer, intent(in), optional :: iflag
    real(8) :: sqrtn
    
    if(present(iflag)) then
        if(iflag.eq. 0 ) imv%c = 0.0
    endif
    
    sqrtn = dsqrt(act)
    if(imv%cdim .eq. 2) then
        imv%c(occ,1) = dcos(ang)*sqrtn
        imv%c(occ,2) = dsin(ang)*sqrtn
    else if(imv%cdim .eq. 4) then
        select case (map4_type)
            case (0) !-- squared, default
                imv%c(occ,1) = dcos(ang/2.0)*sqrtn
                imv%c(occ,2) = dsin(ang/2.0)*sqrtn
                imv%c(occ,3) = dcos(ang/2.0)*sqrtn
                imv%c(occ,4) = dsin(ang/2.0)*sqrtn
            case (1) !-- grouped
                imv%c(occ,1) = dcos(ang)/2.0 * sqrtn
                imv%c(occ,2) = dsin(ang)/2.0 * sqrtn
                imv%c(occ,3) = dcos(ang)/2.0 * sqrtn
                imv%c(occ,4) = dsin(ang)/2.0 * sqrtn
            case (2) !-- crossing
                imv%c(occ,1) = dcos(ang)*sqrtn
                imv%c(occ,2) = dsin(ang)*sqrtn
                imv%c(occ,3) = dcos(ang)*sqrtn
                imv%c(occ,4) = dsin(ang)*sqrtn
            case default
                stop "wrong map4_type"
        endselect
    else
        stop "wrong cdim"
    endif
end subroutine init_mapvar_aa

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Count populations
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine map_population(imv)
    type(mapvar), intent(inout) :: imv
    integer :: k
    
    if(imv%cdim .eq. 2) then
        do k=1, imv%sz
            imv%p(k) = imv%c(k,1)**2 + imv%c(k,2)**2
        enddo
    elseif(imv%cdim .eq. 4) then
        select case(map4_type)
            case (0)
                do k=1, imv%sz
                    imv%p(k) = imv%c(k,1)**2 + imv%c(k,2)**2 &
                    + imv%c(k,3)**2 + imv%c(k,4)**2
                enddo
            case (1)
                do k=1, imv%sz
                    imv%p(k) = ( imv%c(k,1)+imv%c(k,3) )**2 &
                    + ( imv%c(k,2)+imv%c(k,4) )**2
                enddo
            case (2)
                do k=1, imv%sz
                    imv%p(k) = imv%c(k,1) * imv%c(k,3) &
                    + imv%c(k,2) * imv%c(k,4) 
                enddo
            case default
                stop "wrong map4_type"
        endselect
    else
        stop "wrong cdim"
    endif
end subroutine map_population

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- propagators of map variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine map_propagate_x(imv, V, dt)
    type(mapvar), intent(inout) :: imv
    real(8), dimension(:,:), intent(in) :: V
    real(8), intent(in) :: dt
    if(imv%sz .ne. size(V,dim=1)) stop "size error"
    select case (imv%cdim)
        case (2)
            imv%c(:,1) = imv%c(:,1) + matmul(V,imv%c(:,2)) *dt
        case (4)
            imv%c(:,1) = imv%c(:,1) + matmul(V,imv%c(:,2)) *dt
            imv%c(:,4) = imv%c(:,4) - matmul(V,imv%c(:,3)) *dt
        case default
            stop "wrong model"
    endselect
end subroutine map_propagate_x

subroutine map_propagate_p(imv, V, dt)
    type(mapvar), intent(inout) :: imv
    real(8), dimension(:,:), intent(in) :: V
    real(8), intent(in) :: dt
    if(imv%sz .ne. size(V,dim=1)) stop "size error"
    select case (imv%cdim)
        case (2)
            imv%c(:,2) = imv%c(:,2) - matmul(V,imv%c(:,1)) *dt
        case (4)
            imv%c(:,2) = imv%c(:,2) - matmul(V,imv%c(:,1)) *dt
            imv%c(:,3) = imv%c(:,3) + matmul(V,imv%c(:,4)) *dt
        case default
            stop "wrong model"
    endselect
end subroutine map_propagate_p

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- calucate mean-like force
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine map_calculate_fx(fx, imv, dV, zpe)
implicit none
    type(mapvar), intent(in) :: imv
    real(8), dimension(:, :, :), intent(in) :: dV
    real(8), dimension(:), intent(inout) :: fx
    real(8), intent(in) :: zpe
    real(8), dimension(imv%sz) :: c1, c2
    integer :: j,k

    fx = 0
    
    if(imv%cdim .eq. 2) then
        do j=1,size(fx)
            fx(j) = dot_product(imv%c(:,1), matmul( dV(j,:,:), imv%c(:,1)) ) &
            + dot_product(imv%c(:,2), matmul( dV(j,:,:), imv%c(:,2)))
            !-- couting zpe effect
            do k=1, imv%sz
                fx(j) = fx(j) - zpe*dV(j,k,k)
            enddo
        enddo
    elseif(imv%cdim .eq. 4) then
        select case (map4_type)
            case (0)
                do j=1,size(fx)
                    fx(j) = dot_product(imv%c(:,1), matmul( dV(j,:,:), imv%c(:,1)) ) &
                    + dot_product(imv%c(:,2), matmul( dV(j,:,:), imv%c(:,2)) ) &
                    + dot_product(imv%c(:,3), matmul( dV(j,:,:), imv%c(:,3)) ) &
                    + dot_product(imv%c(:,4), matmul( dV(j,:,:), imv%c(:,4)) )
                    do k=1, imv%sz
                        fx(j) = fx(j) - zpe*dV(j,k,k)
                    enddo
                enddo
            case (1)
                do j=1,size(fx)
                    c1 = imv%c(:,1)+imv%c(:,3)
                    c2 = imv%c(:,2)+imv%c(:,4)
                    fx(j) = dot_product(c1, matmul( dV(j,:,:),c1) ) &
                    + dot_product(c2, matmul( dV(j,:,:), c2))
                    do k=1, imv%sz
                        fx(j) = fx(j) - zpe*dV(j,k,k)
                    enddo
                enddo
            case (2)
                do j=1,size(fx)
                    fx(j) = dot_product(imv%c(:,1), matmul( dV(j,:,:), imv%c(:,3)) ) &
                    + dot_product(imv%c(:,2), matmul( dV(j,:,:), imv%c(:,4)))
                    do k=1, imv%sz
                        fx(j) = fx(j) - zpe*dV(j,k,k)
                    enddo
                enddo
            case default
                stop "map4_type error"
        endselect
    else
        stop "cdim error"
    endif
end subroutine map_calculate_fx

end module map_core
