!*
!--- Copyright by --- XShinHe <1500011805@pku.edu.cn>
!------- Date 2018. 10
!--- Acknowledgement to Liu group, of PKU
!*

module map
use base
use linalgebra
use MES_Models
implicit none
    integer, private :: vdim
    integer, private :: cdim
    integer, private :: rdim !-- or said, freedom of nuclues
    integer, private :: nbead
    integer :: map_ipop
    integer :: map_istat
    real(8) :: map_angle = 0.5 
    real(8), dimension(:,:,:), allocatable :: map_c
    real(8), dimension(:,:), allocatable :: map_pop

contains
    subroutine init_map(nb, nc, ipop, istat)
    implicit none
        integer, intent(in) :: nb, nc
        integer, intent(in), optional :: ipop, istat
        if(present(ipop)) then
            map_ipop = ipop
        else
            map_ipop = 0
        endif
        if(present(istat)) then
            map_istat = istat
        else
            map_istat = 1
        endif
        
        rdim = getrdim()
        vdim = getvdim()
        nbead = nb
        cdim  = nc
        allocate(map_c(nbead, vdim, cdim))
        allocate(map_pop(nbead, vdim ))
        call init_stat()
    end subroutine init_map

    subroutine init_stat()
    implicit none
        map_c = 0.0
        if(cdim .eq. 2) then
            map_c(:,map_istat,1) = dcos(map_angle)
            map_c(:,map_istat,2) = dsin(map_angle)
        else if(cdim .eq. 4) then
            select case (map_ipop)
                case (0) !-- squared
                    map_c(:,map_istat,1) = dcos(map_angle/2.0)
                    map_c(:,map_istat,2) = dsin(map_angle/2.0)
                    map_c(:,map_istat,3) = dcos(map_angle/2.0)
                    map_c(:,map_istat,4) = dsin(map_angle/2.0)
                case (1) !-- grouped
                    map_c(:,map_istat,1) = dcos(map_angle)/2.0
                    map_c(:,map_istat,2) = dsin(map_angle)/2.0
                    map_c(:,map_istat,3) = dcos(map_angle)/2.0
                    map_c(:,map_istat,4) = dsin(map_angle)/2.0
                case (2) !-- crossing
                    map_c(:,map_istat,1) = dcos(map_angle)
                    map_c(:,map_istat,2) = dsin(map_angle)
                    map_c(:,map_istat,3) = dcos(map_angle)
                    map_c(:,map_istat,4) = dsin(map_angle)
                case default
                    stop "wrong model"
            endselect
        else
            stop "wrong model"
        endif
    end subroutine init_stat
    
    subroutine map_update_x(Rs, dt)
        real(8), intent(in) :: dt
        real(8), dimension(rdim, nbead), intent(in) :: Rs
        real(8), dimension(vdim, vdim) :: V, Vd, Vnd
        integer :: i,k
        do i=1,nbead
            call Mod_getV(V, Vd, Vnd, Rs(:,i), iflag=0)
            do k=1, vdim
                select case (cdim)
                    case (2)
                        map_c(i,:,1) = map_c(i,:,1) + matmul(V,map_c(i,:,2)) *dt
                    case (4)
                        map_c(i,:,1) = map_c(i,:,1) + matmul(V,map_c(i,:,2)) *dt
                        map_c(i,:,3) = map_c(i,:,3) + matmul(V,map_c(i,:,4)) *dt
                    case default
                        stop "wrong model"
                endselect
            enddo
        enddo
    end subroutine map_update_x
    
    subroutine map_update_p(Rs, dt)
        real(8), intent(in) :: dt
        real(8), dimension(rdim, nbead), intent(in) :: Rs
        real(8), dimension(vdim, vdim) :: V, Vd, Vnd
        integer :: i,k
        do i=1,nbead
            call Mod_getV(V, Vd, Vnd, Rs(:,i), iflag=0)
            do k=1, vdim
                select case (cdim)
                    case (2)
                        map_c(i,:,2) = map_c(i,:,2) - matmul(V,map_c(i,:,1)) *dt
                    case (4)
                        map_c(i,:,2) = map_c(i,:,2) - matmul(V,map_c(i,:,1)) *dt
                        map_c(i,:,4) = map_c(i,:,4) - matmul(V,map_c(i,:,3)) *dt
                    case default
                        stop "wrong model"
                endselect
            enddo
        enddo
    end subroutine map_update_p
    
    !-- QUATERNION SPACE
    subroutine map_stat_pop()
    implicit none
        integer :: i,k
        select case(cdim)
            case (2)
                do i=1, nbead
                    do k=1, vdim
                        map_pop(i,k) = sum( map_c(i,k,:)**2 )
                    enddo
                enddo
            case (4)
                select case(map_ipop)
                    case (0)
                        do i=1, nbead
                            do k=1, vdim
                                map_pop(i,k) = sum( map_c(i,k,:)**2 )
                            enddo
                        enddo
                    case (1)
                        do i=1, nbead
                            do k=1, vdim
                                map_pop(i,k) = ( map_c(i,k,1)+map_c(i,k,3) )**2 + ( map_c(i,k,2)+map_c(i,k,4) )**2
                            enddo
                        enddo
                    case (2)
                        do i=1, nbead
                            do k=1, vdim
                                map_pop(i,k) = map_c(i,k,1) * map_c(i,k,3) + map_c(i,k,2) * map_c(i,k,4) 
                            enddo
                        enddo
                    case default
                        stop "wrong model"
                endselect
            case default
                stop "wrong model"
        endselect
    end subroutine map_stat_pop
    
    subroutine map_fx(fx, Rs)
    implicit none
        real(8), dimension(rdim, nbead), intent(in) :: Rs
        real(8), dimension(rdim, nbead), intent(out) :: fx
        real(8), dimension(rdim, vdim, vdim) :: dV, dVd, dVnd
        real(8), dimension(vdim) :: c1, c2, c3, c4
        integer :: i,j,k
        do i=1, nbead
            call Mod_getdV(dV, dVd, dVnd, Rs(:,i), iflag=0)
            do j=1, rdim
                select case (cdim)
                    case (2)
                        c1 = map_c(i,:,1)
                        c2 = map_c(i,:,2)
                        fx(j,i) = dot_product(c1, matmul( dV(j,:,:),c1) ) + dot_product(c2, matmul( dV(j,:,:), c2))
                    case (4)
                        select case (map_ipop)
                            case (0)
                                c1 = map_c(i,:,1)
                                c2 = map_c(i,:,2)
                                c3 = map_c(i,:,3)
                                c4 = map_c(i,:,4)
                                fx(j,i) = dot_product(c1, matmul( dV(j,:,:),c1) ) + dot_product(c2, matmul( dV(j,:,:), c2)) + &
                                          dot_product(c3, matmul( dV(j,:,:),c3) ) + dot_product(c4, matmul( dV(j,:,:), c4) )
                            case (1)
                                c1 = map_c(i,:,1)+map_c(i,:,3)
                                c2 = map_c(i,:,2)+map_c(i,:,4)
                                fx(j,i) = dot_product(c1, matmul( dV(j,:,:),c1) ) + dot_product(c2, matmul( dV(j,:,:), c2))
                            case (2)
                                c1 = map_c(i,:,1)
                                c2 = map_c(i,:,2)
                                c3 = map_c(i,:,3)
                                c4 = map_c(i,:,4)
                                fx(j,i) = dot_product(c1, matmul( dV(j,:,:),c3) ) + dot_product(c2, matmul( dV(j,:,:), c4))
                            case default
                                stop "wrong model"
                        endselect
                    case default
                        stop "wrong model"
                endselect
            enddo
        enddo
    end subroutine map_fx
    
    
end module map
