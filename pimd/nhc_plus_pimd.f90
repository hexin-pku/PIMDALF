!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- PIMD Simulation procedure
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Nose Hoover Chain Thermostat with RESPA-alg.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module nhc_plus_pimd
use const
implicit none
private
    integer :: ndof   = 0 ! for default
    integer :: nchain = 2 ! for default
    integer :: nrespa = 1 ! for default
    integer, parameter :: nsy = 7
    real(dp), dimension(nsy), parameter :: wsy = & !-- weight of SY
        (/0.784513610477560_dp,0.235573213359357_dp,-1.17767998417887_dp,&
         1.3151863206839063_dp,-1.17767998417887_dp,0.235573213359357_dp,&
         0.784513610477560_dp/)
    real(dp), dimension(nsy) :: wgt !-- weight of SY with RESPA
    real(dp), dimension(:,:,:), allocatable :: nhc_x, nhc_p, nhc_G
    real(dp), dimension(:), allocatable :: nhc_Q

    public :: set_NHC, del_NHC, thermo_NHC
contains

subroutine set_NHC(nchain_in, ndof_in, nrespa_in, gamma_in, T, P)
    integer, intent(in) :: nchain_in, ndof_in, nrespa_in, P
    !-- P: number of beads
    real(dp), intent(in) :: gamma_in, T 
    !-- gamma_in : optimal parameter
    !-- T : temperature
    
    nchain = nchain_in
    ndof = ndof_in
    nrespa = nrespa_in
    
    if( (ndof < 1) .or. ( nchain < 1 ) .or. (nrespa < 1) ) stop 'generate NHC error'

    allocate( nhc_x(nchain, ndof, P), nhc_p(nchain, ndof, P),&
              nhc_G(nchain, ndof, P), nhc_Q(nchain) )
    nhc_x = 0.0_dp; nhc_p = 0.0_dp
    nhc_G = 0.0_dp; nhc_Q = 0.0_dp

    !-- optimal gamma_in = 20*dtime
    !-- optimal nhc_Q
    nhc_Q = T * gamma_in**2
    nhc_Q(1) = ndof * T * gamma_in**2

    wgt = wsy / nrespa !-- weight of each composed-SY-RESPA time steps

end subroutine set_NHC

subroutine del_NHC()
    deallocate(nhc_x, nhc_p, nhc_G, nhc_Q)
end subroutine del_NHC

!-- introduce variable nhc_G
subroutine thermo_NHC(sys_p, sys_m, dt, T, P)
    real(dp), intent(inout) :: sys_p(:,:)
    real(dp), intent(in) :: sys_m(:,:)
    real(dp), intent(in) :: dt, T
    integer, intent(in) :: P
    integer :: h, i, j, k, s
    do i=1,ndof
        do j=1,P
            do k=1,nrespa
                do s=1,nsy
                    !-- update auxiliary variable nhc_G
                    nhc_G(1,i,j) = sys_p(i,j)**2/sys_m(i,j) - T
                    do h=2,nchain
                        nhc_G(h,i,j) = nhc_p(h-1,i,j)*nhc_p(h-1,i,j)/nhc_Q(h-1) - T
                    end do
                    !-- update nhc_p from tail to head
                    nhc_p(nchain,i,j) = nhc_p(nchain,i,j) &
                                + nhc_G(nchain,i,j) * wgt(s) * dt / 2
                    do h=nchain-1,1,-1
                        nhc_p(h,i,j) = nhc_p(h,i,j) &
                                * EXP( -nhc_p(h+1,i,j)/nhc_Q(h+1) * wgt(s) * dt / 4 )
                        nhc_p(h,i,j) = nhc_p(h,i,j) &
                                + nhc_G(h,i,j) * wgt(s) * dt / 2
                        nhc_p(h,i,j) = nhc_p(h,i,j) &
                                * EXP( -nhc_p(h+1,i,j)/nhc_Q(h+1) * wgt(s) * dt / 4 )
                    end do
                    !-- update nhc_x and system_p
                    do h=1,nchain
                        nhc_x(h,i,j) = nhc_x(h,i,j) &
                                + nhc_p(h,i,j) / nhc_Q(h)* wgt(s) * dt
                    end do
                    sys_p(i,j) = sys_p(i,j) &
                                * EXP( -nhc_p(1,i,j)/nhc_Q(1)*wgt(s) * dt )
                    print *,'i j k a p p1', i, j, k, s, sys_p(i,j), nhc_p(1,i,j)
                    !-- update nhc_p from head to tail
                    do h=1,nchain-1,1
                        nhc_p(h,i,j) = nhc_p(h,i,j) &
                                * EXP( -nhc_p(h+1,i,j)/nhc_Q(h+1) * wgt(s) * dt / 4 )
                        nhc_p(h,i,j) = nhc_p(h,i,j) &
                                + nhc_G(h,i,j) * wgt(s) * dt / 2
                        nhc_p(h,i,j) = nhc_p(h,i,j) &
                                * EXP( -nhc_p(h+1,i,j)/nhc_Q(h+1) * wgt(s) * dt / 4 )
                    end do
                    nhc_p(nchain,i,j) = nhc_p(nchain,i,j) + nhc_G(nchain,i,j) *wgt(s) * dt / 2
                end do
            end do
        end do
    end do
end subroutine thermo_NHC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Deprecated codes
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subroutine thermo_NHC1(sys_p, sys_m, dt, T, P)
!    real(dp), intent(inout) :: sys_p(:,:)
!    real(dp), intent(in) :: sys_m(:,:)
!    real(dp), intent(in) :: dt, T
!    integer, intent(in) :: P
!    integer :: h, i, j, k, s
!    do i=1,ndof
!        do j=1,P
!            do k=1,nrespa
!                do s=1,nsy
!                    nhc_p(nchain,i,j) = nhc_p(nchain,i,j) + ( nhc_p(nchain-1,i,j)*nhc_p(nchain-1,i,j)/nhc_Q(nchain-1) - T ) &
!                                    *(wgt(s) * dt / 2)
!                    do h=nchain-1,2,-1
!                        nhc_p(h,i,j) = nhc_p(h,i,j) * EXP( -nhc_p(h+1,i,j)/nhc_Q(h+1) * wgt(s) * dt / 4 )
!                        nhc_p(h,i,j) = nhc_p(h,i,j) + ( nhc_p(h-1,i,j)*nhc_p(h-1,i,j)/nhc_Q(h-1) - T ) * (wgt(s) * dt / 2)
!                        nhc_p(h,i,j) = nhc_p(h,i,j) * EXP( -nhc_p(h+1,i,j)/nhc_Q(h+1) * wgt(s) * dt / 4 )
!                    end do
!                    nhc_p(1,i,j) = nhc_p(1,i,j) * EXP( -nhc_p(2,i,j)/nhc_Q(2) * wgt(s) * dt / 4 )
!                    nhc_p(1,i,j) = nhc_p(1,i,j) + ( sys_p(i,j)**2/sys_m(i,j) - T ) &
!                                * (wgt(s) * dt / 2)
!                    nhc_p(1,i,j) = nhc_p(1,i,j) * EXP( -nhc_p(2,i,j)/nhc_Q(2) * wgt(s) * dt / 4 )
!                    !!
!                    do h=1,nchain
!                        nhc_x(h,i,j) = nhc_x(h,i,j) + nhc_p(h,i,j) / nhc_Q(h)* wgt(s) * dt
!                    end do
!                    sys_p(i,j) = sys_p(i,j) * EXP( -nhc_p(1,i,j)/nhc_Q(1)*wgt(s) * dt )
!                    !!
!                    nhc_p(1,i,j) = nhc_p(1,i,j) * EXP( -nhc_p(2,i,j)/nhc_Q(2) * wgt(s) * dt / 4 )
!                    nhc_p(1,i,j) = nhc_p(1,i,j) + ( sys_p(i,j)**2/sys_m(i,j) - T ) & 
!                                * (wgt(s) * dt / 2)
!                    nhc_p(1,i,j) = nhc_p(1,i,j) * EXP( -nhc_p(2,i,j)/nhc_Q(2) * wgt(s) * dt / 4 )
!                    do h=2,nchain-1,1
!                        nhc_p(h,i,j) = nhc_p(h,i,j) * EXP( -nhc_p(h+1,i,j)/nhc_Q(h+1) * wgt(s) * dt / 4 )
!                        nhc_p(h,i,j) = nhc_p(h,i,j) + ( nhc_p(h-1,i,j)*nhc_p(h-1,i,j)/nhc_Q(h-1) - T ) * (wgt(s) * dt / 2)
!                        nhc_p(h,i,j) = nhc_p(h,i,j) * EXP( -nhc_p(h+1,i,j)/nhc_Q(h+1) * wgt(s) * dt / 4 )
!                    end do
!                    nhc_p(nchain,i,j) = nhc_p(nchain,i,j) + ( nhc_p(nchain-1,i,j)*nhc_p(nchain-1,i,j)/nhc_Q(nchain-1) - T ) &
!                                    *(wgt(s) * dt / 2)
!                end do
!            end do
!        end do
!    end do
!end subroutine thermo_NHC1
!
!subroutine thermo_NHC2(sys_p, sys_m, dt, T, P) !-- polynomial version
!    real(dp), intent(inout) :: sys_p(:,:)
!    real(dp), intent(in) :: sys_m(:,:)
!    real(dp), intent(in) :: dt, T
!    integer, intent(in) :: P
!    integer :: h, i, j, k, s
!    do i=1,ndof
!        do j=1,P
!            do k=1,nrespa
!                do s=1,nsy
!                    nhc_p(nchain,i,j) = nhc_p(nchain,i,j) + ( nhc_p(nchain-1,i,j)*nhc_p(nchain-1,i,j)/nhc_Q(nchain-1) - T ) &
!                                    *( wgt(s)*dt/2 )
!                    do h=nchain-1,2,-1
!                        nhc_p(h,i,j) = nhc_p(h,i,j) - nhc_p(h+1,i,j)/nhc_Q(h+1) * wgt(s) * dt / 4
!                        nhc_p(h,i,j) = nhc_p(h,i,j) + ( nhc_p(h-1,i,j)*nhc_p(h-1,i,j)/nhc_Q(h-1) - T ) * (wgt(s) * dt / 2)
!                        nhc_p(h,i,j) = nhc_p(h,i,j) - nhc_p(h+1,i,j)/nhc_Q(h+1) * wgt(s) * dt / 4
!                    end do
!                    nhc_p(1,i,j) = nhc_p(1,i,j) - nhc_p(2,i,j)/nhc_Q(2) * wgt(s) * dt / 4 
!                    nhc_p(1,i,j) = nhc_p(1,i,j) + ( sys_p(i,j)**2/sys_m(i,j) - T ) &
!                                 * (wgt(s) * dt / 2)
!                    nhc_p(1,i,j) = nhc_p(1,i,j) - nhc_p(2,i,j)/nhc_Q(2) * wgt(s) * dt / 4 
!                    !!
!                    do h=1,nchain
!                        nhc_x(h,i,j) = nhc_x(h,i,j) + nhc_p(h,i,j) / nhc_Q(h)* wgt(s) *dt
!                    end do
!                    sys_p(i,j) = sys_p(i,j) - nhc_p(1,i,j)/nhc_Q(1)* wgt(s) * dt
!                    print *,'i j k a p p1', i, j, k, s, sys_p(i,j), nhc_p(1,i,j)
!                    !!
!                    nhc_p(1,i,j) = nhc_p(1,i,j) - nhc_p(2,i,j)/nhc_Q(2) * wgt(s) * dt / 4
!                    nhc_p(1,i,j) = nhc_p(1,i,j) + ( sys_p(i,j)**2/ sys_m(i,j) - T ) &
!                                * (wgt(s) * dt / 2)
!                    nhc_p(1,i,j) = nhc_p(1,i,j) - nhc_p(2,i,j)/nhc_Q(2) * wgt(s) * dt / 4
!                    do h=2,nchain-1,1
!                        nhc_p(h,i,j) = nhc_p(h,i,j) - nhc_p(h+1,i,j)/nhc_Q(h+1) * wgt(s) * dt / 4
!                        nhc_p(h,i,j) = nhc_p(h,i,j) + ( nhc_p(h-1,i,j)*nhc_p(h-1,i,j)/nhc_Q(h-1) - T ) * (wgt(s) * dt / 2)
!                        nhc_p(h,i,j) = nhc_p(h,i,j) - nhc_p(h+1,i,j)/nhc_Q(h+1) * wgt(s) * dt / 4
!                    end do
!                    nhc_p(nchain,i,j) = nhc_p(nchain,i,j) + ( nhc_p(nchain-1,i,j)*nhc_p(nchain-1,i,j)/nhc_Q(nchain-1) - T ) &
!                                    *(wgt(s) * dt / 2)
!                end do
!            end do
!        end do
!    end do
!end subroutine thermo_NHC2

end module nhc_plus_pimd



