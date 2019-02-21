!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- PIMD Simulation procedure
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Nose Hoover Chain Thermostat with RESPA-alg.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module nhc_plus_md
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
    real(dp), dimension(:,:), allocatable :: nhc_x, nhc_p, nhc_G
    real(dp), dimension(:), allocatable :: nhc_Q

    public :: set_NHC, del_NHC, thermo_NHC
contains

subroutine set_NHC(nchain_in, ndof_in, nrespa_in, gamma_in, T)
    integer, intent(in) :: nchain_in, ndof_in, nrespa_in
    !-- P: number of beads
    real(dp), intent(in) :: gamma_in, T 
    !-- gamma_in : optimal parameter
    !-- T : temperature
    
    nchain = nchain_in
    ndof = ndof_in
    nrespa = nrespa_in
    
    if( (ndof < 1) .or. ( nchain < 1 ) .or. (nrespa < 1) ) stop 'generate NHC error'

    allocate( nhc_x(nchain, ndof), nhc_p(nchain, ndof),&
              nhc_G(nchain, ndof), nhc_Q(nchain) )
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
subroutine thermo_NHC(sys_p, sys_m, dt, T)
    real(dp), intent(inout) :: sys_p(:)
    real(dp), intent(in) :: sys_m(:)
    real(dp), intent(in) :: dt, T

    integer :: h, i, k, s
    do i=1,ndof
        do k=1,nrespa
            do s=1,nsy
                !-- update auxiliary variable nhc_G
                nhc_G(1,i) = sys_p(i)**2/sys_m(i) - T
                do h=2,nchain
                    nhc_G(h,i) = nhc_p(h-1,i)*nhc_p(h-1,i)/nhc_Q(h-1) - T
                end do
                !-- update nhc_p from tail to head
                nhc_p(nchain,i) = nhc_p(nchain,i) &
                            + nhc_G(nchain,i) * wgt(s) * dt / 2
                do h=nchain-1,1,-1
                    nhc_p(h,i) = nhc_p(h,i) &
                            * EXP( -nhc_p(h+1,i)/nhc_Q(h+1) * wgt(s) * dt / 4 )
                    nhc_p(h,i) = nhc_p(h,i) &
                            + nhc_G(h,i) * wgt(s) * dt / 2
                    nhc_p(h,i) = nhc_p(h,i) &
                            * EXP( -nhc_p(h+1,i)/nhc_Q(h+1) * wgt(s) * dt / 4 )
                end do
                !-- update nhc_x and system_p
                do h=1,nchain
                    nhc_x(h,i) = nhc_x(h,i) &
                            + nhc_p(h,i) / nhc_Q(h)* wgt(s) * dt
                end do
                sys_p(i) = sys_p(i) &
                            * EXP( -nhc_p(1,i)/nhc_Q(1)*wgt(s) * dt )
                !-- update nhc_p from head to tail
                do h=1,nchain-1,1
                    nhc_p(h,i) = nhc_p(h,i) &
                            * EXP( -nhc_p(h+1,i)/nhc_Q(h+1) * wgt(s) * dt / 4 )
                    nhc_p(h,i) = nhc_p(h,i) &
                            + nhc_G(h,i) * wgt(s) * dt / 2
                    nhc_p(h,i) = nhc_p(h,i) &
                            * EXP( -nhc_p(h+1,i)/nhc_Q(h+1) * wgt(s) * dt / 4 )
                end do
                nhc_p(nchain,i) = nhc_p(nchain,i) + nhc_G(nchain,i) *wgt(s) * dt / 2
            end do
        end do
    end do
end subroutine thermo_NHC

end module nhc_plus_md




