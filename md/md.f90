!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- PIMD Simulation procedure
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pimd
use const
use simul, dtime=>parmdt, nstep=>parmN, sstep=>parmNs, beta=>parmB
use model, only: ndof, M ! you can revise here #TODO
implicit none
    real(8), dimension(:), pointer :: ptr_x1
    real(8), dimension(:), pointer :: ptr_p1
    real(8), dimension(:), pointer :: ptr_f1
    real(8), dimension(:), pointer :: ptr_m1
    
    integer, parameter :: max_estimators = 20
    real(8), dimension(max_estimators) :: estimators

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- force field loader
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine update_fx()
use model, only: V1 => model_V1
implicit none
    real(dp) :: v
    
    call V1(v, ptr_f1, ptr_x1, 1)

end subroutine update_fx


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- x updating procedure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine update_x(dt)
implicit none
    real(8), intent(in) :: dt
    integer :: j
    !-- propagation on x1 (or said ks)
    do j=1, ndof
        ptr_x1(j) = ptr_x1(j) + ptr_p1(j) * dt / ( ptr_m1(j) )
    enddo
    !-- relative updating--------------------- #TODO

    call update_fx()

end subroutine update_x


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- p updating procedure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine update_p(dt)
implicit none
    real(8), intent(in) :: dt
    integer :: i, j
    real(8) :: bf2

    do j=1, ndof
        ptr_p1(j) = ptr_p1(j) - ptr_f1(j) * dt
    enddo
    
end subroutine update_p


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- thermostat procedure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- flag = 0 : NVE
!-- flag = 1 : langevin thermostat
!-- flag = 2 : andersen thermostat
!-- flag = 3 : nose-hoover chain
subroutine update_t(dt)
use random
use nhc_plus_md
implicit none
    real(8), intent(in) :: dt
    integer :: j
    real(8) :: et, c1, c2p, tmp_rand
    
    !-- NVE ensemble
    if(thermo_flg .eq. 0) return
    
    c1 = exp( - thermo_gamma * dt )
    c2p = sqrt( 1.0 - c1**2 )
    do j=1, ndof
        select case (thermo_flg)
        !-- NVE
        case (0)
            return
        !-- Langevin
        case (1)
            call random_normal(et)
            ptr_p1(j) = c1 * ptr_p1(j) + c2p * sqrt( ptr_m1(j) / beta ) * et
        !-- Andersen
        case (2)
            call random_number(tmp_rand)
            if ( tmp_rand < 1. - c1 ) then
                call random_normal(et)
                ptr_p1(j) = sqrt( ptr_m1(j) / beta ) * et
            endif
        !-- NHC
        case (3)
            !--stop "refer to NHC module"
            call thermo_NHC(ptr_p1, ptr_m1, dt, 1.0/beta )
        endselect
    enddo

end subroutine update_t


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- instructor of pimd
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine instructor(ibox1, ifile, iflag)
use random
use myobj
implicit none
    type(box1d), intent(inout) :: ibox1
    character(*), intent(in), optional :: ifile
    integer, intent(in), optional :: iflag
    integer :: i, j
    real(8), dimension(8) :: readline8

    !-- creat box object and assign pointer to box object
    !-- box1: trans_plus  system, box2: primtive system
    call init_box1d(ibox1, ndof)
    
    ptr_x1 => ibox1%x
    ptr_p1 => ibox1%p
    ptr_f1 => ibox1%f
    ptr_m1 => ibox1%m
    
    !-----------------------------------
    !-- assign mass of each nuclear dof
    ptr_m1 = M
 
    if(.not. present(iflag)) stop "loss of flag"
    select case (iflag)
        !-- with x=0 and p=0
        case(0)
            ptr_x1 = 0.0; ptr_p1 = 0.0
            ptr_f1 = 0.0;
        !-- with x=0 and p dist
        case(1)
            ptr_x1 = 0.0
            ptr_f1 = 0.0; 
            do j=1, ndof
                call random_normal(ptr_p1(j))
                ptr_p1(j) = ptr_p1(j) * sqrt( ptr_m1(j) / beta )
            enddo
        !-- start from restart file "final.rst"/or others
        case(2)
            if(.not. present(ifile)) then
            open(unit=10,file='final.rst')
            do j=1, ndof
                read(10,*) readline8
                ptr_x1(j) = readline8(4)
                ptr_p1(j) = readline8(5)
                ptr_f1(j) = readline8(7)
            enddo
            close(unit=10)
            else
            open(unit=10,file=trim(ifile))
            do j=1, ndof
                read(10,*) readline8
                ptr_x1(j) = readline8(4)
                ptr_p1(j) = readline8(5)
                ptr_f1(j) = readline8(7)
            enddo
            close(unit=10)
            endif
        case default
            stop "argument dismarch"
    endselect
    
    !-- if initialization from guess, complete the info
    if(iflag < 2) then
        call update_fx()
    endif
    return
end subroutine instructor


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- simulation procedure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine pimd_process(traj_file, anlys_file)
use nhc_plus_md
implicit none
    character(*), intent(in) :: traj_file, anlys_file
    real(8) :: halfdt
    integer :: npass
    integer :: i,j
    
    
    open(unit=trj_unit, file=traj_file, status='replace')
    open(unit=ana_unit, file=anlys_file, status='replace')
    
    halfdt = dtime/2.0
    npass = 0
    estimators = 0.0
    
    call set_NHC(5, ndof, 5, thermo_gamma, 1.0/beta ) !-- bead = 1
    
    do i = 1, nstep
        !-- using middle scheme
        call update_p(halfdt)
        call update_x(halfdt)
        call update_t(dtime)
        call update_x(halfdt) !---------------------- #TODO
        call update_p(halfdt)
        npass = npass + 1
        if(mod(npass, sstep) .eq. 0) then
            call estimator()
            call sampler(npass)
        endif
    enddo
    close(unit=trj_unit)
    close(unit=ana_unit)
    
    !-- write to a restart file
    open(unit=trj_unit, file='final.rst',status='replace')
    do j=1, ndof
            !-- only position space is meaningful
        write(trj_unit,*) npass, j, ptr_x1(j), ptr_p1(j), ptr_f1(j)
    enddo
    close(unit=trj_unit)
    
    call del_NHC()

end subroutine pimd_process


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- estimator procedure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine estimator()
use model, V0 => model_V0
implicit none !#TODO
    real(dp), allocatable :: v, dv(:), ddv(:,:)
    real(dp) :: Vtot, Kprim, Kvir
    integer :: i
    
    return
end subroutine estimator


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- sampler procedure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine sampler(npass)
implicit none
    integer, intent(in) :: npass
    integer :: i,j
    write (ana_unit,*) npass, estimators
    do j=1, ndof
        write (trj_unit,*) npass, j, ptr_x1(j), ptr_p1(j), ptr_f1(j)
    enddo
end subroutine sampler

end module pimd




