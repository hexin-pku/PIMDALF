!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- PIMD Simulation procedure
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pimd
use const
use pisimul, dtime=>parmdt, nstep=>parmN, sstep=>parmNs, nbead=>parmP, beta=>parmB
use model, only: ndof, M ! you can revise here #TODO
implicit none
    real(8), dimension(:,:), pointer :: ptr_x1
    real(8), dimension(:,:), pointer :: ptr_p1
    real(8), dimension(:,:), pointer :: ptr_f1
    real(8), dimension(:,:), pointer :: ptr_m1
    real(8), dimension(:,:), pointer :: ptr_x2
    real(8), dimension(:,:), pointer :: ptr_p2
    real(8), dimension(:,:), pointer :: ptr_f2
    real(8), dimension(:,:), pointer :: ptr_m2
    
    integer, parameter :: max_estimators = 20
    real(8), dimension(max_estimators) :: pimd_estimators

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- force field loader
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine update_fx()
use model, only: V2 => model_V2
implicit none
    real(dp) :: v, ddv(ndof, ndof)
    integer :: i
    do i=1,nbead
        !-- give mes-force of a configuration
        call V2(v, ptr_f2(:,i), ddv, ptr_x2(:,i), 1)
    enddo
end subroutine update_fx


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- x updating procedure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine update_x(dt)
use trans_plus
implicit none
    real(8), intent(in) :: dt
    integer :: i, j
    !-- propagation on x1 (or said ks)
    do j=1, ndof
        do i=1, nbead
            ptr_x1(j,i) = ptr_x1(j,i) + ptr_p1(j,i) * dt / ( ptr_m1(j,i) )
        enddo
    enddo
    !-- relative updating--------------------- #TODO
    do j=1, ndof
        call ks2x(ptr_x2(j,:), ptr_x1(j,:), nbead)
    enddo
    call update_fx()
    do j=1, ndof
        call fx2fks(ptr_f1(j,:), ptr_f2(j,:), nbead)
    enddo
end subroutine update_x


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- p updating procedure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine update_p(dt)
implicit none
    real(8), intent(in) :: dt
    integer :: i, j
    real(8) :: bf2
    
    bf2 = (nbead/beta**2)
    do j=1, ndof
        ptr_p1(j,1) = ptr_p1(j,1) - ptr_f1(j,1) * dt
        do i=2, nbead
            ptr_p1(j,i) = ptr_p1(j,i) - ( ptr_f1(j,i) + ptr_m1(j,i) * bf2 * ptr_x1(j,i) ) * dt
        enddo
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
use nhc_plus_pimd
implicit none
    real(8), intent(in) :: dt
    integer :: i, j
    real(8) :: et, c1, c2p, tmp_rand
    
    !-- NVE ensemble
    if(thermo_flg .eq. 0) return
    
    c1 = exp( - thermo_gamma * dt )
    c2p = sqrt( 1.0 - c1**2 )
    do j=1, ndof
        do i=1, nbead
            select case (thermo_flg)
            !-- NVE
            case (0)
                return
            !-- Langevin
            case (1)
                call random_normal(et)
                ptr_p1(j,i) = c1 * ptr_p1(j,i) + c2p * sqrt( ptr_m1(j,i) / beta ) * et
            !-- Andersen
            case (2)
                call random_number(tmp_rand)
                if ( tmp_rand < 1. - c1 ) then
                    call random_normal(et)
                    ptr_p1(j,i) = sqrt( ptr_m1(j,i) / beta ) * et
                endif
            !-- NHC
            case (3)
                !--stop "refer to NHC module"
                call thermo_NHC(ptr_p1, ptr_m1, dt, 1.0/beta, nbead)
            endselect
        enddo
    enddo
end subroutine update_t


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- instructor of pimd
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine instructor(ibox1, ibox2, ifile, iflag)
use random
use trans_plus
use myobj
implicit none
    type(box2d), intent(inout) :: ibox1, ibox2
    character(*), intent(in), optional :: ifile
    integer, intent(in), optional :: iflag
    integer :: i, j
    real(8), dimension(8) :: readline8

    !-- creat box object and assign pointer to box object
    !-- box1: trans_plus  system, box2: primtive system
    call init_box2d(ibox1, ndof, nbead)
    call init_box2d(ibox2, ndof, nbead)
    ptr_x1 => ibox1%x
    ptr_p1 => ibox1%p
    ptr_f1 => ibox1%f
    ptr_m1 => ibox1%m
    ptr_x2 => ibox2%x
    ptr_p2 => ibox2%p
    ptr_f2 => ibox2%f
    ptr_m2 => ibox2%m
    
    !-----------------------------------
    !-- assign mass of each nuclear dof
    do i=1,ndof
        ptr_m2(i,:) = M(i)
    enddo
    call m2mb(ptr_m1, ptr_m2, nbead)
    
    !-- initialization x(or said ks) and p
    
    if(.not. present(iflag)) stop "loss of flag"
    select case (iflag)
        !-- with x=0 and p=0
        case(0)
            ptr_x1 = 0.0; ptr_p1 = 0.0
            ptr_x2 = 0.0
            ptr_f1 = 0.0; ptr_f2 = 0.0
        !-- with x=0 and p dist
        case(1)
            ptr_x1 = 0.0
            ptr_x2 = 0.0
            ptr_f1 = 0.0; ptr_f2 = 0.0
            do j=1, ndof
                do i=1, nbead
                    call random_normal(ptr_p1(j,i))
                    ptr_p1(j,i) = ptr_p1(j,i) * sqrt( ptr_m1(j,i) / beta )
                enddo
            enddo
        !-- start from restart file "final.rst"/or others
        case(2)
            if(.not. present(ifile)) then
            open(unit=10,file='final.rst')
            do j=1, ndof
                do i=1, nbead
                    read(10,*) readline8
                    ptr_x1(j,i) = readline8(4)
                    ptr_p1(j,i) = readline8(5)
                    ptr_x2(j,i) = readline8(6)
                    ptr_f1(j,i) = readline8(7)
                    ptr_f2(j,i) = readline8(8)
                enddo
            enddo
            close(unit=10)
            else
            open(unit=10,file=trim(ifile))
            do j=1, ndof
                do i=1, nbead
                    read(10,*) readline8
                    ptr_x1(j,i) = readline8(4)
                    ptr_p1(j,i) = readline8(5)
                    ptr_x2(j,i) = readline8(6)
                    ptr_f1(j,i) = readline8(7)
                    ptr_f2(j,i) = readline8(8)
                enddo
            enddo
            close(unit=10)
            endif
        case default
            stop "argument dismarch"
    endselect
    
    !-- if initialization from guess, complete the info
    if(iflag < 2) then
        do j=1, ndof
            call ks2x(ptr_x2(j,:), ptr_x1(j,:), nbead)
        enddo
        call update_fx()
        do j=1, ndof
            call fx2fks(ptr_f1(j,:), ptr_f2(j,:), nbead)
        enddo
    endif
    return
end subroutine instructor


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- simulation procedure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine pimd_process(traj_file, anlys_file)
use nhc_plus_pimd
use trans_plus
implicit none
    character(*), intent(in) :: traj_file, anlys_file
    real(8) :: halfdt
    integer :: npass
    integer :: i,j
    
    open(unit=trj_unit, file=traj_file, status='replace')
    open(unit=ana_unit, file=anlys_file, status='replace')
    
    halfdt = dtime/2.0
    npass = 0
    pimd_estimators = 0.0
    
    call set_NHC(5, ndof, 5, thermo_gamma, 1.0/beta, nbead)
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
        do i=1, nbead
            !-- only position space is meaningful
            write(trj_unit,*) npass, j, i, ptr_x1(j,i), ptr_p1(j,i),&
            ptr_x2(j,i), ptr_f1(j,i), ptr_f2(j,i)
        enddo
    enddo
    close(unit=trj_unit)
    
    call del_NHC()
    call del_trans()
end subroutine pimd_process


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- estimator procedure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine estimator()
use model, V2 => model_V2
implicit none !#TODO
    real(dp), allocatable :: dv(:), ddv(:,:), Xc(:)
    real(dp) :: v, Vtot, Kprim, Kvir
    real(dp) :: bf2
    integer :: i,j
    
    allocate(dv(ndof), ddv(ndof,ndof), Xc(ndof))
    
    Vtot = 0.0_8
    Kprim = 0.5_8 * ndof * nbead / beta
    Kvir = 0.5_8 * ndof / beta
    Xc = sum( ptr_x2, dim=2) / nbead
    bf2 = (nbead/beta**2)
    !-- print *, Xc
    !-- stop 'debug'
    
    do i=1,nbead
        call V2(v, dv, ddv, ptr_x2(:,i), 0)
        Vtot = Vtot + v
    enddo
    do i=1,ndof
        Kprim = Kprim - 0.5_8 * sum( ptr_m1(i,2:) * bf2 * ptr_x1(i,2:)**2 )
        Kvir = Kvir + 0.5_8 * dot_product( ptr_x2(i,:) - Xc(i), ptr_f2(i,:) ) / nbead
    enddo
    Vtot = Vtot / nbead
    pimd_estimators(1) = Vtot
    pimd_estimators(2) = Kprim
    pimd_estimators(3) = Kvir
    pimd_estimators(4) = Vtot + Kprim
    pimd_estimators(5) = Vtot + Kvir
    return
end subroutine estimator


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- sampler procedure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine sampler(npass)
implicit none
    integer, intent(in) :: npass
    integer :: i,j
    write (ana_unit,*) npass, pimd_estimators
    do j=1, ndof
        do i=1, nbead
            write (trj_unit,*) npass, j, i, ptr_x1(j,i), ptr_p1(j,i),&
            ptr_x2(j,i), ptr_f1(j,i), ptr_f2(j,i)
        enddo
    enddo
end subroutine sampler

end module pimd




