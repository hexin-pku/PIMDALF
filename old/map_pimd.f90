!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- PIMD Simulation procedure (Copyright Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!-- Module for MD simulations
module map_pimd
use base
implicit none
    !-- finite & numerical
    real(8) :: md_dtime
    integer :: md_nstep
    integer :: md_sstep = 1000
    integer :: md_npass
    
    !-- freedom and bead
    integer :: md_nfree !-- should consistent with the models
    integer :: md_nbead
    real(8) :: md_bfreq
    real(8) :: md_bf2
    
    !-- for thermostat process
    real(8) :: md_temp
    real(8) :: md_beta
    real(8) :: md_coeff
    real(8) :: md_cfreq_optimval
    integer :: md_scheme_flg
    integer :: md_thermo_flg
    
    !-- save for rand
    real(8) :: md_rand
    
    
    !-- IO settings
    integer, parameter :: trj_unit = 10, ana_unit=11, frm_unit=12, elc_unit=20
    integer :: md_start_flg = 1
    
    !-- MD varibles
    real(8), dimension(:,:), allocatable :: box_x
    real(8), dimension(:,:), allocatable :: box_ks
    real(8), dimension(:,:), allocatable :: box_p
    real(8), dimension(:,:), allocatable :: box_fx
    real(8), dimension(:,:), allocatable :: box_fks
    real(8), dimension(:,:), allocatable :: box_m
    real(8), dimension(:), allocatable :: md_ana_vals

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- for Staging transformation
subroutine x2ks(ks, x, P)
implicit none
    integer, intent(in) :: P
    real(8), dimension(P), intent(in) :: x
    real(8), dimension(P), intent(out) :: ks
    integer :: i
    ks(1) = x(1)
    if(P.eq.1) return
    ks(P) = x(P) - x(1)
    do i=2, P-1
        ks(i) = x(i) - ( x(i+1) * real(i - 1) + x(1) ) / real(i)
    enddo
end subroutine x2ks

subroutine ks2x(x, ks, P)
implicit none
    integer, intent(in) :: P
    real(8), dimension(P), intent(in) :: ks
    real(8), dimension(P), intent(out) :: x
    integer :: i
    x(1) = ks(1)
    if(P .eq.1) return
    x(P) = ks(P) + ks(1)
    do i=P-1,2,-1
        x(i) = ks(i) + ( real( i - 1 ) * x(i+1) + ks(1) ) / real(i)
    enddo
end subroutine ks2x

subroutine fx2fks(fks, fx, P)
implicit none
    integer, intent(in) :: P
    real(8), dimension(P), intent(in) :: fx
    real(8), dimension(P), intent(out) :: fks
    integer :: i

    fks(1) = 0.0
    do i=1,P
        fks(1) = fks(1) + fx(i)
    enddo
    
    do i=2,P,1
        fks(i) = fx(i) + fks(i-1) * real(i-2)/real(i-1)
    enddo
    !-- note here all fks is without Harmonic Oscillitor term!
    fks = fks /real(P)
end subroutine fx2fks
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- initialization of molecule dynamics numerical settings
subroutine init_md_parameters(par_file)
implicit none
    character(*), intent(in) :: par_file
    character(len=32) :: tmp_string1, tmp_string2
    integer :: iostat_here
    logical :: exist_here
    logical, dimension(14) :: is_completed = .false.
    integer :: i
    
    inquire(file=trim(par_file), exist=exist_here)
    if (exist_here .eqv. .true.) then
        open(unit=10, file='./'//trim(par_file))
        do while (.true.)
            read(10,*,iostat=iostat_here) tmp_string1, tmp_string2
            if (iostat_here < 0) exit
            select case (tmp_string1)
                case ('beta')
                    read(tmp_string2,*) md_beta
                    md_temp = 1.0/md_beta
                    is_completed(1) = .true.
                case ('cfreq')
                    read(tmp_string2,*) md_cfreq_optimval
                    is_completed(2) = .true.
                case ('dtime')
                    read(tmp_string2,*) md_dtime
                    is_completed(3) = .true.
                case ('coeff')
                    read(tmp_string2,*) md_coeff
                    is_completed(4) = .true.
                case ('nstep')
                    read(tmp_string2,*) md_nstep
                    is_completed(5) = .true.
                case ('ischeme')
                    read(tmp_string2,*) md_scheme_flg
                    is_completed(6) = .true.
                case ('ithermo')
                    read(tmp_string2,*) md_thermo_flg
                    is_completed(7) = .true.
                case ('nfree')
                    read(tmp_string2,*) md_nfree
                    is_completed(8) = .true.
                case ('sstep')
                    read(tmp_string2,*) md_sstep
                    is_completed(9) = .true.
                case ('nbead')
                    read(tmp_string2,*) md_nbead
                    is_completed(10) = .true.
                case default
                    cycle
            end select
        enddo
    else
        call log_message("E",'parameter file is not existing, please check it!')
    endif
    
    !-- if character frequency is given a positive number, it will do optimization 
    !------ according to md_cfreq_optimval 
    if (md_cfreq_optimval > 0) then
        md_dtime = 0.2 / md_cfreq_optimval
        !md_nstep = 10000 / md_cfreq_optimval
        select case (md_thermo_flg)
            !-- Langevin
            case (0)
                md_coeff = md_cfreq_optimval
            !-- Anderson
            case (1)
                md_coeff = sqrt(2.) * md_cfreq_optimval
            !-- NH Chain
            case (2)
                md_coeff = 20 * md_dtime
        end select
    endif
    
    !-- initial settings at time t=0
    md_npass = 0
    md_bfreq = sqrt( real(md_nbead) ) * md_temp
    md_bf2 = real(md_nbead) * md_temp * md_temp
    
    !-- check initialization for parameters
    do i=1,10
        if (is_completed(i) .eqv. .false.) then
            call log_message("E", 'Initialization is uncompleled!')
        endif
    enddo
end subroutine init_md_parameters


subroutine calc_sys_fx()
use MES_Models
use map
implicit none
    !-- MAP-PIMD real dynamics
    call map_fx(box_fx, box_x)
    !-- MES-PIMD for statistic
    !call Mod_esti_dV(box_fx(:,:), box_x(:,:))
    !-- HO Model
    !box_fx = box_x
end subroutine calc_sys_fx

subroutine calc_sys_ana()
use MES_Models
implicit none
    real(8) :: temp1, temp2, temp3
    
    !-- FOR HO MODEL
    !md_ana_vals(3) = 0.5*sum(box_x*box_x) /real(md_nbead)
    !return
    !----------------------------------------------------
    
    md_ana_vals(1) = update_pfd        !-- estimate pfd, since dV calculation, it has been evaluated
    !--
    call Mod_esti_PFF(temp1, box_x)
    md_ana_vals(2) = update_pff        !-- estimate pff
    !--
    call Mod_esti_V(temp1, box_x)
    md_ana_vals(3) = temp1                  !-- estimate V
    !--
    call Mod_esti_K(temp2, temp3, box_x, iflag=1)
    md_ana_vals(4) = temp2                  !-- estimate Kprim
    md_ana_vals(5) = temp3                  !-- estimate Kvir
    md_ana_vals(6) = md_ana_vals(3) + temp2 !-- estimate Eprim
    md_ana_vals(7) = md_ana_vals(3) + temp3 !-- estmate Evir
    !--
    call Mod_esti_coh(temp1)           
    md_ana_vals(8) = update_coh        !-- estimate Lcohenren
    !--
    call Mod_esti_theta(box_x)
    md_ana_vals(9:12) = update_theta(1:4)
end subroutine calc_sys_ana


end module map_pimd
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- where fx and analysis varibales should be specific for models


!-- initialization system at x=0
subroutine init_sys_parameters(sta_file)
use map_pimd
use MES_Models
use map
implicit none
    character(*), intent(in) :: sta_file
    integer :: i, j
    real(8), dimension(8) :: linevals
    real(8), dimension(:), allocatable :: sys_m
    
    !-- Interface to MES_Models
    !------------------------------------------
    !-- re-update md_nfree & mass from models--
    !------ Here take double ES model as example ----
    call init_Model(md_nbead, md_beta, 1 )
    md_nfree = getrdim()
    allocate(sys_m(md_nfree))
    sys_m = Mod_M
    !------ Or test for Harmonic Oscillitor models
    !md_nfree = 1
    !allocate(sys_m(md_nfree))
    !sys_m = 1.0
    !!-- Note what force should call
    !------------------------------------------
    
    !-- Interface to MAP
    call init_map(nb=md_nbead, nc=4, ipop=2, istat=1)
    
    allocate(box_x(md_nfree, md_nbead))
    allocate(box_ks(md_nfree, md_nbead))
    allocate(box_p(md_nfree, md_nbead))
    allocate(box_fx(md_nfree, md_nbead))
    allocate(box_fks(md_nfree, md_nbead))
    allocate(box_m(md_nfree, md_nbead))
    allocate(md_ana_vals(12))
    
    do j=1, md_nfree
        box_m(j,1) = sys_m(j)
        do i=2,md_nbead
            box_m(j,i) = sys_m(j) * real(i) / real( i - 1 )   
        enddo
    enddo
    
    !-- initialization x(or said ks) and p
    select case (md_start_flg)
        !-- with x=0 and p=0
        case(0)
            box_ks = 0.0
            box_p = 0.0
            box_x = 0.0
            box_fks = 0.0
            box_fx = 0.0
        !-- with x=0 and p distribution
        case(1)
            box_ks = 0.0
            do j=1, md_nfree
                do i=1, md_nbead
                    call random_norm(box_p(j,i), md_rand, sqrt( box_m(j,i) * md_temp ) )
                enddo
            enddo
        !-- with x and p distribution
        case(2)
            do j=1, md_nfree
                do i=1, md_nbead
                    call random_norm(box_p(j,i), md_rand, sqrt( md_temp/box_m(j,i) ) / md_cfreq_optimval )
                    call random_norm(box_p(j,i), md_rand, sqrt( box_m(j,i) * md_temp ) )
                enddo
            enddo
        !-- start from restart file, with fixed name "final.rst"
        case(3)
            open(unit=10,file='final.rst')
            do j=1, md_nfree
                do i=1, md_nbead
                    read(10,*) linevals
                    box_ks(j,i) = linevals(4)
                    box_p(j,i) = linevals(5)
                    box_x(j,i) = linevals(6)
                    box_fks(j,i) = linevals(7)
                    box_fx(j,i) = linevals(8)
                enddo
            enddo
            close(unit=10)
        !-- start from a start file
        case(4)
            open(unit=10,file=trim(sta_file))
            do j=1, md_nfree
                do i=1, md_nbead
                    read(10,*) linevals
                    box_ks(j,i) = linevals(4)
                    box_p(j,i) = linevals(5)
                    box_x(j,i) = linevals(6)
                    box_fks(j,i) = linevals(7)
                    box_fx(j,i) = linevals(8)
                enddo
            enddo
            close(unit=10)
        !-- use md_start_flg
        case default
            call log_message("E","arguments mismatch!")
    end select
    
    if(md_start_flg .le. 3) then
        do j=1, md_nfree
            call ks2x(box_x(j,:), box_ks(j,:), md_nbead)
        enddo
        call calc_sys_fx()
        do j=1, md_nfree
            call fx2fks(box_fks(j,:), box_fx(j,:), md_nbead)
        enddo
    endif
    call calc_sys_ana()

    deallocate(sys_m)
end subroutine init_sys_parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- the MD process with x,p,T propagators
subroutine md_process(trj_file, ana_file)
use map_pimd
use map
implicit none
    interface
        subroutine md_update_t(dt)
        use map_pimd
        implicit none
            real(8), intent(in) :: dt
            integer :: i, j
            real(8) :: et
        end subroutine md_update_t
        
        subroutine md_update_x(dt)
        use map_pimd
        implicit none
            real(8), intent(in) :: dt
            integer :: i, j
        end subroutine md_update_x
        
        subroutine md_update_p(dt)
        use map_pimd
        implicit none
            real(8), intent(in) :: dt
            integer :: i, j
        end subroutine md_update_p
        
        subroutine md_sample()
        use map_pimd
        implicit none
            integer :: i, j
        end subroutine
        
    end interface

    character(*), intent(in) :: trj_file, ana_file
    integer :: i,j
    
    open(unit=trj_unit, file=trj_file, status='replace')
    !-- without header
    !write(trj_unit,*) 'nstep    ', 'nfree', 'nbead    ','ks   ','p    ','x    ','fks  ', 'fx   '
    open(unit=ana_unit, file=ana_file, status='replace')
    !-- without header
    !write(ana_unit,*) 'nstep    ', 'anas'
    open(unit=elc_unit, file='elec.dat', status='replace')
    
    !-- middle scheme
    call md_sample()
    do i = 1, md_nstep
        call md_update_p(md_dtime/2.0)
        call md_update_x(md_dtime/2.0)
        !call md_update_t(md_dtime)
        call md_update_x(md_dtime/2.0)
        call md_update_p(md_dtime/2.0)
        md_npass = md_npass + 1
        !-- ADD MAP
        call map_stat_pop()
        if(mod(md_npass, md_sstep) .eq. 0) then
            call md_sample()
        endif
    enddo
    
    close(unit=trj_unit)
    close(unit=ana_unit)
    
    !-- write a restart file
    open(unit=10, file='final.rst',status='replace')
    do j=1, md_nfree
        do i=1, md_nbead
            write(10,*) md_npass, j, i, box_ks(j,i), box_p(j,i), box_x(j,i), box_fks(j,i), box_fx(j,i)
        enddo
    enddo
    close(unit=10)
end subroutine md_process


!-- sampling the trajectory and statistic quantities
subroutine md_sample()
use map_pimd
use map
implicit none
    integer :: i,j
    !-- write analysis file (contain averaged position, all kinds of energy, some other thermodynamics properties)
    write (ana_unit,*) md_npass, md_ana_vals
    do j=1, md_nfree
        do i=1, md_nbead
            write (trj_unit,*) md_npass, j, i, box_ks(j,i), box_p(j,i), box_x(j,i), box_fks(j,i), box_fx(j,i)
        enddo
    enddo
    do i=1, md_nbead
        write (elc_unit,*) md_npass, i, map_pop(i,:), map_c(i,:,:)
    enddo
end subroutine md_sample


!-- an update of thermostat
!-- flag = 1 'lang' : langevin thermostat
!-- flag = 2 'ads' : andersen thermostat
!-- flag = 3 'nhc' : nose-hoover chain
subroutine md_update_t(dt)
use map_pimd
use map
implicit none
    real(8), intent(in) :: dt
    integer :: i, j
    real(8) :: et, c1, c2prime
    
    c1 = exp( - md_coeff * dt )
    c2prime = sqrt( 1 - c1**2 )
    
    do j=1, md_nfree
        do i=1, md_nbead
            select case (md_thermo_flg)
                !-- NVE ensembles
                case (0)
                    return
                !-- Langevin
                case (1)
                    call random_norm(et, md_rand, 1.0_8)
                    box_p(j,i) = c1 * box_p(j,i) + sqrt( box_m(j,i) * md_temp ) * c2prime * et        !--------------------------------- Note
                !-- Anderson
                case (2)
                    call random_number(md_rand)
                    if ( md_rand < 1. - c1 ) then
                        call random_norm(et, md_rand, 1.0_8)
                        box_p(j,i) = sqrt( box_m(j,i) * md_temp ) * et
                    endif
                !-- NHC
                case (3)
                    call log_message("E", 'NHC is invalid for now')
            end select
        enddo
    enddo
end subroutine md_update_t


!-- an update of position
subroutine md_update_x(dt)
use map_pimd
use map
implicit none
    real(8), intent(in) :: dt
    integer :: i, j
    !-- propagation on x( or said ks)
    do j=1, md_nfree
        do i=1, md_nbead
            box_ks(j,i) = box_ks(j,i) + box_p(j,i) * dt / ( box_m(j,i) )
        enddo
    enddo
    call map_update_x(box_x, dt)
    !-- relative updating
    do j=1, md_nfree
        call ks2x(box_x(j,:), box_ks(j,:), md_nbead)
    enddo
    call calc_sys_fx()
    !call calc_sys_ana()
    do j=1, md_nfree
        call fx2fks(box_fks(j,:), box_fx(j,:), md_nbead)
    enddo
end subroutine md_update_x


!-- an update of momentum
subroutine md_update_p(dt)
use map_pimd
use map
implicit none
    real(8), intent(in) :: dt
    integer :: i, j
    do j=1, md_nfree
        box_p(j,1) = box_p(j,1) - box_fks(j,1) * dt
        do i=2, md_nbead
            box_p(j,i) = box_p(j,i) - box_fks(j,i) * dt & !------------------------------------------------------------ bugs 
                - box_m(j,i) * md_bf2 * box_ks(j,i) * dt
        enddo
    enddo
    call map_update_p(box_x, dt)
end subroutine md_update_p
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- main pragram
program main
use map_pimd
implicit none
    interface
        subroutine init_sys_parameters(sta_file)
        use map_pimd
        implicit none
            character(*), intent(in) :: sta_file
            integer :: i, j
            real(8), dimension(7) :: linevals
        end subroutine init_sys_parameters
        
        subroutine md_process( trj_file, ana_file )
        use map_pimd
        implicit none
            character(*), intent(in) :: trj_file, ana_file
            integer :: i, j
        end subroutine md_process
    end interface

    character(len=32) :: trj_file, ana_file, sta_file, par_file
    character(len=32) :: arg1, arg2
    integer :: n, i
    real(8) :: time1, time2
    
    !-- default values
    sta_file = trim('example.smp')
    par_file = trim('put.rc')
    
    n = command_argument_count()
    if ((mod(n,2) .ne. 0)) then
        call log_message("E ","arguments mismatch 1")
    endif
    !-- if nessessary check for some arguments
    do i=1,n/2
        call get_command_argument(2*i-1,arg1)
        call get_command_argument(2*i, arg2)
        select case (trim(arg1))
            !-- output name
            case ("-o")
                trj_file = trim(arg2)//'.trj'
                ana_file = trim(arg2)//'.ana'
            !-- start mod
            case ("-s")
                select case (arg2)
                    case ("0")
                        md_start_flg = 0
                    case ("1")
                        md_start_flg = 1
                    case ("2")
                        md_start_flg = 2
                    case ("r")
                        md_start_flg = 3
                    case ("f")
                        md_start_flg = 4
                    case default
                        call log_message("E","arguments mismarch 2")
                end select
            !-- parameter file location
            case ("-p")
                par_file = trim(arg2)
            case ("-f")
                sta_file = trim(arg2)
            case default
                call log_message("E","arguments mismarch 3")
        end select
    enddo

    
    call init_md_parameters(trim(par_file))
    call init_seed()
    call init_sys_parameters(trim(sta_file))
    call cpu_time(time1)
    call md_process(trim(trj_file), trim(ana_file))
    call cpu_time(time2)
    print *, "Using CPU time:  ", time2-time1
    
    deallocate(box_x, box_ks, box_p, box_fks, box_fx, md_ana_vals)
end program main
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- |  /|  / T ``T`` |  | |``| |  | ``T`` |``> |  | r``\
!-- | / | /  |   |   |--| |  | |  |   |   |-<  |  | |  
!-- |/  |/   |   |   |  | L..| L..L   |   L..> L..L L..T
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



