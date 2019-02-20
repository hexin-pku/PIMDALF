!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- PIMD Simulation procedure
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Harmonic Oscillitor Simulation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


program HO_PIMD
use random
use myobj
use pisimul
use pimd
implicit none
    character(len=32) :: traj_file, anlys_file, start_file, parms_file
    character(len=32) :: arg1, arg2
    integer :: n, i
    integer :: start_flg
    logical :: set_equilibrium = .false.
    real(8) :: time1, time2
    type(box2d) :: box1, box2
    
    !-- default values
    start_file = ''
    parms_file = trim('put.par')
    
    n = command_argument_count()
    if ((mod(n,2) .ne. 0)) then
        stop "args number error"
    endif
    !-- if nessessary check for some arguments
    do i=1,n/2
        call get_command_argument(2*i-1,arg1)
        call get_command_argument(2*i, arg2)
        select case (trim(arg1))
            !-- output name
            case ("-o")
                traj_file = trim(arg2)//'.trj'
                anlys_file = trim(arg2)//'.ana'
            !-- start mod
            case ("-s")
                select case (arg2)
                    case ("0")
                        start_flg = 0
                    case ("1")
                        start_flg = 1
                    case ("e")
                        start_flg = 1
                        set_equilibrium = .true.
                    case ("r")
                        start_flg = 2
                    case ("f")
                        start_flg = 2
                    case default
                        stop "args mismatch at 2"
                endselect
            !-- parameter file location
            case ("-p")
                parms_file = trim(arg2)
            case ("-f")
                start_file = trim(arg2)
            case default
                stop "args mismatch at 3"
        endselect
    enddo
    
    call simul_readparms(trim(parms_file))
    if( set_equilibrium ) then
        parmN = 500
    endif
    
    call init_seed()
    
    if(start_file .eq. '') then
        call instructor(box1, box2, iflag=start_flg)
    else
        call instructor(box1, box2, ifile=trim(start_file), iflag=start_flg)
    endif
    
    call cpu_time(time1)
    print *, 'debug'
    call pimd_process(trim(traj_file), trim(anlys_file))
    call cpu_time(time2)
    
    print *, "Using CPU time:  ", time2-time1
    
    call del_box2d(box1)
    call del_box2d(box2)

end program HO_PIMD



