!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- MES-PIMD with mapping variables for Tully's Model
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program mvpimd
use pisimul, dtime=>parmdt, nstep=>parmN, sstep=>parmNs, nbead=>parmP, beta=>parmB
use md_pimd
use tullymodel
use map_core
use nucl_sample
implicit none
    type(mapvar), allocatable :: imv(:)
    call simul_readparms('s.par')
    call tully_read('tully1.parms')
    
    allocate(imv(nbead))
    
    call init_mes7_smors()
    call simul_readparms(trim(parms_file))
    if( set_equilibrium ) then
        parmN = 500
    endif
    
    call init_seed()
    call init_mes()
    
    if(start_file .eq. '') then
        call instructor(box1, box2, iflag=start_flg)
    else
        call instructor(box1, box2, ifile=trim(start_file), iflag=start_flg)
    endif
    call mes_setlevel(2)
    
    
    call cpu_time(time1)
    call pimd_process(trim(traj_file), trim(anlys_file))
    call cpu_time(time2)
    
    print *, "Using CPU time:  ", time2-time1
    
    call del_box2d(box1)
    call del_box2d(box2)
    call del_mes()

end program mvpimd

