!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! this is special for MMT simulations
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module pisimul
implicit none
private
    !-- simulation steps, sampling steps, dt
    public :: parmN, parmNs, parmdt
    !-- beta & beads number
    public :: parmB, parmP
    !-- thermostat parameter (friction)
    public :: thermo_gamma
    !-- flags
    public :: scheme_flg, thermo_flg
    !-- units
    public :: trj_unit, ana_unit, ele_unit
    
    integer :: parmN, parmNs, parmP
    real(8) :: parmdt, parmB
    real(8) :: thermo_gamma  !-- Langevin/Andersen parameter
    integer :: scheme_flg, thermo_flg
    integer, parameter :: trj_unit = 10, ana_unit=11, ele_unit = 12
    
    public :: simul_readparms
contains
subroutine simul_readparms(parms_file)
    character(*), intent(in) :: parms_file
    open(unit=22, file=trim(parms_file), status='old')
    read(22,*)
    read(22,*) parmN, parmNs, parmdt
    read(22,*)
    read(22,*) parmB, parmP
    read(22,*)
    read(22,*) thermo_gamma
    read(22,*)
    read(22,*) scheme_flg, thermo_flg
    close(unit=22)
end subroutine simul_readparms
end module pisimul
