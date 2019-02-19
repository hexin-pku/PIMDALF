!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- 7-Elec-States coupled to Single Morse potential
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module mes7_smors
implicit none
private
    integer, parameter :: ndof = 1
    integer, parameter :: vdim = 7
    real(8), dimension(ndof) :: M
    real(8) :: a1, a2
    real(8), dimension(vdim) ::  V0, D, W, Req
    real(8), dimension(vdim,vdim) :: C, Rc
public :: M, ndof, vdim
public :: mes7_smors_read, init_mes7_smors, getV, getdV
contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- readparms
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine mes7_smors_read(parmsfile)
character(*), intent(in) :: parmsfile
logical :: my_exist
    inquire(file=parmsfile, exist=my_exist ) 
    if(.not. my_exist ) then
        print *, "warning: parms loss, using default settings"
        call init_mes7_smors()
        return
    endif
    open(unit=20, file=parmsfile, status="old")
    read(20,*)
    read(20,*) V0
    read(20,*)
    read(20,*) D
    read(20,*)
    read(20,*) W
    read(20,*)
    read(20,*) Req
    read(20,*)
    read(20,*) C
    read(20,*)
    read(20,*) Rc
    read(20,*)
    read(20,*) a1, a2
    close(unit=20)
    !print *, M, V0, D, W, Req, C, Rc, a1,a2
end subroutine mes7_smors_read


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- initialization
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine init_mes7_smors()
implicit none
    integer :: i,j
    real(8), dimension(vdim,vdim) :: HFMO, RFMO
    real(8), parameter :: wn_2_hetree = 219474.6313702
    
    !-- setting the system parameters
    !-- data ref. [66] M. H. Cho, H. M. Vaswani, T. Brixner, J. stenger, and G. R. Fleming
    !------ J. Phys. Chem. B 109(21), 10542-10556(2005)
    data ((HFMO(i,j),j=1,vdim),i=1,vdim) &
    / 0._8, -62._8, 17._8, 8._8, -1._8, -9._8, 28._8, &
    -62._8, 175._8, -57._8, -5._8, -70._8, -19._8, 6._8, &
    17._8, -57._8, 260._8, -4._8, -2._8, 32._8, 1._8, &
    8._8, -5._8, -4._8, 280._8, 6._8, -8._8, -106._8, &
    -1._8, -70._8, -2._8, 6._8, 320._8, 40._8, 2._8, &
    -9._8, -19._8, 32._8, -8._8, 40._8, 360._8, 13._8, &
    28._8, 6._8, 1._8, -106._8, 2._8, 13._8, 420._8 /
    !-- convert to atom unit
    HFMO = HFMO / wn_2_hetree

    data ((RFMO(i,j),j=1,vdim),i=1,vdim) &
    / 0._8, 224._8, 157._8, 122._8, 99._8, 93._8, 94._8, &
    224._8, 0._8, 102._8, 67._8, 65._8, 66._8, 72._8, &
    157._8, 102._8, 0._8, 34._8, 48._8, 54._8, 65._8, &
    112._8, 67._8, 34._8, 0._8, 62._8, 65._8, 76._8, &
    99._8, 65._8, 48._8, 62._8, 0._8, 69._8, 84._8, &
    93._8, 66._8, 54._8, 65._8, 69._8, 0._8, 100._8, &
    94._8, 72._8, 65._8, 76._8, 84._8, 100._8, 0._8 /
    !-- initialize parameters
    M = 1.0
    do i=1,vdim
        V0(i) = HFMO(i,i)
        Rc(i,i) = 0.0
        do j=i+1,vdim
            C(i,j) = HFMO(i,j)
            C(j,i) = HFMO(i,j)
            Rc(i,j) = RFMO(i,j)
            Rc(j,i) = RFMO(i,j)
        enddo
    enddo
    D(:)   = (/7.28e-2_8, 7.17e-2_8, 7.06e-2_8, 6.95e-2_8, 6.84e-2_8, 6.73e-2_8, 6.62e-2_8/)
    W(:)   = (/212.2_8, 209.0_8, 205.8_8, 202.7_8, 199.5_8, 196.3_8, 193.1_8/) / wn_2_hetree
    Req(:) = (/0.0_8, 5.0_8, 10.0_8, 15.0_8, 20.0_8, 25.0_8, 30.0_8/)
    a1 = 5.0e-5_8
    a2 = 0.02_8
end subroutine init_mes7_smors


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- getV
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getV(V, R)
implicit none
    real(8), dimension(ndof), intent(in) :: R
    real(8), dimension(vdim,vdim), intent(out) :: V
    integer :: i,j

    !-- calculation MES-PES
    V = 0.0_8
    do i=1,vdim
        print *, D(i), M(1), W(i), R(1), Req(i), V0(i)
        stop 'debug'
        
        V(i,i) = D(i)*(1.0_8-dexp(-dsqrt(0.5_8*M(1)*W(i)**2/D(i))*(R(1)-Req(i))) )**2 + V0(i)
        do j=i+1,vdim
            V(i,j) = C(i,j) * dexp( - a1*(R(1)-Rc(i,j))**2 ) * dcos(a2*(R(1)-Rc(i,j)))
            V(j,i) = V(i,j)
        enddo
    enddo
end subroutine getV


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- getdV
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getdV(dV, R)
implicit none
    real(8), dimension(ndof), intent(in) :: R
    real(8), dimension(ndof,vdim,vdim), intent(out) :: dV
    integer :: i,j
    real(8) :: eterm

    dV = 0.0
    do i=1,vdim
        eterm = dexp(-dsqrt(0.5*M(1)*W(i)**2/D(i))*(R(1)-Req(i)))
        dV(1,i,i) = (1.0-eterm) * eterm * dsqrt(2.0 *M(1)*W(i)**2 * D(i))
        do j=i+1,vdim
            dV(1,i,j) = - C(i,j) * dexp( - a1*(R(1)-Rc(i,j))**2 ) * &
               ( ( 2.0*a1*(R(1)-Rc(i,j)) )* dcos(a2*(R(1)-Rc(i,j))) & 
                + dsin(a2*(R(1)-Rc(i,j))) * a2 )
            dV(1,j,i) = dV(1,i,j)
        enddo
    enddo
end subroutine getdV

end module mes7_smors


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Test
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!program test
!use mes7_smors
!implicit none
!    real(8), dimension(vdim,vdim) :: H, Hd, Hnd, dH, dHd, dHnd 
!    real(8) :: R
!    integer :: i
!    R = 1
!    call init_mes7_smors()
!    call mes7_smors_read('mes7smors.parms')
!    do i=1,3
!        R = 0.05*i
!        call getdV(dH,dHd,dHnd,R)
!        print *, dH
!        !print *, Hd
!    enddo
!end program test

