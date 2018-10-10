!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- MES Models (Copyright Xin He <1500011805@pku.edu.cn>)
!-- It contains:
!------ 0) Double electronic states coulping to a single Morse potential
!------ 1) Seven electronic states coupling to a single Morse potential
!------ 2) Double electronic state -- Spin Boson Model coupling to a set 
!------------ of oscillitor potentials
!------ 3) Seven electronic states coupling to a set of Morse potentials
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


module MES_Models
use base
use linalgebra
implicit none
    !-- setting related to MD parameters
    real(8), private :: betap ! betap = beta / P
    integer, private :: nbead ! <also noted as P>
    
    !-- general settings
    integer, private :: mod_flag
    integer, private :: rdim, vdim

    !-- parameters of Morse potentials & Harmonic oscillitor
    real(8), dimension(:), allocatable      :: Mod_M, Mod_V0
    real(8), dimension(:,:), allocatable      :: Mod_D, Mod_W, Mod_Req, Mod_C, Mod_R
    !-- for the same coupling strength, the following works as: <or to vectorize them>
    real(8) :: Mod_a1, Mod_a2, Mod_delta
    
    !-- Recorders related with model
    real(8), dimension(:,:,:), allocatable   :: Mod_V, Mod_Vd, Mod_Vnd
    real(8), dimension(:,:,:,:), allocatable :: Mod_dVd, Mod_dVnd
    real(8), dimension(:,:), allocatable :: update_rho
    
    !-- Recorders irrelated to model
    real(8) :: update_pfd, update_pff, update_v, update_kprim, update_kvir, update_eprim, update_evir, update_coh
    real(8),dimension(5) :: update_theta

contains
    function getrdim() result (d)
    implicit none
        integer :: d
        d = rdim
    end function getrdim
    
    function getvdim() result (d)
    implicit none
        integer :: d
        d = vdim
    end function getvdim
    
    subroutine init_Model(mduse_bead, mduse_beta, imodel)
    implicit none
        integer, intent(in) :: mduse_bead, imodel
        real(8), intent(in) :: mduse_beta
        nbead = mduse_bead
        betap = mduse_beta / real(nbead)
        mod_flag = imodel
        select case (mod_flag)
            case (0)    !-- Double electronic states coulping to a single Morse potential
                rdim = 1
                vdim = 2
            case (1)    !-- Seven electronic states coupling to a single Morse potential
                rdim = 1
                vdim = 7
            case (2)    !-- Spin Boson Model coupling to a set of HO potentials
                rdim = 100
                vdim = 2
            case (3)    !-- Seven electronic states coupling to a set of Morse potentials
                rdim = 50
                vdim = 7
            case default    !-- 1-dimensional harmonic oscillitor potential
                rdim = 1
                vdim = 1
                stop "no model match!"
        endselect
        allocate(Mod_M(rdim))
        allocate(Mod_V0(vdim))
        allocate(Mod_D(rdim, vdim))
        allocate(Mod_W(rdim, vdim))
        allocate(Mod_Req(rdim, vdim))
        allocate(Mod_C(vdim, vdim))
        allocate(Mod_R(vdim, vdim))
        !-- Recorders as:
        allocate(Mod_Vd(nbead, vdim, vdim))
        allocate(Mod_Vnd(nbead, vdim, vdim))
        allocate(Mod_dVd(rdim, nbead, vdim, vdim))
        allocate(Mod_dVnd(rdim, nbead, vdim, vdim))
        select case (mod_flag)
            case (0)    !-- Double electronic states coulping to a single Morse potential
                call init_ES2_morse()
            case (1)    !-- Seven electronic states coupling to a single Morse potential
                call init_ES7_morse()
            case (2)    !-- Spin Boson Model coupling to a set of HO potentials
                call init_SB2_multiharmonic()
            case (3)    !-- Seven electronic states coupling to a set of Morse potentials
                call init_ES7_multimorse()
            case default    !-- 1-dimensional harmonic oscillitor potential
                stop "no model match!"
        endselect
    end subroutine init_Model
    
    
    subroutine Mod_getV(V, Vd, Vnd, R, iflag)
    implicit none
        real(8), dimension(rdim), intent(in) :: R
        integer, intent(in), optional :: iflag
        real(8), dimension(vdim,vdim), intent(inout) :: V, Vd, Vnd
        select case (mod_flag)
            case (0)
                call Mod0_getV(V, Vd, Vnd, R, iflag)
            case (1)
                call Mod1_getV(V, Vd, Vnd, R, iflag)
            case (2)
                call Mod2_getV(V, Vd, Vnd, R, iflag)
            case (3)
                call Mod3_getV(V, Vd, Vnd, R, iflag)
            case default
                stop "wrong modeling"
        endselect
    end subroutine Mod_getV
    
    
    subroutine Mod_getdV(dV, dVd, dVnd, R, iflag)
    implicit none
        real(8), dimension(rdim), intent(in) :: R
        integer, intent(in), optional :: iflag
        real(8), dimension(vdim,vdim), intent(inout) :: dV, dVd, dVnd
        select case (mod_flag)
            case (0)
                call Mod0_getdV(dV, dVd, dVnd, R, iflag)
            case (1)
                call Mod1_getdV(dV, dVd, dVnd, R, iflag)
            case (2)
                call Mod2_getdV(dV, dVd, dVnd, R, iflag)
            case (3)
                call Mod3_getdV(dV, dVd, dVnd, R, iflag)
            case default
                stop "wrong modeling"
        endselect
    end subroutine Mod_getdV
    
    
    !-----------------------------------------------------------------------------------------------------------------------------------
    subroutine Mod_getPFD(Y, Rs)
    implicit none
        real(8), dimension(rdim, nbead), intent(in) :: Rs
        real(8), dimension(vdim, vdim), intent(out) :: Y
        real(8), dimension(vdim, vdim) :: V,Vd,Vnd
        integer :: i
        Y = lin_eye(vdim)
        do i=1, size(Rs,dim=2)
            call Mod_getV(V,Vd,Vnd,Rs(:,i),iflag=0)        !-- flag==1, only calculate the diagonal force, not here
            Y = matmul(Y, lin_expdiag(-betap*Vd(:,:) ))
            !-- we save all the V matrix ---------------------------------------------------------------- [IMPORTANT]
            Mod_Vd(i,:,:) = Vd
            Mod_Vnd(i,:,:) = Vnd !----------------------------------------------------------------------------------------
        enddo
        !-- after calculate PFD, then we save its trace -------------------------------------------------- [IMPORTANT]
        update_pfd = lin_trace(Y)
    end subroutine Mod_getPFD
    
    
    !-- update the evaluated force, for each step of md running
    subroutine Mod_esti_dV(esval, Rs)
    implicit none
        real(8), dimension(rdim, nbead), intent(in) :: Rs
        real(8), dimension(vdim,vdim) :: PFD_matrix
        real(8), dimension(rdim,vdim,vdim) :: dV, dVd, dVnd
        real(8), dimension(rdim, nbead), intent(out) :: esval
        integer :: i,j
        call Mod_getPFD(PFD_matrix, Rs)     
        do i=1, size(Rs,dim=2)
            !-- for each indenpent frame of system (a set of bead), we calculate a set of force
            call Mod_getdV(dV, dVd, dVnd, Rs(:,i), iflag=0) !-- flag==1, only calculate the diagonal force, not here
            !-- save------------------         ------------------------------------------------------------ [IMPORTANT]
            Mod_dVd(:,i,:,:) = dVd
            Mod_dVnd(:,i,:,:) = dVnd
            !-------------------------
            do j=1, size(Rs,dim=1)
                esval(j,i) = lin_trace( matmul(dVd(j,:,:), PFD_matrix) ) / update_pfd
            enddo
        enddo
        !-- so by different  definition, here needn't to be divided by bead number, ref. < md_pimd: fx2fks, line 94 >
    end subroutine Mod_esti_dV
    
    
    subroutine Mod_esti_PFF(esval, Rs)
    implicit none
        real(8), dimension(rdim, nbead), intent(in) :: Rs
        real(8), dimension(vdim,vdim) :: PFD_matrix, V, Vd, Vnd, EVd, EVnd , O, Y
        real(8), intent(out) :: esval
        integer :: i
        !-- we needn't re-call PFD
        !--call Mod_getPFD(PFD_matrix, Rs)
        Y = lin_eye(vdim)
        do i=1, size(Rs,dim=2)
            !-- SAVEFORIME
            !call Mod_getV(V,Vd,Vnd,Rs(:,i),iflag=0)
            Vd = Mod_Vd(i,:,:)
            Vnd = Mod_Vnd(i,:,:)
            !----------------------------------------
            EVd = lin_expdiag(-0.5*betap*Vd)
            EVnd = (lin_eye(vdim) - 0.5*betap*Vnd)
            O = matmul( EVnd, EVd )
            Y = matmul(Y,  matmul(transpose(O), O) )
        enddo
        !esval = lin_trace(Y) / lin_trace(PFD_matrix)
        esval = lin_trace(Y) / update_pfd
        !-- save data pff, its important in other properties evaluation ------------------------------ [IMPORTANT]
        update_pff = esval
    end subroutine Mod_esti_PFF
    
    !-- The evaluation method of V can be generalize to all nuclear properties --- relating all of frames of R
    subroutine Mod_esti_V(esval, Rs)
    implicit none
        real(8), dimension(rdim, nbead), intent(in) :: Rs
        real(8), dimension(vdim, vdim) :: Y, V, Vd, Vnd, EVd, EVnd , O, OO
        real(8), intent(out) :: esval
        integer :: i,j

        esval = 0.0
        do i=1, size(Rs,dim=2)
            Y = lin_eye(vdim)
            do j=1, size(Rs,dim=2)
                !-- May not re-call
                !call Mod_getV(V,Vd,Vnd,Rs(:,j),iflag=0)
                !-- Or by SAVEFORTIME
                Vd = Mod_Vd(j,:,:)
                Vnd = Mod_Vnd(j,:,:)
                V = Vd + Vnd
                !----------------------------------------
                EVd = lin_expdiag(-0.5*betap*Vd)
                EVnd = (lin_eye(vdim)-0.5*betap*Vnd)
                O = matmul( EVnd, EVd )
                OO = matmul(transpose(O), O)
                if(j/= i) then
                    Y = matmul(Y,  OO)
                else
                    !call Mod_esti_B(B, Rs(:,j))
                    !-- estimate of V
                    Y = matmul(Y, matmul( transpose(O), matmul(V, O) ) )
                endif
            enddo
            esval = esval  + lin_trace(Y)
        enddo
        esval = esval /(real(nbead) * update_pfd) !-- HAVE BEEN SAVED
    end subroutine Mod_esti_V
    
    
    subroutine Mod_getOpR(OpR, R, iframe, iflag)             !-- flag=0 (first order, default), flag=1(hyperbolic), flag=2 diagonalize
        real(8), dimension(rdim), intent(in) :: R  !-- must at one frame!
        real(8), dimension(rdim,vdim,vdim) :: dV,dVd,dVnd
        real(8), dimension(rdim,vdim,vdim), intent(out) :: OpR
        real(8), dimension(vdim,vdim) :: V,Vd,Vnd
        integer, intent(in), optional :: iframe, iflag
        integer :: i,j,frame,flag
        if(present(iflag)) then
            flag = iflag
        else
            flag = 0 !-- default value, first order method
        endif
        if(present(iframe)) then
            frame = iframe
        else
            frame = 0 !-- default value, for this occusion, should re-calc the V and dV
        endif
        
        if(flag .eq. 0) then
            if(frame .eq. 0) then
                call Mod_getV(V, Vd, Vnd, R, iflag=0)
                call Mod_getdV(dV, dVd, dVnd, R, iflag=0)
            else
                !-- SAVEFORTIME
                Vd = Mod_Vd(frame,:,:)
                Vnd = Mod_Vnd(frame,:,:)
                dVd = Mod_dVd(:,frame,:,:)
                dVnd = Mod_dVnd(:,frame,:,:)
            endif
            do j=1,size(R)  !---------------------------------------------------------------------------------------------------------------- BUGS
                OpR(j,:,:) = -0.5*betap*dVnd(j,:,:)*lin_expdiag(-0.5*betap*Vd)  + &
                            ( lin_eye(vdim) - 0.5*betap*Vnd)*lin_expdiag(-0.5*betap* Vd(:,:)) * (-0.5*betap* dVd(j,:,:))
            enddo
            return
        elseif(flag .eq. 1) then
            call log_message('E', "flag error")
        elseif(flag .eq. 2) then
            call log_message('E', "flag error")
        else
            call log_message('E', "flag error")
        endif
    end subroutine Mod_getOpR
    
    subroutine Mod_esti_K(kprim, kvir, Rs, iflag)
    implicit none
        real(8), dimension(rdim, nbead), intent(in) :: Rs
        real(8), dimension(vdim,vdim) :: PFD_matrix, V, Vd, EVd, EVnd, Vnd, dV, dVd, dVnd, O, OO, tmp
        real(8), dimension(rdim,vdim,vdim) :: dOdR, Y
        real(8), intent(out) :: kprim, kvir
        real(8), dimension(rdim) :: Rc, pR
        real(8) :: S, term, PFF, pfd_trace
        integer, intent(in), optional :: iflag
        integer :: i,j,k,flag
        if(present(iflag)) then
            flag = iflag
        else
            flag = 1    !-- default additional virial kenatics
        endif
        
        !-- SAVEFORTIME
        !--call Mod_esti_PFF(PFF, Rs)

        S = 0.5 * size(Rs,dim=1) / betap
        do i=1, size(Rs,dim=2)
            term = 0.0 !------------------------------------------------------------------------------------------------------------bugs
            do j=1, size(Rs, dim=1)   
                term = term + Mod_M(j) * ( Rs(j,i) - Rs(j, mod(i,size(Rs,dim=2)) + 1 ))**2
            enddo
            S = S - 0.5/(nbead*betap*betap) *term !---------------------------------------------------------- om2 ???
        enddo
        kprim = S*update_pff 
        !--
        if(flag > 0) then !--- additional virial estimator
            Rc = sum(Rs,dim=2) / real(nbead)
            !call Mod_getPFD(PFD_matrix, Rs)
            !pfd_trace = lin_trace(PFD_matrix)
            term = 0.0
            do i=1, size(Rs,dim=2)
                !-- for reuse date, we use call Mod_getOpR(dOdR, Rs(:,i), [optional]iframe =i ,iflag=0), 
                !-------------- instanding call Mod_getOpR(dOdR, Rs(:,i),iflag=0)
                call Mod_getOpR(dOdR, Rs(:,i), iframe=i, iflag=0)  
                do j=1, size(Rs, dim=1)
                    Y(j,:,:) = lin_eye(vdim)
                enddo
                do k=1, size(Rs,dim=2)
                    !-- SAVEFORTIME
                    !call Mod_getV(V,Vd,Vnd,Rs(:,k),iflag=0)
                    Vd = Mod_Vd(k,:,:)
                    Vnd= Mod_Vnd(k,:,:)
                    !------------------------------------------
                    EVd = lin_expdiag(-0.5*betap*Vd)
                    EVnd = (lin_eye(vdim) - 0.5*betap*Vnd)
                    O = matmul( EVnd, EVd )
                    OO = matmul( transpose(O), O)
                    
                    if(i.eq. k) then
                        do j=1,size(Rs,dim=1)
                            tmp = matmul( transpose(dOdR(j,:,:)), O) + matmul(transpose(O), dOdR(j,:,:))
                            Y(j,:,:) = matmul(Y(j,:,:), tmp )
                        enddo                  
                    else
                        do j=1,size(Rs,dim=1)
                            Y(j,:,:) = matmul(Y(j,:,:), OO )
                        enddo
                    endif
                enddo
                
                do j=1,size(Rs,dim=1)
                    pR(j) = lin_trace(Y(j,:,:))
                    term = term + (Rs(j,i) - Rc(j)) * pR(j)
                enddo
                
            enddo
            kvir = 0.5_8/(nbead*betap) * ( real(size(Rs,dim=1)) * update_pff - term/update_pfd )
        endif
        return
    end subroutine Mod_esti_K
    
    
    subroutine Mod_esti_theta(Rs)
    implicit none
        real(8), dimension(rdim, vdim), intent(in) :: Rs
        real(8), dimension(vdim,vdim) :: O, OO, Y1, Y2, Y3, Y4, RpR, Vd, Vnd, EVd, EVnd, V_all
        real(8), dimension(rdim,vdim,vdim) :: dV_all, dOdR, dTdR
        real(8) :: T1, T2, T3, T4
        real(8), dimension(rdim) :: Rc, DeltaR
        integer :: i,j,k,m
        
        T1 = 0.0
        T2 = 0.0
        T3 = 0.0
        T4 = 0.0
        Rc = sum(Rs,dim=2) / nbead
        do j=1,nbead
            do k=1,nbead                    !-- from k=j, if k==j, calc Y2, else calc Y1, then note Y1 we count double times !!
                Y1 = lin_eye(vdim)
                Y2 = lin_eye(vdim)
                Y3 = lin_eye(vdim)     !-- note we use Y3,T3 evaluate all theta3, theta4, theta5 term ??
                Y4 = lin_eye(vdim)
                do i=1,nbead
                    Vd = Mod_Vd(i,:,:)
                    Vnd = Mod_Vnd(i,:,:)
                    EVd = lin_expdiag(-0.5*betap*Vd)
                    EVnd = (lin_eye(vdim) - 0.5*betap*Vnd)
                    O = matmul( EVnd, EVd )
                    OO = matmul(transpose(O), O)
                    
                    if(i/=j .and. i/=k .and. j/=k) then
                        Y1 = matmul(Y1, OO)
                        Y3 = matmul(Y3, OO)
                    elseif(i.eq.j .and. j/=k) then
                        V_all = Vd + Vnd
                        Y1 = matmul(Y1, matmul(transpose(O), matmul(V_all, O) ) )
                        !--
                        call Mod_getOpR(dOdR, Rs(:,i), iframe=i, iflag=0)
                        do m=1,size(Rs,dim=1)
                            dTdR(m,:,:) = matmul(transpose(dOdR(m,:,:)), O) + matmul(transpose(O), dOdR(m,:,:))
                        enddo
                        DeltaR = Rs(:,i)-Rc
                        RpR = 0.0
                        do m =1, size(Rs,dim=1)
                            RpR = RpR + DeltaR(m)*dTdR(m,:,:)       !-- size(rdim) (dot) size(rdim,vdim,vdim)
                        enddo
                        Y3 = matmul(Y3, RpR)
                    elseif(i.eq.k .and. j/=k) then
                        V_all = Vd + Vnd
                        !-- for Y1, here we note we count it doubles times, so it will without times 2 in the final result
                        Y1 = matmul(Y1, matmul(transpose(O), matmul(V_all, O) ) )
                        !--
                        Y3 = matmul(Y3, matmul(transpose(O), matmul(V_all, O) ) )
                    elseif(j.eq.k .and. j/=i) then
                        Y2 = matmul(Y2, OO )
                        Y4 = matmul(Y4, OO )                      
                    else !--(i.eq.j .and. i.eq.k)
                        V_all = Vd + Vnd
                        dV_all = Mod_dVd(:,i,:,:) + Mod_dVnd(:,i,:,:)
                        Y2 = matmul(Y2, matmul( transpose(O), matmul(V_all*V_all, O) ) )
                        !--
                        call Mod_getOpR(dOdR, Rs(:,i), iframe=i, iflag=0)
                        do m=1,size(Rs,dim=1)
                            dTdR(m,:,:) = matmul(transpose(dOdR(m,:,:)), matmul(V_all,O))  &
                                        + matmul(transpose(O), matmul(V_all, dOdR(m,:,:)) ) &
                                        + matmul(transpose(O), matmul(dV_all(m,:,:), O ) )
                        enddo
                        DeltaR = Rs(:,i)-Rc
                        RpR = 0.0
                        do m =1, size(Rs,dim=1)
                            RpR = RpR + DeltaR(m)*dTdR(m,:,:)       !-- size(rdim) (dot) size(rdim,vdim,vdim)
                        enddo
                        Y4 = matmul(Y4, RpR)
                    endif
                enddo
                !--
                if(j .eq. k) then
                    T2 = T2 + lin_trace(Y2)
                    T4 = T4 + lin_trace(Y4)
                else
                    T1 = T1 + lin_trace(Y1)
                    T3 = T3 + lin_trace(Y3) 
                endif
            enddo
        enddo
        update_theta(1) = T1 / real(nbead**2) !-------------------- NO FRACTOR 2
        update_theta(2) = T2 / real(nbead**2)
        update_theta(3) = - 0.5 * T3 / real(nbead**2 * betap)
        update_theta(4) = - 0.5 * T4 / real(nbead**2 * betap)
    end subroutine Mod_esti_theta
    
    
    subroutine Mod_esti_coh(esval)
    implicit none
        real(8), dimension(vdim, vdim) :: Rho, Rij, Y, O, OO, Vd, Vnd, EVd, EVnd
        real(8), intent(out) :: esval
        real(8) :: s1, s2
        integer :: i,j,k,l
        esval = 0.0
        Rho = 0.0       !----------------------------------------------------------------------------------------------------------- BUGs
        do i=1, vdim
            do j=i, vdim
                Rij = 0.0
                Rij(i,j) = Rij(i,j) + 0.5
                Rij(j,i) = Rij(j,i) + 0.5
                do k=1, nbead
                    Y = lin_eye(vdim)
                    do l=1, nbead
                        !-- May not re-call
                        !call Mod_getV(V,Vd,Vnd,Rs(:,j),iflag=0)
                        !-- Or by SAVEFORTIME
                        Vd = Mod_Vd(l,:,:)
                        Vnd = Mod_Vnd(l,:,:)
                        !----------------------------------------
                        EVd = lin_expdiag(-0.5*betap*Vd)
                        EVnd = (lin_eye(vdim) - 0.5*betap*Vnd)
                        O = matmul( EVnd, EVd )
                        OO = matmul(transpose(O), O )
                        if(l/= k) then
                            Y = matmul(Y,  OO )
                        else
                            !-- estimate of Rhoij
                            Y = matmul(Y, matmul( transpose(O), matmul(Rij, O) ) )
                        endif
                    enddo
                    Rho(i,j) = Rho(i,j) + lin_trace(Y)
                    if(i/=j) Rho(j,i) = Rho(j,i) + lin_trace(Y) 
                enddo
            enddo
        enddo
        Rho = Rho / real(nbead)
        s1 = 0.0
        s2 = 0.0
        do i=1, vdim
            do j=1, vdim
                s1 = s1 + Rho(i,j)**2
                s2 = s2 + abs( Rho(i,j) )
            enddo
        enddo
        esval = s2**2/(vdim*s1)
        update_coh = esval !------------------------------------------------- SAVE
        return
    end subroutine Mod_esti_coh
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !-- MODEL 0 : 2 ELECTRONIC STATES WITH MORSE POTENTIAL 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !-- When mod_flag == 0
    subroutine init_ES2_morse()
    implicit none
        
        Mod_M(1) = 1.0
        Mod_V0 = (/0.0, 9.8e-5/)
        Mod_D(1,:) = (/4.71e-3, 4.71e-3/)
        Mod_W(1,:) = (/9.43e-5, 9.0e-5/)
        Mod_Req(1,:) = (/-1.75, 0.9226/)
        Mod_C = 0.0
        Mod_C(1,2) = 6.11e-5
        Mod_C(2,1) = 6.11e-5
        Mod_R = 0.0
        Mod_R(1,2) = 0.923
        Mod_R(2,1) = 0.923
        Mod_a1 = 0.05
        Mod_a2 = 0.05
        
    end subroutine init_ES2_morse
    
    subroutine Mod0_getV(V, Vd, Vnd, R, iflag)
    implicit none
        real(8), dimension(rdim), intent(in) :: R
        integer, intent(in), optional :: iflag
        real(8), dimension(vdim,vdim), intent(inout) :: V,Vd,Vnd
        integer :: i,j, flag
        !-- for 3N=1, here is only one freedom of nuclues
        if(present(iflag)) then
            flag = iflag
        else
            flag = 0 !-- default value
        endif
        !-- calculation MES-PES
        Vd =0.0
        Vnd=0.0
        do i=1,vdim
            V(i,i) = Mod_D(1,i)*(1.0-dexp(-dsqrt(0.5*Mod_M(1)*Mod_W(1,i)**2/Mod_D(1,i))*(R(1)-Mod_Req(1,i))) )**2 + & !----------------------- 1/=1
                Mod_V0(i)
            Vd(i,i)=V(i,i)
            if(flag .eq. 1) cycle !-- flag==1, then only calculate diagonal term
            do j=i+1,vdim
                V(i,j) = Mod_C(i,j) * dexp( - Mod_a1*(R(1)-Mod_R(i,j))**2 ) * dcos(Mod_a2*(R(1)-Mod_R(i,j)))
                V(j,i) = V(i,j)
                Vnd(i,j) = V(i,j)
                Vnd(j,i) = V(i,j)
            enddo
        enddo
    end subroutine Mod0_getV
    
    subroutine Mod0_getdV(dV, dVd, dVnd, R, iflag)
    implicit none
        real(8), dimension(rdim), intent(in) :: R
        real(8), dimension(rdim, vdim, vdim), intent(out) :: dV, dVd, dVnd
        integer, intent(in), optional :: iflag
        integer :: i,j,k,flag
        real(8) :: eterm
        if(present(iflag)) then
            flag = iflag
        else
            flag = 0 !-- default value
        endif
        dV = 0.0
        dVd = 0.0
        dVnd = 0.0
        do j=1,rdim
            do i=1,vdim
                eterm = dexp(-dsqrt(0.5*Mod_M(1)*Mod_W(1,i)**2/Mod_D(1,i))*(R(1)-Mod_Req(1,i)))
                dV(j,i,i) = (1.0-eterm) * eterm * dsqrt(2.0 *Mod_M(1)*Mod_W(1,i)**2 * Mod_D(1,i)) !--------------------------------------- 1/=1
                dVd(j,i,i) = dV(j,i,i) 
            enddo
            if (flag.eq. 1) cycle!-- flag==1, only give the dv_diangonal term
            do k=i+1, vdim
                dV(j,i,k) = - Mod_C(i,k) * dexp( - Mod_a1*(R(1)-Mod_R(i,k))**2 ) * &           !------------------------------------- bugs
                   ( ( 2.0*Mod_a1*(R(1)-Mod_R(i,k)) )* dcos(Mod_a2*(R(1)-Mod_R(i,k)))  &
                    + dsin(Mod_a2*(R(1)-Mod_R(i,k))) * Mod_a2 )
                dV(j,k,i) = dV(j,i,k)
                dVnd(j,i,k) = dV(j,i,k)
                dVnd(j,k,i) = dV(j,i,k)
            enddo
        enddo
    end subroutine Mod0_getdV
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !-- MODEL 1 : 7 ELECTRONIC STATES WITH MORSE POTENTIAL 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine init_ES7_morse()
    implicit none
        integer :: i,j
        real(8), dimension(7,7) :: HFMO, RFMO
        
        !-- setting the system parameters
        !-- data ref. [66] M. H. Cho, H. M. Vaswani, T. Brixner, J. stenger, and G. R. Fleming
        !------ J. Phys. Chem. B 109(21), 10542-10556(2005)
        data ((HFMO(i,j),j=1,7),i=1,7) / 0., -62., 17., 8., -1., -9., 28., &
                                        -62., 175., -57., -5., -70., -19., 6., &
                                        17., -57., 260., -4., -2., 32., 1., &
                                        8., -5., -4., 280., 6., -8., -106., &
                                        -1., -70., -2., 6., 320., 40., 2., &
                                        -9., -19., 32., -8., 40., 360., 13., &
                                        28., 6., 1., -106., 2., 13., 420. /
        !-- convert to atom unit
        HFMO = HFMO / wn_2_hetree

        data ((RFMO(i,j),j=1,7),i=1,7) / 0., 224., 157., 122., 99., 93., 94., &
                                        224., 0., 102., 67., 65., 66., 72., &
                                        157., 102., 0., 34., 48., 54., 65., &
                                        112., 67., 34., 0., 62., 65., 76., &
                                        99., 65., 48., 62., 0., 69., 84., &
                                        93., 66., 54., 65., 69., 0., 100., &
                                        94., 72., 65., 76., 84., 100., 0. /
        !-- initialize parameters
        Mod_M(1) = 1.0
        do i=1,vdim
            Mod_V0(i) = HFMO(i,i)
            Mod_R(i,i) = 0.0
            do j=i+1,vdim
                Mod_C(i,j) = HFMO(i,j)
                Mod_C(j,i) = HFMO(i,j)
                Mod_R(i,j) = RFMO(i,j)
                Mod_R(j,i) = RFMO(i,j)
            enddo
        enddo
        Mod_D(1,:)   = (/7.28e-2, 7.17e-2, 7.06e-2, 6.95e-2, 6.84e-2, 6.73e-2, 6.62e-2/)
        Mod_W(1,:)   = (/212.2, 209.0, 205.8, 202.7, 199.5, 196.3, 193.1/) / wn_2_hetree
        Mod_Req(1,:) = (/0., 5., 10., 15., 20., 25., 30./)
        Mod_a1 = 5.0e-5
        Mod_a2 = 0.02
    end subroutine init_ES7_morse

    subroutine Mod1_getV(V, Vd, Vnd, R, iflag)
    implicit none
        real(8), dimension(rdim), intent(in) :: R
        integer, intent(in), optional :: iflag
        real(8), dimension(vdim,vdim), intent(inout) :: V, Vd, Vnd
        integer :: i,j, flag
        !-- for 3N=1
        if(present(iflag)) then
            flag = iflag
        else
            flag = 0 !-- default value
        endif
        !-- calculation MES-PES
        Vd  = 0.0
        Vnd = 0.0
        do i=1,vdim
            V(i,i) = Mod_D(1,i)*(1.0-dexp(-dsqrt(0.5*Mod_M(1)*Mod_W(1,i)**2/Mod_D(1,i))*(R(1)-Mod_Req(1,i))) )**2 &
                + Mod_V0(i)
            Vd(i,i)=V(i,i)
            if(flag .eq. 1) cycle !-- flag==1, then only calculate diagonal term
            do j=i+1,vdim
                V(i,j) = Mod_C(i,j) * dexp( - Mod_a1*(R(1)-Mod_R(i,j))**2 ) * dcos(Mod_a2*(R(1)-Mod_R(i,j)))
                V(j,i) = V(i,j)
                Vnd(i,j) = V(i,j)
                Vnd(j,i) = V(i,j)
            enddo
        enddo
    end subroutine Mod1_getV

    subroutine Mod1_getdV(dV, dVd, dVnd, R, iflag)
    implicit none
        real(8), dimension(rdim), intent(in) :: R
        real(8), dimension(rdim, vdim, vdim), intent(out) :: dV, dVd, dVnd
        integer, intent(in), optional :: iflag
        integer :: i,j,k,flag
        real(8) :: eterm
        if(present(iflag)) then
            flag = iflag
        else
            flag = 0 !-- default value
        endif
        dV = 0.0
        dVd = 0.0
        dVnd = 0.0
        do j=1,rdim
            do i=1,vdim
                eterm = dexp(-dsqrt(0.5*Mod_M(1)*Mod_W(1,i)**2/Mod_D(1,i))*(R(1)-Mod_Req(1,i)))
                dV(j,i,i) = (1.0-eterm) * eterm * dsqrt(2.0 *Mod_M(1)*Mod_W(1,i)**2 * Mod_D(1,i))
                dVd(j,i,i) = dV(j,i,i) 
            enddo
            if (flag.eq. 1) cycle!-- flag==1, only give the dv_diangonal term
            do k=i+1, vdim
                dV(j,i,k) = - Mod_C(i,k) * dexp( - Mod_a1*(R(1)-Mod_R(i,k))**2 ) * &
                    ( ( 2.0*Mod_a1*(R(1)-Mod_R(i,k)) )* dcos(Mod_a2*(R(1)-Mod_R(i,k))) & !----------------------------------------
                    + dsin(Mod_a2*(R(1)-Mod_R(i,k))) * Mod_a2 )
                dV(j,k,i) = dV(j,i,k)
                dVnd(j,i,k) = dV(j,i,k)
                dVnd(j,k,i) = dV(j,i,k)
            enddo
        enddo
    end subroutine Mod1_getdV


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !-- MODEL 2 : SPIN BOSON MODEL
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !-- When mod_flag == 2
    subroutine init_SB2_multiharmonic()
        integer :: i,j
        !-- prev initialize the parameters
        real(8) :: delta, epsilonau, beta, omegac, alpha
        real(8), dimension(rdim) :: omegaj, cj

        delta = 0.1
        epsilonau = 0.0
        beta = 1
        omegac = 0.25
        alpha = 0.09
        !-- for Model I :   / 0.1, 0, 1,  0.25 , 0.09/
        !-- for Model II :  / 0.1, 0, 50, 0.25 , 0.09/
        !-- for Model III : / 0.1, 1, 50, 0.25 , 0.1 /
        !-- for Model IV :  / 0.1, 1, 2.5,0.1 ,  0.1 /
        !-- for Model V :   / 0.1, 0, 10, 0.1 ,  2   /
        !-- for Model VI :  / 0.1, 5, 1,  0.2 ,  4   /

        Mod_delta = delta
        do i=1,rdim
            Mod_M(i) = 1.0
            omegaj(i) = - omegac*dlog(1._8-real(i)/(rdim+1))
            cj(i) = dsqrt(alpha*omegaj(i)**2*omegac/real(rdim+1))
        enddo
        do i=1,vdim
            Mod_V0(i) = (-1)**(i+1)*epsilonau
            do j=1,rdim
                Mod_W(j,i) = omegaj(j)
                if(mod(i,2).eq. 0) then
                    Mod_Req(j,i) = - cj(j)/(Mod_M(j)*omegaj(j)**2)
                else
                    Mod_Req(j,i) = cj(j)/(Mod_M(j)*omegaj(j)**2)
                endif
                Mod_V0(i) = Mod_V0(i) - 0.5 * cj(j)/(Mod_M(j)*omegaj(j)**2)
            enddo
        enddo
    end subroutine init_SB2_multiharmonic


    subroutine  Mod2_getV(V, Vd, Vnd, R, iflag)
        real(8), dimension(rdim), intent(in) :: R
        integer, intent(in), optional :: iflag
        real(8), dimension(vdim,vdim), intent(inout) :: V, Vd, Vnd
        integer :: i,j,k, flag
        if(present(iflag)) then
            flag = iflag
        else
            flag = 0
        endif
        Vd =0.0
        Vnd=0.0
        do i=1,vdim
            V(i,i) = Mod_V0(i)
            do j=1,rdim
                V(i,i) = V(i,i) + 0.5*Mod_M(j)*( Mod_W(j,i) * (R(j)-Mod_Req(j,i)) )**2
            enddo
            Vd(i,i)=V(i,i)
            if(flag .eq. 1) cycle !-- flag==1, then only calculate diagonal term
            do k=i+1,vdim
                V(i,k) = Mod_delta
                V(k,i) = V(i,k)
                Vnd(i,k) = V(i,k)
                Vnd(k,i) = V(i,k)
            enddo
        enddo
    end subroutine Mod2_getV
    
    subroutine  Mod2_getdV(V, Vd, Vnd, R, iflag)
    implicit none
        real(8), dimension(rdim), intent(in) :: R
        integer, intent(in), optional :: iflag
        real(8), dimension(vdim,vdim), intent(inout) :: V, Vd, Vnd
    end subroutine Mod2_getdV


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !-- MODEL 3 : 7 ELECTRONIC STATES WITH MULTI-MORSE POTENTIAL 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine init_ES7_multimorse()
    implicit none
        real(8) :: omegac, eta
        real(8), dimension(7,7) :: HFMO, RFMO
        integer :: i,j
        
        omegac = 212.2 / wn_2_hetree
        eta = 35.0 / wn_2_hetree
        !-- data ref. [66] M. H. Cho, H. M. Vaswani, T. Brixner, J. stenger, and G. R. Fleming
        !------ J. Phys. Chem. B 109(21), 10542-10556(2005)
        data ((HFMO(i,j),j=1,7),i=1,7) / 0., -62., 17., 8., -1., -9., 28., &
                                        -62., 175., -57., -5., -70., -19., 6., &
                                        17., -57., 260., -4., -2., 32., 1., &
                                        8., -5., -4., 280., 6., -8., -106., &
                                        -1., -70., -2., 6., 320., 40., 2., &
                                        -9., -19., 32., -8., 40., 360., 13., &
                                        28., 6., 1., -106., 2., 13., 420. /
        !-- convert to atom unit
        HFMO = HFMO / wn_2_hetree

        data ((RFMO(i,j),j=1,7),i=1,7) / 0., 224., 157., 122., 99., 93., 94., &
                                        224., 0., 102., 67., 65., 66., 72., &
                                        157., 102., 0., 34., 48., 54., 65., &
                                        112., 67., 34., 0., 62., 65., 76., &
                                        99., 65., 48., 62., 0., 69., 84., &
                                        93., 66., 54., 65., 69., 0., 100., &
                                        94., 72., 65., 76., 84., 100., 0. /
        !-- initialize parameters
        Mod_M = 1.0
        do i=1,vdim
            Mod_V0(i) = HFMO(i,i)
            Mod_R(i,i) = 0.0
            do j=i+1,vdim
                Mod_C(i,j) = HFMO(i,j)
                Mod_C(j,i) = HFMO(i,j)
                Mod_R(i,j) = RFMO(i,j)
                Mod_R(j,i) = RFMO(i,j)
            enddo
        enddo

        do i=1,vdim
            do j=1,rdim
                Mod_W(j,i) = - omegac* (1.-0.015*(i-1)) * dlog(1._8-real(j)/real(rdim+1))
                Mod_Req(j,i) = 1./(Mod_M(j)*Mod_W(j,1)) * dsqrt(2*eta/real(rdim+1)) + 5.0*(i-1)
                Mod_D(j,i) = 0.0728*( 1.-0.015*real(i-1))
            enddo
        enddo
        Mod_a1 = 5.e-5
        Mod_a2 = 0.02
    end subroutine init_ES7_multimorse

    subroutine Mod3_getV(V, Vd, Vnd, R, iflag)
    implicit none
        real(8), dimension(vdim, vdim), intent(inout) :: V, Vd, Vnd
        real(8), dimension(rdim), intent(in) :: R
        integer, intent(in), optional :: iflag
        integer :: i,j,k, flag
        if(present(iflag)) then
            flag = 1
        else
            flag = 0
        endif
        Vd =0.0
        Vnd=0.0
        do i=1,vdim
            V(i,i) = Mod_V0(i)
            do k=1,rdim
                V(i,i) = V(i,i) + Mod_D(k,i)* &
                (1.0-dexp(-dsqrt(0.5*Mod_M(k)*Mod_W(k,i)**2/Mod_D(k,i))*(R(k)-Mod_Req(k,i))) )**2
            enddo
            Vd(i,i)=V(i,i)
            if(flag .eq. 1) cycle !-- flag==1, only calculate the diagonal term
            do j=i+1,vdim
                V(i,j) = 0.0
                do k=1,rdim
                    V(i,j) = V(i,j) + Mod_C(i,j) * &
                    dexp( - Mod_a1*(R(k)-Mod_R(i,j))**2 ) * dcos(Mod_a2*(R(k)-Mod_R(i,j)))
                enddo
                V(j,i) = V(i,j)
                Vnd(i,j) = V(i,j)
                Vnd(j,i) = V(i,j)
            enddo
        enddo
    end subroutine Mod3_getV
    
    subroutine Mod3_getdV(V, Vd, Vnd, R, iflag)
    implicit none
        real(8), dimension(rdim), intent(in) :: R
        integer, intent(in), optional :: iflag
        real(8), dimension(vdim,vdim), intent(inout) :: V, Vd, Vnd
    end subroutine Mod3_getdV
    
end module MES_Models












    
    

