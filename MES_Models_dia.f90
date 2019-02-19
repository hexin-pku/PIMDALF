!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- MES  routines
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


module MES_Models
use const
use linalgebra
!-- take parms from pisimul
use pisimul, nbead=>parmP, beta=>parmB
!-- XXX Interface XXX
use mes7_smors, link_V=>getV, link_dV=>getdV
implicit none
private
    real(8), private :: betap ! betap = beta / P
    real(8), private :: wp2   ! wp2   = P / beta^2
    integer, private :: savelevel = 0
    !-- Recorders related with model
    real(8), dimension(:,:,:), allocatable   :: save_Vds, save_Vnds
    real(8), dimension(:,:,:), allocatable   :: save_EVds, save_EVnds, save_Os, save_OOs
    real(8), dimension(:,:,:,:), allocatable :: save_dVds, save_dVnds
    real(8), dimension(:,:), allocatable :: save_rho

    !-- Recorders irrelated to model
    real(8) :: save_pfd, save_pfr, save_v, save_kprim, save_kvir, save_coh
    real(8), dimension(5) :: save_theta
    
    public :: save_pfd, save_pfr, save_v, save_kprim, &
              save_kvir, save_coh, save_theta, &
              init_mes, del_mes, mes_setlevel, &
              esti_pfd, esti_pfr, esti_dV, esti_V, esti_K, esti_theta, esti_coh 
    
contains


subroutine init_mes()
    allocate(save_Vds(nbead,vdim,vdim))
    allocate(save_Vnds(nbead,vdim,vdim))
    allocate(save_EVds(nbead,vdim,vdim))
    allocate(save_EVnds(nbead,vdim,vdim))
    allocate(save_Os(nbead,vdim,vdim))
    allocate(save_OOs(nbead,vdim,vdim))
    allocate(save_dVds(ndof,nbead,vdim,vdim))
    allocate(save_dVnds(ndof,nbead,vdim,vdim))
    allocate(save_rho(vdim,vdim))
    betap = beta / nbead
    wp2 = real(nbead) / beta**2
end subroutine init_mes


subroutine del_mes()
    deallocate(save_Vds)
    deallocate(save_Vnds)
    deallocate(save_EVds)
    deallocate(save_EVnds)
    deallocate(save_Os)
    deallocate(save_OOs)
    deallocate(save_dVds)
    deallocate(save_dVnds)
    deallocate(save_rho)
end subroutine del_mes

subroutine mes_setlevel(n)
    integer, intent(in) :: n
    savelevel = n
end subroutine mes_setlevel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- calc prefactor diagonal (PFD)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine esti_pfd(Y, Rs)
implicit none
    real(8), dimension(ndof, nbead), intent(in) :: Rs
    real(8), dimension(vdim, vdim), intent(out) :: Y
    real(8), dimension(vdim, vdim) :: V,Vd,Vnd
    integer :: i
    
    Y = lin_eye(vdim)
    do i=1, size(Rs,dim=2)
        call link_V(V,Rs(:,i))
        Vd = lin_diag(V)
        Vnd = V - Vd
        !-- we also save all V matrix this step
        !-- XXX
        save_Vds(i,:,:) = Vd
        save_Vnds(i,:,:) = Vnd 
        Y = matmul(Y, lin_expdiag(-betap*Vd(:,:) ))
    enddo
    save_pfd = lin_trace(Y)
    savelevel = 1
end subroutine esti_pfd


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- evaluate force, provide for md simulation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine esti_dV(esval, Rs)
implicit none
    real(8), dimension(ndof, nbead), intent(in) :: Rs
    real(8), dimension(vdim,vdim) :: PFD_matrix
    real(8), dimension(ndof,vdim,vdim) :: dV, dVd, dVnd
    real(8), dimension(ndof, nbead), intent(out) :: esval
    integer :: i,j
    
    call esti_pfd(PFD_matrix, Rs)   
      
    do i=1, size(Rs,dim=2)
        !-- for each indenpent frame of system (a set of bead)
        call link_dV(dV, Rs(:,i)) 
        do j=1, size(Rs,dim=1)
            dVd(j,:,:) = lin_diag(dV(j,:,:))
            dVnd(j,:,:) = dV(j,:,:) - dVd(j,:,:)
        enddo
        !-- we also save all dV matrix this step XXX
        save_dVds(:,i,:,:) = dVd
        save_dVnds(:,i,:,:) = dVnd     
        
        do j=1, size(Rs,dim=1)
            esval(j,i) = lin_trace( matmul(dVd(j,:,:), PFD_matrix) ) / save_pfd
        enddo
    enddo
    savelevel = 2
end subroutine esti_dV


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- calc prefactor factor
!-- * you should have called esti_pfd before
!-- * PFR = tr(PFA)/tr(PFD)
!--   PFA = e^{-b*V(all)}, PFD=e^{-b*V(diagonal)}
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine esti_pfr(esval, Rs)
implicit none
    real(8), dimension(ndof, nbead), intent(in) :: Rs
    real(8), dimension(vdim,vdim) :: PFD_matrix, V, Vd, Vnd, EVd, EVnd, O, OO, Y
    real(8), intent(out) :: esval
    integer :: i
    
    if(savelevel < 2) stop "savelevel error"
    Y = lin_eye(vdim)
    do i=1, size(Rs,dim=2)
        Vd = save_Vds(i,:,:)
        Vnd = save_Vnds(i,:,:)
        EVd = lin_expdiag(-0.5*betap*Vd)
        EVnd = (lin_eye(vdim) - 0.5*betap*Vnd)
        !-- save XXX
        save_EVds(i,:,:) = EVd
        save_EVnds(i,:,:) = EVnd
        O = matmul( EVnd, EVd )

        OO = matmul(transpose(O), O)
        
        !-- save XXX
        save_Os(i,:,:) = O
        save_OOs(i,:,:) = OO
        Y = matmul(Y, OO)
    enddo
    
    esval = lin_trace(Y) / save_pfd
    save_pfr = esval
    savelevel = 3
    
end subroutine esti_pfr


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- estimate the potential
!-- * you should have called esti_pfd before
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine esti_V(esval, Rs)
implicit none
    real(8), dimension(ndof, nbead), intent(in) :: Rs
    real(8), dimension(vdim, vdim) :: Y, O, OO, V
    real(8), intent(out) :: esval
    integer :: i,j

    if(savelevel < 3) stop "savelevel error"

    esval = 0.0
    do i=1, size(Rs,dim=2)
        Y = lin_eye(vdim)
        do j=1, size(Rs,dim=2)
            if(j/= i) then
                OO = save_OOs(j,:,:)
                Y = matmul(Y,  OO)
            else
                O = save_Os(j,:,:)
                V = save_Vds(j,:,:) + save_Vnds(j,:,:)
                Y = matmul(Y, matmul( transpose(O), matmul(V, O) ) )
            endif
        enddo
        esval = esval  + lin_trace(Y)
    enddo
    esval = esval /(real(nbead) * save_pfd)
    save_v = esval
end subroutine esti_V


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Return \partial O / \partial R
!-- * you should have called esti_pfd before
!-- * you should have called esti_dV before
!-- * flag=0 (1st order), 1(hyperbolic), 2(diagonalize)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getOpR(OpR, R, frame, iflag)             
    real(8), dimension(ndof), intent(in) :: R  !-- only one of the frames XXX
    real(8), dimension(ndof,vdim,vdim), intent(out) :: OpR
    real(8), dimension(vdim,vdim) :: L,T,dV
    integer, intent(in) :: frame
    integer, intent(in), optional :: iflag
    integer :: i,j,flag
    if(present(iflag)) then
        flag = iflag
    else
        flag = 2 !-- default value
    endif
    
    if(savelevel < 2) stop "savelevel error"
    
    if(flag .eq. 0) then
        do j=1,size(R) !#XXX bugs --> dot confused with multiply
            OpR(j,:,:) = -0.5*betap*( matmul( save_dVnds(j,frame,:,:), save_EVds(frame,:,:) ) + &
            matmul( save_EVnds(frame,:,:), matmul( save_EVds(frame,:,:), save_dVds(j,frame,:,:) ) ) )
        enddo
        return
    else
        stop "flag failed"
    endif
end subroutine getOpR



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- estimate kinetic energy
!-- * you should have called esti_pfd before
!-- * you should have called esti_dV before
!-- * flag=0, only return primitive value
!-- * flag=1, both primitive and virial value
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine esti_K(kprim, kvir, Rs, iflag)
implicit none
    real(8), dimension(ndof, nbead), intent(in) :: Rs
    real(8), dimension(vdim,vdim) :: PFD_matrix, O, OO, tmp
    real(8), dimension(ndof,vdim,vdim) :: dOdR, Y
    real(8), intent(out) :: kprim, kvir
    real(8), dimension(ndof) :: Rc, pR
    real(8) :: S, term
    integer, intent(in), optional :: iflag
    integer :: i,j,k,flag
    if(present(iflag)) then
        flag = iflag
    else
        flag = 1    !-- both primitive and virial
    endif

    if(savelevel < 3) stop "savelevel error"

    !-- primitive estimator
    S = 0.5 * size(Rs,dim=1) / betap
    do i=1, size(Rs,dim=2)-1
        term = 0.0
        do j=1, size(Rs, dim=1)
            term = term + M(j) * ( Rs(j,i) - Rs(j,i+1) )**2
        enddo
        S = S - 0.5*wp2*term
    enddo
    term = 0.0
    do j=1, size(Rs, dim=1)
        term = term + M(j) * ( Rs(j,1) - Rs(j,size(Rs,dim=2)) )**2
    enddo
    S = S - 0.5*wp2*term
    kprim = S*save_pfr
    save_kprim = kprim
        
    if(flag < 1) return 
    
    !-- additional virial estimator
    Rc = sum(Rs,dim=2) / real(nbead)
    term = 0.0
    do i=1, size(Rs,dim=2)
        call getOpR(dOdR, Rs(:,i), frame=i, iflag=0)  
        do j=1, size(Rs, dim=1)
            Y(j,:,:) = lin_eye(vdim)
        enddo
        do k=1, size(Rs,dim=2)
            if(i.eq. k) then
                O = save_Os(k,:,:)
                do j=1,size(Rs,dim=1)
                    tmp = matmul( transpose(dOdR(j,:,:)), O) + matmul(transpose(O), dOdR(j,:,:))
                    Y(j,:,:) = matmul(Y(j,:,:), tmp )
                enddo            
            else
                OO = save_OOs(k,:,:)
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
    kvir = 0.5_8/(beta) * ( size(Rs,dim=1) * save_pfr - term/save_pfd )
    save_kvir = kvir
    
end subroutine esti_K


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- estimate heat capacity's theta
!-- * you should have called esti_pfd before
!-- * you should have called esti_dV before
!-- actually, the T1 T2 is calulate for
!-- $$ \frac{\partial^2}{\partial \beta^2} F_A $$
!-- and T3 T4 is calculate
!-- $$ (R-R_c) \cdot \frac{\partial^2}{\partial \beta \parial R} F_A  $$
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine esti_theta(Rs)
implicit none
    real(8), dimension(ndof, nbead), intent(in) :: Rs !#XXX fixed bug --> rank definition error
    real(8), dimension(vdim,vdim) :: O, OO, Y1, Y2, Y3, Y4, RpR, V_all
    real(8), dimension(ndof,vdim,vdim) :: dV_all, dOdR, dTdR
    real(8) :: T1, T2, T3, T4
    real(8), dimension(ndof) :: Rc, DeltaR
    integer :: i,j,k,m
    
    if(savelevel < 3) stop "savelevel error"
    
    T1 = 0.0
    T2 = 0.0
    T3 = 0.0
    T4 = 0.0
    Rc = sum(Rs,dim=2) / nbead
    
    do j=1,nbead
        do k=1,nbead               !-- from k=j, if k==j, calc Y2, else calc Y1, then note Y1 we count double times !!
            Y1 = lin_eye(vdim)
            Y2 = lin_eye(vdim)
            Y3 = lin_eye(vdim)     !-- note we use Y3,T3 evaluate all theta3, theta4, theta5 term ??
            Y4 = lin_eye(vdim)
            do i=1,nbead
                O = save_Os(i,:,:)
                OO = save_OOs(i,:,:)
                V_all = save_Vds(i,:,:) + save_Vnds(i,:,:)
                
                if(i/=j .and. i/=k .and. j/=k) then
                    Y1 = matmul(Y1, OO)
                    Y3 = matmul(Y3, OO)
                elseif(i.eq.j .and. j/=k) then
                    Y1 = matmul(Y1, matmul(transpose(O), matmul(V_all, O) ) )

                    call getOpR(dOdR, Rs(:,i), frame=i, iflag=0)
                    do m=1,size(Rs,dim=1)
                        dTdR(m,:,:) = matmul(transpose(dOdR(m,:,:)), O) + matmul(transpose(O), dOdR(m,:,:))
                    enddo
                    DeltaR = Rs(:,i)-Rc
                    RpR = 0.0
                    do m =1, size(Rs,dim=1)
                        RpR = RpR + DeltaR(m)*dTdR(m,:,:)       !-- size(ndof) (dot) size(ndof,vdim,vdim)
                    enddo
                    Y3 = matmul(Y3, RpR)
                elseif(i.eq.k .and. j/=k) then
                    !-- for Y1, here we note we count it doubles times, so it will without times 2 in the final result
                    Y1 = matmul(Y1, matmul(transpose(O), matmul(V_all, O) ) )
                    Y3 = matmul(Y3, matmul(transpose(O), matmul(V_all, O) ) )
                elseif(j.eq.k .and. j/=i) then
                    Y2 = matmul(Y2, OO )
                    Y4 = matmul(Y4, OO )                 
                elseif(j.eq.k .and. j.eq.i) then!--(i.eq.j .and. i.eq.k)
                    Y2 = matmul(Y2, matmul( transpose(O), matmul(matmul(V_all,V_all),O) ) )
                    
                    dV_all = save_dVds(:,i,:,:) + save_dVnds(:,i,:,:)
                    call getOpR(dOdR, Rs(:,i), frame=i, iflag=0)
                    do m=1,size(Rs,dim=1)
                        dTdR(m,:,:) = matmul(transpose(dOdR(m,:,:)), matmul(V_all,O))  &
                                    + matmul(transpose(O), matmul(V_all, dOdR(m,:,:)) ) &
                                    + matmul(transpose(O), matmul(dV_all(m,:,:), O ) )
                    enddo
                    DeltaR = Rs(:,i)-Rc
                    RpR = 0.0
                    do m =1, size(Rs,dim=1)
                        RpR = RpR + DeltaR(m)*dTdR(m,:,:)       !-- size(ndof) (dot) size(ndof,vdim,vdim)
                    enddo
                    Y4 = matmul(Y4, RpR)
                else
                    stop 'error'
                endif
            enddo
            !--
            if(j /= k) then
                T1 = T1 + lin_trace(Y1)
                T3 = T3 + lin_trace(Y3)
            else          
                T2 = T2 + lin_trace(Y2)
                T4 = T4 + lin_trace(Y4)
            endif
        enddo
    enddo
    save_theta(1) = T1 / real(nbead**2) !-------------------- NO FRACTOR 2
    save_theta(2) = T2 / real(nbead**2)
    save_theta(3) = - 0.5 * T3 / real(nbead * beta)
    save_theta(4) = - 0.5 * T4 / real(nbead * beta)
end subroutine esti_theta


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- estimate coherence length
!-- * you should have called esti_pfd before
!-- * you should have called esti_dV before
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine esti_coh(esval)
implicit none
    real(8), dimension(vdim, vdim) :: Rho, Rij, Y, O, OO
    real(8), intent(out) :: esval
    real(8) :: s1, s2
    integer :: i,j,k,l
    
    if(savelevel < 3) stop "savelevel error"
    
    esval = 0.0
    Rho = 0.0
    do i=1, vdim
        do j=i, vdim
            Rij = 0.0
            Rij(i,j) = Rij(i,j) + 0.5
            Rij(j,i) = Rij(j,i) + 0.5
            do k=1, nbead
                Y = lin_eye(vdim)
                do l=1, nbead                   
                    if(l/= k) then
                        OO = save_OOs(l,:,:)
                        Y = matmul(Y,  OO )
                    else
                        O = save_Os(l,:,:)
                        Y = matmul(Y, matmul( transpose(O), matmul(Rij, O) ) )
                    endif
                enddo
                Rho(i,j) = Rho(i,j) + lin_trace(Y)
            enddo
            if(i/=j) Rho(j,i) = Rho(i,j)
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
    save_coh = esval
    return
end subroutine esti_coh

end module MES_Models















