! ##################################################
! ##################################################
!       This subroutine propagates classical
!       tracjectories under the MMST 
!       Hamiltonian with the MInt algorithm.
! ##################################################
! ##################################################
! Last edited : 12/12/17
subroutine PropagateFwd(Initialq,Initialp,coord,Sf,Monodromy,flag)
use parameters, only: TimeStep, Ntime, Mqq, Mqp, Mpq, Mpp, EnergyTolerance, Ndof, Mass
implicit none

real*8, intent(in)     :: Initialq(Ndof), Initialp(Ndof)
real*8, intent (out)   :: coord(2,Ndof,0:Ntime), Sf(0:Ntime), Monodromy(Ndof,Ndof,4,0:Ntime)
integer, intent(inout) :: flag   
                         
integer             :: i, j, n, nd2, k
real*8              :: InitialEnergy, dtp(3), Energy, s
real*8              :: maxdiff, h, dt, hdtm
real*8              :: dv11, dv12, dv22,d2v11, d2v12, d2v22
real*8              :: ham, dv, d2v, v11, v12, v22, V(2,2)
real*8, allocatable :: Um(:,:), q(:),p(:),mono(:,:), qs(:), ps(:), bigmono(:,:), tmono(:,:)

nd2 = 2*Ndof

allocate(Um(Ndof,Ndof),mono(Ndof,Ndof),q(Ndof),p(Ndof), qs(2), ps(2))
allocate(bigmono(nd2,nd2),tmono(nd2,nd2))

s = 0.d0

do i = 1, Ndof
    q(i) = initialq(1+mod(i+1,Ndof)) ! q(1) = initialq(3), q(2) = initialq(1), q(3) = initialq(2)
    p(i) = initialp(1+mod(i+1,Ndof))
end do

h  = TimeStep
dt = TimeStep

! Set initial monodromy matrix to the identity
Mqq  = 0.d0
Mqp  = 0.d0
Mpq  = 0.d0
Mpp  = 0.d0
Um   = 0.d0
mono = 0.d0
do i = 1, Ndof
    Mqq(i,i) = 1.d0
    Mpp(i,i) = 1.d0
    Um(i,i) = 1.d0
enddo

! Compute initial energy
call potbitsM1(q(1),v11,v22,v12,dv11,dv22,dv12,d2v11,d2v22,d2v12) !--#XXX fixed bug
V(1,1) = v11
V(1,2) = v12
V(2,1) = v12
V(2,2) = v22
call hamiltonian(q,p,V,ham) !--#XXX fixed bugs
InitialEnergy = ham
! test print *, q,p,InitialEnergy

! Initial conditions
Sf(0)        = s
coord(1,:,0) = Initialq
coord(2,:,0) = Initialp

! Initialize monodromy matrix to the identity
bigmono = 0.0
do k = 1, nd2
    bigmono(k,k) = 1.d0
enddo

Monodromy(:,:,1,0) = Mqq
Monodromy(:,:,2,0) = Mqp
Monodromy(:,:,3,0) = Mpq
Monodromy(:,:,4,0) = Mpp

hdtm = 0.5d0*dt/Mass

do i = 1, Ntime         ! Loop over timesteps

    ! Evolve nuclear position under H1 for half timestep
    q(1) = q(1) + p(1)*hdtm

    ! Evolve monodromy matrix elem    dRt/dx0  x={all 6 variables} 
    bigmono(4,:) = bigmono(4,:) + bigmono(1,:)*hdtm
    
    ! Compute potential and its derivatives
    call potbitsM1(q(1),v11,v22,v12,dv11,dv22,dv12,d2v11,d2v22,d2v12)
    
    !stop 'debug'
    ! Evolve electronic variables and nuclear momentum under H2
    call elpmat(p(1),q(2:3),p(2:3),dt,v11,v22,v12,dv11,dv22,dv12,d2v11,d2v22,d2v12,tmono)
    
    ! Update monodromy matrix after evolution under H2
    bigmono = matmul(tmono,bigmono)
    
    ! Evolve nuclear position again for half timestep
    q(1) = q(1) + p(1)*hdtm

    ! Evolve monodromy matrix elements again
    bigmono(4,:) = bigmono(4,:) + bigmono(1,:)*hdtm

    ! Translate bigmono to Mqq etc
    Mpp(:,:) = bigmono(1:Ndof,1:Ndof)
    Mpq(:,:) = bigmono(1:Ndof,Ndof+1:nd2)
    Mqp(:,:) = bigmono(Ndof+1:nd2,1:Ndof)
    Mqq(:,:) = bigmono(Ndof+1:nd2,Ndof+1:nd2)

    ! Define symplecticity criteria
    Mono    = matmul(transpose(Mqq),Mpp)-matmul(transpose(Mpq),Mqp)
    MaxDiff = maxval(dabs(Mono-Um))

    ! Compute  energy
    call potbitsM1(q(1),v11,v22,v12,dv11,dv22,dv12,d2v11,d2v22,d2v12)
    V(1,1) = v11
    V(1,2) = v12       
    V(2,1) = v12
    V(2,2) = v22
    call hamiltonian(q,p,V,ham)
    Energy = ham !--#XXX
    !-- test print *, q,p,Energy

    ! Update output arrays
    do j = 1, Ndof
        ! Trajectory
        coord(1,j,i) = q(1+mod(j,Ndof))
        coord(2,j,i) = p(1+mod(j,Ndof))
        ! Monodromy matrix
        Monodromy(j,:,1,i) = Mqq(1+mod(j,Ndof),:)
        Monodromy(j,:,2,i) = Mqp(1+mod(j,Ndof),:)
        Monodromy(j,:,3,i) = Mpq(1+mod(j,Ndof),:)
        Monodromy(j,:,4,i) = Mpp(1+mod(j,Ndof),:)
    end do

    ! Compute action
    dtp(2) = p(2)*v11 + p(3)*v12
    dtp(3) = p(3)*v22 + p(2)*v12
    dtp(1) = p(1)/mass
    Sf(i) = Sf(i-1) + TimeStep*(dot_product(p,dtp) - Energy)

    !Check for energy conservation and symplecticity
    !print *, dabs(1.d0-Energy/InitialEnergy), MaxDiff
    ! stop 'debug'
    if (dabs(1.d0-Energy/InitialEnergy).ge.EnergyTolerance.or.MaxDiff.ge.EnergyTolerance) then
        flag = 1
        goto 112 
    endif

end do
112 continue
end subroutine PropagateFwd


! ##################################################
! ##################################################
!       Propagation of nuclear momentum and
!       electronic variables under H2.
! ##################################################
! ##################################################
subroutine elpmat(pp,qs,ps,dt,a,b,c,dv11,dv22,dv12,d2v11,d2v22,d2v12,tmono)
implicit none
real*8, dimension (2,2) :: Vmat, Smat, VPmat, Wmat, Gammat, Ximat
real*8, dimension (2,2) :: Emat, Fmat, dSmat, Cmat, dCmat, Dmat
real*8, dimension (2,2) :: dGammat, dXimat, dDmat, tempmat, Vppmat, dWmat
real*8, dimension (2)   :: lam, qs, ps, qc, pc, co, si, dlam
real*8                  :: dt, pp, a, b, c, sinterm
real*8                  :: dv11,dv22,dv12,lam12, delp, ar, br, cr 
real*8                  :: hmr, hpr, d2v11,d2v22,d2v12,dlam12
real*8                  :: dlaminv, dcosterm, costerm, temp
integer                 :: k
real*8, dimension (6,6), intent(out) :: tmono

    ! Compute the eigenvalues (lam) and eigenvectors (Smat), as well
    ! as their derivatives, of the matrix V
    call matderiv(a,b,c,dv11,dv22,dv12,lam,dlam,Smat,dSmat)

    do k = 1, 2
       co(k) = dcos(lam(k)*dt)
       si(k) = dsin(-lam(k)*dt)
    enddo

    ! Copy original coordinates
    qc = qs
    pc = ps
    call Threm(Smat,co,transpose(Smat),Cmat)
    call Threm(Smat,si,transpose(Smat),Dmat)

    ! Calculate the time-evolved electronic positions and momenta
    qs = matmul(Cmat,qc) - matmul(Dmat,pc)
    ps = matmul(Cmat,pc) + matmul(Dmat,qc)

    ! Initialize monodromy matrix to identity
    tmono(:,:) = 0.d0
    do k = 1, 6
       tmono(k,k) = 1.d0 
    enddo

    ! Mono. mat. elements from electronic-only evolution
    tmono(2:3,2:3) =  Cmat(:,:)
    tmono(5:6,5:6) =  Cmat(:,:)
    tmono(2:3,5:6) =  Dmat(:,:)
    tmono(5:6,2:3) = -Dmat(:,:)

    ! Nuclear momentum motion requires V' and V''
    VPmat(1,1)  = dv11
    VPmat(2,2)  = dv22
    VPmat(1,2)  = dv12
    VPmat(2,1)  = dv12
    VPPmat(1,1) = d2v11
    VPPmat(1,2) = d2v12
    VPPmat(2,1) = d2v12
    VPPmat(2,2) = d2v22

    ! Compute potential derivative matrix in adiabatic basis
    Wmat       = matmul(transpose(Smat),matmul(VPmat,Smat))
    lam12      = lam(1) - lam(2)
    dlam12     = dlam(1) - dlam(2)
    Ximat(:,:) = 0.d0
    Gammat     = Wmat*dt

    if (lam12.ne.0.d0) then
        sinterm     =  sin(lam12*dt)
        Gammat(1,2) =  sinterm*Wmat(1,2)/lam12
        Gammat(2,1) =  Gammat(1,2)
        costerm     =  1.d0-cos(lam12*dt)
        Ximat(1,2)  =  Wmat(1,2)*costerm/lam12
        Ximat(2,1)  = -Ximat(1,2) ! skew-symmetric
    endif

    Emat = matmul(matmul(Smat,Gammat),transpose(Smat))
    Fmat = matmul(matmul(Smat,Ximat),transpose(Smat))
    delp = dot_product(qc,matmul(Emat,qc)) + dot_product(pc,matmul(Emat,pc)) & 
         - 2.d0*dot_product(qc,matmul(Fmat,pc))
    ! Evolution of nuclear momentum
    PP = PP - 0.5d0*delp + 0.5d0*(VPmat(1,1) + VPmat(2,2))*dt

    ! Mono. mat. elements dP/dp and dP/dq
    tmono(1,2:3) = -matmul(pc,Emat) + matmul(qc,Fmat)
    tmono(1,5:6) = -matmul(qc,Emat) - matmul(pc,Fmat)

    ! Mono. mat. elements dp/dR and dq/dR
    call threm(dSmat,co,transpose(Smat),dCmat)
    call threm(Smat,si*dlam*dt,transpose(Smat),tempmat)
    dCmat = dCmat + tempmat + transpose(dCmat)
    call threm(dSmat,si,transpose(Smat),dDmat)
    call threm(Smat,-co*dlam*dt,transpose(Smat),tempmat)
    dDmat = dDmat + tempmat + transpose(dDmat)
    tmono(2:3,4) = matmul(dCmat,pc) + matmul(dDmat,qc)
    tmono(5:6,4) = matmul(dCmat,qc) - matmul(dDmat,pc)

    ! Compute mono. mat. elements dP/dR
    dWmat = matmul(transpose(dSmat),matmul(vpmat,smat)) + matmul(transpose(Smat),matmul(vppmat,smat)) &
         + matmul(transpose(Smat),matmul(vpmat,dSmat))
    dGammat = dWmat*dt
    dXimat = 0.d0
    if (lam12 .ne. 0.d0) then
        dGammat(1,2) = sinterm*dWmat(1,2)/lam12 + &
            Wmat(1,2)*(dlam12/lam12)*(-sinterm/lam12 + dcos(lam12*dt)*dt)
        dGammat(2,1) = dGammat(1,2)
        dXimat(1,2)  = -(dlam12/lam12)*ximat(1,2) & 
            + (dlam12/lam12)*sinterm*dt*wmat(1,2) &
            + dWmat(1,2)*costerm/lam12
        dXimat(2,1)  = -dXimat(1,2)
    endif
    tempmat    = matmul(matmul(dSmat,Gammat),transpose(Smat))
    tempmat    = tempmat + matmul(matmul(Smat,dGammat),transpose(Smat)) + transpose(tempmat)
    temp       = dot_product(qc,matmul(tempmat,qc)) + dot_product(pc,matmul(tempmat,pc))
    tempmat    = matmul(matmul(Smat,dXimat),transpose(Smat)) 
    temp       = temp - 2.d0*dot_product(qc,matmul(tempmat,pc))
    tmono(1,4) = -0.5d0*temp + 0.5d0*(vppmat(1,1)+vppmat(2,2))*dt
    
    return
end subroutine elpmat

! ##################################################
! ##################################################
!       Calculate the potential energy matrix
!       as well as its first and second
!       derivative for models 1 and 2.
!       provide:                                            
!                Tully's model I
! ##################################################
! ##################################################
subroutine potbitsM1(R,v11,v22,v12,dv11,dv22,dv12,d2v11,d2v22,d2v12)
use parameters, only : a_tanh, gaus_pref, b_exp, v0
implicit none

real*8   :: R, v11, v22, v12, dv11, dv22, dv12, d2v11, d2v22, d2v12

v11   =  v0*(1.d0+dtanh(a_tanh*R))
v22   =  v0*(1.d0-dtanh(a_tanh*R))
v12   =  gaus_pref*dexp(-b_exp*R**2) 
dv11  =  v0*a_tanh*(1.d0-dtanh(a_tanh*R)**2)
dv22  = -v0*a_tanh*(1.d0-dtanh(a_tanh*R)**2)
dv12  = -2.d0*b_exp*R*v12
d2v11 = -2.d0*dtanh(a_tanh*R)*a_tanh*dv11
d2v22 = -2.d0*dtanh(a_tanh*R)*a_tanh*dv22
d2v12 = (4.d0*(b_exp*R)**2 - 2.d0*b_exp)*v12

return
end subroutine

subroutine matderiv(a,b,c,ar,br,cr,lam,dlam,Smat,dSmat)
! ##################################################
! ##################################################
!       Finds the eigenvalues and eigenvectors 
!       of a 2x2 real symmetric matrix. Also finds 
!       the derivatives of eigenvalues and derivatives 
!       of eigenvectors. Matrices are of the form
!           M = (a c)  ,  dM/dR = (ar cr)
!               (c b)  ,          (cr br)
! ##################################################
! ##################################################
    implicit none
    double precision :: a,b,c,ar,br,cr,pi,hp,hm,hpr,hmr,theta,disc,dtdr
    double precision, dimension (2) :: lam, dlam
    double precision, dimension (2,2) :: Smat, dSmat

    pi   = dacos(-1.d0)
    hp   = 0.5d0*(a+b)
    hm   = 0.5d0*(a-b)
    disc = dsqrt(hm**2 + c**2)

    lam(1) = hp + disc
    lam(2) = hp - disc ! lam(1) > lam(2)

    hpr = 0.5d0*(ar+br)
    hmr = 0.5d0*(ar-br)

    !-- derivatives of eigenvalues
    dlam(1) = hpr + (hmr*hm+c*cr)/disc
    dlam(2) = hpr - (hmr*hm+c*cr)/disc

    !-- derivatives of eigenvectors
    if (c.eq.0.d0) then
       if (a.ge.b) then ! includes trivial case c=0, a=b
          theta = 0.d0
       else 
          theta = pi*0.5d0
       endif
    elseif (a.eq.b) then ! implies c .ne. 0
       if (c.gt.0) then
          theta = pi*0.25d0
       else
          theta = pi*0.75d0
       endif
    else
       theta = datan(2*c/(a-b))
       if (theta.le.0) theta = theta + pi
       theta = theta*0.5d0
       if (c.le.0) theta = theta + pi*0.5d0
    end if

    Smat(1,1)  = dcos(theta)
    Smat(2,2)  = Smat(1,1)
    Smat(1,2)  = -dsin(theta)
    Smat(2,1)  = -Smat(1,2)

    dtdr       = 0.5d0*(cr*hm - c*hmr)/(hm**2+c**2) ! d theta/d R

    dSmat(1,1) = Smat(1,2)*dtdr
    dSmat(1,2) = -Smat(1,1)*dtdr
    dSmat(2,1) = -dSmat(1,2)
    dSmat(2,2) = dSmat(1,1)

    return
end subroutine matderiv


! ##################################################
! ##################################################
!             Compute MMST Hamiltonian
! ##################################################
! ##################################################
subroutine hamiltonian(q,p,V,ham) !-- q = qN, q1, q2
use parameters, only : Ndof, Mass
implicit none

real*8, intent(in) :: q(Ndof), p(Ndof), V(2,2)
real*8, intent(out):: ham

ham = 0.5d0*p(1)**2.d0/Mass &
          + 0.5d0*dot_product(q(2:3),matmul(V,q(2:3))) &
          + 0.5d0*dot_product(p(2:3),matmul(V,p(2:3))) - (V(1,1)+V(2,2))/2.d0 !--#XXX bugs fixed

return
end subroutine


! ##################################################
! ##################################################
!       Calculates Tmat = Amat*lam*Bmat where 
!       lam is diagonal (given as 2x2 vector)
! ##################################################
! ##################################################
subroutine threm(Amat,lam,Bmat,Tmat)
implicit none

double precision, dimension (2,2) :: Amat, Bmat, Tmat
double precision, dimension (2)   :: lam
integer                           :: j, k

do j = 1, 2
    do k = 1, 2
      Tmat(j,k) = Amat(j,1)*lam(1)*Bmat(1,k) + Amat(j,2)*lam(2)*Bmat(2,k)
    enddo
enddo

return
end subroutine
