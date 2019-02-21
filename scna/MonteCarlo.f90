! ##################################################
! ##################################################
!       This subroutine samples the initial
!       phase space points for both electronic
!       and nuclear degrees of freedom, for
!       both pairs of trajectories.
! ##################################################
! ##################################################
! Last edited : 11/13/17
subroutine montecarlosamp(PScoord,icount1,icount2q,icount2p,icountNq,icountNp)
use parameters, only : q0, p0, Iu, Width, Ndof, Tuningq, Tuningp, &
jump1, jump2, jumpNq, jumpNp
implicit none

integer                 ::      i, j, k, r2
integer, intent(inout)  ::      icount1(4), icount2q(2), icount2p(2)
integer, intent(inout)  ::      icountNq(2), icountNp(2)
real*8, intent(inout)   ::      PScoord(4,Ndof)
real*8                  ::      z01(4), z02q(2), z02p(2)
real*8                  ::      z0Nq(2), z0Np(2)
real*8                  ::      zt1(4), zt2q(2), zt2p(2)
real*8                  ::      ztNq(2), ztNp(2)
real*8                  ::      r1, W01, W02q, W02p
real*8                  ::      Wt1, Wt2q, Wt2p, gauss
real*8                  ::      W0Nq, W0Np, WtNq, WtNp
real*8                  ::      accept1(4), accept2q(2), accept2p(2)
real*8                  ::      acceptNq(2), acceptNp(2)

! Read input configurations
!-- first electronic freedom, z01 --- q0,q0',p0,p0'
do i = 1, 4
  z01(i) = PScoord(i,1)
enddo

!-- second electronic freedom, z02q --- q0,q0'; z02p --- p0,p0'; 
z02q(1) = PScoord(1,2)
z02q(2) = PScoord(3,2)
z02p(1) = PScoord(2,2)
z02p(2) = PScoord(4,2)

!-- the nuclear freedom, z0Nq --- R0,R0'; z0Np --- P0,P0'; 
z0Nq(1) = PScoord(1,3)
z0Nq(2) = PScoord(3,3)
z0Np(1) = PScoord(2,3)
z0Np(2) = PScoord(4,3)

!-- randomly move one of 1-electonic mapping vals, (x10,p10,x10',p10')
do i = 1, 4 ! Do it four times so that, on average, 
            ! each gets moved once each time.
    call random_number(r1)
    call onefour(r2)

    zt1(r2) = z01(r2) + (r1-0.5d0)*jump1

    call probs1(z01,W01)
    call probs1(zt1,Wt1)

    call random_number(r1)
    if (Wt1/W01.gt.r1) then
        W01         = Wt1
        z01(r2)     = zt1(r2)
        icount1(r2) = icount1(r2) + 1
    endif
enddo

!-- randomly move one of 2-electonic mapping vals
do i = 1, 2
    ! Randomly move one of (x20,x20')
    call random_number(r1)
    call onetwo(r2)
    zt2q(r2) = z02q(r2) + (r1-0.5d0)*jump2

    call probs2q(z02q,W02q)
    call probs2q(zt2q,Wt2q)

    call random_number(r1)
    if (Wt2q/W02q.gt.r1) then
        W02q        = Wt2q
        z02q(r2)     = zt2q(r2)
        icount2q(r2) = icount2q(r2) + 1
    endif
    ! Randomly move one of (p20,p20')
    call random_number(r1)
    call onetwo(r2)
    zt2p(r2) = z02p(r2) + (r1-0.5d0)*jump2

    call probs2p(z02p,W02p)
    call probs2p(zt2p,Wt2p)

    call random_number(r1)
    if (Wt2p/W02p.gt.r1) then
        W02p        = Wt2p
        z02p(r2)     = zt2p(r2)
        icount2p(r2) = icount2p(r2) + 1
    endif
enddo


!-- randomly move one of nuclear vals
do i = 1, 2
    ! Randomly move one of (R0,R0')
    call random_number(r1)
    call onetwo(r2)
    ztNq(r2) = z0Nq(r2) + (r1-0.5d0)*jumpNq

    call probsNq(z0Nq,W0Nq)
    call probsNq(ztNq,WtNq)

    call random_number(r1)
    if (WtNq/W0Nq.gt.r1) then
      W0Nq        = WtNq
      z0Nq(r2)     = ztNq(r2)
      icountNq(r2) = icountNq(r2) + 1
    endif
    ! Randomly move one of (PR0,PR0')
    call random_number(r1)
    call onetwo(r2)
    ztNp(r2) = z0Np(r2) + (r1-0.5d0)*jumpNp

    call probsNp(z0Np,W0Np)
    call probsNp(ztNp,WtNp)

    call random_number(r1)
    if (WtNp/W0Np.gt.r1) then
      W0Np         = WtNp
      z0Np(r2)     = ztNp(r2)
      icountNp(r2) = icountNp(r2) + 1
    endif
enddo

!-- collect the random walk result
do i = 1, 4
    PScoord(i,1) = z01(i)
enddo

PScoord(1,2) = z02q(1)
PScoord(2,2) = z02p(1)
PScoord(3,2) = z02q(2)
PScoord(4,2) = z02p(2)

PScoord(1,3) = z0Nq(1)
PScoord(2,3) = z0Np(1)
PScoord(3,3) = z0Nq(2)
PScoord(4,3) = z0Np(2)

end subroutine

! Sampling distribution for initially excited oscillator,
! i.e. electronic phase space variables for initially
! occupied electronic state
!-- z0 --- q,p, q',p'
subroutine probs1(z0,w)
use parameters, only : Ndof, Width, InverseWidth, Tuningq, Tuningp, pi
implicit none

real*8, intent(in)   :: z0(4)
real*8, intent(out)  :: w

w =   dexp(-0.5d0*Tuningq(1,1)*(z0(3)-z0(1))**2.d0)*&
      dexp(-0.5d0*Tuningp(1,1)*(z0(4)-z0(2))**2.d0)

w = w*dsqrt(z0(1)**2.d0+z0(2)**2.d0)*dexp(-0.25d0*z0(1)**2.d0-0.25d0*z0(2)**2.d0)*&
      dsqrt(z0(3)**2.d0+z0(4)**2.d0)*dexp(-0.25d0*z0(3)**2.d0-0.25d0*z0(4)**2.d0)
w = w**2 !--#XXX
end subroutine

! Sampling distribution for oscillator initially in its
! ground state configuration, i.e. the phase space variables
! of the initially unoccupied electronic state
! 1) Positions
subroutine probs2q(z0,w)
use parameters, only : Ndof, Width, InverseWidth, Tuningq, Tuningp, pi
implicit none

real*8, intent(in)   :: z0(2)
real*8, intent(out)  :: w

w =   dexp(-0.5d0*Tuningq(2,2)*(z0(2)-z0(1))**2.d0)

w = w*dexp(-0.25d0*z0(1)**2.d0-0.25d0*z0(2)**2.d0)
w = w**2 !--#XXX
end subroutine
! 2) Momenta
subroutine probs2p(z0,w)
use parameters, only : Ndof, Width, InverseWidth, Tuningq, Tuningp, pi
implicit none

real*8, intent(in)   :: z0(2)
real*8, intent(out)  :: w

w =   dexp(-0.5d0*Tuningp(2,2)*(z0(2)-z0(1))**2.d0)

w = w*dexp(-0.25d0*z0(1)**2.d0-0.25d0*z0(2)**2.d0)
w = w**2 !--#XXX
end subroutine

! Sampling distribution for initial nuclear phase
! variables.
! 1) Positions
subroutine probsNq(z0,w)
use parameters, only : Ndof, Width, InverseWidth, Tuningq, Tuningp, pi,&
q0, p0
implicit none

real*8, intent(in)   :: z0(2)
real*8, intent(out)  :: w

w =   dexp(-0.5d0*Tuningq(3,3)*(z0(2)-z0(1))**2.d0)

w = w*dexp(-0.25d0*Width(3,3)*(z0(1)-q0(3))**2.d0)*&
      dexp(-0.25d0*Width(3,3)*(z0(2)-q0(3))**2.d0)
w = w**2 !--#XXX
end subroutine
! Momenta
subroutine probsNp(z0,w)
use parameters, only : Ndof, Width, InverseWidth, Tuningq, Tuningp, pi,&
q0, p0
implicit none

real*8, intent(in)   :: z0(2)
real*8, intent(out)  :: w

w =   dexp(-0.5d0*Tuningp(3,3)*(z0(2)-z0(1))**2.d0)

w = w*dexp(-0.25d0*InverseWidth(3,3)*(z0(1)-p0(3))**2.d0)*&
      dexp(-0.25d0*InverseWidth(3,3)*(z0(2)-p0(3))**2.d0)
w = w**2 !--#XXX
end subroutine

! Gaussian random number generator
double precision function gauss(sigma)

implicit none
real*8, intent(in):: sigma
real*8:: w, v1, v2, l
real*8:: s1, s2

w = 2.d0

do

call random_number(s1)
call random_number(s2)

v1 = 2.d0*s1 - 1.d0
v2 = 2.d0*s2 - 1.d0
w = v1*v1 + v2*v2

if (w.lt.1.d0) exit

end do

l = v1*sqrt(-2.d0*log(w)/(w))
l = sigma*l
gauss = l

end function gauss

! Arbitrary starting points for sampling
subroutine InitMC(PScoord)
use parameters, only : Ndof
implicit none

real*8, intent(inout)   ::      PScoord(4,Ndof)

PScoord(1,1) = 2.5d0          
PScoord(2,1) = 2.15d0         
PScoord(3,1) = 2.35d0         
PScoord(4,1) = 2.25d0         

PScoord(1,2) = 1.45d0
PScoord(2,2) = 1.5d0
PScoord(3,2) = 1.5d0
PScoord(4,2) = 1.55d0

PScoord(1,3) = -6.d0
PScoord(2,3) = 18.1d0
PScoord(3,3) = -4.d0
PScoord(4,3) = 20.1d0

end subroutine

! Generate random integer from 1 to 4
subroutine onefour(r)
implicit none

integer, intent(out)    ::      r
real*8                  ::      r1

call random_number(r1)

if (r1.gt.0.75d0) then
        r = 4
elseif (r1.lt.0.25d0.or.r1.eq.0.25d0) then
        r = 3
elseif (r1.lt.0.5d0) then
        r = 2
else
        r = 1
endif

end subroutine

! Generate random integer from 1 to 2
subroutine onetwo(r)
implicit none

integer, intent(out)    ::      r
real*8                  ::      r1

call random_number(r1)

if (r1.lt.0.5d0.or.r1.eq.0.5d0) then
  r = 1
else
  r = 2
endif

end subroutine
