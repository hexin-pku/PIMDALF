! ##################################################
! ##################################################
!       The matrix element of operator  
!       B = delta(Pf-P) in the basis of 
!       coherent states
! ##################################################
! ##################################################
subroutine operB(pf,nuc,nucp,wf)
use parameters, only : Width, q0, p0, Iu, Ndof, pi
implicit none

real*8, intent(in)      :: pf, nuc(2), nucp(2)
complex*16, intent(out) :: wf
real*8                  :: g

g = Width(3,3)

wf = dexp(-(nuc(2)-pf)**2.d0/2.d0/g-(nucp(2)-pf)**2/2.d0/g)*&
    cdexp(Iu*pf*(nucp(1)-nuc(1)))/dsqrt(g*pi)

end subroutine


! ##################################################
! ##################################################
!      Overlap of electronic coherent states
! ##################################################
! ##################################################
subroutine elecoverlap(q,p,qprime,pprime,overlap)
use parameters, only : WidthT, InverseWidthT, Iu, Ndof
implicit none

real*8, intent(in) :: q, p, qprime, pprime
complex*16 :: overlap

overlap = dexp(-0.25d0*WidthT(1,1)*(qprime-q)**2-0.25d0*InverseWidthT(1,1)*(pprime-p)**2)*cdexp(0.5d0*Iu*(pprime+p)*(qprime-q))

end subroutine elecoverlap

! ##################################################
! ##################################################
!               MQC-IVR Prefactor
! ##################################################
! ##################################################

subroutine MQCprefactor(Mfwd,Mbck,coeff)
use parameters, only : Width, WidthT, Tuningq, Tuningp, InverseWidth, InverseWidthT, &
Iu, Ndof
implicit none

real*8, intent(in) :: Mfwd(Ndof,Ndof,4), Mbck(Ndof,Ndof,4)
real*8 ::  G(Ndof,Ndof), Ginv(Ndof,Ndof), Un(Ndof,Ndof), tuningmult(Ndof,Ndof)
complex*16, intent(out) :: coeff

complex*16 :: detRR, detSS
complex*16 :: coeff2(Ndof,Ndof), coeff3(Ndof,Ndof), coeff4(Ndof,Ndof)
complex*16 :: AA(Ndof,Ndof), BB(Ndof,Ndof), DD(Ndof,Ndof), EE(Ndof,Ndof), GGp(Ndof,Ndof), HHp(Ndof,Ndof)
complex*16 :: HH(Ndof,Ndof), GG(Ndof,Ndof), FF(Ndof,Ndof), JJ(Ndof,Ndof), KK(Ndof,Ndof), LL(Ndof,Ndof)
complex*16 :: OO(Ndof,Ndof), PP(Ndof,Ndof), QQ(Ndof,Ndof)

integer:: i, info, sgn
integer, allocatable :: ipiv(:)

Un=0.d0
Ginv=0.d0
G=matmul(Tuningq,InverseWidth)+matmul(Tuningp,Width)+2.d0*matmul(Tuningq,Tuningp)

do i=1, Ndof
        Ginv(i,i)=1.d0/G(i,i)
        Un(i,i)=1.d0
end do

AA = Mfwd(:,:,4)-Iu*matmul(WidthT,Mfwd(:,:,2))
BB = matmul(Mbck(:,:,4),WidthT)+Iu*Mbck(:,:,3)
JJ = Mbck(:,:,1)-Iu*matmul(Mbck(:,:,2),WidthT)
FF = matmul(WidthT,Mfwd(:,:,1))+Iu*Mfwd(:,:,3)

DD = matmul((Ginv+Un),BB)
EE = 0.5d0*matmul(AA,DD)

HH = matmul(Ginv,BB)
GGp = matmul((0.5d0*InverseWidth+Tuningp),HH)
HHp = matmul(FF,GGp)

KK = matmul((Ginv+Un),JJ)
LL = 0.5d0*matmul(FF,KK)

OO = matmul(Ginv,JJ)
PP = matmul((0.5d0*Width+Tuningq),OO)
QQ = matmul(AA,PP)

coeff3 = matmul(InverseWidthT/2.d0,(EE+HHp+LL+QQ))
coeff4 = matmul(G,coeff3)

coeff2 = coeff4
allocate(ipiv(Ndof))

call zgetrf(Ndof,Ndof,coeff2,Ndof,ipiv,info)

sgn=1
coeff=1.d0

do i=1, Ndof
   coeff=coeff*coeff2(i,i)
   if (ipiv(i).ne.i) sgn=-sgn
enddo

coeff=coeff*sgn

end subroutine MQCprefactor

! ##################################################
! ##################################################
!       Divide date between processes
! ##################################################
! ##################################################
subroutine para_range(n1,n2,nprocs,irank,ista,iend)
implicit none
integer, intent(in) :: n1,n2,nprocs,irank
integer, intent(out) :: ista,iend
integer :: iwork1,iwork2
iwork1 = (n2 - n1 + 1)/nprocs
iwork2 = mod(n2 - n1 + 1, nprocs)
ista = irank*iwork1 + n1 + min(irank, iwork2)
iend = ista + iwork1 -1
if (iwork2 .gt. irank) then
iend = iend + 1
endif
return
end subroutine para_range

! ##################################################
! ##################################################
!               Initiate random seed
! ##################################################
! ##################################################
subroutine init_random_seed()
integer :: i, n, clock
integer, allocatable :: seed(:)

call random_seed(size = n)
allocate(seed(n))

call system_clock(count=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /) + node

call random_seed(put = seed)

deallocate(seed)
end subroutine init_random_seed


