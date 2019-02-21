subroutine MQCprefactor(Mfwd,Mbck,coeff)
use parameters, only : Width, WidthT, Tuningq, Tuningp, InverseWidth, InverseWidthT, &
Iu, Ndof
implicit none

real(8), intent(in) :: Mfwd(Ndof,Ndof,4), Mbck(Ndof,Ndof,4)
real(8) ::  G(Ndof,Ndof), Ginv(Ndof,Ndof), Un(Ndof,Ndof)
complex(16), intent(out) :: coeff

complex(16) :: Mcoeff(Ndof,Ndof), TMP(Ndof, Ndof)
complex(16) :: F1(Ndof,Ndof), F2(Ndof,Ndof), B1(Ndof,Ndof), B2(Ndof,Ndof)
complex(16) :: O1(Ndof,Ndof), O2(Ndof,Ndof), O3(Ndof,Ndof), O4(Ndof,Ndof)

integer :: i, info, sgn
integer :: ipiv(Ndof)

!-- G matrix
G=matmul(Tuningq,InverseWidth)+matmul(Tuningp,Width)+2.d0*matmul(Tuningq,Tuningp)

Ginv=0.d0
Un=0.d0
do i=1, Ndof
        Ginv(i,i)=1.d0/G(i,i) !-- Ginv, what G is diagonal ?
        print *, G
        stop 'debug check G is whether diagnal! if not, error'
        Un(i,i)=1.d0 !-- unit matrix
end do

!-- for example
!-- Mfwd(:,:,1)  Mf_qq
!-- Mfwd(:,:,2)  Mf_qp
!-- Mfwd(:,:,3)  Mf_pq
!-- Mfwd(:,:,4)  Mf_pp

F1 = Mfwd(:,:,4)-Iu*matmul(WidthT,Mfwd(:,:,2))
B1 = matmul(Mbck(:,:,4),WidthT)+Iu*Mbck(:,:,3)
B2 = Mbck(:,:,1)-Iu*matmul(Mbck(:,:,2),WidthT)
F2 = matmul(WidthT,Mfwd(:,:,1))+Iu*Mfwd(:,:,3)

TMP = matmul((Ginv+Un),B1)
O1 = 0.5d0*matmul(F1,TMP)

TMP = matmul(Ginv,B1)
TMP = matmul((0.5d0*InverseWidth+Tuningp),TMP)
O2 = matmul(F2,TMP)

TMP = matmul((Ginv+Un),B2)
O3 = 0.5d0*matmul(F2,TMP)

TMP = matmul(Ginv,B2)
TMP = matmul((0.5d0*Width+Tuningq),TMP)
O4 = matmul(F1,TMP)

!-- the order of multiplication, differ from original
print *, 'warning, change the multiplication order of Mcoeff'

TMP = matmul(G,(O1+O2+O3+O4))
TMP = matmul(InverseWidthT/2.d0,TMP)
Mcoeff = TMP !- save a copy, not nesesarry


!-- LU decomposition, from lapack
!-- IPIV    (output) INTEGER array, the pivot indices; for 1 <= i <= min(M,N), row i of
!-- the matrix was interchanged with row IPIV(i).
call zgetrf(Ndof,Ndof,Mcoeff,Ndof,ipiv,info)

sgn=1
coeff=1.d0

do i=1, Ndof
   coeff=coeff*coeff2(i,i)
   if (ipiv(i).ne.i) sgn=-sgn
enddo

coeff=coeff*sgn
return

end subroutine MQCprefactor



