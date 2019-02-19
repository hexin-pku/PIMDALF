module nucl_sample
use random
implicit none
private
public :: wignerHO_0, wignerHOs_0 
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- groud state of wigner distribution (nuclear)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine wignerHO_0(R, P, m, w, b, Re)
    real(8), intent(out) :: R, P
    real(8), intent(in) :: m, w, b      !-- mass. frequency, beta
    real(8), intent(in), optional :: Re ! equilibrium position
    real(8) :: Q
    
    !-- quantum correct factor Q
    if(b*w<1.0e-7) then
        Q = 1.0_dp
    else
        Q = 0.5*b*w/dtanh(0.5*b*w)
    endif

    R = dsqrt(Q/(b*m*w*w)) * gaussian_rn()
    if(present(Re)) R = R + Re
    
    P = dsqrt(Q*m/b) * gaussian_rn()
end subroutine wignerHO_0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- groud states of wigner distribution (nuclear)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine wignerHOs_0(R, P, m, w, b, Re)
    real(8), intent(in) :: m, w(:), b      !-- mass. frequency, beta
    real(8), intent(out) :: R(size(w)), P(size(w))
    real(8), intent(in), optional :: Re(size(w)) ! equilibrium position
    real(8) :: Q(size(w))
    
    !-- quantum correct factor Q
    Q = 0.5*b*w/dtanh(0.5*b*w)
    
    R = dsqrt(Q/(b*m*w*w)) * gaussian_rns(size(w))
    if(present(Re)) R = R + Re
    
    P = dsqrt(Q*m/b) * gaussian_rns(size(w))
end subroutine wignerHOs_0

end module nucl_sample


