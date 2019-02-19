program testla
use linalgebra
implicit none
    real(8), dimension(8,8) :: X,Y
    real(8), dimension(8) :: Z
    real(8) :: eps = 1e-15
    integer :: i,j
    X = reshape((/ ((1/real(i+j),i=1,8), j=1,8) /), (/ 8, 8 /))
    Z = (/(real(i),i=1,8)/)
    print *, X
    print *, ' ### '
    Y = lin_random_norm(8,8)
    
    
    call lin_ev(Z,Y,X)
    print *, X
    print *, ' ### '
    print *, Y
    print *, ' ### '
    print *, Z
    print *, ' ### '
    X = lin_eps( matmul( transpose(Y), matmul(X, Y)), eps)
    print *, X
end program testla
