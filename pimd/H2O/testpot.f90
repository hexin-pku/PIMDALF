program test
use model
implicit none
    real(8) :: v, dv(9), ddv(9,9), r(9), cart(3,3), pole(3)
    real(8) :: geom(3) = (/1.81_8,1.81_8,1.8239_8/)
    integer :: i,j,k,N
    
    N = 10
    
    open(unit=20, file='_testpot.out', status='replace')
    do i=-N,N
    pole(1) = geom(1) + i*0.005_8
    !do j=-N,N
    pole(2) = geom(2) + i*0.005_8
    do k=-N,N
    pole(3) = geom(3) + k*0.01_8
        call convert_geom_to_cart(cart, pole)
        call vec_rearrange(r, cart)
        call model_V2(v, dv, ddv, r, 2)
        write(20,*) pole, v, dv
    enddo
    !enddo
    enddo
    close(unit=20)

end program test
