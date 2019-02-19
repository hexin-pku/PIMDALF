!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- random module for sci-computations 
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module random
use const
implicit none
    logical, private :: has_initialed = .false.
    
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- initial seed for random number's generation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine init_seed(my_pe)
	integer :: n, ival(8), v(3), i
	integer, allocatable :: seed(:)
	integer, optional :: my_pe
	
	call date_and_time(values=ival)
	v(1) = 101 * ival(8) + 256*ival(7) + mod(ival(8), 103)
    v(2) = ival(6) + 64*ival(5) + ( mod(1993, 2+ival(7) ) + 3 ) * 97 * ival(8) 
    v(3) = ival(3) + 4*ival(2) + 16*ival(1) + ( mod(997, 1+ival(7)) + 5 ) * 101* ival(8)
    
 	call random_seed(size=n)
    allocate(seed(n))
    
    !-- Give the seed an implementation-dependent kick
	call random_seed()
	call random_seed(get=seed)
	
	!-- first Bias
	do i=1, n
    	seed(i) = seed(i) + v(mod(i-1, 3) + 1) + ival(8)
  	enddo
  	  	
  	!-- second Bias
  	if ( n >= size(ival) ) then
        seed( 1 : size(ival) ) = seed( 1 : size(ival) ) + ival(:)
    else
        seed(:) = seed(:) + ival( 1 : n )
    endif
    
    !-- if parallel
    if ( present(my_pe) ) seed(:) = seed(:) + 5 * my_pe
    
  	call random_seed(put=seed)
  	deallocate(seed)
  	has_initialed = .true.
end subroutine init_seed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- uniform random number function (array)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uniform_rn() result(res)
    real(dp) :: res
    call random_number(harvest=res)
end function uniform_rn

function uniform_rns(arraysize) result(res)
    integer, intent(in) :: arraysize
    real(dp) :: res(max(1,arraysize))
    call random_number(harvest=res)
end function uniform_rns

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- catalog random number function (array)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function catalog_rn(a) result(res)
    integer, intent(in) :: a
    integer :: res
    real(dp) :: tmp
    call random_number(harvest=tmp)
    res = int(tmp*a) + 1
end function catalog_rn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- normal (gaussian) random number subroutine
!-- scalar subroutine: random_normal
!-- array  subroutine: random_normals
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--#XXX, #TODO how to fix with scalar input harvest
subroutine random_normal(harvest)
use const
    real(dp), intent(inout) :: harvest
    real(dp) :: un(2)
    !-- Box-Muller algorithm (cos side)
    call random_number(harvest = un)
    harvest = dsqrt( -2.0_dp * dlog(un(1)) ) * dcos( 2.0_dp * pi * un(2))
    return
end subroutine random_normal

subroutine random_normals(harvest)
use const
    real(dp), intent(inout) :: harvest(:)
    real(dp) :: un(2, size(harvest))
    !-- Box-Muller algorithm (cos side)
    call random_number(harvest = un)
    harvest = dsqrt( -2.0_dp * dlog(un(1,:)) ) * dcos( 2.0_dp * pi * un(2,:))
    return
end subroutine random_normals

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- guassian random number function (array)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Box-Muller algorithm (cos side)
function gaussian_rn(sigma) result(res)
use const
    real(dp), intent(in), optional :: sigma
    real(dp) :: res
    real(dp) :: un(2)
    call random_number(harvest=un)
    res = dsqrt( -2.0_dp * dlog(un(1)) ) * dcos( 2.0_dp * pi * un(2))
    if(present(sigma)) res = res* sigma
end function gaussian_rn

!-- Metropolis algorithm
function gaussian_rn1(sigma) result(res)
implicit none
real(dp), intent(in), optional :: sigma
real(dp) :: res
real(dp):: w, v1, v2, l
real(dp):: s1, s2
    w = 2.d0
    do
    call random_number(s1)
    call random_number(s2)
    v1 = 2.d0*s1 - 1.d0
    v2 = 2.d0*s2 - 1.d0
    w = v1*v1 + v2*v2
    if (w.lt. 1.d0) exit
    end do
    
    l = v1*sqrt(-2.d0*log(w)/(w))
    if(present(sigma)) l = sigma*l
    res = l
end function gaussian_rn1

!-- Box-Muller algorithm (cos side)
function gaussian_rns(arraysize, sigma) result(res)
use const
    integer, intent(in) :: arraysize
    real(dp), intent(in), optional :: sigma
    real(dp) :: res(arraysize)
    real(dp) :: un(2, arraysize)
    call random_number(harvest=un)
    res(:) = dsqrt( -2.0_dp * dlog(un(1,:)) ) * dcos( 2.0_dp * pi * un(2,:))
    if(present(sigma)) res = res * sigma
end function gaussian_rns

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- exponential distribution
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exponent_rn(lambda) result(res)
implicit none
    real(dp), intent(in), optional :: lambda
    real(dp) :: res
    call random_number(res)
    res = - dlog(1.0_dp - res)
    if(present(lambda)) res = res / lambda
end function exponent_rn

function exponent_rns(arraysize, lambda) result(res)
implicit none
    integer, intent(in) :: arraysize
    real(dp), intent(in), optional :: lambda
    real(dp) :: res(arraysize)
    call random_number(res)
    res = - dlog(1.0_dp - res)
    if(present(lambda)) res = res / lambda
end function exponent_rns

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Possion distribution
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Knuth algorithm
function possion_rn(lambda) result(res)
implicit none
    real(dp), intent(in), optional :: lambda
    real(dp) :: L, p, randu
    integer :: res 
    res = 0
    p = 1.0_dp
    
    if(lambda <=0 ) stop 'error'
    if(present(lambda)) then
        L = dexp(-lambda)
    else
        L = dexp(-1.0_dp)
    endif
    do while( p > L)
        call random_number(randu)
        p = p* randu
        res = res + 1
    enddo
    res = res - 1
end function possion_rn

end module random

!program test
!use random
!implicit none
!    real(dp) :: rn(3)
!    call init_seed()
!    call random_normals(rn)
!    print *, rn
!    rn = gaussian_rns(size(rn))
!    print *, rn
!    print *, catalog_rn(10)
!end program test

