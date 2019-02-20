! --- Copyright by Shin He <hx0824916@pku.edu.cn> ---

module elim_plus
use MyDef
use AM_script
use md_info
implicit none
    real(dp), dimension(3) :: elim_rc, elim_pc  !! position and momenta of centroid
    real(dp), dimension(3,3) :: elim_inertia    !! the intertia
    real(dp), dimension(3) :: elim_omega        !! the angular vecolity
    real(dp), dimension(3) :: elim_AMomenta     !! the angular momenta
    real(dp) :: ww !! total mass-factor
    !!
    integer :: elim_succ !! note for lapack complication
    logical :: elim_set0 !! note whether centroid is on zero point
    !!
contains
    subroutine set_elim(mobj)
        type(mole), intent(inout) :: mobj
        integer :: i
        elim_succ = 0
        elim_set0 = .false.
        ww=0._dp
        do i=1,mobj%nb
            ww=ww+real(md_bead)*mobj%a(i)%m
        end do
    end subroutine set_elim
    !!    
    subroutine outer(y,x1,x2) !! here return y=x1 outer x2
        real(dp), dimension(3), intent(inout) :: y, x1, x2
        y(1)=x1(2)*x2(3) - x1(3)*x2(2)
        y(2)=x1(3)*x2(1) - x1(1)*x2(3)
        y(3)=x1(1)*x2(2) - x1(2)*x2(1)
        return 
    end subroutine outer
    !!
    subroutine elim_trs(mobj)
    implicit none
        type(mole), intent(inout) :: mobj
        integer :: i,j
        !! real(dp), dimension(3) :: tmpvec
        !! calculate elim_rc and elim_pc
        elim_rc=0._dp
        elim_pc=0._dp
        do i=1,mobj%nb
            do j=1,3
                elim_rc(j)=elim_rc(j)+sum( mobj%a(i)%x(j::3) ) * mobj%a(i)%m
                elim_pc(j)=elim_pc(j)+sum( mobj%a(i)%p(j::3) )
            end do
        end do
        elim_rc=elim_rc/ww
        !! elim_pc is needn't changed
        !! then we should eliminate the transistion motion (rotation elimination need!)
        do i=1,mobj%nb
            mobj%a(i)%x(1::3) = mobj%a(i)%x(1::3) - elim_rc(1)
            mobj%a(i)%x(2::3) = mobj%a(i)%x(1::3) - elim_rc(2)
            mobj%a(i)%x(3::3) = mobj%a(i)%x(1::3) - elim_rc(3)
        end do
        elim_set0 = .true.
    end subroutine elim_trs
    !!
    subroutine elim_rot(mobj)
        type(mole), intent(inout) :: mobj
        integer, dimension(3) :: ipiv=1
        real(dp), dimension(3) :: tmpvec
        real(dp), dimension(3,3) :: tmpvv, workplace
        integer :: nwork = 3
        integer :: i,j,k
        !!
        if(elim_set0.eqv..false.) then
            call elim_trs(mobj)
        else
            elim_set0 = .false.
        end if
        !!
        elim_inertia=0._dp
        elim_AMomenta=0._dp
        do i=1,mobj%nb
            !! Angular momenta
            do j=1,md_bead
                call outer(tmpvec, mobj%a(i)%x(3*j-2:3*j:1), mobj%a(i)%p(3*j-2:3*j:1))
                elim_AMomenta = elim_AMomenta + tmpvec
            end do
            !!
            !! Inertia
            elim_inertia(1,1) = elim_inertia(1,1) + mobj%a(i)%m * dot_product(mobj%a(i)%x(2::3),mobj%a(i)%x(2::3)) + &
                mobj%a(i)%m * dot_product(mobj%a(i)%x(3::3),mobj%a(i)%x(3::3))
            elim_inertia(2,2) = elim_inertia(2,2) + mobj%a(i)%m * dot_product(mobj%a(i)%x(1::3),mobj%a(i)%x(1::3)) + &
                mobj%a(i)%m * dot_product(mobj%a(i)%x(3::3),mobj%a(i)%x(3::3))
            elim_inertia(3,3) = elim_inertia(3,3) + mobj%a(i)%m * dot_product(mobj%a(i)%x(1::3),mobj%a(i)%x(1::3)) + &
                mobj%a(i)%m * dot_product(mobj%a(i)%x(2::3),mobj%a(i)%x(2::3))
            do j=1,3
                do k=j+1,3
                    elim_inertia(j,k) = elim_inertia(j,k) - mobj%a(i)%m* dot_product(mobj%a(i)%x(j::3),mobj%a(i)%x(k::3))
                    elim_inertia(k,j) = elim_inertia(j,k)
                end do
            end do
        end do
        !!
        !! calc elim_omega
        elim_omega = elim_AMomenta
        tmpvv = elim_inertia ! maybe we needn't, and just use elim_inertia directly 
        call DSYSV('U',3,1,tmpvv,3,ipiv,elim_omega,3,workplace,nwork,elim_succ)
        if(elim_succ.ne.0) then
            print *,'Lapack error'
            elim_succ = 1
            stop
        end if
        
        !!
        !! rotation elimnation
        do i=1,mobj%nb
            do j=1,md_bead
                call outer(tmpvec, elim_omega, mobj%a(i)%x(3*j-2:3*j:1) )
                mobj%a(i)%p(3*j-2:3*j:1) = mobj%a(i)%p(3*j-2:3*j:1) - mobj%a(i)%m * tmpvec
            end do
        end do
    end subroutine elim_rot
end module elim_plus
