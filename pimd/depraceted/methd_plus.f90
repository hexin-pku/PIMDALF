! --- Copyright by Shin He <hx0824916@pku.edu.cn> ---

module trans_plus
use const
use pisimul
implicit none
    real(dp), dimension(:,:), allocatable, private :: OM

contains
subroutine m2mb(mb, m, P, iflag)
    integer, intent(in) :: P
    real(8), dimension(P), intent(in) :: m
    real(8), dimension(P), intent(out) :: mb
    integer :: i,j
    
    if (md_methd.eq.1) then
        masscoeff(1) = 1.0_dp
        do i=2,bead
            masscoeff(i) = real(i)/(real(i) - 1.0_dp)
        end do
    else
        if( mod(bead,2) .ne. 0) then
            call send_err('norm mode pimd, beads should be even!')
        end if
        !!
        allocate( OM(bead, bead) )
        masscoeff(1) = real(bead,kind=dp)
        OM(1,:) = 1.0_dp/ DSQRT( real(bead,kind=dp) )
        do i=2,bead/2
            masscoeff(2*i-2) = 2.0_dp*( 1 - DCOS(twopi * real(i-1)/real(bead,kind=dp) ) )*real(bead,kind=dp)
            masscoeff(2*i-1) = 2.0_dp*( 1 - DCOS(twopi * real(i-1)/real(bead,kind=dp) ) )*real(bead,kind=dp)
            do j=1,bead
                OM(2*i-2,j) = DCOS( twopi* real(i-1)*real(j-1)/real(bead,kind=dp) ) * DSQRT(2.0_dp/real(bead,kind=dp) )
                OM(2*i-1,j) = -DSIN( twopi* real(i-1)*real(j-1)/real(bead,kind=dp) ) * DSQRT(2.0_dp/ real(bead,kind=dp) )
            end do
        end do
        masscoeff(bead) = 4.0_dp*real(bead,kind=dp)
        OM(bead,1::2) = 1.0_dp/ DSQRT( real(bead,kind=dp) )
        OM(bead,2::2) = -1.0_dp/ DSQRT( real(bead,kind=dp) )
    end if
    !!
end subroutine set_trsfrm

	!-- normal mode transformation
subroutine calc_ks_norm(mobj)
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i,j,k
    do i=1,mobj%nb
        do j=1,md_bead
            mobj%a(i)%ks(j) = 0.0_dp
            do k=1,md_bead
                mobj%a(i)%ks(j) = mobj%a(i)%ks(j) + mobj%a(i)%x(k)*OM(j,k)
            end do
            mobj%a(i)%ks(j) = mobj%a(i)%ks(j) / DSQRT(real(md_bead,kind=dp))
        end do
    end do
end subroutine calc_ks_norm

subroutine calc_x_norm(mobj)
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i, j, k
    !!
    do i=1,mobj%nb
        do j=1,md_bead
            mobj%a(i)%x(j) = 0.0_dp
            do k=1, md_bead
                mobj%a(i)%x(j) = mobj%a(i)%x(j) + mobj%a(i)%ks(k) * OM(k,j)
            end do
            mobj%a(i)%x(j) = mobj%a(i)%x(j) * DSQRT(real(md_bead,kind=dp))
        end do
    end do
end subroutine calc_x_norm

subroutine calc_fks_norm(mobj)
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i, j, k
    !!
    do i=1,mobj%nb
        do j=1,md_bead,1
            mobj%a(i)%fks(j) = 0.0_dp
            do k=1,md_bead
                mobj%a(i)%fks(j) = mobj%a(i)%fks(j) + mobj%a(i)%fx(k) * OM(j,k)
            end do
            mobj%a(i)%fks(j) = mobj%a(i)%fks(j)
        end do
        mobj%a(i)%fks(1) = mobj%a(i)%fks(1) / DSQRT( real(md_bead,kind=dp) )
        do j=2,md_bead
            mobj%a(i)%fks(j) = mobj%a(i)%fks(j)/ DSQRT(real(md_bead,kind=dp)) + masscoeff(j) &
                * mobj%a(i)%m * md_bfreq2 * mobj%a(i)%ks(j)
        end do
    end do
end subroutine calc_fks_norm
	
	!-- staging transformation
subroutine calc_ks_stag(mobj)
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i,j
    do i=1,mobj%nb
        mobj%a(i)%ks(1) = mobj%a(i)%x(1)
        if(md_bead.eq.1) cycle
        mobj%a(i)%ks(md_bead) = mobj%a(i)%x(md_bead) - mobj%a(i)%x(1)
        do j=2,md_bead-1,1
            mobj%a(i)%ks(j) = mobj%a(i)%x(j) - ( mobj%a(i)%x(j+1) * real(j - 1) + mobj%a(i)%x(1) ) / real(j)
        end do
    end do
end subroutine calc_ks_stag

subroutine calc_x_stag(mobj)
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i, j
    !!
    do i=1,mobj%nb
        mobj%a(i)%x(1) = mobj%a(i)%ks(1)
        if(md_bead.eq.1) cycle
        mobj%a(i)%x( md_bead ) = mobj%a(i)%ks( md_bead ) + mobj%a(i)%x(1)
        do j=md_bead-1,2,-1
            mobj%a(i)%x(j) = mobj%a(i)%ks(j) + ( real( j - 1 ) * mobj%a(i)%x(j+1) + mobj%a(i)%ks(1) ) / real(j)
        end do
    end do
end subroutine calc_x_stag

subroutine calc_fks_stag(mobj)
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i, j
    !!
    do i=1,mobj%nb
        mobj%a(i)%fks(1) = 0.0_dp
        do j=1,md_bead
            mobj%a(i)%fks(1) = mobj%a(i)%fks(1) + mobj%a(i)%fx(j)
        end do
        do j=2,md_bead,1
            mobj%a(i)%fks(j) = mobj%a(i)%fx(j) + mobj%a(i)%fks(j-1) * real(j-2)/real(j-1) 
        end do
        mobj%a(i)%fks(1) = mobj%a(i)%fks(1)/real(md_bead,kind=dp)
        do j=2,md_bead
            mobj%a(i)%fks(j) = mobj%a(i)%fks(j)/real(md_bead,kind=dp) + masscoeff(j) * mobj%a(i)%m * md_bfreq2 * mobj%a(i)%ks(j)
        end do
    end do
end subroutine calc_fks_stag

	!-- public function
subroutine calc_fx_ofmethd(mobj)
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i, j
    do i=1,mobj%nb
        do j=1, md_bead
            mobj%a(i)%fx(j) = Fpn2( 0.5_dp*mobj%a(i)%m,  mobj%a(i)%x(j) )
        end do
    end do
end subroutine calc_fx_ofmethd


!--------------------------------------------------------------------------------------
!! ------------------------ unit normal mode transformation ---------------------------
subroutine x2ks_norm(ks, x)
implicit none
    real(dp), dimension(:), intent(inout) :: ks,x
    integer :: j,k
    do j=1,md_bead
        ks(j) = 0.0_dp
        do k=1,md_bead
            ks(j) = ks(j) + x(k) * OM(j,k)
        end do
        ks(j) = ks(j) / DSQRT(real(md_bead,kind=dp))
    end do
end subroutine x2ks_norm

subroutine ks2x_norm(x,ks)
implicit none
    real(dp), dimension(:), intent(inout) :: ks,x
    integer :: j, k
    !!
    do j=1,md_bead
        x(j) = 0.0_dp
        do k=1, md_bead
            x(j) = x(j) + ks(k) * OM(k,j)
        end do
        x(j) = x(j) * DSQRT(real(md_bead,kind=dp))
    end do
end subroutine ks2x_norm

subroutine fx2fks_norm(fks, fx, ks, mi)
implicit none
    real(dp), dimension(:), intent(inout) :: fks,fx,ks
    real(dp) :: mi
    integer :: j, k
    !!
    do j=1,md_bead,1
        fks(j) = 0.0_dp
        do k=1,md_bead
            fks(j) = fks(j) + fx(k) * OM(j,k)
        end do
        fks(j) = fks(j) * DSQRT( real(md_bead,kind=dp) )
    end do
    fks(1) = fks(1)/DSQRT(real(md_bead,kind=dp))
    do j=2,md_bead
        fks(j) = fks(j)/DSQRT(real(md_bead,kind=dp)) + masscoeff(j) &
            * mi * md_bfreq2 * ks(j)
    end do
end subroutine fx2fks_norm
	!-- staging transformation 
subroutine x2ks_stag(ks,x)
implicit none
    real(dp), dimension(:), intent(inout) :: ks,x
    integer :: j
    ks(1) = x(1)
    if(md_bead.eq.1) return
    ks(md_bead) = x(md_bead) - x(1)
    do j=2,md_bead-1,1
        ks(j) = x(j) - ( x(j+1) * real(j - 1) + x(1) ) / real(j)
    end do
end subroutine x2ks_stag

subroutine ks2x_stag(x,ks)
implicit none
    real(dp), dimension(:), intent(inout) :: ks,x
    integer :: j
    !!
    x(1) = ks(1)
    if(md_bead.eq.1) return
    x( md_bead ) = ks( md_bead ) + x(1)
    do j=md_bead-1,2,-1
        x(j) = ks(j) + ( real( j - 1 ) * x(j+1) + ks(1) ) / real(j)
    end do
end subroutine ks2x_stag

subroutine fx2fks_stag(fks, fx, ks, mi)
implicit none
    real(dp), dimension(:), intent(inout) :: fks,fx,ks
    real(dp) :: mi
    integer :: j
    !!
    fks(1) = 0.0_dp
    do j=1,md_bead
        fks(1) = fks(1) + fx(j)
    end do
    do j=2,md_bead,1
        fks(j) = fx(j) + fks(j-1) * real(j-2)/real(j-1) 
    end do
    fks(1) = fks(1)/real(md_bead,kind=dp)
    do j=2,md_bead
        fks(j) = fks(j)/real(md_bead,kind=dp) + masscoeff(j) * mi * md_bfreq2 * ks(j)
    end do
end subroutine fx2fks_stag    

end module trans_plus


