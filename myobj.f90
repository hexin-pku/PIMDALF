!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- myobj module (provide box type)
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module myobj
implicit none
!-- box1d for md
type box1d
    real(8), dimension(:), pointer :: x
    real(8), dimension(:), pointer :: p
    real(8), dimension(:), pointer :: f
    real(8), dimension(:), pointer :: m
end type box1d
!-- box2d for pimd
type box2d
    real(8), dimension(:,:), pointer :: x
    real(8), dimension(:,:), pointer :: p
    real(8), dimension(:,:), pointer :: f
    real(8), dimension(:,:), pointer :: m
end type box2d

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- initializations
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine init_box1d(ibox, ifdom)
    type(box1d), intent(inout) :: ibox
    integer, intent(in) :: ifdom
    allocate( ibox%x( ifdom ) )
    allocate( ibox%p( ifdom ) )
    allocate( ibox%f( ifdom ) )
    allocate( ibox%m( ifdom ) )
end subroutine init_box1d

subroutine init_box2d(ibox, ifdom, ibead)
    type(box2d), intent(inout) :: ibox
    integer, intent(in) :: ifdom, ibead
    allocate( ibox%x( ifdom, ibead ) )
    allocate( ibox%p( ifdom, ibead ) )
    allocate( ibox%f( ifdom, ibead ) )
    allocate( ibox%m( ifdom, ibead ) )
end subroutine init_box2d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- destructions
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine del_box1d(ibox)
    type(box1d), intent(inout) :: ibox
    deallocate( ibox%x )
    deallocate( ibox%p )
    deallocate( ibox%f )
    deallocate( ibox%m )
end subroutine del_box1d

subroutine del_box2d(ibox)
    type(box2d), intent(inout) :: ibox
    deallocate( ibox%x )
    deallocate( ibox%p )
    deallocate( ibox%f )
    deallocate( ibox%m )
end subroutine del_box2d

end module



