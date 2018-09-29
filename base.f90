!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Base module for sci-computations (Copyright Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module base
implicit none
    !-- constant numbers
    real(8), parameter :: pi = 3.14159265358979323846
    real(8), parameter :: e = 2.718281828459045
    real(8), parameter :: hb = 1.054571726E-19
    real(8), parameter :: kb = 1.38064852E-23
    !-- atom
    real(8), parameter :: au_m = 9.10938291E-31
    real(8), parameter :: au_e = 1.602176565E-19
    real(8), parameter :: au_hb = 1.054571726E-34
    real(8), parameter :: au_ke = 8.9875517873681E+9
    real(8), parameter :: au_c = 137.0359991
    real(8), parameter :: au_a0 = 5.2917721092E-11
    real(8), parameter :: au_eh = 4.35974417E-18
    real(8), parameter :: au_t = 2.418884326505E-17
    real(8), parameter :: au_temp = 3.1577464E+5
    real(8), parameter :: au_kb = 1.00000000000
    real(8), parameter :: au_beta = 3.166815367E-6
    real(8), parameter :: wn_2_hetree = 219474.6313702
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- for random number's generation
subroutine init_seed()
	integer :: n, ival(8), v(3), i
	integer, allocatable :: seed(:)
	call date_and_time(values=ival)
	v(1) = ival(8) + 2048*ival(7)
	v(2) = ival(6) + 64*ival(5)     !-- value(4) isn't real(8)ly 'random'
	v(3) = ival(3) + 32*ival(2) + 32*8*ival(1)
 	call random_seed(size=n)
	allocate(seed(n))
	call random_seed()   !-- Give the seed an implementation-dependent kick
	call random_seed(get=seed)
	do i=1, n
    	seed(i) = seed(i) + v(mod(i-1, 3) + 1)
  	enddo
  	call random_seed(put=seed)
  	deallocate(seed)
end subroutine init_seed

subroutine random_norm(n_norm, n_u, n_sigma, flag)
implicit none
    real(8), intent(inout) :: n_norm, n_u
    real(8), intent(in) :: n_sigma
    character, optional :: flag
    real(8) :: tmp_c

    call random_number(n_u)
    tmp_c = sqrt(-2*log(n_u))*n_sigma
    call random_number(n_u)
    
    if (present(flag) .and. flag == 'c') then
        n_norm = tmp_c*cos(2*pi*n_u)
    else
        n_norm = tmp_c*sin(2*pi*n_u)
    endif
end subroutine random_norm
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Log: report errors, or warnings
recursive subroutine log_message(tag, message)
implicit none
    character(*), intent(in) :: tag
    character(*), intent(in) :: message
    select case (trim(tag))
        case ("E")
            print *, "Error: ", message
            stop
        case ("A")
            print *, "Alert: ", message
        case ("W")
            print *, "Warning: ", message
        case ("D")
            print *, "Debug: ", message
        case default
            call log_message("E", "it gives an invalid tag in (subroutine::log_message)")
    end select
end subroutine log_message
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module
