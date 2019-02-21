! ##################################################
! ##################################################
!       This program computes a 1D particle's
!       distribution of nuclear momentum after
!       transiting a 2-state curve crossing within
!       a 2-level system.
! ##################################################
! ##################################################
! Last edited : 12/12/17
program CF
use parameters
implicit none

include 'mpif.h'

integer              :: i, j, k, l
integer              :: ierr, myrank, nprocs, istart, iend, length, imc
integer              :: bad, badtot, flag, Maslov
integer, allocatable :: icount1(:), icount2q(:), icount2p(:), icountNq(:), icountNp(:)
integer, allocatable :: icount12(:), icount2q2(:), icount2p2(:), icountNq2(:), icountNp2(:) 

real*8, allocatable  :: PScoord(:,:), qpfwd(:,:,:), qpbck(:,:,:)
real*8, allocatable  :: Mfwd(:,:,:,:), Mfprime(:,:,:,:), Sfwd(:), Sbck(:), grid(:)
real*8, allocatable  :: realtcf(:), imagtcf(:)
real*8, allocatable  :: Initialq(:), Initialp(:), q0p(:), p0p(:), coord(:,:)
real*8, allocatable  :: MonodromyFwd(:,:,:), MonodromyBck(:,:,:)
real*8               :: dpf, summ

complex*16              :: coeff, prev, TCF, overlap1, overlap2, overlap12
complex*16              :: elemB, oR, o1, o2
complex*16, allocatable :: probs(:)

! Read input file
call input

allocate(PScoord(4,Ndof),qpfwd(2,Ndof,0:Ntime),qpbck(2,Ndof,0:Ntime))
allocate(Mfwd(Ndof,Ndof,4,0:Ntime),Mfprime(Ndof,Ndof,4,0:Ntime),Sfwd(0:Ntime),Sbck(0:Ntime))
allocate(coord(2,Ndof),MonodromyFwd(Ndof,Ndof,4),MonodromyBck(Ndof,Ndof,4))
allocate(Initialq(Ndof),Initialp(Ndof),q0p(Ndof),p0p(Ndof))
allocate(realtcf(0:Npf),imagtcf(0:Npf),grid(0:Npf))
allocate(probs(0:Npf))
allocate(icount1(4),icount2q(2),icount2p(2),icountNq(2),icountNp(2))
allocate(icount12(4),icount2q2(2),icount2p2(2),icountNq2(2),icountNp2(2))

bad     = 0
badtot  = 0
TCF     = 0.d0
realtcf = 0.d0
imagtcf = 0.d0
probs   = 0.d0

! Grid of Pf values
dpf = (pfmax-pfmin)/Npf
do i = 0, Npf
   grid(i) = pfmin + i*dpf
enddo

icount1  = 0
icount2q = 0
icount2p = 0
icountNq = 0
icountNp = 0

! Initialize phase space coordinates before random walk
!! Arbitrary starting points for sampling
PScoord = 0.d0
call InitMC(PScoord)

! Setting up the parallelization

call mpi_init(ierr) !<--- Create child processes
call mpi_comm_size(mpi_comm_world,nprocs,ierr) !<--- Find the total number of processes
call mpi_comm_rank(mpi_comm_world,myrank,ierr) !<--- Find the rank of each process
call para_range(1,NumberMC,nprocs,myrank,istart,iend) !<--- Divide configurations between processes
call mpi_barrier(mpi_comm_world,ierr)

call init_random_seed()

do imc = istart, iend  ! Loop over initial conditions
  ! Random walk obtians new initial phase space point.
  ! The counters track MC acceptance ratios.
  do i = 1, 4  ! Keep every four random moves. Four is arbitrary
   call montecarlosamp(PScoord,icount1,icount2q,icount2p,icountNq,icountNp)
  enddo

  flag   = 0           ! Reset bad trajectory tracker
  Maslov = 0           ! Reset Maslov index
  prev   = 1.d0        ! Reset prefactor tracker #XXX

  Initialq = PScoord(1,:) ! q0
  Initialp = PScoord(2,:) ! p0
  q0p      = PScoord(3,:) ! q0'
  p0p      = PScoord(4,:) ! p0'

  ! Forward propagation outputs trajectory, action, M
  call PropagateFwd(Initialq,Initialp,qpfwd,Sfwd,Mfwd,flag)
  
  ! Check for broken trajectory
  if (flag.eq.1) then
     bad = bad + 1
     goto 222
  endif
   
  ! Reset for second trajectory of pair
  flag = 0
   
  ! Propagate second trajectory of pair
  call PropagateFwd(q0p,p0p,qpbck,Sbck,Mfprime,flag)
  
  ! Check for broken trajectory
  if (flag.eq.1) then
     bad = bad + 1
     goto 222
  endif

  ! Remaining terms of <z0|A|z0'> after division by
  ! sampling distribution, most for imatginary part
  oR = cdexp(0.5d0*Iu*(Initialp(3)+p0(3))*(Initialq(3)-q0(3)))*&
       cdexp(-0.5d0*Iu*(p0p(3)+p0(3))*(q0p(3)-q0(3)))

  o1 = 0.5d0*(Initialq(1)-Iu*Initialp(1))*(q0p(1)+Iu*p0p(1))*&
       cdexp(0.5d0*Iu*Initialp(1)*Initialq(1)-0.5d0*Iu*q0p(1)*p0p(1))/&
       dsqrt(Initialq(1)**2.d0+Initialp(1)**2.d0)/&
       dsqrt(q0p(1)**2.d0+p0p(1)**2.d0)

  o2 = cdexp(0.5d0*Iu*Initialp(2)*Initialq(2)-0.5d0*Iu*p0p(2)*q0p(2))


  ! Prefactor and phase tracking
  do j = 0, Ntime

      Monodromybck(:,:,1) =  transpose(Mfprime(:,:,4,j))
      Monodromybck(:,:,2) = -transpose(Mfprime(:,:,2,j))
      Monodromybck(:,:,3) = -transpose(Mfprime(:,:,3,j))
      Monodromybck(:,:,4) =  transpose(Mfprime(:,:,1,j))
      
      MonodromyFwd(:,:,:) = Mfwd(:,:,:,j)

      call MQCprefactor(MonodromyFwd,MonodromyBck,coeff)

      ! Track the Maslov index 
      !-- cross the -x axis
      if (dble(coeff).lt.0.d0 .and. dimag(coeff)*dimag(prev).lt.0.d0) Maslov = Maslov + 1
      prev = coeff

  enddo

  ! Overlap of electronic coherent states at final time step
  !-- usage: elecoverlap(q,p,qprime,pprime,overlap)
  call elecoverlap(qpfwd(1,1,Ntime),qpfwd(2,1,Ntime),qpbck(1,1,Ntime),qpbck(2,1,Ntime),overlap1)
  call elecoverlap(qpfwd(1,2,Ntime),qpfwd(2,2,Ntime),qpbck(1,2,Ntime),qpbck(2,2,Ntime),overlap2)
  overlap12 = overlap1*overlap2

  ! Compute correlation function
  TCF = (-1)**Maslov*cdsqrt(coeff)*cdexp(Iu*(Sfwd(Ntime)-Sbck(Ntime)))*&
        overlap12*oR*o1*o2/((2.d0*pi)**(2.d0*Ndof))

  ! Inclue matrix elements of B
  do i = 0, Npf

        ! <Rt' PRt'|delta(Pf-P)|Rt PRt>
        call operB(grid(i),qpfwd(:,3,Ntime),qpbck(:,3,Ntime),elemB)

        probs(i) = probs(i) + TCF*elemB

  enddo
222 continue
enddo

! Collect data after paralellization 
call mpi_barrier(mpi_comm_world,ierr)
call mpi_reduce(dble(probs),realtcf,Npf+1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(dimag(probs),imagtcf,Npf+1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(icount1,icount12,4,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(icount2q,icount2q2,2,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(icount2p,icount2p2,2,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(icountNq,icountNq2,2,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(icountNp,icountNp2,2,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(bad,badtot,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)

if (myrank.eq.0) then

do i = 0, Npf
   realtcf(i) = realtcf(i)/dble(NumberMC-badtot)
enddo

! Normalize statistics to obtain probability distribution
summ = sum(realtcf)
do i = 0, Npf
   realtcf(i) = realtcf(i)/dpf/summ
enddo

open(100,file='Broken_Traj.out',status='unknown')
open(101,file='MC_Acceptance.out',status='unknown')
open(102,file='distribution.out',status='unknown')

   write(100,*) 'Number of broken trajectories: ', badtot 
   do i = 0, Npf
        write(102,*) grid(i), realtcf(i)
   enddo
   
   write(*,*) 'MC Acceptance Rates'

   write(101,*) 'Electronic State 1:',   (dble(icount12(j))/dble(4*NumberMC)*100,j=1,4)
   write(101,*) 'Electronic State 2(p)', (dble(icount2q2(j))/dble(4*NumberMC)*100,j=1,2)
   write(101,*) 'Electronic State 2(q)', (dble(icount2p2(j))/dble(4*NumberMC)*100,j=1,2)
   write(101,*) 'Nuclear coordinate(q)', (dble(icountNq2(j))/dble(4*NumberMC)*100,j=1,2)
   write(101,*) 'Nuclear coordinate(p)', (dble(icountNp2(j))/dble(4*NumberMC)*100,j=1,2)

   write(*,*) 'Electronic State 1:',   (dble(icount12(j))/dble(4*NumberMC)*100,j=1,4)
   write(*,*) 'Electronic State 2(p)', (dble(icount2q2(j))/dble(4*NumberMC)*100,j=1,2)
   write(*,*) 'Electronic State 2(q)', (dble(icount2p2(j))/dble(4*NumberMC)*100,j=1,2)
   write(*,*) 'Nuclear coordinate(q)', (dble(icountNq2(j))/dble(4*NumberMC)*100,j=1,2)
   write(*,*) 'Nuclear coordinate(p)', (dble(icountNp2(j))/dble(4*NumberMC)*100,j=1,2)

close(100)
close(101)
close(102)

endif 

call mpi_finalize(ierr)
print *, 'succesflluy done'

end program CF

