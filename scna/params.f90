! ##################################################
! ##################################################
!       This file contains the global module
!       as well as the routine that reads the
!       input file and initializes parameters
!       and arrays.
! ##################################################
! ##################################################
! Last edited : 12/12/17
module parameters
implicit none

! Define pi and the imaginary unit
real*8, parameter     :: pi = dacos(-1.d0)
complex*16, parameter :: Iu = (0.d0,1.d0)

! Number of trajectory pairs
integer :: NumberMC

! MC jump parameters (state 1, state2, nuclear pos, nuclear mom)
real*8  :: jump1, jump2, jumpNq, jumpNp

! Propagator parameters (number timesteps, timestep, energy tolerance)
integer :: Ntime
real*8  :: TimeStep, EnergyTolerance

! Full system dimensionality ( #nuc + #el )
integer :: Ndof  
integer :: Ndof_n
integer :: Ndof_e

! System parameters
real*8  :: Mass
real*8  :: v0, a_tanh, gaus_pref, b_exp

! Grid of final nuclear momenta (min, max, number of grid points)
real*8  :: pfmin, pfmax
integer :: Npf

! Monodromy matrix elements
real*8, allocatable :: Mpp(:,:), Mqq(:,:), Mpq(:,:), Mqp(:,:)

! Coherent state widths, tuning parameters, center of initial wavepacket 
real*8, allocatable :: Width(:,:), WidthT(:,:), Tuningq(:,:), Tuningp(:,:), InverseWidth(:,:)
real*8, allocatable :: InverseWidthT(:,:), InverseTuningq(:,:), InverseTuningp(:,:)
real*8, allocatable :: p0(:), q0(:)

end module

! Read the input file
subroutine input
use parameters
implicit none

integer :: i, j, k
character*75 infostr

open(555,file='input',status='old')

read(555,'(a75)') infostr
read(555,*) Ndof

allocate(p0(Ndof),q0(Ndof))
allocate(Width(Ndof,Ndof),Tuningq(Ndof,Ndof),Tuningp(Ndof,Ndof),InverseWidth(Ndof,Ndof))
allocate(WidthT(Ndof,Ndof),InverseWidthT(Ndof,Ndof))
allocate(InverseTuningq(Ndof,Ndof),InverseTuningp(Ndof,Ndof))
allocate(Mpp(Ndof,Ndof),Mqq(Ndof,Ndof),Mpq(Ndof,Ndof),Mqp(Ndof,Ndof))

read(555,'(a75)') infostr
read(555,*) Mass

Width           = 0.d0
WidthT          = 0.d0
Tuningq         = 0.d0
Tuningp         = 0.d0
InverseWidth    = 0.d0
InverseWidthT   = 0.d0
InverseTuningq  = 0.d0
InverseTuningp  = 0.d0

read(555,'(a75)') infostr
do k = 1, Ndof
read(555,*) q0(k), p0(k), Width(k,k), WidthT(k,k), Tuningq(k,k), Tuningp(k,k)
enddo

read(555,'(a75)') infostr
read(555,*) TimeStep, Ntime, EnergyTolerance

read(555,'(a75)') infostr
read(555,*) NumberMC

read(555,'(a75)') infostr
read(555,*) pfmin, pfmax, Npf

read(555,'(a75)') infostr
read(555,*) v0, a_tanh, gaus_pref, b_exp

! MonteCarlo jump parameters
open(556,file='JumpValues',status='old')

read(556,*) jump1
read(556,*) jump2
read(556,*) jumpNq
read(556,*) jumpNp

close(555)
close(556)

! Inverse width and tuning matrices
do i = 1, Ndof
   InverseWidth(i,i)    = 1.d0/Width(i,i)
   InverseWidthT(i,i)   = 1.d0/WidthT(i,i)
   InverseTuningq(i,i)  = 1.d0/Tuningq(i,i)
   InverseTuningp(i,i)  = 1.d0/Tuningp(i,i)
enddo

end subroutine input
