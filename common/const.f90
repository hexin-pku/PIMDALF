!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- base module for sci-computations 
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module const
implicit none
    !-- precision of float number and length of strings
    integer, parameter :: sp=kind(0.0)
    integer, parameter :: dp=kind(0.d0)
    integer, parameter :: len0=4
    integer, parameter :: len1=20
    integer, parameter :: len2=100
    integer, parameter :: len3=500

    !-- mathematical constant numbers
    real(8), parameter :: pi = 3.14159265358979323846_8
    real(8), parameter :: pi2 = pi**2
    real(8), parameter :: twopi = 2._8 * pi
    real(8), parameter :: sqrtpi = dsqrt(pi)
    real(8), parameter :: eu = 2.71828182845904523536_8
    real(8), parameter :: ec = 0.57721566490153286060_8

    !-- physics constant
    real(8), parameter :: hb = 1.0545718001D-34
    real(8), parameter :: kb = 1.3806490351D-23
    
    !-- atom unit
    !-- quantity(a.u.) == au2si_quantity * quantity(si)
    real(8), parameter :: au2si_mass   = 9.10938291D-31       ! (kg)
    real(8), parameter :: au2si_charge = 1.602176565D-19      ! (C)
    real(8), parameter :: au2si_electric_const = 8.9875517873681D+9 ! ( kg*m^3*s^-2*C^-2 )
    real(8), parameter :: au2si_action = 1.054571726D-34      ! (J*s)
    real(8), parameter :: au2si_length = 5.2917721092D-11     ! (m)
    real(8), parameter :: au2si_momentum = 1.99285188224D-24  ! (kg*m/s)
    real(8), parameter :: au2si_energy = 4.35974417D-18       ! (J)
    real(8), parameter :: au2si_time = 2.418884326505D-17     ! (s)
    real(8), parameter :: au2si_temperatue = 3.157746455D+5   ! (K)
    real(8), parameter :: au2si_beta = 3.166815367D-6         ! (1/K)
    
    !-- others:
    !----- a_in_b : 1(a) == a_in_b (b) 
    !-- energy transfrom
    real(8), parameter :: au_in_wn = 219474.6313702_8
    
    !-- length transfrom
    real(8), parameter :: au_in_ai = 5.2917721092D-1
    
    !-- time transfrom
    real(8), parameter :: au_in_fs = 2.418884326505D-2
    real(8), parameter :: au_in_ps = 2.418884326505D-5
    
end module



