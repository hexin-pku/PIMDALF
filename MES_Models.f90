!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- MES Models (Copyright Xin He <1500011805@pku.edu.cn>)
!-- It contains:
!------ 1) Seven electronic states coupling to Morse potential
!------ 2) Double electronic state -- Spin Boson Model couplinf to a set of bath
!------ 3) Seven electronic states coupling to a set of Morse potential bath
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


module MES_Models
use base
implicit none
    !-- (1) First model
    !-- for Seven electronic states coupling to Morse potential
    integer, parameter :: SES_vdim = 7
    integer, parameter :: SES_rdim = 1
    real(8) :: SES_M = 1.0
    real(8), dimension(SES_vdim) :: SES_D, SES_W, SES_Req, SES_V0
    real(8), dimension(SES_vdim,SES_vdim) :: SES_C, SES_R
    !-- alpha1, alpha2 all same for every coupling
    real(8) :: SES_a1, SES_a2

    !-- (2) Second model
    !-- for double-ES Spin Boson model
    integer, parameter :: DSB_vdim = 2
    integer, parameter :: DSB_rdim = 100
    real(8), dimension(DSB_vdim) :: DSB_V0
    real(8), dimension(DSB_rdim) :: DSB_M
    real(8), dimension(DSB_vdim,DSB_rdim) :: DSB_Req, DSB_W
    real(8) :: DSB_Delta

    !-- (3) third model
    !-- for Seven electronic states coupling to multi Morse potential
    integer, parameter :: mSES_vdim = 7
    integer, parameter :: mSES_rdim = 50
    real(8), dimension(mSES_vdim) :: mSES_V0
    real(8), dimension(mSES_rdim) :: mSES_M
    real(8), dimension(mSES_vdim, mSES_rdim) :: mSES_D, mSES_Req, mSES_W
    real(8), dimension(mSES_vdim, mSES_vdim) :: mSES_C, mSES_R
    !-- alpha1, alpha2 all same for every coupling
    real(8) :: mSES_a1, mSES_a2

contains
    subroutine init_SevenES_cMors()
    implicit none
        integer :: i,j
        real(8), dimension(7,7) :: HFMO, RFMO

        !-- data ref. [66] M. H. Cho, H. M. Vaswani, T. Brixner, J. stenger, and G. R. Fleming
        !------ J. Phys. Chem. B 109(21), 10542-10556(2005)
        data ((HFMO(i,j),j=1,7),i=1,7) / 0., -62., 17., 8., -1., -9., 28., &
                                        -62., 175., -57., -5., -70., -19., 6., &
                                        17., -57., 260., -4., -2., 32., 1., &
                                        8., -5., -4., 280., 6., -8., -106., &
                                        -1., -70., -2., 6., 320., 40., 2., &
                                        -9., -19., 32., -8., 40., 360., 13., &
                                        28., 6., 1., -106., 2., 13., 420. /
        !-- convert to atom unit
        HFMO = HFMO / wn_2_hetree

        data ((RFMO(i,j),j=1,7),i=1,7) / 0., 224., 157., 122., 99., 93., 94., &
                                        224., 0., 102., 67., 65., 66., 72., &
                                        157., 102., 0., 34., 48., 54., 65., &
                                        112., 67., 34., 0., 62., 65., 76., &
                                        99., 65., 48., 62., 0., 69., 84., &
                                        93., 66., 54., 65., 69., 0., 100., &
                                        94., 72., 65., 76., 84., 100., 0. /
        !-- initialize parameters
        do i=1,SES_vdim
            SES_V0(i) = HFMO(i,i)
            SES_R(i,i) = 0.0
            do j=i+1,SES_vdim
                SES_C(i,j) = HFMO(i,j)
                SES_C(j,i) = HFMO(i,j)
                SES_R(i,j) = RFMO(i,j)
                SES_R(j,i) = RFMO(i,j)
            enddo
        enddo

        SES_D = (/7.28e-2, 7.17e-2, 7.06e-2, 6.95e-2, 6.84e-2, 6.73e-2, 6.62e-2/)
        SES_W = (/212.2, 209.0, 205.8, 202.7, 199.5, 196.3, 193.1/)/ wn_2_hetree
        SES_Req = (/0., 5., 10., 15., 20., 25., 30./)
        SES_a1 = 5.e-5
        SES_a2 = 0.02
    end subroutine init_SevenES_cMors

    subroutine SevenES_cMors_V(V, Vd, Vnd, R)
    implicit none
        real(8), dimension(SES_rdim), intent(in) :: R
        real(8), dimension(SES_rdim,SES_vdim,SES_vdim), intent(out) :: V,Vd,Vnd
        !-- for 3N=1, here is only one freedom of nuclues
        Vd =0.0
        Vnd=0.0
        integer :: i,j
        do i=1,SES_vdim
            V(:,i,i) = SES_D(i)*(1.0-dexp(-dsqrt(0.5*SES_M*SES_W(i)**2/SES_D(i))*(R-SES_Req(i))) )**2 + &
                SES_V0(i)
            Vd(:,i,i)=V(:,i,i)
            do j=i+1,SES_vdim
                V(:,i,j) = SES_C(i,j) * dexp( - SES_a1*(R-SES_R(i,j))**2 ) * dcos(SES_a2*(R-SES_R(i,j)))
                V(:,j,i) = V(:,i,j)
                Vnd(:,i,j) = V(:,i,j)
                Vnd(:,j,i) = V(:,i,j)
            enddo
        enddo
    end subroutine SevenES_cMors_V

    subroutine SES_dV_diag(dV_diag, R)
    implicit none
        real(8), dimension(SES_rdim), intent(in) :: R
        real(8), dimension(SES_rdim) :: eterm
        real(8), dimension(SES_rdim, SES_vdim, SES_vdim), intent(out) :: dV_diag
        dv_diag = 0.0
        do i=1,SES_vdim
            eterm = dexp(-dsqrt(0.5*SES_M*SES_W(i)**2/SES_D(i))*(R-SES_Req(i)))
            dV_diag(:,i,i) = 2.0*(1.0-eterm) * eterm * dsqrt(0.5*SES_M*SES_W(i)**2 * SES_D(i))
        enddo
    end subroutine SES_dV_diag

    subroutine SES_trace_Z(Z, R)
    implicit none
        real(8), dimension(SES_rdim), intent(in) :: R
        real(8), intent(out) :: Z
        real(8), dimension(SES_vdim, SES_vdim) :: V,Vd,Vnd
        call SevenES_cMors_V(V,Vd,Vnd,R)

    end subroutine SES_trace_Z

    subroutine init_DoubleSB()
        integer :: i,j
        !-- prev initialize the parameters
        real(8) :: delta, epsilonau, beta, omegac, alpha
        real(8), dimension(DSB_rdim) :: omegaj, cj

        delta = 0.1
        epsilonau = 0.0
        beta = 1
        omegac = 0.25
        alpha = 0.09
        !-- for Model I :   / 0.1, 0, 1,  0.25 , 0.09/
        !-- for Model II :  / 0.1, 0, 50, 0.25 , 0.09/
        !-- for Model III : / 0.1, 1, 50, 0.25 , 0.1 /
        !-- for Model IV :  / 0.1, 1, 2.5,0.1 ,  0.1 /
        !-- for Model V :   / 0.1, 0, 10, 0.1 ,  2   /
        !-- for Model VI :  / 0.1, 5, 1,  0.2 ,  4   /

        DSB_Delta = delta
        do i=1,DSB_rdim
            DSB_M(i) = 1.0
            omegaj(i) = - omegac*dlog(1._8-real(i)/(DSB_rdim+1))
            cj(i) = dsqrt(alpha*omegaj(i)**2*omegac/real(DSB_rdim+1))
        enddo
        do i=1,DSB_vdim
            DSB_V0(i) = (-1)**(i+1)*epsilonau
            do j=1,DSB_rdim
                DSB_W(i,j) = omegaj(j)
                if(mod(i,2).eq. 0) then
                    DSB_Req(i,j) = - cj(j)/(DSB_M(j)*omegaj(j)**2)
                else
                    DSB_Req(i,j) = cj(j)/(DSB_M(j)*omegaj(j)**2)
                endif
                DSB_V0(i) = DSB_V0(i) - 0.5 * cj(j)/(DSB_M(j)*omegaj(j)**2)
            enddo
        enddo
    end subroutine init_DoubleSB

    subroutine  DoubleSB_V(V, Vd, Vnd, R)
        real(8), dimension(DSB_rdim), intent(in) :: R
        real(8), dimension(DSB_vdim,DSB_vdim), intent(out) :: V,Vd,Vnd
        integer :: i,j,k
        Vd =0.0
        Vnd=0.0
        do i=1,DSB_vdim
            V(i,i) = DSB_V0(i)
            do j=1,DSB_rdim
                V(i,i) = V(i,i) + 0.5*DSB_M(j)*( DSB_W(i,j) * (R(j)-DSB_Req(i,j)) )**2
            enddo
            Vd(i,i)=V(i,i)
            do k=i+1,DSB_vdim
                V(i,k) = DSB_Delta
                V(k,i) = V(i,k)
                Vnd(i,k) = V(i,k)
                Vnd(k,i) = V(i,k)
            enddo
        enddo
    end subroutine DoubleSB_V

    subroutine init_SevenES_mulMors()
    implicit none
        integer :: i,j
        real(8), dimension(7,7) :: HFMO, RFMO
        real(8) :: omegac, eta

        omegac = 212.2 / wn_2_hetree
        eta = 35.0 / wn_2_hetree
        !-- data ref. [66] M. H. Cho, H. M. Vaswani, T. Brixner, J. stenger, and G. R. Fleming
        !------ J. Phys. Chem. B 109(21), 10542-10556(2005)
        data ((HFMO(i,j),j=1,7),i=1,7) / 0., -62., 17., 8., -1., -9., 28., &
                                        -62., 175., -57., -5., -70., -19., 6., &
                                        17., -57., 260., -4., -2., 32., 1., &
                                        8., -5., -4., 280., 6., -8., -106., &
                                        -1., -70., -2., 6., 320., 40., 2., &
                                        -9., -19., 32., -8., 40., 360., 13., &
                                        28., 6., 1., -106., 2., 13., 420. /
        !-- convert to atom unit
        HFMO = HFMO / wn_2_hetree

        data ((RFMO(i,j),j=1,7),i=1,7) / 0., 224., 157., 122., 99., 93., 94., &
                                        224., 0., 102., 67., 65., 66., 72., &
                                        157., 102., 0., 34., 48., 54., 65., &
                                        112., 67., 34., 0., 62., 65., 76., &
                                        99., 65., 48., 62., 0., 69., 84., &
                                        93., 66., 54., 65., 69., 0., 100., &
                                        94., 72., 65., 76., 84., 100., 0. /
        !-- initialize parameters
        mSES_M = 1.0

        do i=1,mSES_vdim
            mSES_V0(i) = HFMO(i,i)
            mSES_R(i,i) = 0.0
            do j=i+1,mSES_vdim
                mSES_C(i,j) = HFMO(i,j)
                mSES_C(j,i) = HFMO(i,j)
                mSES_R(i,j) = RFMO(i,j)
                mSES_R(j,i) = RFMO(i,j)
            enddo
        enddo

        do i=1,mSES_vdim
            do j=1,mSES_rdim
                mSES_W(i,j) = - omegac* (1.-0.015*(i-1)) * dlog(1._8-real(j)/real(mSES_rdim+1))
                mSES_Req(i,j) = 1./(mSES_M(j)*mSES_W(1,j)) * dsqrt(2*eta/real(mSES_rdim+1)) + 5.0*(i-1)
                mSES_D(i,j) = 0.0728*( 1.-0.015*real(i-1))
            enddo
        enddo
        mSES_a1 = 5.e-5
        mSES_a2 = 0.02
    end subroutine init_SevenES_mulMors

    subroutine SevenES_mulMors_V(V, Vd, Vnd, R)
    implicit none
        real(8), dimension(mSES_vdim, mSES_vdim), intent(out) :: V, Vd, Vnd
        real(8), dimension(mSES_rdim), intent(in) :: R
        integer :: i,j,k

        Vd =0.0
        Vnd=0.0
        do i=1,mSES_vdim
            V(i,i) = mSES_V0(i)
            do k=1,mSES_rdim
                V(i,i) = V(i,i) + mSES_D(i,k)* &
                (1.0-dexp(-dsqrt(0.5*mSES_M(k)*mSES_W(i,k)**2/mSES_D(i,k))*(R(k)-mSES_Req(i,k))) )**2
            enddo
            Vd(i,i)=V(i,i)
            do j=i+1,mSES_vdim
                V(i,j) = 0.0
                do k=1,mSES_rdim
                    V(i,j) = V(i,j) + SES_C(i,j) * &
                    dexp( - SES_a1*(R(k)-SES_R(i,j))**2 ) * dcos(SES_a2*(R(k)-SES_R(i,j)))
                enddo
                V(j,i) = V(i,j)
                Vnd(i,j) = V(i,j)
                Vnd(j,i) = V(i,j)
            enddo
        enddo
    end subroutine SevenES_mulMors_V

    subroutine
end module MES_Models
