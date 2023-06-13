!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute pressure (updated 20/30/16 from croco)
!! for NONLIN_EOS only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# define NONLIN_EOS
# define CROCO_QH

      subroutine get_pressure_croco_QH(Lm,Mm,N,T,S,
     &         z_r,z_w,rho0,pm,pn,rmask,
     &         u,v,e,cosa,sina,
     &         P)

      implicit none

      integer Lm,Mm,N, istrR,iendR,jstrR,jendR, i,j,k,
     &      istr,iend,jstr,jend,istrU,jstrV

      real*8 T(0:Lm+1,0:Mm+1,N),  S(0:Lm+1,0:Mm+1,N),
     &      rho1(0:Lm+1,0:Mm+1,N), qp1(0:Lm+1,0:Mm+1,N),
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      u(1:Lm+1,0:Mm+1,N), v(0:Lm+1,1:Mm+1,N),
     &      e(0:Lm+1,0:Mm+1),
     &      cosa(0:Lm+1,0:Mm+1), sina(0:Lm+1,0:Mm+1),
     &      rmask(0:Lm+1,0:Mm+1),
     &      Tt,Ts,sqrtTs, rho0,  K0, dr00,
     &      cff ,cff1, cff2, cfr, HalfGRho, GRho,
     &      var1, var2,var3, var4

      real*8 P(0:Lm+1,0:Mm+1,N),
     &      rho(0:Lm+1,0:Mm+1,N), dpth,
     &      dR(0:Lm+1,0:N), dZ(0:Lm+1,0:N)

      real*8, parameter :: OneFifth=0.2, OneTwelfth=1./12., epsil=0.

      real*8, parameter :: r00=999.842594,   r01=6.793952E-2,
     &          r02=-9.095290E-3, r03=1.001685E-4,  r04=-1.120083E-6,
     &                                              r05=6.536332E-9,
     &          r10=0.824493,     r11=-4.08990E-3,  r12=7.64380E-5,
     &                            r13=-8.24670E-7,  r14=5.38750E-9,
     &          rS0=-5.72466E-3,  rS1=1.02270E-4,   rS2=-1.65460E-6,
     &          r20=4.8314E-4,
     &          K00=19092.56,     K01=209.8925,     K02=-3.041638,
     &                            K03=-1.852732e-3, K04=-1.361629e-5,
     &          K10=104.4077,     K11=-6.500517,    K12=0.1553190,
     &                                              K13=2.326469e-4,
     &          KS0=-5.587545,    KS1=+0.7390729,   KS2=-1.909078e-2,
     &          qp2=0.0000172

#  ifdef DUKO_2001
      real rho1_0, K0_Duk
#  endif
# include "scalars.h"



Cf2py intent(in) Lm,Mm,N,T,S,z_r,z_w,rho0,pm,pn,rmask,u,v,e,cosa,sina
Cf2py intent(out) P


!
!======================================================================
! Compute density anomaly via Equation Of State (EOS) for seawater.
! Following Jackett and McDougall, 1995, physical EOS is assumed to
! have form
!
!                           rho0 + rho1(T,S)
!           rho(T,S,z) = ------------------------                 (1)
!                         1 - 0.1*|z|/K(T,S,|z|)
!
! where rho1(T,S) is sea-water density perturbation [kg/m^3] at
! standard pressure of 1 Atm (sea surface); |z| is absolute depth,
! i.e. distance from free-surface to the point at which density is
! computed, and
!
!     K(T,S,|z|) = K00 + K01(T,S) + K1(T,S)*|z| + K2(T,S)*|z|^2.  (2)
!
! To reduce errors of pressure-gradient scheme associated with
! nonlinearity of compressibility effects, as well as to reduce
! roundoff errors, the dominant part of density profile,
!
!                           rho0
!                     ----------------                            (3)
!                      1 - 0.1|z|/K00
!
! is removed from from (1). [Since (3) is purely a function of z,
! it does not contribute to pressure gradient.]  This results in
!
!                   rho1 - rho0*[K01+K1*|z|+K2*|z|^2]/[K00-0.1|z|]
!    rho1 + 0.1|z| -----------------------------------------------
!                        K00 + K01 + (K1-0.1)*|z| + K2*|z|^2
!                                                                 (4)
! which is suitable for pressure-gradient calculation.
!
! Optionally, if CPP-switch SPLIT_EOS is defined, term proportional
! to |z| is linearized using smallness 0.1|z|/[K00 + K01] << 1 and
! the resultant EOS has form
!
!              rho(T,S.z) = rho1(T,S) + qp1(T,S)*|z|               (5)
!
! where
!                            rho1 - rho0*K01(T,S)/K00
!              qp1(T,S)= 0.1 --------------------------            (6)
!                                K00 + K01(T,S)
!
! is stored in a special array.
!
! This splitting allows representation of spatial derivatives (and
! also differences) of density as sum of adiabatic derivatives and
! compressible part according to
!
!         d rho     d rho1           d qp1             d |z|
!        ------- = -------- + |z| * -------  +  qp1 * -------      (7)
!         d x,s      d x,s           d x,s            d x,s
!
!                  |<----- adiabatic ----->|   |<- compress ->|
!
! so that constraining of adiabatic derivative for monotonicity is
! equivalent to enforcement of physically stable stratification.
! [This separation and constraining algorithm is subsequently used
! in computation of pressure gradient within prsgrd32ACx-family
! schemes.]
!
! If so prescribed compute the Brunt-Väisäla frequency [1/s^2] at
! horizontal RHO-points and vertical W-points,
!
!                          g    d rho  |
!             bvf^2 = - ------ ------- |                          (8)
!                        rho0    d z   | adiabatic
!
! where density anomaly difference is computed by adiabatically
! rising/lowering the water parcel from RHO point above/below to
! the W-point depth at "z_w".
!
! References:
! ----------
!  Shchepetkin, A.F., McWilliams, J.C., 2003: A method for computing
!  horizontal pressure-gradient force in an oceanic model with a
!  non-aligned vertical coordinate. J. Geophys. Res. 108 (C3), 3090.
!
!======================================================================
!

        istr=1
        istrU=1
        iend=Lm
        jstr=1
        jstrV=1
        jend=Mm

        istrR=istr
        iendR=iend
        jstrR=jstr
        jendR=jend


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute density (copy-pasted from rho_eos.F CROCO 20/03/16)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#  ifdef DUKO_2001
      Tt=3.8D0
      Ts=34.5D0
      sqrtTs=sqrt(Ts)
      K0_Duk= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &       +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &            +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
#  endif


!
!  compute rho as a perturbation to rho0 (at the surface)
!
      dr00=r00-rho0
!

!---------------------------------------------------------------------------------------

      do j=jstrR,jendR

!---------------------------------------------------------------------------------------
        do k=1,N
          do i=istrR,iendR
            Tt=T(i,j,k)
            Ts=max(S(i,j,k), 0.)
            sqrtTs=sqrt(Ts)
            rho1(i,j,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     &                                           r04+Tt*r05 ))))
     &                         +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     &                                            r13+Tt*r14 )))
     &                              +sqrtTs*(rS0+Tt*(
     &                                   rS1+Tt*rS2 ))+Ts*r20 ))
#  ifdef MASKING
     &                                                *rmask(i,j)
#  endif

            K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &         +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &              +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
 


#  ifdef SPLIT_EOS
#   ifdef DUKO_2001
            qp1(i,j,k)= 0.1D0*(rho0+rho1(i,j,k))*(K0_Duk-K0)
     &                               /((K00+K0)*(K00+K0_Duk))
#   else
            qp1(i,j,k)=0.1*(K00*rho1(i,j,k)-rho0*K0)/(K00*(K00+K0))
#   endif
#   ifdef MASKING
     &                                                  *rmask(i,j)
#   endif
            dpth=z_w(i,j,N)-z_r(i,j,k)
            rho(i,j,k)=rho1(i,j,k) +qp1(i,j,k)*dpth*(1.-qp2*dpth)

#  else /* ! SPLIT_EOS */

            K1=B00+Tt*(B01+Tt*(B02+Tt*B03)) +Ts*( B10+Tt*( B11
     &                                    +Tt*B12 )+sqrtTs*BS1 )

            K2=E00+Tt*(E01+Tt*E02) +Ts*(E10+Tt*(E11+Tt*E12))

            dpth=z_w(i,j,N)-z_r(i,j,k)
            cff=K00-0.1*dpth
            cff1=K0+dpth*(K1+K2*dpth)
            rho(i,j,k)=( rho1(i,j,k)*cff*(K00+cff1)
     &                           -0.1*dpth*rho0*cff1
     &                            )/(cff*(cff+cff1))
#  endif /* SPLIT_EOS */
#  ifdef MASKING
            rho(i,j,k)=rho(i,j,k)*rmask(i,j)
#  endif

          enddo
        enddo
      enddo    ! <-- j


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ifdef CROCO_QH
!======================================================================
!  Quasi-hydrostatique correction for non-traditional Coriolis force
!======================================================================
!         dR = -rho0/g* e (U cos(a) - V sin(a) )
!         with e = 2 Omega cos(Phi)
!              a = angle between North and meridional grid axis
!  --> QH pressure gradient is DPdz=-(rho+dR)*g/rho0
!-----------------------------------------------------------------------

      cff=0.5*rho0/g
      do j=jstr,jend
        do k=1,N
          do i=istr,iend
            rho1(i,j,k)=rho1(i,j,k)
     &                         - cff*e(i,j)* (
     &          cosa(i,j)*(u(i,j,k)+u(i+1,j,k))
     &        - sina(i,j)*(v(i,j,k)+v(i,j+1,k)) )
#  ifdef MASKING
            rho1(i,j,k)=rho1(i,j,k)*rmask(i,j)
#  endif
          enddo
        enddo
      enddo


      do j=jstrR,jendR ! resume j loop
        do k=1,N
          do i=istrR,iendR
#   ifdef SPLIT_EOS
            dpth=z_w(i,j,N)-z_r(i,j,k)
            rho(i,j,k)=rho1(i,j,k) +qp1(i,j,k)*dpth*(1.-qp2*dpth)
#   else
            K1=B00+Tt*(B01+Tt*(B02+Tt*B03)) +Ts*( B10+Tt*( B11
     &                                    +Tt*B12 )+sqrtTs*BS1 )
            K2=E00+Tt*(E01+Tt*E02) +Ts*(E10+Tt*(E11+Tt*E12))

            dpth=z_w(i,j,N)-z_r(i,j,k)
            cff=K00-0.1*dpth
            cff1=K0+dpth*(K1+K2*dpth)
            rho(i,j,k)=( rho1(i,j,k)*cff*(K00+cff1)
     &                           -0.1*dpth*rho0*cff1
     &                            )/(cff*(cff+cff1))
#   endif /* SPLIT_EOS */
#  ifdef MASKING
            rho(i,j,k)=rho(i,j,k)*rmask(i,j)
#  endif
          enddo
        enddo
      enddo    ! <-- j
# endif /* CROCO_QH */


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Preliminary step (same for XI- and ETA-components:
!------------ ---- ----- --- --- --- ---------------
!
      GRho=g/rho0
      HalfGRho=0.5*GRho

      do j=jstrV-1,jend
        do k=1,N-1
          do i=istrU-1,iend
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
# ifdef SPLIT_EOS 
            dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))

            dR(i,k)=rho1(i,j,k+1)-rho1(i,j,k)            ! Elementary
     &              +(qp1(i,j,k+1)-qp1(i,j,k))           ! adiabatic
     &                     *dpth*(1.-qp2*dpth)           ! difference
# else
            dR(i,k)=rho(i,j,k+1)-rho(i,j,k)
# endif
          enddo
        enddo
        do i=istrU-1,iend
          dR(i,N)=dR(i,N-1)
          dR(i,0)=dR(i,1)
          dZ(i,N)=dZ(i,N-1)
          dZ(i,0)=dZ(i,1)
        enddo
        do k=N,1,-1               !--> irreversible
          do i=istrU-1,iend
            cff=2.*dZ(i,k)*dZ(i,k-1)
            dZ(i,k)=cff/(dZ(i,k)+dZ(i,k-1))

            cfr=2.*dR(i,k)*dR(i,k-1)
            if (cfr.gt.epsil) then
              dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0.
            endif
# ifdef SPLIT_EOS
            dpth=z_w(i,j,N)-z_r(i,j,k)
            dR(i,k)=dR(i,k)  -qp1(i,j,k)*dZ(i,k)*(1.-2.*qp2*dpth)
# endif
          enddo
        enddo

        do i=istrU-1,iend
          P(i,j,N)=g*z_w(i,j,N) + GRho*( rho(i,j,N)
     &            +0.5*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))
     &              /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
        do k=N-1,1,-1
          do i=istrU-1,iend
            P(i,j,k)=P(i,j,k+1)+HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))
     &                                    *(z_r(i,j,k+1)-z_r(i,j,k))
     &     -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     &                             -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )
     &                -(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)
     &                             -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
     &                                                            ))
          enddo
        enddo
      enddo   !<-- j








      return
      end

 
