! $Id: bulk_flux.F 1531 2014-04-17 15:05:26Z marchesiello $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!

      subroutine bulk_stress (Lm, Mm, sst, tair, rhum,
     &                        uwnd, vwnd, wspd, rho0,
     &                        sustr, svstr)


!=======================================================================
!                                                                      !
!  This routine computes the bulk parameterization of surface wind     !
!  stress.                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Fairall, C.W., E.F. Bradley, D.P. Rogers, J.B. Edson and G.S.     !
!      Young, 1996:  Bulk parameterization of air-sea fluxes for       !
!      tropical ocean-global atmosphere Coupled-Ocean Atmosphere       !
!      Response Experiment, JGR, 101, 3747-3764.                       !
!                                                                      !
!    Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B.        !
!      Edson, and G.S. Young, 1996:  Cool-skin and warm-layer          !
!      effects on sea surface temperature, JGR, 101, 1295-1308.        !
!                                                                      !
!    Liu, W.T., K.B. Katsaros, and J.A. Businger, 1979:  Bulk          !
!        parameterization of the air-sea exchange of heat and          !
!        water vapor including the molecular constraints at            !
!        the interface, J. Atmos. Sci, 36, 1722-1735.                  !
!                                                                      !
!  Adapted from COARE code written originally by David Rutgers and     !
!  Frank Bradley.                                                      !
!                                                                      !
!  EMINUSP option for equivalent salt fluxes added by Paul Goodman     !
!  (10/2004).                                                          !
!                                                                      !
!  Modified by Kate Hedstrom for COARE version 3.0 (03/2005).          !
!  Modified by Jim Edson to correct specific hunidities.               !
!                                                                      !
!  References:                                                          !
!                                                                      !
!     Fairall et al., 2003: J. Climate, 16, 571-591.                   !
!                                                                      !
!     Taylor, P. K., and M. A. Yelland, 2001: The dependence of sea    !
!     surface roughness on the height and steepness of the waves.      !
!     J. Phys. Oceanogr., 31, 572-590.                                 !
!                                                                      !
!     Oost, W. A., G. J. Komen, C. M. J. Jacobs, and C. van Oort, 2002:!
!     New evidence for a relation between wind stress and wave age     !
!     from measurements during ASGAMAGE. Bound.-Layer Meteor., 103,    !
!     409-438.                                                         !
!
!  Modified for CROCO by Fred marin (2008)
!  Modified for CROCO by P. Marchesiello (2009)
!                                                                      !
!=======================================================================
!
      implicit none
#include "scalars.h"
      integer Lm, Mm, i,j, Istr,Iend,Jstr,Jend,
     &        imin,imax,jmin,jmax
      integer IterMax,Iter
      real a,cff
      real rho0i,cpi
      real TseaC,TseaK,Qsea
      real TairC,TairK,rhoAir,Qair
      real Q,RH,VisAir,Hlv
      real delW,delT,delQ
      real u10,Zo10,Cd10,Ch10,Ct10,Cd
      real Ct,CC,Ri,Ribcu,Zetu,L10
      real Wstar,Tstar,Qstar
      real ZoW,ZoT,ZoT10,ZoQ,ZoL,L,Rr,Bf
      real Wpsi,Tpsi,Qpsi
      real wspd0,Wgus,charn

      real bulk_psiu,bulk_psit

      real hfsen,hflat,hflw,upvel,evap

      real blk_Rgas,blk_ZW,blk_ZT,blk_ZQ,blk_beta,blk_Cpa,blk_Zabl
      parameter (blk_Rgas=287.1)
      parameter (blk_ZW=10.0)
      parameter (blk_ZT=10.0)
      parameter (blk_ZQ=10.0)
      parameter (blk_Zabl=600.0)
      parameter (blk_beta=1.2)
      parameter (blk_Cpa=1004.67)
      real rhow,patm
      parameter (rhow=1000.0)
      parameter (patm=1010.0)
      real eps
      parameter (eps=1.d-20)
      real(kind=8), parameter :: r3=1.0/3.0

      real aer(0:Lm+1,0:Mm+1),
     &     cer(0:Lm+1,0:Mm+1)

      ! INPUT
      real*8  sst(0:Lm+1,0:Mm+1),  tair(0:Lm+1,0:Mm+1)
     &      ,rhum(0:Lm+1,0:Mm+1),  wspd(0:Lm+1,0:Mm+1)
     &      ,uwnd(1:Lm+1,0:Mm+1),  vwnd(0:Lm+1,1:Mm+1)
     &      ,rho0
      ! OUTPUTS
      real*8  sustr(1:Lm+1,0:Mm+1), svstr(0:Lm+1,1:Mm+1)
!
!----------------------------------------------------------------------

Cf2py intent(in) Lm, Mm, sst, tair, rhum, uwnd, vwnd, wspd, rho0
Cf2py intent(out) sustr, svstr

!----------------------------------------------------------------------

!
      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1
!
!----------------------------------------------------------------------
!  Input variables:
!
!      - radsw/srflx
!      - radlw         #ifdef BULK_LW: downward longwave radiation) 
!                      #else: net outgoing longwave radiation)
!      - prate
!      - sst           t(N,itemp) or sstsk
!      - sss           t(N,isalt)
!      - tair
!      - rhum
!      - wspd
!      - uwnd
!      - vwnd
!
!----------------------------------------------------------------------
!
!  define inverse seawater density, use mean value for seawater density.
      rho0i=1.0/rho0
!
!  set inverse of specific heat for seawater (kg-degC/Joule).
!  cp is defined in scalars.h
      cpi=1.0/cp
!
!
      do j=jmin,jmax
        do i=imin,imax
!
!
!  Input bulk parameterization fields
!
          wspd0=wspd(i,j)
          TairC=tair(i,j)
          TairK=TairC+273.16
          TseaC=sst(i,j)
          TseaK=TseaC+273.16
          RH=rhum(i,j)

!
!-----------------------------------------------------------------------
!  Compute specific humidities (kg/kg).
!
!    note that Qair is the saturation specific humidity at Tair
!                 Q is the actual specific humidity
!              Qsea is the saturation specific humidity at Tsea

!          Saturation vapor pressure in mb is first computed and then
!          converted to specific humidity in kg/kg
!
!          The saturation vapor pressure is computed from Teten formula
!          using the approach of Buck (1981):
!
!          Esat(mb) = (1.0007+3.46E-6*patm)*6.1121*
!                  EXP(17.502*TairC(C)/(240.97+TairC(C)))
!
!          The ambient vapor is found from the definition of the
!          Relative humidity:
!
!          RH = W/Ws*100 ~ E/Esat*100   E = RH/100*Esat if RH is in %
!                                       E = RH*Esat     if RH fractional
!
!          The specific humidity is then found using the relationship:
!
!          Q = 0.622 E/(P + (0.622-1)e)
!
!          Q(kg/kg) = 0.62197*(E(mb)/(patm(mb)-0.378*E(mb)))
!
!-----------------------------------------------------------------------
!
!  Compute air saturation vapor pressure (mb), using Teten formula.
!
          cff=(1.0007+3.46e-6*patm)*6.1121*
     &        exp(17.502*TairC/(240.97+TairC))

!
!  Compute specific humidity, Q (kg/kg).
!
          if (RH.lt.2.0) then                          !RH fraction
            cff=cff*RH                                 !Vapor pres (mb)
            Q=0.62197*(cff/(patm-0.378*cff))           !Spec hum (kg/kg)
          else          !RH input was actually specific humidity in g/kg
            Q=RH/1000.0                                !Spec Hum (kg/kg)
          endif
!
!  Compute water saturation vapor pressure (mb), using Teten formula.
!
          cff=(1.0007+3.46e-6*patm)*6.1121*
     &            exp(17.502*TseaC/(240.97+TseaC))
!
!  Vapor Pressure reduced for salinity (Kraus & Businger, 1994, pp 42).
!
          cff=cff*0.98
!
!  Compute Qsea (kg/kg) from vapor pressure.
!
          Qsea=0.62197*(cff/(patm-0.378*cff))
!
!-----------------------------------------------------------------------
!  Compute Monin-Obukhov similarity parameters for wind (Wstar),
!  heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!-----------------------------------------------------------------------
!
!  Moist air density (kg/m3).
!
          rhoAir=patm*100.0/(blk_Rgas*TairK*
     &                              (1.0+0.61*Q))
!
!  Kinematic viscosity of dry air (m2/s), Andreas (1989).
!
          VisAir=1.326E-5*(1.0+TairC*(6.542E-3+TairC*
     &               (8.301e-6-4.84e-9*TairC)))

!
!  Assume that wind is measured relative to sea surface and include
!  gustiness.
!
          Wgus=0.5
          delW=SQRT(wspd0*wspd0+Wgus*Wgus)
          delQ=Qsea-Q
          delT=TseaC-TairC
!
!  Neutral coefficients.
!
          ZoW=0.0001
          u10=delW*LOG(10.0/ZoW)/LOG(blk_ZW/ZoW)
          Wstar=0.035*u10
          Zo10=0.011*Wstar*Wstar/g+0.11*VisAir/Wstar
          Cd10=(vonKar/LOG(10.0/Zo10))**2
          Ch10=0.00115
          Ct10=Ch10/sqrt(Cd10)
          ZoT10=10.0/exp(vonKar/Ct10)
          Cd=(vonKar/LOG(blk_ZW/Zo10))**2
!
!  Compute Richardson number.
!
          Ct=vonKar/LOG(blk_ZT/ZoT10)  ! T transfer coefficient
          CC=vonKar*Ct/Cd
          Ribcu=-blk_ZW/(blk_Zabl*0.004*blk_beta**3)
          Ri=-g*blk_ZW*(delT+0.61*TairK*delQ)/
     &          (TairK*delW*delW)
          if (Ri.lt.0.0) then
            Zetu=CC*Ri/(1.0+Ri/Ribcu)       ! Unstable
          else
            Zetu=CC*Ri/(1.0+3.0*Ri/CC)      ! Stable
          endif
          L10=blk_ZW/Zetu
          if (Zetu.gt.50.0) then
            IterMax=1
          else
            IterMax=3
          endif
!
!  First guesses for Monin-Obukhov similarity scales.
!
          Wstar=delW*vonKar/(LOG(blk_ZW/Zo10)-
     &                             bulk_psiu(blk_ZW/L10,pi))
          Tstar=-delT*vonKar/(LOG(blk_ZT/ZoT10)-
     &              bulk_psit(blk_ZT/L10,pi))
          Qstar=-delQ*vonKar/(LOG(blk_ZQ/ZoT10)-
     &              bulk_psit(blk_ZQ/L10,pi))
!
!  Modify Charnock for high wind speeds. The 0.125 factor below is for
!  1.0/(18.0-10.0).
!
          if (delW.gt.18.0) then
            charn=0.018
          elseif ((10.0.lt.delW).and.(delW.le.18.0)) then
            charn=0.011+0.125*(0.018-0.011)*(delW-10.)
          else
            charn=0.011
          endif
!
!  Iterate until convergence. It usually converges within four
!  iterations.
!
        do Iter=1,IterMax
          ZoW=charn*Wstar*Wstar/g+0.11*VisAir/(Wstar+eps)
          Rr=ZoW*Wstar/VisAir
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
          ZoQ=MIN(1.15e-4,5.5e-5/Rr**0.6)
          ZoT=ZoQ
          ZoL=vonKar*g*blk_ZW*
     &             (Tstar*(1.0+0.61*Q)+0.61*TairK*Qstar)/
     &             (TairK*Wstar*Wstar*(1.0+0.61*Q)+eps)
          L=blk_ZW/(ZoL+eps)
!
!  Evaluate stability functions at Z/L.
!
          Wpsi=bulk_psiu(ZoL,pi)
          Tpsi=bulk_psit(blk_ZT/L,pi)
          Qpsi=bulk_psit(blk_ZQ/L,pi)
!
!  Compute wind scaling parameters, Wstar.
!
          Wstar=MAX(eps,delW*vonKar/(LOG(blk_ZW/ZoW)-Wpsi))
          Tstar=-delT*vonKar/(LOG(blk_ZT/ZoT)-Tpsi)
          Qstar=-delQ*vonKar/(LOG(blk_ZQ/ZoQ)-Qpsi)
!
!  Compute gustiness in wind speed.
!
          Bf=-g/TairK*Wstar*(Tstar+0.61*TairK*Qstar)
          if (Bf.gt.0.0) then
            Wgus=blk_beta*(Bf*blk_Zabl)**r3
          else
            Wgus=0.2
          endif
        enddo
!
!-----------------------------------------------------------------------
!  Compute Atmosphere/Ocean fluxes.
!-----------------------------------------------------------------------
!
!
!  Compute transfer coefficients for momentum (Cd).
!
          wspd0=SQRT(wspd0*wspd0+Wgus*Wgus)
          Cd=Wstar*Wstar/(wspd0*wspd0+eps)

!
!  Save for the calculation of wind stress
!
          aer(i,j)=rhoAir*wspd0*rho0i
          cer(i,j)=Cd
          !if ( i.eq.26 .and. j.eq.27 ) then
          !    write(*,*) '            '
          !    write(*,*) 'in bulk_stress'
          !    write(*,*) 'rhoAir,wspd0,rho0i', rhoAir,wspd0,rho0i
          !    write(*,*) 'aer(i,j)', aer(i,j)
          !    write(*,*) 'cer(i,j)', cer(i,j)
          !endif

        enddo
      enddo


!
!--------------------------------------------------------------
!  Compute kinematic, surface wind stress (m2/s2).
!--------------------------------------------------------------
!
# ifdef BULK_SM_UPDATE
!
!     sustr = rhoA/rhoW * Cd * wspd*uwnd  --> u point
!
      do j=jmin,jmax
        do i=imin+1,imax 
          a=0.5*(aer(i-1,j)+aer(i,j))
          cff=0.5*(cer(i-1,j)+cer(i,j))
          sustr(i,j)=a*cff*uwnd(i,j)
          !if ( i.eq.26 .and. j.eq.27 ) then
          !    write(*,*) '            '
          !    write(*,*) 'in bulk_stress'
          !    write(*,*) 'a,cff,uwnd(i,j)', a,cff,uwnd(i,j)
          !    write(*,*) 'sustr(i,j)', sustr(i,j)
          !endif
        enddo
      enddo
!
!     svstr = rhoA/rhoW * Cd * wspd*vwnd  --> v point
!
      do j=jmin+1,jmax
        do i=imin,imax
          a=0.5*(aer(i,j-1)+aer(i,j))
          cff=0.5*(cer(i,j-1)+cer(i,j))
          svstr(i,j)=a*cff*vwnd(i,j)
          !if ( i.eq.26 .and. j.eq.27 ) then
          !    write(*,*) 'a,cff,vwnd(i,j)', a,cff,vwnd(i,j)
          !    write(*,*) 'svstr(i,j)', svstr(i,j)
          !endif
        enddo
      enddo
# endif

      return
      end

      FUNCTION bulk_psiu (ZoL, pi)
!
!=======================================================================
!                                                                      !
!  This function evaluates the stability function for  wind speed      !
!  by matching Kansas  and free convection forms.  The convective      !
!  form follows Fairall et al. (1996) with profile constants from      !
!  Grachev et al. (2000) BLM.  The  stable  form is from Beljaars      !
!  and Holtslag (1991).                                                !
!                                                                      !
!=======================================================================
!
!
!  Function result
!
      real bulk_psiu
!
!  Imported variable declarations.
!
      real ZoL, pi
!
!  Local variable declarations.
!
      real r3
      parameter (r3=1.0/3.0)
      real Fw, cff, psic, psik, x, y
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
      if (ZoL.lt.0.0) then
        x=(1.0-15.0*ZoL)**0.25
        psik=2.0*LOG(0.5*(1.0+x))+LOG(0.5*(1.0+x*x))-
     &       2.0*ATAN(x)+0.5*pi
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0)
        y=(1.0-10.15*ZoL)**r3
        psic=1.5*LOG(r3*(1.0+y+y*y))-
     &       cff*ATAN((1.0+2.0*y)/cff)+pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        cff=ZoL*ZoL
        Fw=cff/(1.0+cff)
        bulk_psiu=(1.0-Fw)*psik+Fw*psic
!
!  Stable conditions.
!
      else
        cff=MIN(50.0,0.35*ZoL)
        bulk_psiu=-((1.0+ZoL)+0.6667*(ZoL-14.28)/
     &            EXP(cff)+8.525)
      endif
      return
      END FUNCTION bulk_psiu

      FUNCTION bulk_psit (ZoL, pi)
!
!=======================================================================
!                                                                      !
!  This function evaluates the  stability function  for moisture and   !
!  heat by matching Kansas and free convection forms. The convective   !
!  form follows Fairall et al. (1996) with  profile  constants  from   !
!  Grachev et al. (2000) BLM.  The stable form is from  Beljaars and   !
!  and Holtslag (1991).                                                !
!
!=======================================================================
!
!  Function result
!
      real bulk_psit
!
!  Imported variable declarations.
!
      real ZoL, pi
!
!  Local variable declarations.
!
      real r3
      parameter (r3=1.0/3.0)
      real Fw, cff, psic, psik, x, y
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
      if (ZoL.lt.0.0) then
        x=(1.0-15.0*ZoL)**0.5
        psik=2.0*LOG(0.5*(1.0+x))
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0)
        y=(1.0-34.15*ZoL)**r3
        psic=1.5*LOG(r3*(1.0+y+y*y))-
     &       cff*ATAN((1.0+2.0*y)/cff)+pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        cff=ZoL*ZoL
        Fw=cff/(1.0+cff)
        bulk_psit=(1.0-Fw)*psik+Fw*psic
!
!  Stable conditions.
!
      else
        cff=MIN(50.0,0.35*ZoL)
        bulk_psit=-((1.0+2.0*ZoL)**1.5+
     &            0.6667*(ZoL-14.28)/EXP(cff)+8.525)
      endif
      return
      END FUNCTION bulk_psit


