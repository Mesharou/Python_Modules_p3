
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define TS_HADV_C6
#define TS_VADV_SPLINES



      subroutine get_tracer_advection_croco_c6 (Lm,Mm,N,u,v, z_r,z_w
     &  ,pm,pn,dt
     &  ,t,We,rmask
     &  ,TXadv,TYadv,TVadv)    


      integer Lm,Mm,N,NT, i,j,k
     &      ,Istr,Iend,Jstr,Jend
     &      ,imin,imax,jmin,jmax
     &      ,itemp,isalt,dt,itrc

      parameter (NT=2)  
      parameter (itemp=1,isalt=2)  
      parameter (epsil=1.E-16)

      ! INPUTS
      REAL t(0:Lm+1,0:Mm+1,N,NT)

      REAL u(1:Lm+1,0:Mm+1,N),  v(0:Lm+1,1:Mm+1,N)
     &      ,z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N)
     &      ,pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1)
     &      ,We(0:Lm+1,0:Mm+1,0:N), rmask(0:Lm+1,0:Mm+1)

      ! OUTPUTS
      REAL TXadv(0:Lm+1,0:Mm+1,N,NT)
     &      ,TYadv(0:Lm+1,0:Mm+1,N,NT)
     &      ,TVadv(0:Lm+1,0:Mm+1,N,NT)


      ! LOCAL
      REAL FX(0:Lm+1,0:Mm+1), FE(0:Lm+1,0:Mm+1)
     &      ,FC(0:Lm+1,0:N), DC(0:Lm+1,0:N)
     &      ,CF(0:Lm+1,0:N)
     &      ,HZR(0:Lm+1,0:Mm+1,N)
     &      ,dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1)
     &      ,Huon(1:Lm+1,0:Mm+1,N),  Hvom(0:Lm+1,1:Mm+1,N)


!
!--------------------------------------------------------------------
! Definition of flux operators: 1st, 2nd, 3rd, 4th, 5th or 6th order,
! used in UP5 and C6 advection schemes (and order degradation near
! land masks). cdiff is part of laplacian diffusion in flux1 (used
! near mask):
!    0 --> flux1=flux2 (second order C2 advection scheme)
!    1 --> flux1 gives 1st order monotonic UP1 advection scheme
!--------------------------------------------------------------------
!
      REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2
      REAL    :: ua, vel, cdiff, cdif
      REAL    :: flux1, flux2, flux3, flux4, flux5, flux6
      REAL    :: flx2, flx3, flx4, flx5
      REAL    :: mask0, mask1, mask2, mask3

      flux2(q_im1, q_i, ua, cdiff) = 0.5*( q_i + q_im1 )

      flux1(q_im1, q_i, ua, cdiff) = flux2(q_im1, q_i, ua, cdiff) -
     &      0.5*cdiff*sign(1.,ua)*(q_i-q_im1)

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =
     &      ( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua) =
     &      flux4(q_im2, q_im1, q_i, q_ip1, ua) +
     &      sign(1.,ua)*((q_ip1 -
     &      q_im2)-3.*(q_i-q_im1))/12.0

      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =
     &      ( 37.*(q_i+q_im1) - 8.*(q_ip1+q_im2)
     &      +(q_ip2+q_im3) )/60.0

      flux5(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =
     &      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)
     &      -sign(1.,ua)*(
     &      (q_ip2-q_im3)-5.*(q_ip1-q_im2)+10.*(q_i-q_im1) )/60.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,pm,pn,dt,t,W
Cf2py intent(out) TXadv,TYadv,TVadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        Istr=1
        Iend=Lm
        Jstr=1
        Jend=Mm

        imin=0
        imax=Lm+1
        jmin=0
        jmax=Mm+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=jmin,jmax
        do i=imin,imax
          do k=0,N,+1
              We(i,j,k) = We(i,j,k)/(pm(i,j)*pn(i,j))
            enddo
          enddo
        enddo

      do j=jmin,jmax
        do i=imin,imax
          do k=1,N,+1
           HZR(i,j,k) = z_w(i,j,k) - z_w(i,j,k-1)
          enddo
        enddo
      enddo

      do j=jmin,jmax
        do i=imin+1,imax
            dn_u(i,j) = 2./(pn(i,j)+pn(i-1,j))
            do k=1,N,+1
              Huon(i,j,k) = 0.5*(HZR(i,j,k)+HZR(i-1,j,k))*dn_u(i,j)
     &      * u(i,j,k) 
            enddo
          enddo
      enddo

      do j=jmin+1,jmax
        do i=imin,imax
            dm_v(i,j) = 2./(pm(i,j)+pm(i,j-1))
            do k=1,N,+1
              Hvom(i,j,k) = 0.5*(HZR(i,j,k)+HZR(i,j-1,k))*dm_v(i,j)
     &      * v(i,j,k) 
            enddo
          enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This part of the code is
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Centered scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!======================================================================
!
!  Compute horizontal advection for TS_HADV_C6
!
!======================================================================
!
      do itrc=1,NT
        do k=1,N

        imin=3
        imax=Lm-1
        jmin=3
        jmax=Mm-1

#   define FLUX5 flux6
#   define FLUX3 flux4
#   define FLUX2 flux2
#   undef  UP5_MASKING

          cdif=1.
#  include "t3dadv_order5.h"

#  undef FLUX5
#  undef FLUX3
#  undef FLUX2
#  undef UP5_MASKING

          do j=Jstr,Jend
            do i=Istr,Iend
              TXadv(i,j,k,itrc)=-(FX(i+1,j)-FX(i,j))* rmask(i,j)
              TYadv(i,j,k,itrc)=-(FE(i,j+1)-FE(i,j))* rmask(i,j)
            enddo
          enddo

        enddo    !<-- k
      enddo    !<-- itrc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=Jstr,Jend
        do itrc=1,NT
 

# include "compute_vert_tracer_fluxes_croco.h"
 
            do k=1,N
              do i=Istr,Iend
                TVadv(i,j,k,itrc) = -(FC(i,k)-FC(i,k-1))
     &                                 * rmask(i,j)
              enddo
            enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       enddo      ! <-- itrc 
      enddo      ! <-- j




      do itrc=1,NT

!---------------------------------------------------------------
! Compute the tendency term of tracer diagnostics
! Divide all diagnostic terms by the cell volume
! (HZR(i,j,k,itrc)/(pm(i,j).*pn(i,j)). There after the unit
! of diagnostic terms will be: (unit of tracers)* s-1.
!
! Note: the Horizontal mixing term is computed in t3dmix
! where Trate is updated accordingly
!---------------------------------------------------------------

       do k=1,N
         do j=Jstr,Jend
           do i=Istr,Iend

              cff=pm(i,j)*pn(i,j)/HZR(i,j,k)
              TXadv(i,j,k,itrc)=TXadv(i,j,k,itrc)*cff
              TYadv(i,j,k,itrc)=TYadv(i,j,k,itrc)*cff
              TVadv(i,j,k,itrc)=TVadv(i,j,k,itrc)*cff

           enddo
         enddo

       enddo

      enddo  ! <-- itrc

!---------------------------------------------------------------

      return
      end



#undef TS_HADV_C6
#undef TS_VADV_SPLINES
















