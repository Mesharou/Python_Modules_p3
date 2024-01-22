!
!----------------------------------------------------------
! Compute vertical advective fluxes
! using 5th-order WENO scheme
!----------------------------------------------------------
!
          do k=3,N-3
            do i=Istr,Iend
              FC(i,k)=We(i,j,k)*
# ifdef PREDICTOR
     &              flux6(
# else
     &              flux5_weno(
# endif
     &             t(i,j,k-2,itrc), t(i,j,k-1,itrc),
     &             t(i,j,k  ,itrc), t(i,j,k+1,itrc),
     &             t(i,j,k+2,itrc), t(i,j,k+3,itrc), We(i,j,k))
            enddo
          enddo

          do i=Istr,Iend
            FC(i,2)=We(i,j,2)*
# ifdef PREDICTOR
     &              flux4(
# else
     &              flux3_weno(
# endif
     &             t(i,j,1,itrc), t(i,j,2,itrc),
     &             t(i,j,3,itrc), t(i,j,4,itrc), We(i,j,2))

            FC(i,N-2)=We(i,j,N-2)*
# ifdef PREDICTOR
     &              flux4(
# else
     &              flux3_weno(
# endif
     &             t(i,j,N-3,itrc), t(i,j,N-2,itrc),
     &             t(i,j,N-1,itrc), t(i,j,N  ,itrc), We(i,j,N-2))

            FC(i,1  )=We(i,j,1)*
# ifdef PREDICTOR
     &              flux2(
# else
     &              flux1(
# endif
     &             t(i,j,1  ,itrc),
     &             t(i,j,2  ,itrc), We(i,j,1  ), 1.)

            FC(i,N-1)=We(i,j,N-1)*
# ifdef PREDICTOR
     &              flux2(
# else
     &              flux1(
# endif
     &             t(i,j,N-1,itrc),
     &             t(i,j,N  ,itrc), We(i,j,N-1), 1.)

            FC(i,0)=0.
            FC(i,N )=0.
            CF(i,0)=dt*pm(i,j)*pn(i,j)
          enddo
