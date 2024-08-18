# 1 "R_tools_fort_gigatl.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 31 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 32 "<command-line>" 2
# 1 "R_tools_fort_gigatl.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CROCO ROUTINES (for GIGATL)
!!
!! copied from actual CROCO scripts
!!
!! compile with:
!! "cpp R_tools_fort_gigatl.F R_tools_fort_gigatl.f"
!! "f2py -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m R_tools_fort_gigatl R_tools_fort_gigatl.f" for python use
!!
!! print R_tools_fort_gigatl.rho_eos.__doc__
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 1 "R_tools_fort_routines_gigatl/online_interp.F" 1
! $Id: online_interp.F 1458 2014-02-03 15:01:25Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http:
!======================================================================
!
!
! This is the "interp.F" script (based on the interpolation implemented in the
! Roms Rutgers version)
!------------------------------------------------------------------------------
! This file contains the subfunctions enabling the online interpolation of
! forcing datasets on the simulation domain using linear or cubic approach.
! These functions applied for all discretisations of the domain, MPI or OPENMP.
!------------------------------------------------------------------------------
# 1 "R_tools_fort_routines_gigatl/cppdefs.h" 1
! $Id: cppdefs.h 1628 2015-01-10 13:53:00Z marchesiello $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http:
!======================================================================
!
# 1445 "R_tools_fort_routines_gigatl/cppdefs.h"
# 1 "R_tools_fort_routines_gigatl/cppdefs_dev.h" 1
! $Id: set_global_definitions.h 1616 2014-12-18 14:39:51Z rblod $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http:
!======================================================================
!
# 1446 "R_tools_fort_routines_gigatl/cppdefs.h" 2
# 1 "R_tools_fort_routines_gigatl/set_global_definitions.h" 1
! $Id: set_global_definitions.h 1618 2014-12-18 14:39:51Z rblod $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http:
!======================================================================
!
# 231 "R_tools_fort_routines_gigatl/set_global_definitions.h"
!
# 250 "R_tools_fort_routines_gigatl/set_global_definitions.h"
!! Gc remove because incompatibe with AGRIF
!#elif defined && defined Ifort
!
!# define 16 16
!# define 0.Q0 0.Q0
!! Gc remove because incompatibe with AGRIF
# 303 "R_tools_fort_routines_gigatl/set_global_definitions.h"
!-# define float dfloat
!-# define FLoaT dfloat
!-# define FLOAT dfloat
!-# define sqrt dsqrt
!-# define SQRT dsqrt
!-# define exp dexp
!-# define EXP dexp
!-# define dtanh dtanh
!-# define TANH dtanh
# 1447 "R_tools_fort_routines_gigatl/cppdefs.h" 2
# 21 "R_tools_fort_routines_gigatl/online_interp.F" 2

      SUBROUTINE linterp2d (Lm,Mm,
     & LBx, UBx, LBy, UBy,
     & Xinp, Yinp, Finp,
     & Istr, Iend, Jstr, Jend,
     & Xout, Yout, Fout)
!
!=======================================================================
! !
! Given any gridded 2D field, Finp, this routine linearly interpolate !
! to locations (Xout,Yout). To facilitate the interpolation within !
! any irregularly gridded 2D field, the fractional grid cell indices !
! (Iout,Jout) with respect Finp are needed at input. Notice that the !
! routine "hindices" can be used to compute these indices. !
! !
! On Input: !
! !
! LBx I-dimension lower bound of gridded field, Finp. !
! UBx I-dimension upper bound of gridded field, Finp. !
! LBy J-dimension lower bound of gridded field, Finp. !
! UBy J-dimension upper bound of gridded field, Finp. !
! Xinp X-locations of gridded field, Finp. !
! Yinp Y-locations of gridded field, Finp. !
! Finp 2D field to interpolate from. !
! Istr Starting data I-index to interpolate, Fout. !
! Iend Ending data I-index to interpolate, Fout. !
! Jstr Starting data J-index to interpolate, Fout. !
! Jend Ending data J-index to interpolate, Fout. !
! Xout X-locations to interpolate, Fout. !
! Yout Y-locations to interpolate, Fout. !
! !
! On Output: !
! !
! Fout Interpolated 2D field. !
! !
!=======================================================================
!
!
      implicit none

      integer Lm,Mm
!
! Imported variable declarations.
!
      integer, intent(in) :: LBx, UBx, LBy, UBy
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
      real(kind=8), intent(in) :: Xinp(LBx:UBx)
      real(kind=8), intent(in) :: Yinp(LBy:UBy)
      real(kind=8), intent(in) :: Finp(LBx:UBx,LBy:UBy)
!
      real(kind=8), intent(in) :: Xout(0:Lm+1,0:Mm+1)
      real(kind=8), intent(in) :: Yout(0:Lm+1,0:Mm+1)
!
      real(kind=8), intent(out) :: Fout(0:Lm+1,0:Mm+1)
!
! Local variable declarations.
!
      integer i, i1, i2, j, j1, j2, ii, jj
      real(kind=8) cff, x, x1, x2, y, y1, y2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Cf2py intent(in) Lm,Mm,LBx, UBx, LBy, UBy,Xinp, Yinp, Finp,Istr, Iend, Jstr, Jend,Xout, Yout
Cf2py intent(out) Fout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
! Linearly interpolate requested field
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
! i1=INT(Iout(i,j))
! i2=i1+1
! j1=INT(Jout(i,j))
! j2=j1+1
           DO ii=LBx,(UBx-1)
             if ((Xinp(ii).le.Xout(i,j)).and.
     & (Xinp(ii+1).gt.Xout(i,j))) then
               i1=ii
               i2=ii+1
               goto 10
             endif
           enddo
           print*, 'Did not find i1 and i2'
           goto 100
10 continue
           DO jj=LBy,(UBy-1)
             if ((Yinp(jj).le.Yout(i,j)).and.
     & (Yinp(jj+1).gt.Yout(i,j))) then
               j1=jj
               j2=jj+1
               goto 20
             endif
           enddo
           print*, 'Did not find j1 and j2'
           goto 100
20 continue

          IF (((LBx.le.i1).and.(i1.le.UBx)).and.
     & ((LBy.le.j1).and.(j1.le.UBy))) THEN
            x1=Xinp(i1)
            x2=Xinp(i2)
            y1=Yinp(j1)
            y2=Yinp(j2)
            x=Xout(i,j)
            y=Yout(i,j)

            cff= Finp(i1,j1)*(x2-x )*(y2-y )
     & +Finp(i2,j1)*(x -x1)*(y2-y )
     & +Finp(i1,j2)*(x2-x )*(y -y1)
     & +Finp(i2,j2)*(x -x1)*(y -y1)

            Fout(i,j)=cff/((x2-x1)*(y2-y1))
          END IF
        END DO
      END DO
      RETURN
100 continue
      print*, 'error in linterp2d'
      END SUBROUTINE linterp2d

! ******************************************************************************
      SUBROUTINE cinterp2d (Lm,Mm,
     & LBx, UBx, LBy, UBy,
     & Xinp, Yinp, Finp,
     & Istr, Iend, Jstr, Jend,
     & Xout, Yout, Fout)
!
!=======================================================================
! !
! Given any gridded 2D field, Finp, at locations (Xinp,Yinp) this !
! routine performs bicubic interpolation at locations (Xout,Yout). !
! To facilitate the interpolation within any irregularly gridded !
! field, the fractional grid cell indices (Iout,Jout) with respect !
! Finp are needed at input. Notice that the routine "hindices" can !
! be used to compute these indices. !
! !
! On Input: !
! !
! LBx I-dimension lower bound of gridded field, Finp. !
! UBx I-dimension upper bound of gridded field, Finp. !
! LBy J-dimension lower bound of gridded field, Finp. !
! UBy J-dimension upper bound of gridded field, Finp. !
! Xinp X-locations of gridded field, Finp. !
! Yinp Y-locations of gridded field, Finp. !
! Finp 2D field to interpolate from. !
! Istr Starting data I-index to interpolate, Fout. !
! Iend Ending data I-index to interpolate, Fout. !
! Jstr Starting data J-index to interpolate, Fout. !
! Jend Ending data J-index to interpolate, Fout. !
! Xout X-locations to interpolate, Fout. !
! Yout Y-locations to interpolate, Fout. !
! !
! On Output: !
! !
! Fout Interpolated 2D field. !
! !
!=======================================================================
!
      implicit none

      integer Lm,Mm
!
! Imported variable declarations.
!
      integer, intent(in) :: LBx, UBx, LBy, UBy
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
      real(kind=8), intent(in) :: Xinp(LBx:UBx)
      real(kind=8), intent(in) :: Yinp(LBy:UBy)
      real(kind=8), intent(in) :: Finp(LBx:UBx,LBy:UBy)
!
      real(kind=8), intent(in) :: Xout(0:Lm+1,0:Mm+1)
      real(kind=8), intent(in) :: Yout(0:Lm+1,0:Mm+1)
!
      real(kind=8), intent(out) :: Fout(0:Lm+1,0:Mm+1)
!
! Local variable declarations.
!
      integer i, ic, iter, i1, i2, j, jc, j1, j2, ii, jj

      real(kind=8) :: a11, a12, a21, a22
      real(kind=8) :: e11, e12, e21, e22
      real(kind=8) :: cff, d1, d2, dfc, dx, dy, eta, xi, xy, yx
      real(kind=8) :: f0, fx, fxx, fxxx, fxxy, fxy, fxyy, fy, fyy, fyyy

      real(kind=8), parameter :: C01 = 1.0/48.0
      real(kind=8), parameter :: C02 = 1.0/32.0
      real(kind=8), parameter :: C03 = 0.0625 ! 1/16
      real(kind=8), parameter :: C04 = 1.0/6.0
      real(kind=8), parameter :: C05 = 0.25
      real(kind=8), parameter :: C06 = 0.5
      real(kind=8), parameter :: C07 = 0.3125 ! 5/16
      real(kind=8), parameter :: C08 = 0.625 ! 5/8
      real(kind=8), parameter :: C09 = 1.5
      real(kind=8), parameter :: C10 = 13.0/24.0

      real(kind=8), parameter :: LIMTR = 3.0
      real(kind=8), parameter :: spv = 0.0 ! HGA need work

      real(kind=8), dimension(-1:2,-1:2) :: dfx, dfy, ff
!
Cf2py intent(in) Lm,Mm,LBx, UBx, LBy, UBy,Xinp, Yinp, Finp,Istr, Iend, Jstr, Jend,Xout, Yout
Cf2py intent(out) Fout
!-----------------------------------------------------------------------
! Interpolates requested field locations (Xout,Yout).
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
! i1=INT(Iout(i,j))
! i2=i1+1
! j1=INT(Jout(i,j))
! j2=j1+1
           DO ii=LBx,(UBx-1)
             if ((Xinp(ii).le.Xout(i,j)).and.
     & (Xinp(ii+1).gt.Xout(i,j))) then
               i1=ii
               i2=ii+1
               goto 10
             endif
           enddo
           print*, 'Did not find i1 and i2',
     & Istr,Iend,Jstr,Jend,i,j,Xout(i,j),Xout(i-1,j)
           goto 100
10 continue
           DO jj=LBy,UBy-1
             if ((Yinp(jj).le.Yout(i,j)).and.
     & (Yinp(jj+1).gt.Yout(i,j))) then
               j1=jj
               j2=jj+1
               goto 20
             endif
           enddo
           print*, 'Did not find j1 and j2'
           goto 100
20 continue

          IF (((LBx.le.i1).and.(i1.le.UBx)).and.
     & ((LBy.le.j1).and.(j1.le.UBy))) THEN
!
! Determine local fractional coordinates (xi,eta) corresponding to
! the target point (Xout,Yout) on the grid (Xinp,Yinp). Here, "xi"
! and "eta" are defined, in such a way, that xi=eta=0 corresponds
! to the middle of the cell (i1:i1+1,j1:j1+1), while xi=+/-1/2 and
! eta=+/-1/2 (any combination +/- signs) corresponds to the four
! corner points of the cell. Inside the cell it is assumed that
! (Xout,Yout) are expressed via bi-linear functions of (xi,eta),
! where term proportional to xi*eta does not vanish because
! coordinate transformation may be at least weakly non-orthogonal
! due to discretization errors. The associated non-linear system
! is solved by iterative method of Newton.
!
            xy=Xinp(i2)-Xinp(i1)-Xinp(i2)+Xinp(i1)
            yx=Yinp(j2)-Yinp(j2)-Yinp(j1)+Yinp(j1)
            dx=Xout(i,j)-0.25*(Xinp(i2)+Xinp(i1)+
     & Xinp(i2)+Xinp(i1))
            dy=Yout(i,j)-0.25*(Yinp(j2)+Yinp(j2)+
     & Yinp(j1)+Yinp(j1))
!
! The coordinate transformation matrix:
!
! e11 e12
! e21 e22
!
! contains derivatives of (Xinp,Yinp) with respect to (xi,eta). Because
! the coordinates may be non-orthogonal (at least due to discretization
! errors), the nonlinear system
!
! e11*xi+e12*eta+xy*xi*eta=dx
! e21*xi+e22*eta+yx*xi*eta=dy
!
! needs to be solved in order to retain symmetry.
!
            e11=0.5*(Xinp(i2)-Xinp(i1)+Xinp(i2)-Xinp(i1))
            e12=0.5*(Xinp(i2)+Xinp(i1)-Xinp(i2)-Xinp(i1))
            e21=0.5*(Yinp(j2)-Yinp(j2)+Yinp(j1)-Yinp(j1))
            e22=0.5*(Yinp(j2)+Yinp(j2)-Yinp(j1)-Yinp(j1))
!
            cff=1.0/(e11*e22-e12*e21)
            xi=cff*(e22*dx-e12*dy)
            eta=cff*(e11*dy-e21*dx)
!
            DO iter=1,4
              d1=dx-e11*xi-e12*eta-xy*xi*eta
              d2=dy-e21*xi-e22*eta-yx*xi*eta
              a11=e11+xy*eta
              a12=e12+xy*xi
              a21=e21+yx*eta
              a22=e22+yx*xi
              cff=1.0/(a11*a22-a12*a21)
              xi =xi +cff*(a22*d1-a12*d2)
              eta=eta+cff*(a11*d2-a21*d1)
            END DO


!
! Genuinely two-dimensional, isotropic cubic interpolation scheme
! using 12-point stencil. In the code below the interpolated field,
! Fout, is expanded into two-dimensional Taylor series of local
! fractional coordinates "xi" and "eta", retaining all terms of
! combined power up to third order (that is, xi, eta, xi^2, eta^2,
! xi*eta, xi^3, eta^3, xi^2*eta, and xi*eta^2), with all
! coefficients (i.e, derivatives) computed via x x
! two-dimensional finite difference expressions | |
! of "natural" order of accuracy: 4th-order for x--x--x--x
! the field itself and its first derivatives in | |
! both directions; and 2nd-order for all higher- x--x--x--x
! order derivatives. The permissible range of | |
! of coordinates is -1/2 < xi,eta < +1/2, which x--x
! covers the central cell on the stencil, while
! xi=eta=0 corresponds to its center. This interpolation scheme has
! the property that if xi,eta=+/-1/2 (any combination of +/- signs)
! it reproduces exactly value of the function at the corresponding
! corner of the central "working" cell. However, it does not pass
! exactly through the extreme points of the stencil, where either
! xi=+/-3/2 or eta+/-3/2. And, unlike a split-directional scheme,
! when interpolating along the line eta=+/-1/2 (similarly xi=+/-1/2),
! it has non-zero contribution from points on the side from the line,
! except if xi=-1/2; 0; +1/2 (similarly eta=-1/2; 0; +1/2).
!
            DO jc=-1,2
              DO ic=-1,2
                ff(ic,jc)=Finp(MAX(1,MIN(UBx,i1+ic)),
     & MAX(1,MIN(UBy,j1+jc)))
              END DO
            END DO

            f0=C07*(ff(1,1)+ff(1,0)+ff(0,1)+ff(0,0))-
     & C02*(ff(2,0)+ff(2,1)+ff(1,2)+ff(0,2)+
     & ff(-1,1)+ff(-1,0)+ff(0,-1)+ff(1,-1))

            fx=C08*(ff(1,1)+ff(1,0)-ff(0,1)-ff(0,0))-
     & C01*(ff(2,1)+ff(2,0)-ff(-1,1)-ff(-1,0))-
     & C03*(ff(1,2)-ff(0,2)+ff(1,-1)-ff(0,-1))

            fy=C08*(ff(1,1)-ff(1,0)+ff(0,1)-ff(0,0))-
     & C01*(ff(1,2)+ff(0,2)-ff(1,-1)-ff(0,-1))-
     & C03*(ff(2,1)-ff(2,0)+ff(-1,1)-ff(-1,0))

            fxy=ff(1,1)-ff(1,0)-ff(0,1)+ff(0,0)

            fxx=C05*(ff(2,1)-ff(1,1)-ff(0,1)+ff(-1,1)+
     & ff(2,0)-ff(1,0)-ff(0,0)+ff(-1,0))

            fyy=C05*(ff(1,2)-ff(1,1)-ff(1,0)+ff(1,-1)+
     & ff(0,2)-ff(0,1)-ff(0,0)+ff(0,-1))

            fxxx=C06*(ff(2,1)+ff(2,0)-ff(-1,1)-ff(-1,0))-
     & C09*(ff(1,1)+ff(1,0)-ff(0,1)-ff(0,0))

            fyyy=C06*(ff(1,2)+ff(0,2)-ff(1,-1)-ff(0,-1))-
     & C09*(ff(1,1)-ff(1,0)+ff(0,1)-ff(0,0))

            fxxy=C06*(ff(2,1)-ff(1,1)-ff(0,1)+ff(-1,1)-
     & ff(2,0)+ff(1,0)+ff(0,0)-ff(-1,0))

            fxyy=C06*(ff(1,2)-ff(1,1)-ff(1,0)+ff(1,-1)-
     & ff(0,2)+ff(0,1)+ff(0,0)-ff(0,-1))
# 568 "R_tools_fort_routines_gigatl/online_interp.F"
            Fout(i,j)=f0+
     & fx*xi+
     & fy*eta+
     & C06*fxx*xi*xi+
     & fxy*xi*eta+
     & C06*fyy*eta*eta+
     & C04*fxxx*xi*xi*xi+
     & C06*fxxy*xi*xi*eta+
     & C04*fyyy*eta*eta*eta+
     & C06*fxyy*xi*eta*eta
          END IF
        END DO
      END DO

      RETURN

100 continue
      print*, 'error in cinterp2d'
      END SUBROUTINE cinterp2d
# 17 "R_tools_fort_gigatl.F" 2
