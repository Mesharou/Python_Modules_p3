 
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

# include "R_tools_fort_routines/cppdefs.h"


# include "R_tools_fort_routines_gigatl/online_interp.F"

# include "R_tools_fort_routines_gigatl/bulk_stress.F"









