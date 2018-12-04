MODULE domini
  !!======================================================================
  !!                         *** MODULE domini ***
  !! Domain Initialization : Initialization of array for tmx system 
  !!=====================================================================

  !!---------------------------------------------------------------------
  !!  dom_init   :  initialization of field variables and array   
  !!---------------------------------------------------------------------
  USE par_tide   !  TMX parameters and variables  
  USE phycst     !  Physical constant used by TMX system (Routine phys_cst) 
  USE in_manager !  Input manager for TMX system (Routine name_input)
  USE grn_tide   !  library for green function   (Routine init_gren)
  USE lib_otps   !  library for barotropic tide  (routine init_tide)
  USE tmx_array  !  array manager for TMX system (Routine getdimension, tmx_alloc) 
  USE arayini    !  Initialization of array used by TMX (Routine init_hydro, init_metric &
                 !  & init_deriv, init_fcor, init_uv, comp_cutoff, comp_crt ) 

  IMPLICIT NONE
  PRIVATE
  
  PUBLIC ::  dom_init   !call by tmx.f90

CONTAINS

  SUBROUTINE dom_init
    !!---------------------------------------------------------------------------
    !!                   *** ROUTINE dom_init ***
    !!
    !! ** Purpose :    initialization of dominag by calling different routines
    !!---------------------------------------------------------------------------
    !    
    CALL phys_cst
    !
    CALL name_input
    !
    CALL getdimension 
    !
    CALL tmx_alloc(kx = nx, ky = ny)
    !
    CALL init_hydro
    !
    CALL init_metric
    !
    CALL init_deriv
    !
    CALL init_fcor
    !
    IF ( .NOT. ln_modal ) CALL init_gren
    !
    CALL init_tide 
    !
    CALL init_uv
    !
    CALL comp_cutoff
    !
    CALL comp_crt
    !
  END SUBROUTINE dom_init

  !!====================================================================================
END MODULE domini
