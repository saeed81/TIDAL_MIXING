MODULE phycst
  !!==============================================================================
  !!                      *** MODULE physcst ***
  !! Physical constant :  definition of physical constant for tmx 
  !!==============================================================================

  !!------------------------------------------------------------------------------
  !!  physc_cst  : physical constant used by tmx system 
  !!------------------------------------------------------------------------------
  USE par_tide   !  parameters and variable of TMX system

  IMPLICIT NONE
  PRIVATE

  REAL(wp), PUBLIC                  :: ra = 6371229.0_wp    !: Radius of Earth (m)                                                                  
  REAL(wp), PUBLIC                  :: omega = 7.292e-5_wp !: earth rotation rate  (s^-1)                                                        
  REAL(wp), PUBLIC                  :: rho0 = 1040.0_wp    !: density of water (kg/(m^3))                                                         
  REAL(wp), PUBLIC                  :: rpi                 !: pi variable
  REAL(wp), PUBLIC                  :: rad                 !: radian
  REAL(wp), PUBLIC, DIMENSION(10)   :: tfrq                !: frequecncy of tidal constituents
  
  CHARACTER(len=2), PUBLIC, DIMENSION(10) :: constit       !: tidal consituents

  PUBLIC :: phys_cst             ! called by  domini.f90 

CONTAINS

  SUBROUTINE phys_cst
    !!-----------------------------------------------------------------------------  
    !!                    *** ROUTINE phys_cst ***
    !!
    !! ** Purpose :  calculate some physical constants for tmx model
    !!-----------------------------------------------------------------------------
    !
    rpi = ATAN(1.0) * 4.0_wp
    rad = rpi / 180.0_wp
    !
    tfrq = (/ 1.405189e-04_wp, 1.454441e-04_wp, 1.378797e-04_wp, 1.458423e-04_wp, &
         &   7.292117e-05_wp, 6.759774e-05_wp, 7.252295e-05_wp, 6.495854e-05_wp, &
         &   0.053234e-04_wp, 0.026392e-04_wp /)
    !                                                                                                                          
    constit = (/'m2','s2','n2','k2','k1','o1','p1','q1','mf','mm'/)
    !
    IF(lwp) THEN
       WRITE(numout,*)
       WRITE(numout,*) 'physical constant'
       WRITE(numout,*) '~~~~~~~~~~~~~~~~~'
       WRITE(numout,*) '  the radius of earth is       ', ra 
       WRITE(numout,*) '  the earth rotation rate is   ', omega 
       WRITE(numout,*) '  the density of water is      ', rho0
       WRITE(numout,*) '  the value of pi              ', rpi
       WRITE(numout,*) '  the radian is                ', rad
       WRITE(numout,*)
    END IF
    !
  END SUBROUTINE phys_cst

  !!===============================================================================
END MODULE phycst
