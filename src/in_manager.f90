MODULE in_manager
  !!============================================================================
  !!                    *** MODULE in_manager ***
  !! Input manager : input manager of tmx system 
  !!============================================================================

  !!----------------------------------------------------------------------------
  !!  getarea      :  give the information of name of area and index of that
  !!  name_input   :  give the name of input files
  !!  opn_namelist :  open the namelist file
  !!----------------------------------------------------------------------------
  USE par_tide     !  TMX system parameters

  IMPLICIT NONE
  PRIVATE

  INTEGER       :: nn_isplt = 1
  INTEGER       :: nn_jsplt = 1
  INTEGER       :: narg     = 0
  LOGICAL       :: lmpp = .FALSE.

  PUBLIC :: getarea              ! called by tmx.f90   
  PUBLIC :: name_input           ! called by domini.f90
  PUBLIC :: opn_namelist         ! called by tmx.f90
  PUBLIC :: getunit              ! called by various routines

CONTAINS

  SUBROUTINE getarea
    !!----------------------------------------------------------------------------
    !!                       *** ROUTINE getarea ***
    !!
    !! ** Purpose :  give the information on name of area and number of area 
    !!----------------------------------------------------------------------------    
    INTEGER :: ji, jj
    INTEGER :: ixarea, iyarea
    LOGICAL :: llstop 
    !!
    NAMELIST/nammpp/ nn_isplt, nn_jsplt
    !!----------------------------------------------------------------------------
    !
    REWIND(numnam)
    READ(numnam,nammpp)

    IF( nn_isplt == 1 .AND. nn_jsplt == 1 ) THEN
       narea = 1
       carea = 'glob_reg'
       lmpp = .FALSE.
    END IF

    IF( nn_isplt /= 1 .OR. nn_jsplt /= 1 ) THEN
       narg = iargc()
       IF ( narg /= 2 ) THEN
          PRINT *, ' two arguments are required for parallel run'
          STOP 
       END IF

       lmpp = .TRUE.

       CALL getarg(1, cxarea) 
       CALL getarg(2, cyarea) 

       carea = TRIM(cxarea)//TRIM(cyarea)

       READ(cxarea, '(I4)') ixarea
       READ(cyarea, '(I4)') iyarea

       narea = 0
       llstop = .TRUE.
       jj = 1
       !
       DO WHILE ( jj <= nn_jsplt )
          ji = 1
          DO WHILE ( ji <= nn_isplt .AND. llstop )

             narea = narea + 1
             llstop = ((ji == ixarea) .AND. (jj == iyarea) )
             ji = ji + 1
             llstop = .NOT. llstop
          END DO
          jj = jj + 1
       END DO
    END IF
    !
  END SUBROUTINE getarea


  SUBROUTINE name_input
    !!------------------------------------------------------------------------
    !!                   *** ROUTINE name_input ***
    !!
    !! ** Purpose :  give the name of input files    
    !!------------------------------------------------------------------------  
    CHARACTER(len=3)      :: clsuffix
    CHARACTER(len=3)      :: clmode
    LOGICAL, SAVE         ::   llok
    !!
    NAMELIST/namhydro/ cn_dir
    !!------------------------------------------------------------------------  
    !
    cn_dir = './' 

    REWIND(numnam)
    READ(numnam,namhydro)

    clsuffix = '.nc'

    IF( LEN_TRIM(cn_dir) /= INDEX(cn_dir,"/", back = .TRUE.)) THEN
       cn_dir = TRIM(cn_dir)//"/"   ! check whether name of dir out ends with / or not                                                                                                     
    END IF

    IF ( .NOT. lmpp ) THEN

       IF ( .NOT. ln_modal ) THEN

          cfile_h     = TRIM(cn_dir)//'bathy'//clsuffix
          cvar_h      = 'bathymetry' 
          cfile_bf    = TRIM(cn_dir)//'bottom_frequency'//clsuffix
          cvar_bf     = 'bot_freq'
          cfile_intbf = TRIM(cn_dir)//'int_frequency'//clsuffix
          cvar_intbf  = 'int_freq'
          cvarlon     = 'lon'
          cvarlat     = 'lat'

          INQUIRE( FILE = TRIM(cfile_h), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_h) 
             STOP 
          END IF

          INQUIRE( FILE = TRIM(cfile_bf), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_bf) 
             STOP 
          END IF

          INQUIRE( FILE = TRIM(cfile_intbf), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_intbf) 
             STOP 
          END IF
          !
       ELSE

          WRITE(clmode,'(I2.2)') nn_mode

          cfile_h     = TRIM(cn_dir)//'bathy'//clsuffix
          cvar_h      = 'bathymetry' 
          cfile_speed    = TRIM(cn_dir)//'c'//TRIM(clmode)//clsuffix
          cvar_speed     = 'speed'
          cfile_sn    = TRIM(cn_dir)//'sn'//TRIM(clmode)//clsuffix
          cvar_sn     = 'sn'
          cfile_eps = TRIM(cn_dir)//'eps'//TRIM(clmode)//clsuffix
          cvar_eps  = 'epsilon'
          cvarlon     = 'lon'
          cvarlat     = 'lat'

          INQUIRE( FILE = TRIM(cfile_h), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_h) 
             STOP 
          END IF

          INQUIRE( FILE = TRIM(cfile_sn), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_sn) 
             STOP 
          END IF

          INQUIRE( FILE = TRIM(cfile_eps), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_eps) 
             STOP 
          END IF
          INQUIRE( FILE = TRIM(cfile_speed), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_speed) 
             STOP 
          END IF

          IF(lwp) WRITE(numout,*) 'bathymetry file is :'
          IF(lwp) WRITE(numout,*) TRIM(cfile_h) 
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) 'phase speed file is :'
          IF(lwp) WRITE(numout,*) TRIM(cfile_speed) 
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) 'sn file is :'
          IF(lwp) WRITE(numout,*) TRIM(cfile_sn) 
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) 'epsilon file is :'
          IF(lwp) WRITE(numout,*) TRIM(cfile_eps) 
          IF(lwp) WRITE(numout,*)


       END IF

    END IF

    IF ( lmpp ) THEN
       !
       IF ( .NOT. ln_modal ) THEN

          cfile_h     = TRIM(cn_dir)//'bathy'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//clsuffix
          cvar_h      = 'bathymetry' 
          cfile_bf    = TRIM(cn_dir)//'bottom_frequency'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//clsuffix
          cvar_bf     = 'bot_freq'
          cfile_intbf = TRIM(cn_dir)//'int_frequency'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//clsuffix
          cvar_intbf  = 'int_freq'
          cvarlon     = 'lon'
          cvarlat     = 'lat'

          INQUIRE( FILE = TRIM(cfile_h), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_h) 
             STOP 
          END IF

          INQUIRE( FILE = TRIM(cfile_bf), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_bf) 
             STOP 
          END IF

          INQUIRE( FILE = TRIM(cfile_intbf), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_intbf) 
             STOP 
          END IF

          IF(lwp) WRITE(numout,*) 'bathymetry file is :'
          IF(lwp) WRITE(numout,*) TRIM(cfile_h) 
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) 'bottom frequency file is :'
          IF(lwp) WRITE(numout,*) TRIM(cfile_bf) 
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) 'integral of buoyancy frequency file  is :'
          IF(lwp) WRITE(numout,*) TRIM(cfile_intbf) 
          IF(lwp) WRITE(numout,*)

          !
       ELSE

          WRITE(clmode,'(I2.2)') nn_mode

          cfile_h     = TRIM(cn_dir)//'bathy'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//clsuffix
          cvar_h      = 'bathymetry' 
          cfile_speed    = TRIM(cn_dir)//'c'//TRIM(clmode)//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//clsuffix
          cvar_speed     = 'speed'
          cfile_sn    = TRIM(cn_dir)//'sn'//TRIM(clmode)//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//clsuffix
          cvar_sn     = 'sn'
          cfile_eps = TRIM(cn_dir)//'eps'//TRIM(clmode)//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//clsuffix
          cvar_eps  = 'epsilon'
          cvarlon     = 'lon'
          cvarlat     = 'lat'

          INQUIRE( FILE = TRIM(cfile_h), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_h) 
             STOP 
          END IF

          INQUIRE( FILE = TRIM(cfile_sn), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_sn) 
             STOP 
          END IF

          INQUIRE( FILE = TRIM(cfile_eps), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_eps) 
             STOP 
          END IF
          INQUIRE( FILE = TRIM(cfile_speed), EXIST = llok )

          IF( .NOT. llok ) THEN
             PRINT *, 'the file was not found', TRIM(cfile_speed) 
             STOP 
          END IF

          IF(lwp) WRITE(numout,*) 'bathymetry file is :'
          IF(lwp) WRITE(numout,*) TRIM(cfile_h) 
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) 'phase speed file is :'
          IF(lwp) WRITE(numout,*) TRIM(cfile_speed) 
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) 'sn file is :'
          IF(lwp) WRITE(numout,*) TRIM(cfile_sn) 
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) 'epsilon file is :'
          IF(lwp) WRITE(numout,*) TRIM(cfile_eps) 
          IF(lwp) WRITE(numout,*)



       END IF

    END IF
    !
  END SUBROUTINE name_input


  SUBROUTINE opn_namelist
    !!------------------------------------------------------------------                                                                                                                     
    !!              *** ROUTINE opn_namelist ***                                 
    !!                                                                                                                                                                                       
    !! ** Purpose : Open the namelist and check on its existence                                                                                                                             
    !!                                                                                                                                                                                       
    !!------------------------------------------------------------------                                                                                                                     
    LOGICAL, SAVE         ::   llok      ! check the existence                                                                                                                               
    INTEGER               ::   istat     ! check the status of open file                                                                                                                     
    !!------------------------------------------------------------------
    istat = 0
    INQUIRE( FILE = "namelist", EXIST = llok )

    IF( .NOT. llok ) THEN
       PRINT '(/,a,/)', "the namelist file was not found"
       STOP 
    END IF

    CALL getunit(numnam)

    OPEN( unit=numnam, FILE="namelist", FORM="FORMATTED", STATUS="OLD", iostat = istat)

    IF ( istat /= 0 ) THEN
       PRINT *, 'error in opening the namelist'
       STOP
    END IF
    !                                        
  END SUBROUTINE opn_namelist


  SUBROUTINE getunit (kunit) 
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE  getunit  ***
    !!
    !! ** Purpose :   return the index of an unused logical unit
    !!----------------------------------------------------------------------
    INTEGER, INTENT(out) :: kunit
    LOGICAL              :: llopn
    !!----------------------------------------------------------------------
    !
    kunit = 15     ! choose a unit that is big enough then it is not already used
    llopn = .TRUE.
    !
    DO WHILE( (kunit <= 998) .AND. llopn )
       kunit = kunit + 1
       INQUIRE( unit = kunit, opened = llopn )
    END DO
    !
    IF( (kunit == 999 ) .AND. llopn ) THEN
       STOP 'getunit: All logical units until 998 are used...'
    ELSE
       PRINT *, 'unit ok ', kunit
    ENDIF
    !
  END SUBROUTINE getunit

  !!=====================================================================
END MODULE in_manager
