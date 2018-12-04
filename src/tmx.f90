MODULE tmx
  !!=================================================================================
  !!                           *** MODULE tmx ***
  !! Tidal mixing model : TMX ocean tidal mixing model
  !!=================================================================================   

  !!--------------------------------------------------------------------------------
  !! tmx_model     : compute energy conversion and baroclinic internal tide velocity
  !! tmx_init      : initialization of tmx model
  !! tmx_closefile : close remaing open files
  !!--------------------------------------------------------------------------------
  USE par_tide     ! parameters and flags
  USE in_manager   ! inuput manager (getunit, getarea and opn_namelsit routines)
  USE tmx_array    ! array manager (tmx_dealloc routine) 
  USE domini       ! domain initialization (dom_init routine)   
  USE soltmx       ! tmx solvers (sol_par, sol_energy and sol_velocity routines) 

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: tmx_model        !called by model.f90

CONTAINS 

  SUBROUTINE tmx_model
    !!-----------------------------------------------------------------------------
    !!                   ***  ROUTINE tmx_model ***
    !!
    !! ** Purpose :  calculate energy conversion and baroclinic tidal velocity 
    !!               in the ocean using regular grid on sphere 
    !! 
    !! ** Method  : - initialize the model 
    !!              - using enery or velocity solver (sol_energy and
    !!                sol_velocity routine)
    !!              - finalizing the computation by closing the files
    !! 
    !!------------------------------------------------------------------------------
    !                                    !---------------------!
    !                                    !== Initializations ==! 
    !                                    !---------------------!
    CALL tmx_init               

    IF(lwp) WRITE(numout,*)   
    IF(lwp) WRITE(numout,*) 'TM_TM_TM'   ! control flag  
    IF(lwp) CALL FLUSH(numout)
    !                                    !-----------------------!
    !                                    !== computational part  !
    !                                    !-----------------------!
    IF( ln_energy) THEN
       CALL sol_energy                   ! compute the energy conversion
       PRINT *, 'energy computation was done for area', TRIM(carea), '   ', narea
    END IF
    !
    IF( ln_velocity) THEN
       CALL sol_velocity                 ! compute the baroclinic tidal velocity  
       PRINT *, 'velocity computation was done for area','    ', TRIM(carea), '   ', narea
    END IF
    !                                   !--------------------------------!
    !                                   !== Finalizing the computation ==!
    !                                   !--------------------------------!

    PRINT *, 'TM_TM_TM'                 ! control flag              

    IF(lwp) WRITE(numout,*)   
    IF(lwp) WRITE(numout,*) 'TM_TM_TM'  
    IF(lwp) CALL FLUSH(numout)
    !                                   !--------------------------------!
    !                                   !==   Freeing up the memory    ==!
    !                                   !--------------------------------!
    CALL tmx_dealloc
    !                                   !---------------------------------!
    !                                   !==   Closing the files         ==! 
    !                                   !---------------------------------!
    CALL tmx_closefile
    !
  END SUBROUTINE tmx_model


  SUBROUTINE tmx_init
    !!------------------------------------------------------------------------------
    !!                      *** ROUTINE tmx_init ***
    !!
    !! ** Purpose :   initialization of tmx model
    !!
    !! ** Method  :  -open the namelist and find which variable sould be computed
    !!               -get the information of area of computation
    !!               -call the dom_init routine to initialize the domain
    !!
    !!--------------------------------------------------------------------------------------
    INTEGER :: istat     !temporary integer to check if open status of tide.output is ok
    !!--------------------------------------------------------------------------------------
    !
    CALL opn_namelist    !open the namelist
    !
    CALL sol_par         !reading solver parameters used by soltmx.f90 
    !
    CALL getarea         !get the number and name of area 
    !
    lwp = narea == 1     !flag for printing, only true for processor one
    !
    !
    istat = 0
    !
    !   
    !
    IF (lwp) THEN
       !
       CALL getunit(numout)
       !
       OPEN(UNIT=numout, FILE='tide.output', STATUS='REPLACE', FORM="FORMATTED", ACCESS="SEQUENTIAL", iostat=istat)
       !
       IF ( istat == 0 ) THEN
          WRITE(numout,*)'        Stockholm University           '
          WRITE(numout,*)'     Department of Meteorology         '
          WRITE(numout,*)'               MISU                    '
          WRITE(numout,*)'        Tidal Mixing Model             '
          WRITE(numout,*)'       Version 1.0 (2011)              '
          WRITE(numout,*)'            Authors:                   '
          WRITE(numout,*)'   Saeed Falahat     saeed@misu.su.se  ' 
          WRITE(numout,*)'   Jonas Nycander    jonas@misu.su.se  ' 
          WRITE(numout,*)'   Maxime Ballarotta maxime@misu.su.se '
          WRITE(numout,*)'                                    '
          WRITE(numout,*)   
          WRITE(numout,*) '  solver parameters '
          WRITE(numout,*) '  ~~~~~~~~~~~~~~~~~'
          !
          IF( ln_modal ) THEN
             WRITE(numout,*) ' modal computation will be done for mode = ', nn_mode
             WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~'
             WRITE(numout,*)
          END IF
          !
          IF( ln_energy ) THEN
             WRITE(numout,*) ' energy conversion and wave-drag tensor will be computed'
             WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~'
             WRITE(numout,*)
          END IF
          !
          IF( ln_velocity ) THEN
             IF ( ln_bt ) THEN
                WRITE(numout,*) ' barotropic velocity will be computed'
             ELSE
                WRITE(numout,*) ' baroclinic velocity will be computed'

             END IF
             WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~'
             WRITE(numout,*)
          END IF
          !
          WRITE(numout,'(a,i4)') 'the number of grid points for buffer zone = ', nn_bfz
          IF ( .NOT. ln_modal) THEN
             WRITE(numout,'(a,i2)') 'the method of correction for computation  = ', nn_corr
          END IF
          !
          WRITE(numout,*) 
          WRITE(numout,*) 'TM_TM_TM'
          !
       ELSE 
          PRINT *, 'error in opening tide.output for processor one'
          STOP      
       END IF
       !
    END IF
    !
    CALL dom_init     ! initialize the domain
    !
    !
  END SUBROUTINE tmx_init


  SUBROUTINE tmx_closefile
    !!-----------------------------------------------------------------------------------
    !!                ***  ROUTINE tmx_closefile  ***
    !!
    !!  ** Purpose :   Close the files if they are open                                                                        
    !!-----------------------------------------------------------------------------------                                           
    INTEGER       :: ji           !dummy integer variable for loop                                                                   
    INTEGER       :: inumq        !maximum number of quries for inquire command                                                      
    LOGICAL, SAVE :: llopn, llok  !logical variable for checking existenss and openness                                              
    !!------------------------------------------------------------------------------------                                           
    !
    inumq = 1000       ! maximum number for quering unit
    !
    DO ji = 15, inumq
       INQUIRE( unit = ji, exist = llok, opened = llopn )
       IF( llok .AND. llopn ) THEN
          CLOSE(ji)                ! close the unit that is open
       END IF
    END DO
    !                                                                            
  END SUBROUTINE tmx_closefile

  !!======================================================================================
END MODULE tmx
