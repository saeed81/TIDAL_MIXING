MODULE wrtout
  !!=====================================================================================
  !!                           *** MODULE wrtout ***                     
  !! Output manager :    Writing ouput for tmx system 
  !!=====================================================================================

  !!-------------------------------------------------------------------------------------
  !!  name_outeng        :  output name for energy computation
  !!  output_energy      :  output file for energy computation 
  !!  name_outvel        :  output name for velocity computation
  !!  output_velocity    :  output file for velocity computation 
  !!-------------------------------------------------------------------------------------
  USE par_tide           !  tmx parameters
  USE arraydef           !  definition of array 
  USE io_ezcdf           !  netcdf library 

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=20)  :: cn_varlon  = '' 
  CHARACTER(len=20)  :: cn_varlat  = ''
  CHARACTER(len=20)  :: cn_vartime = ''
  CHARACTER(len=20)  :: cn_varout  = ''
  CHARACTER(len=20)  :: cn_unit    = ''
  CHARACTER(len=300) :: cn_ln      = ''

  PUBLIC :: name_outeng          !called by soltmx.f90
  PUBLIC :: output_energy        !called by soltmx.f90
  PUBLIC :: name_outvel          !called by soltmx.f90
  PUBLIC :: output_bt            !called by soltmx.f90
  PUBLIC :: output_velocity      !called by soltmx.f90
  PUBLIC :: output_cutoff        !called by soltmx.f90
  PUBLIC :: output_steep         !called by soltmx.f90
  PUBLIC :: output_hslope        !called by soltmx.f90

CONTAINS

  SUBROUTINE name_outeng
    !!-----------------------------------------------------------------------------------
    !!                    *** ROUTINE name_outeng ***
    !!
    !! ** Purpose :  Read the output names for energy computation form namelist
    !!-----------------------------------------------------------------------------------
    NAMELIST/namouteng/ cn_dir, cn_varlon, cn_varlat, cn_vartime, cn_varout, &
         & cn_unit, cn_ln, ln_cutoff, ln_steep, ln_hslope, ln_bt
    !!-----------------------------------------------------------------------------------
    !
    cn_dir = './'

    REWIND(numnam)
    READ(numnam,namouteng)

    IF( LEN_TRIM(cn_dir) /= INDEX(cn_dir,"/", back = .TRUE.)) THEN
       cn_dir = TRIM(cn_dir)//"/"   ! check whether name of dir out ends with / or not                                                                                                     
    END IF
    !
  END SUBROUTINE name_outeng


  SUBROUTINE output_energy
    !!---------------------------------------------------------------------------------------
    !!                    *** ROUTINE output_energy ***
    !!
    !! ** Purpose :  Write the output of energy computation 
    !!---------------------------------------------------------------------------------------
    CHARACTER(len=300)       :: clfile, clfile1, clfile2, clfile3
    CHARACTER(len=300)       :: cl_ln
    CHARACTER(len=30)        :: cl_varout, cl_unit
    CHARACTER(len=3)         :: clmode
    INTEGER                  :: idx_f, idx_v, idx_f1, idx_v1, idx_f2, idx_v2, idx_f3, idx_v3 
    !!---------------------------------------------------------------------------------------
    !
    !
    IF( TRIM(carea) == 'glob_reg' ) THEN  

       IF ( ln_modal  ) THEN

          WRITE(clmode,'(I2.2)') nn_mode

          clfile = TRIM(cn_dir)//'energy'//TRIM(clmode)//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'
       ELSE
          
          clfile = TRIM(cn_dir)//'energy'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'

          clfile1 = TRIM(cn_dir)//'ftxx'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'
          
          clfile2 = TRIM(cn_dir)//'ftyy'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'

          clfile3 = TRIM(cn_dir)//'ftxy'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'


       END IF

    ELSE

       IF ( ln_modal  ) THEN

          WRITE(clmode,'(I2.2)') nn_mode

          clfile = TRIM(cn_dir)//'energy'//TRIM(clmode)//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'

       ELSE

          clfile = TRIM(cn_dir)//'energy'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'

          clfile1 = TRIM(cn_dir)//'ftxx'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'

          clfile2 = TRIM(cn_dir)//'ftxy'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'

          clfile3 = TRIM(cn_dir)//'ftyy'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'


       END IF

    END IF

    CALL P2D_T(idx_f, idx_v, nx, ny, 1, 1, REAL(vlons, kind=sp), REAL(vlats, kind=sp), &
         &   REAL(vtime,kind=sp), REAL(encon, kind=sp), TRIM(clfile), &
         &   TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cn_varout), &
         &   TRIM(cn_unit), TRIM(cn_ln), 0.0)
    !                                                                         
  
    CALL P2D_T(idx_f1, idx_v1, nx, ny, 1, 1, REAL(vlons, kind=sp), REAL(vlats, kind=sp), &
         &   REAL(vtime,kind=sp), REAL(ftxx, kind=sp), TRIM(clfile1), &
         &   TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cn_varout), &
         &   TRIM(cn_unit), TRIM(cn_ln), 0.0)
    !                                                                         
    CALL P2D_T(idx_f2, idx_v2, nx, ny, 1, 1, REAL(vlons, kind=sp), REAL(vlats, kind=sp), &
         &   REAL(vtime,kind=sp), REAL(ftxy, kind=sp), TRIM(clfile2), &
         &   TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cn_varout), &
         &   TRIM(cn_unit), TRIM(cn_ln), 0.0)
    !                                                                         
    CALL P2D_T(idx_f3, idx_v3, nx, ny, 1, 1, REAL(vlons, kind=sp), REAL(vlats, kind=sp), &
         &   REAL(vtime,kind=sp), REAL(ftyy, kind=sp), TRIM(clfile3), &
         &   TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cn_varout), &
         &   TRIM(cn_unit), TRIM(cn_ln), 0.0)
    !                                                                         
  
  END SUBROUTINE output_energy


  SUBROUTINE output_cutoff
    !!---------------------------------------------------------------------------------------
    !!                    *** ROUTINE output_energy ***
    !!
    !! ** Purpose :  Write the output of energy computation 
    !!---------------------------------------------------------------------------------------
    CHARACTER(len=300)       :: clcutoff
    CHARACTER(len=300)       :: cl_ln
    CHARACTER(len=30)        :: cl_varout, cl_unit
    CHARACTER(len=3)         :: clmode
    INTEGER                  :: idx_f, idx_v
    !!---------------------------------------------------------------------------------------
    !
    !
    IF( TRIM(carea) == 'glob_reg' ) THEN  

       IF ( ln_modal  ) THEN

          WRITE(clmode,'(I2.2)') nn_mode

          clcutoff = TRIM(cn_dir)//'cutoff'//TRIM(clmode)//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'
       
       ELSE
         
          clcutoff = TRIM(cn_dir)//'cutoff'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'

       END IF

    ELSE

       IF ( ln_modal  ) THEN

          WRITE(clmode,'(I2.2)') nn_mode

          clcutoff = TRIM(cn_dir)//'cutoff'//TRIM(clmode)//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'

       ELSE

          clcutoff = TRIM(cn_dir)//'cutoff'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'
          
       END IF

    END IF

    
       cl_varout = 'cutoff'
       cl_ln = "cutoff length for component"//TRIM(cn_cu)
       cl_unit = 'm'

       CALL P2D_T(idx_f, idx_v, nx, ny, 1, 1, REAL(vlon, kind=sp), REAL(vlat, kind=sp), &
            &   REAL(vtime,kind=sp), REAL(vcut, kind=sp), TRIM(clcutoff), &
            &   TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cl_varout), &
            &   TRIM(cl_unit), TRIM(cl_ln), 0.0)
    
    !                                                                         
  END SUBROUTINE output_cutoff


  SUBROUTINE output_steep
    !!---------------------------------------------------------------------------------------
    !!                    *** ROUTINE output_energy ***
    !!
    !! ** Purpose :  Write the output of energy computation 
    !!---------------------------------------------------------------------------------------
    CHARACTER(len=300)       :: cl_ln, clsteep
    CHARACTER(len=30)        :: cl_varout, cl_unit
    CHARACTER(len=3)         :: clmode
    INTEGER                  :: idx_f, idx_v
    !!---------------------------------------------------------------------------------------
    !
    !
    IF( TRIM(carea) == 'glob_reg' ) THEN  

       IF ( ln_modal  ) THEN

          WRITE(clmode,'(I2.2)') nn_mode

          clsteep = TRIM(cn_dir)//'steepness'//TRIM(clmode)//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'
       
       ELSE
          
          clsteep = TRIM(cn_dir)//'steepness'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'
          
       END IF

    ELSE

       IF ( ln_modal  ) THEN

          WRITE(clmode,'(I2.2)') nn_mode

          clsteep = TRIM(cn_dir)//'steepness'//TRIM(clmode)//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'

       ELSE

          clsteep = TRIM(cn_dir)//'steepness'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'

          
       END IF

    END IF

    
       cl_varout = 'steepness'
       cl_ln = "steepness parameter for component"//TRIM(cn_cu)
       cl_unit = ''

       CALL P2D_T(idx_f, idx_v, nx, ny, 1, 1, REAL(vlon, kind=sp), REAL(vlat, kind=sp), &
            &   REAL(vtime,kind=sp), REAL(vsteep, kind=sp), TRIM(clsteep), &
            &   TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cl_varout), &
            &   TRIM(cl_unit), TRIM(cl_ln), 0.0)

       !                                                                         
  END SUBROUTINE output_steep

  SUBROUTINE output_hslope
    !!---------------------------------------------------------------------------------------
    !!                    *** ROUTINE output_energy ***
    !!
    !! ** Purpose :  Write the output of energy computation 
    !!---------------------------------------------------------------------------------------
    CHARACTER(len=300)       :: cl_ln, clhslope
    CHARACTER(len=30)        :: cl_varout, cl_unit
    CHARACTER(len=3)         :: clmode
    INTEGER                  :: idx_f, idx_v
    !!---------------------------------------------------------------------------------------
    !
    !
    IF( TRIM(carea) == 'glob_reg' ) THEN  

       IF ( ln_modal  ) THEN

          WRITE(clmode,'(I2.2)') nn_mode

          clhslope = TRIM(cn_dir)//'bathyslope'//TRIM(clmode)//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'
       ELSE
          
          clhslope = TRIM(cn_dir)//'bathyslope'//'_'//TRIM(carea)//'.nc'

       END IF

    ELSE

       IF ( ln_modal  ) THEN

          WRITE(clmode,'(I2.2)') nn_mode

          clhslope = TRIM(cn_dir)//'bathyslope'//TRIM(clmode)//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'

       ELSE

          clhslope = TRIM(cn_dir)//'bathyslope'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'.nc'

       END IF

    END IF

    
       cl_varout = 'bathyslope'
       cl_ln = "sqaured slope of the bathymetry" 
       cl_unit = ''

       CALL P2D_T(idx_f, idx_v, nx, ny, 1, 1, REAL(vlon, kind=sp), REAL(vlat, kind=sp), &
            &   REAL(vtime,kind=sp), REAL(vhslope, kind=sp), TRIM(clhslope), &
            &   TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cl_varout), &
            &   TRIM(cl_unit), TRIM(cl_ln), 0.0)
       !                                                                         
  END SUBROUTINE output_hslope


  SUBROUTINE name_outvel
    !!---------------------------------------------------------------------------------------
    !!                         *** ROUTINE name_outvel ***
    !!
    !! ** Purpose  :  Read the output names for velocity computation form namelist
    !!---------------------------------------------------------------------------------------
    NAMELIST/namoutvel/ cn_dir, cn_varlon, cn_varlat, cn_vartime, cn_varout, &
         & cn_unit, cn_ln, ln_bt
    !!---------------------------------------------------------------------------------------                                                                                                
    !
    cn_dir = './'

    REWIND(numnam)
    READ(numnam,namoutvel)
    !
    IF( LEN_TRIM(cn_dir) /= INDEX(cn_dir,"/", back = .TRUE.)) THEN
       cn_dir = TRIM(cn_dir)//"/"   ! check whether name of dir out ends with / or not                                                                                                     
    END IF
    !
  END SUBROUTINE name_outvel


  SUBROUTINE output_bt
    !!--------------------------------------------------------------------------------------
    !!                       *** ROUTINE output_velocity *** 
    !!
    !! ** Purpose  :    Write the output of velocity computation
    !!--------------------------------------------------------------------------------------
    CHARACTER(len=300)      :: clfile, cl_varout, cl_ln, cl_unit
    INTEGER                 :: idx_f1, idx_v1
    INTEGER                 :: idx_f2, idx_v2
    INTEGER                 :: idx_f3, idx_v3
    INTEGER                 :: idx_f4, idx_v4
    !!--------------------------------------------------------------------------------------
    !

    cl_varout = 'velocity' 
    cl_ln = 'computed velocity on topographic points for'//TRIM(cn_cu)
    cl_unit = 'm/s'

    IF( TRIM(carea) == 'glob_reg' ) THEN  
       clfile = TRIM(cn_dir)//'uimag'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'_bt'//'.nc'
    ELSE
       clfile = TRIM(cn_dir)//'uimag'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'_bt'//'.nc'
    END IF

    CALL P2D_T(idx_f1, idx_v1, nx, ny, 1, 1, REAL(vlon, kind=sp), REAL(vlat, kind=sp),&
         &   REAL(vtime, kind=sp), REAL(uimag, kind=sp), TRIM(clfile), &
         &   TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cl_varout), TRIM(cl_unit), &
         &   TRIM(cl_ln), 0.0)
    !                                                                                                                          
    IF( TRIM(carea) == 'glob_reg' ) THEN  
       clfile = TRIM(cn_dir)//'ureal'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'_bt'//'.nc'
    ELSE
       clfile = TRIM(cn_dir)//'ureal'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'_bt'//'.nc'
    END IF
    !
    CALL P2D_T(idx_f2, idx_v2, nx, ny, 1, 1, REAL(vlon, kind=sp), REAL(vlat, kind=sp),&
         &  REAL(vtime, kind=sp), REAL(ureal, kind=sp), TRIM(clfile), &
         &  TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cl_varout), TRIM(cl_unit), &
         &  TRIM(cl_ln), 0.0)
    !                                                                                                                          
    IF( TRIM(carea) == 'glob_reg' ) THEN  
       clfile = TRIM(cn_dir)//'vimag'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'_bt'//'.nc'
    ELSE
       clfile = TRIM(cn_dir)//'vimag'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'_bt'//'.nc'
    END IF
    !
    CALL P2D_T(idx_f3, idx_v3, nx, ny, 1, 1, REAL(vlon,kind=sp), REAL(vlat,kind=sp),&
         &   REAL(vtime, kind=sp), REAL(vimag, kind=sp), TRIM(clfile), &
         &   TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cl_varout), TRIM(cl_unit), &
         &   TRIM(cl_ln), 0.0)
    !

    IF( TRIM(carea) == 'glob_reg' ) THEN  
       clfile = TRIM(cn_dir)//'vreal'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'_bt'//'.nc'
    ELSE
       clfile = TRIM(cn_dir)//'vreal'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'_bt'//'.nc'
    END IF
    !
    CALL P2D_T(idx_f4, idx_v4, nx, ny, 1, 1, REAL(vlon, kind=sp), REAL(vlat, kind=sp),&
         &  REAL(vtime, kind=sp), REAL(vreal, kind=sp), TRIM(clfile), &
         &  TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cl_varout), TRIM(cl_unit), &
         &  TRIM(cl_ln), 0.0)
    !   
  END SUBROUTINE output_bt


  SUBROUTINE output_velocity
    !!--------------------------------------------------------------------------------------
    !!                       *** ROUTINE output_velocity *** 
    !!
    !! ** Purpose  :    Write the output of velocity computation
    !!--------------------------------------------------------------------------------------
    CHARACTER(len=300)      :: clfile
    INTEGER                 :: idx_f1, idx_v1
    INTEGER                 :: idx_f2, idx_v2
    INTEGER                 :: idx_f3, idx_v3
    INTEGER                 :: idx_f4, idx_v4
    !!--------------------------------------------------------------------------------------
    !
    IF( TRIM(carea) == 'glob_reg' ) THEN  
       clfile = TRIM(cn_dir)//'uimag'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'
    ELSE
       clfile = TRIM(cn_dir)//'uimag'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'
    END IF

    CALL P2D_T(idx_f1, idx_v1, nx, ny, 1, 1, REAL(vlon, kind=sp), REAL(vlat, kind=sp),&
         &   REAL(vtime, kind=sp), REAL(uimag, kind=sp), TRIM(clfile), &
         &   TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cn_varout), TRIM(cn_unit), &
         &   TRIM(cn_ln), 0.0)
    !                                                                                                                          
    IF( TRIM(carea) == 'glob_reg' ) THEN  
       clfile = TRIM(cn_dir)//'ureal'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'
    ELSE
       clfile = TRIM(cn_dir)//'ureal'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'
    END IF
    !
    CALL P2D_T(idx_f2, idx_v2, nx, ny, 1, 1, REAL(vlon, kind=sp), REAL(vlat, kind=sp),&
         &  REAL(vtime, kind=sp), REAL(ureal, kind=sp), TRIM(clfile), &
         &  TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cn_varout), TRIM(cn_unit), &
         &  TRIM(cn_ln), 0.0)
    !                                                                                                                          
    IF( TRIM(carea) == 'glob_reg' ) THEN  
       clfile = TRIM(cn_dir)//'vimag'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'
    ELSE
       clfile = TRIM(cn_dir)//'vimag'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'
    END IF
    !
    CALL P2D_T(idx_f3, idx_v3, nx, ny, 1, 1, REAL(vlon,kind=sp), REAL(vlat,kind=sp),&
         &   REAL(vtime, kind=sp), REAL(vimag, kind=sp), TRIM(clfile), &
         &   TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cn_varout), TRIM(cn_unit), &
         &   TRIM(cn_ln), 0.0)
    !

    IF( TRIM(carea) == 'glob_reg' ) THEN  
       clfile = TRIM(cn_dir)//'vreal'//'_'//TRIM(carea)//'_'//TRIM(cn_cu)//'.nc'
    ELSE
       clfile = TRIM(cn_dir)//'vreal'//'_'//'r'//'_'//TRIM(cxarea)//'_'//TRIM(cyarea)//'_'//TRIM(cn_cu)//'.nc'
    END IF
    !
    CALL P2D_T(idx_f4, idx_v4, nx, ny, 1, 1, REAL(vlon, kind=sp), REAL(vlat, kind=sp),&
         &  REAL(vtime, kind=sp), REAL(vreal, kind=sp), TRIM(clfile), &
         &  TRIM(cn_varlon), TRIM(cn_varlat), TRIM(cn_vartime), TRIM(cn_varout), TRIM(cn_unit), &
         &  TRIM(cn_ln), 0.0)
    !   
  END SUBROUTINE output_velocity

  !!========================================================================================
END MODULE wrtout
