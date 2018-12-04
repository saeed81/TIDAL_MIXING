MODULE grn_tide
  !!=========================================================================
  !!                    *** MODULE grn_tide  ***
  !! green function library : green function and derivative of green function
  !!=========================================================================

  !!-------------------------------------------------------------------------
  !! init_gren    : Initialization of green function by reading the files  
  !! gren         : Function to compute green function in arbitrary point
  !! grendr       : Function to compute derivative of green function 
  !!--------------------------------------------------------------------------
  USE par_tide    ! TMX system parametes and variables
  USE phycst      ! Physical constant used by tmx system
  USE in_manager  ! In-put manager for tmx system

  IMPLICIT NONE
  PRIVATE

  INTEGER,  PARAMETER                  :: jpr = 20000  ! the size of green function and its derivative
  REAL(wp), DIMENSION(jpr)             :: gren_d       ! the array ontaing the green function
  REAL(wp), DIMENSION(jpr)             :: gren_dr      ! the array containg the derivative of green function
  REAL(wp)                             :: rdr = 0.0010_wp  ! variable use for interpolation

  PUBLIC :: init_gren            ! called by domini.f90
  PUBLIC :: gren, grendr         ! called by soltmx.f90

CONTAINS

  SUBROUTINE init_gren  
    !!-----------------------------------------------------------------------------
    !!                       *** ROUTINE init_gren ***
    !!
    !! ** Purpose :   Initialization of green function 
    !!
    !! ** Methods :  - choose the method of computing green function
    !!               - open the apropriate files that stored green functions 
    !!                 and derivative of them 
    !!-----------------------------------------------------------------------------
    INTEGER             :: ji       ! integer used for loops    
    INTEGER             :: istog    ! intger used for debugging of open files
    INTEGER             :: istodg   ! integer used for debugging of open files
    INTEGER             :: istrg    ! integer used for debugging when reading the file
    INTEGER             :: istrdg   ! integer used for debugging when reading the file
    CHARACTER(len=300)  :: clfgrn   ! local character for file name of green function
    CHARACTER(len=300)  :: clfdgrn  ! local character for file name of derivative of green function
    !!
    NAMELIST/namgrn/nn_gr, cn_dir, rn_beta
    !!------------------------------------------------------------------------------
    !


    cn_dir = './'            ! initialization of input directory
    !
    REWIND(numnam)
    READ(numnam,namgrn)      !read the green function part

    IF( LEN_TRIM(cn_dir) /= INDEX(cn_dir,"/", back = .TRUE.)) THEN
       cn_dir = TRIM(cn_dir)//"/"   ! check whether name of dir out ends with / or not
    END IF

    IF(lwp) THEN
       WRITE(numout,*)
       WRITE(numout,*) 'green function paramters'
       WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~'
    END IF

    gren_d(:)  = 0.0_wp
    gren_dr(:) = 0.0_wp

    SELECT CASE ( nn_gr )    ! read the funcions from a file and
       !
    CASE ( 1 )               ! Agnesi green function 
       !                      nn_gr = 1: filter 1/(1+k^2) (Agnesi)
       !
       IF(lwp) WRITE(numout,*)
       IF(lwp) WRITE(numout,'(a,a,i2)')   'Agnesi filter is chosen', 'nn_gr = ', nn_gr
       IF(lwp) WRITE(numout,*)   '~~~~~~'
       IF(lwp) WRITE(numout,*)

       clfgrn  = TRIM(cn_dir)//'gren_Agnesi' 
       clfdgrn = TRIM(cn_dir)//'dgren_Agnesi' 

       CALL getunit(numgr)      ! unit for green function file

       OPEN(UNIT=numgr, FILE=TRIM(clfgrn), STATUS='OLD', FORM="FORMATTED", ACCESS="SEQUENTIAL", iostat=istog)

       CALL getunit(numdgr)     ! unit for derivative of green function 

       OPEN(UNIT=numdgr, FILE=TRIM(clfdgrn), STATUS='OLD', FORM="FORMATTED", ACCESS="SEQUENTIAL", iostat=istodg)
       !
       IF ( istog /=0 .OR. istodg /= 0 ) THEN
          IF (lwp) WRITE(numout,*) ' error in opening green function files, nn_gr = ', nn_gr
          STOP 
       END IF
       !
       ji = 0
       istrg = 0
       istrdg = 0
       !
       DO ji = 1, jpr
          READ(numgr, *, iostat = istrg  ) gren_d(ji)
          READ (numdgr, *, iostat = istrdg ) gren_dr(ji)
       END DO
       !
       IF ( istrg /= 0 .OR. istrdg /= 0 ) THEN
          IF(lwp) WRITE(numout,*) 'error in reading files for green function, nn_gr =', nn_gr
          STOP
       END IF
       !
    CASE ( 2 )                  ! Gaussian Green function
       !                        nn_gr = 2: filter exp(-k^2) 
       !
       IF(lwp) WRITE(numout,'(a,a,i2)') 'Gaussian filter is chosen', '    nn_gr = ', nn_gr
       IF(lwp) WRITE(numout,*) '~~~~~~~~'
       IF(lwp) WRITE(numout,*)

       clfgrn  = TRIM(cn_dir)//'gren_Gauss' 
       clfdgrn = TRIM(cn_dir)//'dgren_Gauss' 

       CALL getunit(numgr)      ! unit for green function file

       OPEN(UNIT=numgr, FILE=TRIM(clfgrn), STATUS='OLD', FORM="FORMATTED", ACCESS="SEQUENTIAL", iostat=istog)

       CALL getunit(numdgr)     ! unit for derivative of green function 

       OPEN(UNIT=numdgr, FILE=TRIM(clfdgrn), STATUS='OLD', FORM="FORMATTED", ACCESS="SEQUENTIAL", iostat=istodg)

       IF ( istog /=0 .OR. istodg /= 0 ) THEN
          IF (lwp) WRITE(numout,*) ' error in opening green function files, nn_gr =', nn_gr
          STOP
       END IF
       !
       istrg = 0
       istrdg = 0

       DO ji = 1, jpr
          READ(numgr, *, iostat = istrg  ) gren_d(ji)
          READ (numdgr, *, iostat = istrdg ) gren_dr(ji)
       END DO
       !
       IF ( istrg /= 0 .OR. istrdg /= 0 ) THEN
          IF(lwp) WRITE(numout,*) 'error in reading files for green function, nn_gr =', nn_gr
          STOP
       END IF
       !
    CASE ( 0 )            ! no green function
       !                  ! nn_gr = 0: nonfiltered
       IF(lwp) WRITE(numout,*)
       IF(lwp) WRITE(numout,'(a,a,i2)') 'no filter is chosen', '     nn_gr = ', nn_gr
       IF(lwp) WRITE(numout,*) '~~~~~~~~~'
       IF(lwp) WRITE(numout,*)
       !
    CASE DEFAULT
       STOP ' error in choosing the type of green function; nn_gr in [0,1,2]'
    END SELECT


    !
  END SUBROUTINE init_gren


  FUNCTION gren( pr )
    !!----------------------------------------------------------------------------
    !!                  ***  FUNCTION gren  ***
    !!
    !! ** Purpose : calculate green function for a given value pr 
    !!
    !! ** Method  :
    !!----------------------------------------------------------------------------
    REAL(wp), INTENT(in) :: pr
    REAL(wp)             :: gren
    INTEGER              :: ii
    REAL(wp)             :: zfrac
    !!----------------------------------------------------------------------------
    !
    SELECT CASE (nn_gr) 
       !   
    CASE ( 1 )     ! Agnesi Filter
       !
       IF ( pr < 0.1 ) THEN
          gren = 0.5_wp / ( rpi * pr )- 0.25_wp + 0.5_wp *pr / rpi

       ELSEIF ( pr <= 3.0) THEN
          ii = INT ( pr / rdr )
          zfrac = pr / rdr - ii
          gren = ( 1.0 - zfrac ) * gren_d(ii) + zfrac * gren_d(ii+1) 

       ELSEIF ( pr < rdr * jpr ) THEN
          ii =  NINT(REAL( pr / rdr , wp))
          gren = gren_d(ii)

       ELSE
          gren = ( -1.0 / pr**3.0 - 9.0 / pr**5.0 )/ ( 2.0 * rpi)

       ENDIF
       !
    CASE ( 2 )    ! Guassian filter  
       !
       IF ( pr < 0.2 ) THEN
          gren = 0.5 / ( rpi * pr) + ( -1. + pr**2 / 8.0 ) / ( 4.0 *SQRT(rpi) )

       ELSEIF ( pr <= 3.0) THEN
          ii = INT( pr / rdr )
          zfrac = pr/ rdr - REAL(ii, wp)
          gren = ( 1. - zfrac ) * gren_d(ii) + zfrac * gren_d(ii+1)

       ELSEIF ( pr < rdr * jpr ) THEN
          ii = NINT(REAL(pr / rdr, wp))
          gren = gren_d(ii)

       ELSE
          gren = ( -1. / pr**3.0 - 4.5 / pr**5.0 )/( 2.0 * rpi )

       ENDIF
       !
    CASE ( 0 )     ! no filter
       gren = 1.0 / ( 2.0 * rpi * pr )
       !
    END SELECT
    !
  END FUNCTION gren


  FUNCTION grendr ( pr )
    !!-----------------------------------------------------------------------------------
    !!                      *** FUNCTION dgren ***
    !! ** Purpose : calculate derivative of green function for a given value pr
    !!
    !! ** Method  
    !!-----------------------------------------------------------------------------------
    REAL(wp), INTENT(in) :: pr
    REAL(wp)             :: grendr
    REAL(wp)             :: zfrac
    INTEGER              :: ii
    !!----------------------------------------------------------------------------------
    !
    SELECT CASE (nn_gr) 
       !   
    CASE (1)                       ! Agnessi filter
       IF ( pr < 0.1) THEN
          grendr = ( 0.5 / rpi ) * ( 1.0 -1.0 / pr**2.0 )

       ELSEIF ( pr <= 3.0) THEN
          ii = INT( pr / rdr )
          zfrac = pr / rdr - REAL(ii, wp)
          grendr = ( 1.0 - zfrac) * gren_dr(ii) + zfrac * gren_dr(ii+1)

       ELSEIF ( pr < rdr * jpr ) THEN
          ii = NINT(REAL( (pr / rdr), wp)) 
          grendr = gren_dr(ii)

       ELSE
          grendr = (3.0 / pr**4.0 + 45.0 / pr**6.0 )/ ( 2.0 * rpi ) 

       ENDIF
       !
    CASE ( 2 )   !the Gaussian filter
       !         -1/(2pi pr^2)+r/(16sqrt(pi))*exp(pr*pr/8)*(I_0(pr*pr/8)-I_1(pr*pr/8))

       IF ( pr < 0.20) THEN
          grendr = -0.5_wp / ( rpi * pr **2.0 ) + pr / (16.0* SQRT(rpi) )

       ELSEIF ( pr <= 3.0) THEN
          ii = INT(pr / rdr)
          zfrac = pr / rdr - REAL(ii, wp)
          grendr = ( 1.0 - zfrac ) * gren_dr(ii) + zfrac * gren_dr(ii+1)

       ELSEIF ( pr < rdr * jpr ) THEN
          ii = NINT(REAL( pr / rdr, wp))
          grendr = gren_dr(ii)

       ELSE
          grendr = ( 3.0 / pr**4.0 + 22.50 / pr**6.0 ) / (2.0 * rpi)

       ENDIF
       !
    CASE ( 0 )    !no filter 
       !          ! -1./(2pi*pr**2)
       grendr = -1.0 / ( 2.0 * rpi * pr**2.0 )
       !
    END SELECT
    !
  END FUNCTION grendr

  !!===========================================================================================
END MODULE grn_tide
