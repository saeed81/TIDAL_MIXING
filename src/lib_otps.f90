MODULE lib_otps
  !!============================================================================
  !!                          *** MODULE lib_otps ***
  !! barotropic tide library  : extracting barotropic tide velocity
  !!
  !!============================================================================

  !!----------------------------------------------------------------------------
  !! init_tide        :  initialize tidal field for the requested component
  !! get_uv           :  extract horizontal barotropic tide from input data
  !!----------------------------------------------------------------------------
  USE par_tide        ! TMX parameters and variables
  USE in_manager      ! TMX input managers
  
  IMPLICIT NONE
  PRIVATE 

  INTERFACE getuv
     MODULE PROCEDURE get_uv
  END INTERFACE
  
  INTERFACE intpuv
     MODULE PROCEDURE interp_uv
  END INTERFACE

  INTERFACE mkmask
     MODULE PROCEDURE mk_mask
  END INTERFACE
     

  INTEGER, PARAMETER                        :: jpi = 1440     ! x size of barotropic field  
  INTEGER, PARAMETER                        :: jpj = 721      ! y size of barotropic field
  INTEGER, PARAMETER                        :: jpc = 10       ! the number of tidal components
  INTEGER                                   :: numtpx = -1    ! unit for tidal field
  INTEGER                                   :: numgrd = -1    ! unit for grid field
  INTEGER,  DIMENSION(2,jpi,jpj)            :: muv            ! maks array for uv
  INTEGER,  DIMENSION(jpi,jpj)              :: mz             ! aray 
  REAL(sp), DIMENSION(jpi,jpj)              :: rdepth         ! depth array
  COMPLEX(sp), DIMENSION(2,jpi,jpj)         :: yuv            ! uv barotropic field
  REAL(sp)                                  :: thmin, thmax   !
  REAL(sp)                                  :: rphmin, rphmax !

  PUBLIC :: init_tide !called by domini.f90   
  PUBLIC :: get_uv    !called by arayini.f90
 
CONTAINS 

  SUBROUTINE init_tide 
    !!------------------------------------------------------------------------------------
    !!                    *** ROUTINE init_tide ***
    !!
    !! ** Purpose  : get the tidal data for a given constituent
    !! 
    !! ** Method  : -Reads input files for tidal transport and bathymetry grid
    !!               
    !!------------------------------------------------------------------------------------- 
    INTEGER                            :: istat 
    INTEGER                            :: ij, ii
    INTEGER                            :: jk
    INTEGER                            :: ic
    CHARACTER(len=30)                  :: cl_tpx
    CHARACTER(len=30)                  :: cl_grd
    CHARACTER(len=4), DIMENSION(jpc)   :: cl_id
    CHARACTER(len=300)                 :: clfile
    !!
    NAMELIST/namtide/ cn_dir, cn_cu 
    !!-------------------------------------------------------------------------------------
    cn_dir='./'

    REWIND(numnam)
    READ(numnam,namtide)
    !
    IF( LEN_TRIM(cn_dir) /= INDEX(cn_dir,"/", back = .TRUE.)) THEN
       cn_dir = TRIM(cn_dir)//"/"   ! check whether name of dir out ends with / or not                                                                                                     
    END IF
    !
    SELECT CASE(TRIM(cn_cu))
       !
    CASE('m2')
       ncu = 1
       !
    CASE('s2')
       ncu = 2
       !
    CASE('n2')
       ncu = 3
       !
    CASE('k2')
       ncu = 4
       !
    CASE('k1')
       ncu = 5
       !
    CASE('o1')
       ncu = 6
       !
    CASE('p1')
       ncu = 7
       !
    CASE('q1')
       ncu = 8
       !
    CASE default
       STOP ' at least one component has to be provided by user in namelist'
    END SELECT

    IF(lwp) THEN
       WRITE(numout,*)
       WRITE(numout,*) 'tidal parameter '
       WRITE(numout,*) '~~~~~~~~~~~~~~~~'
       WRITE(numout,'(a,a,a,i2)') ' the tidal constituent ', TRIM(cn_cu),'   ncu = ', ncu 
    END IF


    cl_tpx = 'u_TPXO.6.1'
    clfile = TRIM(cn_dir)//TRIM(cl_tpx)

    CALL getunit(numtpx)

    OPEN( unit = numtpx, FILE = TRIM(clfile), FORM = 'UNFORMATTED', STATUS='OLD', iostat = istat)

    IF( istat /=0 ) THEN
       IF (lwp) THEN
          WRITE(numout,*) ' error in opening u_TPXO6.1 '
       END IF
       STOP
    END IF

    !
    muv(:,:,:) = 0
    !
    mz(:,:) = 0
    !
    rdepth(:,:) = 0.0
    !
    yuv(:,:,:) = (0.0, 0.0)
    !

    READ(numtpx, iostat = istat) ii, ij, ic, thmin, thmax, &
         &   rphmin, rphmax, cl_id

    IF( istat /= 0 ) THEN  
       IF(lwp) THEN
          WRITE(numout,*) ' error in reading u_TPXO6.1 '
       END IF
       STOP
    END IF

    IF( jpi /= ii .OR. jpj /= ij .OR. jpc /= ic ) THEN
       IF(lwp) THEN  
          WRITE(numout,*)'Error in init_tide'
          WRITE(numout,*)'eror in grid size '
          WRITE(numout,*)'jpj = ', ij
          WRITE(numout,*)'jpi = ', ii
          WRITE(numout,*)'jpc = ', ic
       ENDIF
       STOP
    END IF

    IF(lwp) THEN
       WRITE(numout,*)'Model file:'
       WRITE(numout,*) TRIM(clfile)
       WRITE(numout,*)'Grid zise:', jpi, jpj
       WRITE(numout,*)'Limits:', thmin,' - ', thmax
       WRITE(numout,*)'       ', rphmin,' - ', rphmax
       WRITE(numout,*)'Constituents: ', cl_id
       WRITE(numout,*)
    END IF
    !

    DO jk = 1, ncu
       READ(numtpx) yuv
    ENDDO

    CALL getunit(numgrd)

    cl_grd = 'grid'
    clfile = TRIM(cn_dir)//TRIM(cl_grd)

    OPEN( UNIT = numgrd, FILE = TRIM(clfile), FORM = 'UNFORMATTED', STATUS='OLD', iostat=istat  )

    IF( istat /=0 ) THEN
       IF(lwp) THEN
          WRITE(numout,*) ' error in opening grid file '
       END IF
       STOP
    END IF

    READ(numgrd, iostat = istat ) ii, ij

    READ(numgrd)
    !
    IF( ii /= jpi .OR. ij /= jpj .OR. istat /= 0  ) THEN
       IF(lwp) THEN 
          WRITE(numout,*)'Error in init_tide'
          WRITE(numout,*) 'Size or arrays in H_ file', ii, ij
          WRITE(numout,*) 'Wrong grid!'
       END IF
       STOP
    END IF
    !
    READ(numgrd) rdepth
    READ(numgrd) mz

    CALL mk_mask( yuv, muv ) ! make array mask
    !
  END SUBROUTINE init_tide


  SUBROUTINE get_uv(plat, plon, yduv)
    !!---------------------------------------------------------------------------
    !!                   *** ROUTINE get_uv ***
    !!
    !! ** Purpose : extract barotropic tidal velocity in lateral direction
    !!---------------------------------------------------------------------------
    REAL(sp), INTENT(in)                 :: plon, plat
    COMPLEX,  DIMENSION(:), INTENT(out)  :: yduv
    REAL(sp)                             :: zlon, zlat
    REAL(sp)                             :: zdepth
    LOGICAL                              :: llphase
    !!----------------------------------------------------------------------------
    !
    ! convert to single precision
    zlon = REAL(plon, sp)
    zlat = REAL(plat, sp)

    IF( zlon < rphmin) zlon = zlon + 360.0 ! This is necessary so that get_uv can handle negative longitudes.
    IF( zlon > rphmax) zlon = zlon - 360.0
    !
    ! output transports or velocities with our +iwt phase convention
    ! (l_phase = .true. for standard phase -iwt)
    llphase = .FALSE.
    !
    CALL interp_depth( zlat, zlon, zdepth ) !interpolates the array depth to the point lat,lon 


    CALL interp_uv( zlat, zlon, yduv )

    IF( llphase ) THEN
       yduv(1) = CONJG(yduv(1))
       yduv(2) = CONJG(yduv(2))
    ENDIF

    IF( zdepth > 0.0) THEN
       ! transports m**2/s ->  currents cm/s
       !         uv1(1) = 100.*uv1(1)/depth1
       !         uv1(2) = 100.*uv1(2)/depth1
       ! transports m**2/s ->  currents m/s
       yduv(1) = yduv(1) / zdepth
       yduv(2) = yduv(2) / zdepth
    
    ELSEIF( zdepth == 0.0 ) THEN
       yduv(1) = 0.0
       yduv(2) = 0.0
         
    ENDIF
    !
  END SUBROUTINE get_uv


  SUBROUTINE interp_uv ( plat, plon, yduv )
    !!-----------------------------------------------------------------------------
    !!                        *** ROUTINE interp_uv ***
    !!
    !! ** Purpose : interpolates complex array yduv onto point plat, plon
    !!      
    !!-----------------------------------------------------------------------------
    REAL(sp), INTENT(in)    :: plat
    REAL(sp), INTENT(in)    :: plon
    COMPLEX(sp), DIMENSION(:), INTENT(out) :: yduv  
    !!
    REAL(sp), DIMENSION(2) ::  zx, zy
    INTEGER                ::  jk           ! integer used for loops
    REAL(wp)               ::  zdphi, zdtheta
    REAL(sp)               ::  zx1 , zy1
    REAL(sp)               ::  zx2 , zy2
    REAL(sp)               ::  zw12, zw14
    REAL(sp)               ::  zw32, zw34
    REAL(sp)               ::  zwtot
    LOGICAL                ::  llbndi1, llbndi2
    LOGICAL                ::  llbndi3, llbndi4
    INTEGER                ::  i1, i2, i3, i4
    !!------------------------------------------------------------------------------
    !
    zx = (/1.0, 0.5/)
    zy = (/.50, 1.0/)

    zdphi   = ( rphmax - rphmin  ) / jpi
    zdtheta = ( thmax  -  thmin  ) / jpj

    DO jk = 1, SIZE(yduv,1)
       !      U component jk = 1
       !      V component jk = 2
       zx1 = ( plon - rphmin ) / zdphi   + zx(jk)
       zy1 = ( plat - thmin  ) / zdtheta + zy(jk)

       IF(  zx1 < 1.0      )  zx1 = zx1 +  jpi
       IF(  zx1 > jpi + 1  )  zx1 = zx1 - jpi

       i1  = INT(zx1)
       zx2 = zx1 - i1
       i2  = INT(zy1)
       zy2 = zy1 - i2

       i3 = ipshft( i1, 1, jpi )
       i4 = ipshft( i2, 1, jpj )

       llbndi1 = i1 < 1 .OR. i1 > jpi
       llbndi2 = i2 < 1 .OR. i2 > jpj
       llbndi3 = i3 < 1 .OR. i3 > jpi
       llbndi4 = i4 < 1 .OR. i4 > jpj

       IF ( llbndi1 .OR. llbndi2 .OR. llbndi3 .OR. llbndi4 ) THEN
          yduv(jk) = (0.0,0.0)
       ELSE
          zw12  = ( 1.0 - zx2 ) * (1.0 - zy2 ) * muv(jk,i1,i2)
          zw14  = ( 1.0 - zx2 ) * zy2 * muv(jk,i1,i4)
          zw32  = zx2 * ( 1.0 - zy2 ) * muv(jk,i3,i2)
          zw34  = zx2 * zy2 * muv(jk,i3,i4)
          zwtot = zw12 + zw14 + zw32 + zw34

          IF(  zwtot == 0.0 ) THEN

             yduv(jk) = (0.0,0.0)
          ELSE

             yduv(jk) = yuv(jk,i1,i2) * zw12 + yuv(jk,i3,i2) * zw32 + &
                  &   yuv(jk,i1,i4) * zw14 + yuv(jk,i3,i4) * zw34
             yduv(jk) = yduv(jk) / zwtot

          ENDIF
       END IF
       !
    ENDDO

  END SUBROUTINE interp_uv


  SUBROUTINE interp_depth( plat, plon, pdepth )
    !!--------------------------------------------------------------------------------     
    !!                   *** ROUTINE interp_depth ***   
    !!
    !!  ** Purpose :   Interpolates real array depth onto point plat, plon
    !!  
    !!---------------------------------------------------------------------------------  
    REAL(sp), INTENT(in )  :: plat, plon
    REAL(sp), INTENT(out)  :: pdepth
    REAL(sp), DIMENSION(1) :: zx, zy
    REAL(sp)               :: ztheta 
    REAL(sp)               :: zphi
    REAL(sp)               :: zx1 , zy1
    REAL(sp)               :: zx2 , zy2
    REAL(sp)               :: zw12, zw14
    REAL(sp)               :: zw32, zw34
    REAL(sp)               :: zwtot
    LOGICAL                :: llbndi1, llbndi2
    LOGICAL                :: llbndi3, llbndi4
    INTEGER                :: i1, i2, i3, i4
    !!---------------------------------------------------------------------------------
    zx = (/0.5/)
    zy = (/0.5/)

    ztheta = ( thmax -  thmin  ) / jpj
    zphi   = ( rphmax - rphmin ) / jpi

    zx1 = ( plon - rphmin ) / zphi   + zx(1)
    zy1 = ( plat - thmin  ) / ztheta + zy(1)

    IF( zx1 < 1.0     ) zx1 =  jpi + zx1
    IF( zx1 > jpi + 1 ) zx1 =  zx1 - jpi

    i1  = INT(zx1)
    zx2 = zx1 - i1
    i2  = INT(zy1)
    zy2 = zy1 - i2
    i3  = ipshft(i1, 1, jpi)
    i4  = ipshft(i2, 1, jpj)

    llbndi1 = i1 < 1 .OR. i1 > jpi
    llbndi2 = i2 < 1 .OR. i2 > jpj
    llbndi3 = i3 < 1 .OR. i3 > jpi
    llbndi4 = i4 < 1 .OR. i4 > jpj

    IF ( llbndi1 .OR. llbndi2 .OR. llbndi3 .OR. llbndi4 ) THEN
       pdepth = 0.0
    ELSE
       zw12  = ( 1.0 - zx2 ) * (1.0 - zy2) * mz(i1,i2)
       zw14  = ( 1.0 - zx2 ) * zy2 * mz(i1,i4)
       zw32  = zx2 * ( 1.0 - zy2 ) * mz(i3,i2)
       zw34  = zx2 * zy2 * mz(i3,i4)
       zwtot = zw12 + zw14 + zw32 + zw34

       IF( zwtot == 0.0 ) THEN

          pdepth = 0.0
       ELSE

          pdepth = rdepth(i1,i2) * zw12 + rdepth(i3,i2) * zw32 + &
               &   rdepth(i1,i4) * zw14 + rdepth(i3,i4) * zw34
          pdepth = pdepth / zwtot
       ENDIF
    END IF
    !
  END SUBROUTINE interp_depth


  FUNCTION ipshft( ki, kish, kn ) 
    !!-----------------------------------------------------------------
    !!                *** ROUTINE ipshift
    !! 
    !!  ** Purpose :  periodic shift maps ki to ki + kish, mod kn
    !!                always between 1 and kn;  never 0 
    !!-----------------------------------------------------------------
    INTEGER, INTENT(in) :: ki, kn, kish
    INTEGER             :: ipshft
    !!------------------------------------------------------------------
    !
    ipshft = MOD(ki + kish + kn -1 , kn) + 1 
    !
  END FUNCTION ipshft


  SUBROUTINE mk_mask ( yduv, kuv )
    !!---------------------------------------------------------------------
    !!                     *** ROUTINE mk_mask ***
    !!
    !! ** Purpose :   create the integer mask file from  magnitude of yduv  
    !!---------------------------------------------------------------------
    COMPLEX(sp), DIMENSION(:,:,:), INTENT(in )  :: yduv
    INTEGER,     DIMENSION(:,:,:), INTENT(out)  :: kuv
    !!---------------------------------------------------------------------
    !
    WHERE( ABS( yduv(:,:,:)) > 0.0 ) 
       kuv(:,:,:) = 1
    ELSEWHERE
       kuv(:,:,:) = 0
    END WHERE
    !
  END SUBROUTINE mk_mask

  !!==========================================================================================
END MODULE lib_otps
