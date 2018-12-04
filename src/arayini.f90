MODULE arayini
  !!==================================================================================
  !!                            *** MODULE arayini ***   
  !! Array Initialization   : Initialization of arrays for tmx system
  !!==================================================================================

  !!----------------------------------------------------------------------------------
  !!   getdimension         :
  !!   init_hydro           :
  !!   init_metric          :
  !!   init_deriv           :
  !!   init_fcor            :
  !!   comp_cutoff          :
  !!   comp_crt             :
  !!   init_uv              :
  !!----------------------------------------------------------------------------------
  USE par_tide              !
  USE phycst                !
  USE io_ezcdf              !
  USE lib_otps              !
  USE arraydef              !

  IMPLICIT NONE
  PRIVATE

  PUBLIC ::  getdimension   !called by domini.f90
  PUBLIC ::  init_hydro     !called by domini.f90
  PUBLIC ::  init_metric    !called by domini.f90
  PUBLIC ::  init_deriv     !called by domini.f90
  PUBLIC ::  init_fcor      !called by domini.f90
  PUBLIC ::  comp_cutoff    !called by domini.f90
  PUBLIC ::  comp_crt       !called by domini.f90
  PUBLIC ::  init_uv        !called by domini.f90

CONTAINS

  SUBROUTINE getdimension
    !!---------------------------------------------------------------------------------
    !!                       *** ROUTINE getdimension ***
    !!
    !! ** Purpose : get the dimension of input field
    !!
    !!---------------------------------------------------------------------------------
    INTEGER :: ix1, ix2, ix3, ix4 
    INTEGER :: iy1, iy2, iy3, iy4
    INTEGER :: iz1, iz2, iz3, iz4
    INTEGER :: it1, it2, it3, it4
    !!---------------------------------------------------------------------------------
    IF ( .NOT. ln_modal ) THEN

       CALL DIMS( cfile_h, cvar_h, ix1, iy1, iz1, it1 )  
       CALL DIMS( cfile_bf, cvar_bf, ix2, iy2, iz2, it2 ) 
       CALL DIMS( cfile_intbf, cvar_intbf, ix3, iy3, iz3, it3 )

       IF ( (ix1 /= ix2) .OR. (ix1 /= ix3) .OR. (ix2 /= ix3) ) THEN
          PRINT *, 'the size of input variables in x direction are not the same'
          !STOP
       END IF

       IF ( (iy1 /= iy2) .OR. (iy1 /= iy3) .OR. (iy2 /= iy3) ) THEN
          PRINT *, 'the size of input variables in y direction are not the same'
          !STOP
       END IF

       IF ( (iz1 /= iz2) .OR. (iz1 /= iz3) .OR. (iz2 /= iz3) ) THEN
          PRINT *, 'the size of input variables in z direction are not the same'
          !STOP
       END IF

       IF ( (it1 /= it2) .OR. (it1 /= it3) .OR. (it2 /= it3) ) THEN
          PRINT *, 'the size of input variables in t direction are not the same'
          !STOP

       END IF

    ELSE
       CALL DIMS( cfile_h, cvar_h, ix1, iy1, iz1, it1 )  
       CALL DIMS( cfile_sn, cvar_sn, ix2, iy2, iz2, it2 ) 
       CALL DIMS( cfile_eps, cvar_eps, ix3, iy3, iz3, it3 )
       CALL DIMS( cfile_speed, cvar_speed, ix4, iy4, iz4, it4 )

       IF ( (ix1 /= ix2) .OR. (ix1 /= ix3) .OR. (ix1 /= ix4) .OR. (ix2 /= ix3) .OR. (ix2 /= ix4) .OR. (ix3 /= ix4)) THEN
          PRINT *, 'the size of input variables in x direction are not the same'
          STOP
       END IF

       IF ( (iy1 /= iy2) .OR. (iy1 /= iy3) .OR. (iy1 /= iy4) .OR. (iy2 /= iy3) .OR. (iy2 /= iy4) .OR. (iy3 /= iy4)) THEN
          PRINT *, 'the size of input variables in y direction are not the same'
          STOP
       END IF

       IF ( (iz1 /= iz2) .OR. (iz1 /= iz3) .OR. (iz1 /= iz4) .OR. (iz2 /= iz3) .OR. (iz2 /= iz4) .OR. (iz3 /= iz4) ) THEN
          PRINT *, 'the size of input variables in z direction are not the same'
          STOP
       END IF

       IF ( (it1 /= it2) .OR. (it1 /= it3) .OR. (it1 /= it4) .OR. (it2 /= it3) .OR. (it2 /= it4) .OR. (it3 /= it4) ) THEN
          PRINT *, 'the size of input variables in t direction are not the same'
          STOP
       END IF

    END IF
    !


    nx = ix1
    ny = iy1
    nz = iz1
    nt = it1
    !
    IF(lwp) THEN
       WRITE(numout,*) 'Dimension of input field is obtained'
       WRITE(numout,*) '~~~~~~~~~'
       WRITE(numout,*) ' the number of points in x direction', nx
       WRITE(numout,*) ' the number of points in y direction', ny
       WRITE(numout,*) ' the number of points in z direction', nz
       WRITE(numout,*) ' the number of points in t direction', nt
    END IF
    !
  END SUBROUTINE getdimension


  SUBROUTINE init_hydro
    !!--------------------------------------------------------------------------------
    !!                  *** ROUTINE init_hydro ***
    !!
    !! ** Purpose : Initialization of hydrography for tmx model
    !!--------------------------------------------------------------------------------
    INTEGER :: idx_f1, idx_f2, idx_f3, idx_f4
    INTEGER :: idx_v1, idx_v2, idx_v3, idx_v4
    !!--------------------------------------------------------------------------------
    !
    CALL GETVAR_1D (cfile_h, cvarlon, nx, vlon)
    CALL GETVAR_1D (cfile_h, cvarlat, ny, vlat)

    
    CALL GETVAR_2D (idx_f1, idx_v1, cfile_h, cvar_h, nx, ny, nt, 1, 1, h)

    IF ( .NOT. ln_modal ) THEN
       CALL GETVAR_2D (idx_f2, idx_v2, cfile_bf, cvar_bf, nx, ny, nt, 1, 1, bfbt)
       CALL GETVAR_2D (idx_f3, idx_v3, cfile_intbf, cvar_intbf, nx, ny, nt, 1, 1, bfint)
    ELSE

       CALL GETVAR_2D (idx_f2, idx_v2, cfile_sn, cvar_sn, nx, ny, nt, 1, 1, sn)
       CALL GETVAR_2D (idx_f3, idx_v3, cfile_eps, cvar_eps, nx, ny, nt, 1, 1, eps)
       CALL GETVAR_2D (idx_f4, idx_v4, cfile_speed, cvar_speed, nx, ny, nt, 1, 1, speed)

    END IF

    h(:,:) =  h(:,:)
    !
  END SUBROUTINE init_hydro


  SUBROUTINE init_metric
    !!------------------------------------------------------------
    !!                 *** ROUTINE init_metric ***
    !!
    !! ** Purpose : initialize metric for tmx model
    !!------------------------------------------------------------
    REAL(wp) :: zdxdeg 
    REAL(wp) :: zdydeg 
    INTEGER  :: ji, jj
    !!------------------------------------------------------------
    !
    zdxdeg = ( vlon(nx) - vlon(1) ) / REAL(nx - 1, wp)

    dx(:) = rad * zdxdeg * ra * COS( rad * vlat(:) )

    DO jj = 1, ny - 1
       zdydeg = vlat(jj+1) - vlat(jj)
       dy(jj) = rad * zdydeg * ra
    ENDDO
    !
    dy(ny) = dy(ny-1)

    DO ji = 1, nx - 1                                                                                                         
       vlons(ji) = 0.5 * ( vlon(ji) + vlon(ji+1) )                                                                             
    ENDDO

    vlons(nx) = vlon(nx)                                                                                                      
    DO jj = 1, ny - 1                                                                                                         
       vlats(jj) = 0.5 * ( vlat(jj) + vlat(jj+1) )                                                                          
    ENDDO

    vlats(ny) = vlat(ny)

    DO ji = 1, nx

       DA(ji,:) = dx(:) * dy(:)

    END DO
    !
  END SUBROUTINE init_metric


  SUBROUTINE init_deriv
    !!--------------------------------------------------------------------
    !!                      *** ROUTINE init_deriv ***
    !!
    !! ** Purpose : initialize derivative arrays  
    !!--------------------------------------------------------------------
    INTEGER :: ji, jj
    !!--------------------------------------------------------------------
    !
    DO jj = 1, ny - 1
       DO ji = 1, nx - 1
          dhdx(ji,jj) = 0.5*(h(ji+1,jj) + h(ji+1,jj+1) - h(ji,jj) - h(ji,jj+1) ) / dx(jj)
          dhdy(ji,jj) = 0.5*(h(ji,jj+1) + h(ji+1,jj+1) - h(ji,jj) - h(ji+1,jj) ) / dy(jj)
       ENDDO
    ENDDO

    DO ji = 1, nx - 1
       dhdxs(ji,:) = ( h(ji+1,:) - h(ji,:) ) / dx(:)
    ENDDO

    DO jj = 1, ny - 1
       dhdys(:,jj) = ( h(:,jj+1) - h(:,jj) ) / dy(jj)
    ENDDO

    dqdx(:,:)    = 0.0
    dqdy(:,:)    = 0.0
    d2qdx2(:,:)  = 0.0
    d2qdy2(:,:)  = 0.0
    d2qdxdy(:,:) = 0.0
    !
  END SUBROUTINE init_deriv


  SUBROUTINE init_fcor
    !!----------------------------------------------------------------------------
    !!                    *** ROUTINE init_fcor ***
    !!
    !! ** Purpose : compute Coriolis parameter with the use of latitude array
    !!-----------------------------------------------------------------------------

    fcor(:)  = 2.0 * omega * SIN( rad * vlat(:) )

    fcors(:) = 2.0 * omega * SIN( rad * vlats(:) )                                                                           
    !
  END SUBROUTINE init_fcor


  SUBROUTINE comp_cutoff
    !!-------------------------------------------------------------
    !!                    *** ROUTINE comp_cutoff
    !!
    !! ** Purpose : compute cutoff length for tmx system
    !!-------------------------------------------------------------
    REAL(wp) :: zcut, znomega, zslope
    INTEGER  :: ji, jj 
    !!-------------------------------------------------------------
    vcut(:,:) = -999.00
    !                                                                                                                          
    IF ( ln_energy .AND. ( .NOT. ln_modal )) THEN
       DO jj = 1, ny
          zcut = ( tfrq(ncu) * tfrq(ncu) ) - ( fcors(jj) * fcors(jj) )
          IF (zcut > 0.0 ) THEN
             DO ji = 1, nx
                
                IF ( h(ji,jj) > 0.0 ) THEN

                vcut(ji,jj) = ( rn_beta * bfint(ji,jj) ) / ( rpi * SQRT(zcut) )
                nxinterval(ji,jj) = INT( rn_trunc * vcut(ji,jj) / dx(jj) )
                nyinterval(ji,jj) = INT( rn_trunc * vcut(ji,jj) / dy(jj) )
             
             END IF

             END DO
          END IF
       END DO

    END IF

    IF ( ln_energy .AND. ln_modal ) THEN

       DO ji = 1 , nx
          WHERE( ABS(tfrq(ncu)) >= ABS(fcors(:)) .AND. speed(ji,:) > 0.0 ) 
             vcut(ji,:) = SQRT(tfrq(ncu) * tfrq(ncu) - fcors(:) * fcors(:)) / speed(ji,:)
             nxinterval(ji,:) = INT( rn_trunc / ( vcut(ji,:) * dx(:) ) ) 
             nyinterval(ji,:) = INT( rn_trunc / ( vcut(ji,:) * dy(:) ) )

          END WHERE
       END DO

    END IF



    IF ( ln_velocity ) THEN
       DO jj = 1, ny
          zcut = ( tfrq(ncu) * tfrq(ncu) ) - ( fcor(jj) * fcor(jj) )
          IF (zcut > 0.0 ) THEN
             DO ji = 1, nx
                vcut(ji,jj) = ( rn_beta * bfint(ji,jj) ) / ( rpi * SQRT(zcut) )
                nxinterval(ji,jj) = INT( rn_trunc * vcut(ji,jj) / dx(jj) )
                nyinterval(ji,jj) = INT( rn_trunc * vcut(ji,jj) / dy(jj) )
             END DO
          END IF
       END DO
    END IF



    tmask(:,:) = 1.0_wp

    WHERE ( vcut(:,:) == -999.0 ) tmask(:,:) = 0.0_wp


    IF ( .NOT. ln_modal ) THEN


       vsteep(:,:) = -999.0

       vhslope(:,:) = -999.0

    

       DO jj = 1, ny

          zcut = ( tfrq(ncu) * tfrq(ncu) ) - ( fcors(jj) * fcors(jj) )

          IF (zcut > 0.0 ) THEN

             DO ji = 1, nx

                znomega = bfbt(ji,jj) * bfbt(ji,jj)- tfrq(ncu) * tfrq(ncu)

                IF ( h(ji,jj) /= 0.0 .AND. znomega >= 0.0 ) THEN

                   zslope = SQRT ( znomega / zcut )

                   vsteep(ji,jj)  = SQRT ( dhdxs(ji,jj) * dhdxs(ji,jj) + dhdys(ji,jj) * dhdys(ji,jj) ) * zslope   

                
                END IF

             END DO

          END IF

       END DO

       vhslope(:,:) = dhdxs(:,:) * dhdxs(:,:) + dhdys(:,:) * dhdys(:,:) 
    END IF
       !
  END SUBROUTINE comp_cutoff


  SUBROUTINE comp_crt
    !!------------------------------------------------------------------------
    !!                    ***  ROUTINE comp_crt ***
    !!
    !! ** Purpose : compute the correction array used for energy and velocity 
    !!              computation 
    !!------------------------------------------------------------------------
    !
    asx(:) = dx(:) / dy(:)
    asy(:) = dy(:) / dx(:)
    !
    IF ( ln_energy ) THEN  
       SELECT CASE(nn_corr) 

       CASE ( 1 )
          vcorx(:) = -4.0 * asx(:)**2.0 / ( 1.0 + asx(:)**2)**1.5 + 2.0 * LOG( asx(:) + &
               &   SQRT( 1.0 + asx(:)**2.0 )) / asx(:)
          vcory(:) = -4.0 * asy(:)**2.0 / ( 1.0 + asy(:)**2)**1.5 + 2.0 * LOG( asy(:) + &
               &   SQRT( 1.0 + asy(:)**2.0) ) / asy(:)
          !
       CASE ( 2 ) 
          vcorx(:) = -asx(:)**2.0 / ( 1.0 + asx(:)**2.0 )**1.5 + LOG( asx(:) + &
               &   SQRT( 1.0 + asx(:)**2.0 ) ) / asx(:)
          vcory(:) = -asy(:)**2.0 / ( 1.0 + asy(:)**2 )**1.5 + LOG( asy(:) + &
               &   SQRT( 1.0 + asy(:)**2.0 ) ) / asy(:)
          !
       CASE ( 0 )
          vcorx(:) = 0.0_wp
          vcory(:) = 0.0_wp
          !
       END SELECT
       !
    END IF
    !
    IF ( ln_velocity ) THEN  

       SELECT CASE(nn_corr) 
          !
       CASE ( 1, 2 )
          !
          vcorx(:) =  (LOG ( 2.0 * asx(:) + SQRT( 1.0 + 4.0 * asx(:)**2.0_wp ) ) - 2.0 ) * dy(:) / dx(:)**2.0

          vcory(:) =  (LOG( 2.0 * asy(:) + SQRT(1.0 + 4.0 * asy(:)**2.0) ) - 2.0_wp) * dx(:) / dy(:)**2.0 


          vcorxy(:) = 0.5 * ( LOG( asx(:) + SQRT( 1.0 + asx(:)**2.0 )) / dx(:)+ &
               &   LOG( asy(:) + SQRT( 1.0_wp + asy(:)**2.0 )) / dy(:) - 2.0 / SQRT( dx(:)**2.0 &
               & + dy(:)**2.0 ))
          !
       CASE ( 0 )
          vcorx(:)  = 0.0
          vcory(:)  = 0.0
          vcorxy(:) = 0.0 
          !
       END SELECT
       !
    END IF
    !
  END SUBROUTINE comp_crt


  SUBROUTINE init_uv
    !!----------------------------------------------------------------
    !!                *** ROUTINE init_uv ***
    !!
    !! ** Purpose : initialize basic barotropic tidal currrent
    !!
    !!----------------------------------------------------------------
    INTEGER :: ji, jj
    COMPLEX, DIMENSION(2) :: yluv
    !!----------------------------------------------------------------
    yluv(:) = (0.0,0.0) 
    !
    ureal(:,:) = 0.0
    uimag(:,:) = 0.0
    vreal(:,:) = 0.0
    vimag(:,:) = 0.0

    IF ( ln_energy ) THEN
       DO jj = 1, ny
          DO ji = 1, nx
             CALL get_uv( vlats(jj), vlons(ji), yluv )
             ureal(ji,jj) = REAL(yluv(1))
             uimag(ji,jj) = AIMAG(yluv(1))
             vreal(ji,jj) = REAL(yluv(2))
             vimag(ji,jj) = AIMAG(yluv(2))
          END DO
       END DO
    END IF
    !
    IF ( ln_velocity ) THEN
       DO jj = 1, ny
          DO ji = 1, nx
             CALL get_uv( vlat(jj), vlon(ji), yluv )
             ureal(ji,jj) = REAL(yluv(1))
             uimag(ji,jj) = AIMAG(yluv(1))
             vreal(ji,jj) = REAL(yluv(2))
             vimag(ji,jj) = AIMAG(yluv(2))
          END DO
       END DO
    END IF

    WHERE( h(:,:) == 0.0 ) 

       ureal(:,:) = 0.0
       uimag(:,:) = 0.0
       vreal(:,:) = 0.0
       vimag(:,:) = 0.0

    END WHERE

    !
  END SUBROUTINE init_uv

  !!=================================================================================================
END MODULE arayini
