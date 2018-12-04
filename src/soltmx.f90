MODULE soltmx
  !!===============================================================
  !!                         *** MOULE soltmx ***
  !! Tidal mixing solver :  energy and velocity solvers
  !!===============================================================

  !!----------------------------------------------------------------
  !! sol_par         :   reading solver parametes from namelist
  !! sol_energy      :   compute energy conversion rate
  !! sol_velocity    :   compute baroclinic tidal velocity
  !!----------------------------------------------------------------
  USE par_tide       ! TMX system parameters and variables
  USE arraydef       ! Definition of arrays used in tmx model
  USE phycst         ! Physical constant used in the model 
  USE grn_tide       ! Definition of green functions
  USE wrtout         ! Out-put manager for tmx system
  USE bessel_jy0     ! library for bessel function of zero order 
  USE bessel_jy1     ! library for bessel function of first order

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: sol_par             ! called by tmx.f90
  PUBLIC :: sol_energy          ! called by tmx.f90
  PUBLIC :: sol_velocity        ! called by tmx.f90

CONTAINS 

  SUBROUTINE sol_par
    !!----------------------------------------------------------------------------------------
    !!                   *** ROUTINE par_sol ***
    !!
    !! ** Purpose :  Reading solver parameters from namelist 
    !!----------------------------------------------------------------------------------------
    NAMELIST/namsol/ ln_modal, ln_energy, ln_velocity, nn_bfz, rn_trunc, nn_corr, nn_mode
    !!----------------------------------------------------------------------------------------
    READ(numnam, namsol)  ! read which field variable should be computed; enery and velocity        

    IF( ( .NOT. ln_energy ) .AND. ( .NOT. ln_velocity ) ) THEN
       PRINT *, " At least one variable should be chosen in namelist to be computed "
       STOP     " change the flag for energy or velocity in namelist "
    END IF
    !
  END SUBROUTINE sol_par


  SUBROUTINE sol_energy
    !!-------------------------------------------------------------------
    !!                      *** ROUTINE sol_energy  ***
    !! 
    !! ** Purpose :  Calculate enery conversion rate from the interaction 
    !!               of barotropic tide with bottom topography
    !!
    !! ** Method  : -The energy is calculated on the staggered grid, with the regular grid                                                     
    !!               defined as the points with topographic data h(ji,jj).
    !!              -An improved finite difference scheme is implemented, with two versions                                                    
    !!               for the energy calculation:                                                                                               
    !!              -nn_corr=1: the convolution integrals over the four cells closest to                                                      
    !!               the field point are calculated analytically, using the taylor                                                          
    !!               expansion of the topography around the field point.                                                                    
    !!              -nn_corr=2: analytic integration is done between the four h-point                                                         
    !!               closest to the staggered field point, using the linearly                                                               
    !!               interpolated topography between these four points. this is the                                                         
    !!               same as nn_corr, except that the analytic integration is done only 
    !!               over 1/4 of the four closest cells.                 
    !!                                                                                                           
    !! 
    !! References : References : Nycander, J. 2005 Generation of internal waves in the deep
    !!              ocean by tides. J. Geophys. Res., C10028, 110.   
    !!--------------------------------------------------------------------
    INTEGER  :: ji, jj, jk, jl     ! dummy integers used for loops
    INTEGER  :: is, iu, iw, ie     ! local integer used for loops
    INTEGER  :: iis, iiu, iiw, iie ! local integer used for loops
    REAL(wp) :: zdist              ! local real variable for calculating distance
    REAL(wp) :: zda                ! local real variable for calcualting area
    REAL(wp) :: zxdif, zydif       ! local real varaibles for calculating distance in x &
    ! and y direction
    REAL(wp) :: zfcrm, zprefac
    LOGICAL  :: llbathy = .TRUE.
    LOGICAL  :: llbsl0  = .TRUE.   
    REAL(wp) :: zin1, zin2, zin3, zint
    !!--------------------------------------------------------------------
    !
    CALL name_outeng
    
    IF ( ln_bt     ) CALL output_bt

    IF ( ln_cutoff ) CALL output_cutoff

    IF ( ln_steep  ) CALL output_steep

    IF ( ln_hslope ) CALL output_hslope

!
    encon(:,:) = 0.0_wp

    is = nn_bfz + 2
    iu = SIZE(h,2) - nn_bfz
    iw = nn_bfz + 2
    ie = SIZE(h,1) - nn_bfz
    !

    IF ( .NOT. ln_modal )   THEN

       DO jj = is, iu

          !!PRINT *, jj, narea       !monitoring computation for each lattitude

          IF ( lwp ) WRITE(numout,*)  jj, narea       !monitoring computation for each lattitude

          DO ji = iw, ie

             IF ( vcut(ji,jj) > 0.0 ) THEN

                !vcut(ji,jj) = ABS(vcut(ji,jj))

                iiw = MAX( 1, ji - nxinterval(ji,jj) )
                !                                                                                                                                                                             
                iie = MIN( nx, ji + nxinterval(ji,jj) )
                !                                                                                                                                                                             
                iis = MAX( 1, jj - nyinterval(ji,jj) )
                !                                                                                                                                                                              
                iiu = MIN( ny, jj + nyinterval(ji,jj) )
                !                                                                                                                                                                              

                DO jl = iis, iiu
                   DO jk = iiw, iie

                      zda = dx(jl) * dy(jl)

                      zxdif = dx(INT(( jj + jl ) / 2.0)) * REAL( ( ji - jk + 0.5 ), wp)
                      zydif = dy(INT(( jj + jl ) / 2.0)) * REAL( ( jj - jl + 0.5 ), wp)
                      !
                      zdist = SQRT( zxdif * zxdif + zydif * zydif)
                      zdist = MIN( zdist, rn_trunc * vcut(ji,jj))

                      IF ( llbathy ) THEN  
                         dqdx(ji,jj) = dqdx(ji,jj) + zda * ( h(jk,jl) - h(ji,jj) )  * zxdif * ( 1.0 / vcut(ji,jj) )**2.0 &
                              &   * grendr( zdist / vcut(ji,jj) ) / zdist
                         dqdy(ji,jj) = dqdy(ji,jj) + zda * ( h(jk,jl) - h(ji,jj) )* zydif * ( 1.0 / vcut(ji,jj) )**2.0 &
                              &   * grendr( zdist / vcut(ji,jj) ) / zdist
                      ELSE
                         dqdx(ji,jj) = dqdx(ji,jj) + zda * h(jk,jl) * zxdif * ( 1.0 / vcut(ji,jj) )**2.0 &
                              &   * grendr( zdist / vcut(ji,jj) ) / zdist
                         dqdy(ji,jj) = dqdy(ji,jj) + zda * h(jk,jl)* zydif * ( 1.0 / vcut(ji,jj) )**2.0 &
                              &   * grendr( zdist / vcut(ji,jj) ) / zdist
                      END IF
                   ENDDO
                ENDDO

                zfcrm = fcors(jj) / tfrq(ncu)
                zprefac = SQRT( 1.0 - zfcrm**2.0 )
                !
                dqdx(ji,jj) = dqdx(ji,jj) + dhdx(ji,jj) * dx(jj) * vcorx(jj) / rpi

                dqdy(ji,jj) = dqdy(ji,jj) + dhdy(ji,jj) * dy(jj) * vcory(jj) / rpi
                !
                encon(ji,jj) = ( dhdx(ji,jj) * dqdx(ji,jj) * ( ureal(ji,jj)**2.0 + uimag(ji,jj)**2.0 ) &
                     &  + dhdy(ji,jj) * dqdy(ji,jj) * ( vreal(ji,jj)**2.0 + vimag(ji,jj)**2.0 ) &
                     &  + ( dhdx(ji,jj) * dqdy(ji,jj) + dhdy(ji,jj) * dqdx(ji,jj) ) &
                     &   * ( ureal(ji,jj) * vreal(ji,jj) + uimag(ji,jj) * vimag(ji,jj) ) ) * &
                     & 0.5_wp *rho0 * bfbt(ji,jj) * zprefac


                ftxx(ji,jj) = rho0 * bfbt(ji,jj) * zprefac * dhdx(ji,jj) * dqdx(ji,jj)       !!! Green and Nycander

                ftyy(ji,jj) = rho0 * bfbt(ji,jj) * zprefac * dhdy(ji,jj) * dqdy(ji,jj)       !!! Green and Nycander 

                ftxy(ji,jj) = 0.5_wp *rho0 * bfbt(ji,jj) * zprefac * ( dhdx(ji,jj) * dqdy(ji,jj) + dhdy(ji,jj) * dqdx(ji,jj) ) !!! Green and Nycander



             END IF
          END DO

       END DO

       !encon(:,:) = encon(:,:) !* tmask(:,:)

    END IF
    ! mdal computation 

     OPEN(6,FORM='FORMATTED',RECL=256)


    IF ( ln_modal ) THEN

       DO jj = is, iu

          PRINT *, jj, narea       !monitoring computation for each lattitude

          DO ji = iw, ie

             IF ( vcut(ji,jj) > 0.0 ) THEN

                !vcut(ji,jj) = ABS(vcut(ji,jj))

                iiw = MAX( 1, ji - nxinterval(ji,jj) )
                !                                                                                                                                                                             
                iie = MIN( nx, ji + nxinterval(ji,jj) )
                !                                                                                                                                                                             
                iis = MAX( 1, jj - nyinterval(ji,jj) )
                !                                                                                                                                                                              
                iiu = MIN( ny, jj + nyinterval(ji,jj) )
                !                                                                                                                                                                              

                DO jl = iis, iiu
                   DO jk = iiw, iie

                      zda = dx(jl) * dy(jl)

                      IF ( llbsl0 ) THEN
                         zxdif = dx(INT(( jj + jl ) / 2.0)) * REAL( ji - jk , wp)
                         zydif = dy(INT(( jj + jl ) / 2.0)) * REAL( jj - jl , wp)
                         !
                         zdist = SQRT( zxdif * zxdif + zydif * zydif)
                         zdist = MIN( zdist, rn_trunc / vcut(ji,jj) )

                         !PRINT *, 'distnace', rn_trunc / vcut(ji,jj)

                         dqdx(ji,jj) = dqdx(ji,jj) + zda * dhdx(jk,jl) *  &
                              &    besj0( REAL( zdist * vcut(ji,jj), wp ) ) 
                         dqdy(ji,jj) = dqdy(ji,jj) + zda * dhdy(jk,jl) *  &
                              &    besj0( REAL( zdist * vcut(ji,jj), wp )  )

                         !!PRINT *, besj0( REAL( zdist * vcut(ji,jj), wp ) )

                         !!PRINT *, zdist * vcut(ji,jj), jj, besj0( REAL( zdist * vcut(ji,jj), wp )  )


                      ELSE

                         zxdif = dx(INT(( jj + jl ) / 2.0)) * REAL( ( ji - jk + 0.5 ), wp)
                         zydif = dy(INT(( jj + jl ) / 2.0)) * REAL( ( jj - jl + 0.5 ), wp)
                         !
                         zdist = SQRT( zxdif * zxdif + zydif * zydif)
                         zdist = MIN( zdist, rn_trunc / vcut(ji,jj))

                         IF ( llbathy ) THEN  
                            dqdx(ji,jj) = dqdx(ji,jj) + zda * ( h(jk,jl) - h(ji,jj) )  * zxdif * ( vcut(ji,jj) ) &
                                 &   * ( - besj1( REAL(zdist * vcut(ji,jj), wp) ) ) / zdist  
                            dqdy(ji,jj) = dqdy(ji,jj) + zda * ( h(jk,jl) - h(ji,jj) )* zydif * ( vcut(ji,jj) ) &
                                 &   * ( - besj1( REAL(zdist * vcut(ji,jj), wp) ) ) / zdist
                         ELSE
                            dqdx(ji,jj) = dqdx(ji,jj) + zda * h(jk,jl) * zxdif * ( vcut(ji,jj) ) &
                                 &   * ( - besj1( REAL(zdist * vcut(ji,jj), wp) ) ) / zdist
                            dqdy(ji,jj) = dqdy(ji,jj) + zda * h(jk,jl)* zydif * ( vcut(ji,jj) ) &
                                 &   * ( - besj1( REAL(zdist * vcut(ji,jj), wp) ) ) / zdist
                         END IF

                      END IF


                   ENDDO
                ENDDO

                !dqdx(ji,jj) = dqdx(ji,jj) + dhdx(ji,jj) * dx(jj) * vcorx(jj) / rpi

                !dqdy(ji,jj) = dqdy(ji,jj) + dhdy(ji,jj) * dy(jj) * vcory(jj) / rpi
                !

                zfcrm = fcors(jj) / tfrq(ncu)
                zprefac = SQRT( 1.0 - zfcrm**2.0 ) !* fcors(jj) 

                !encon(ji,jj) = ( dhdx(ji,jj) * dqdx(ji,jj) * ( ureal(ji,jj)**2.0 + uimag(ji,jj)**2.0 ) &
                !     &  + dhdy(ji,jj) * dqdy(ji,jj) * ( vreal(ji,jj)**2.0 + vimag(ji,jj)**2.0 ) &
                !     &  + ( dhdx(ji,jj) * dqdy(ji,jj) + dhdy(ji,jj) * dqdx(ji,jj) ) &
                !     &   * ( ureal(ji,jj) * vreal(ji,jj) + uimag(ji,jj) * vimag(ji,jj) ) ) * &
                !     & 0.1250_wp *rho0 * vcut(ji,jj) * eps(ji,jj) * sn(ji,jj) * sn(ji,jj) * zprefac

                encon(ji,jj) = ( dhdx(ji,jj) * dqdx(ji,jj) * ( ureal(ji,jj)**2.0 + uimag(ji,jj)**2.0 ) &
                     &  + dhdy(ji,jj) * dqdy(ji,jj) * ( vreal(ji,jj)**2.0 + vimag(ji,jj)**2.0 ) &
                     &  + ( dhdx(ji,jj) * dqdy(ji,jj) + dhdy(ji,jj) * dqdx(ji,jj) ) &
                     &   * ( ureal(ji,jj) * vreal(ji,jj) + uimag(ji,jj) * vimag(ji,jj) ) ) & 
                     &   * 0.1250_wp *rho0 *  eps(ji,jj) * vcut(ji,jj) * sn(ji,jj) * zprefac !* sn(ji,jj)

                

                zin1 = ( dhdx(ji,jj) * dqdx(ji,jj) * ( ureal(ji,jj)**2.0 + uimag(ji,jj)**2.0 ))
                zin2 = dhdy(ji,jj) * dqdy(ji,jj) * ( vreal(ji,jj)**2.0 + vimag(ji,jj)**2.0 )
                zin3 =  ( dhdx(ji,jj) * dqdy(ji,jj) + dhdy(ji,jj) * dqdx(ji,jj) ) &
                     &   * ( ureal(ji,jj) * vreal(ji,jj) + uimag(ji,jj) * vimag(ji,jj) )

                zint = zin1 + zin2 + zin3 
 
                !!encon(ji,jj) = zint
                !WRITE(6,*) zin1, zin2, zin3, zint



             END IF

          END DO

       END DO

       !WHERE ( encon(:,:) < 0.0 ) encon(:,:) = -1.0e-3 
    END IF
    !
    CALL output_energy 

    
    

    
    !!WHERE ( encon(:,:) < 0.0 ) encon(:,:) = -9999.0 
    

    !WHERE ( encon(:,:) > 0.01 ) encon(:,:) = 0.0 
    

   !!IF ( lwp ) WRITE(numout,*) COUNT ( encon(:,:) /= -9999.0 ), 'positive'


   !!IF ( lwp ) WRITE(numout,*) COUNT ( encon(:,:) == -9999.0 ) , 'negative'


    
  !!IF ( lwp) WRITE(numout,*)  SIZE ( encon(:,:) ),               'size'


  IF ( lwp) WRITE(numout,*)  REAL ( REAL(COUNT ( encon(:,:) < 0.0 ),wp) / SIZE ( encon(:,:) ), wp ) * 100.0,"%",' negative points' 



    
  IF ( lwp ) WRITE(numout,*) REAL ( REAL(COUNT ( encon(:,:) >= 0.0 ),wp) / SIZE ( encon(:,:) ), wp ) * 100.0,"%",' positive points'
  
   
  PRINT *, SUM(encon(:,:) , mask = encon(:,:) >= 0.0 ) / 1.e6_wp
 
  PRINT *, SUM(encon(:,:) , mask = encon(:,:) < 0.0 ) / 1.e6_wp
    
  IF ( lwp ) WRITE(numout,*) REAL( SUM(encon(:,:) ) / SIZE ( encon(:,:) ) / 1.e6_wp, 8)
 


  !IF ( lwp ) WRITE(numout,*) REAL ( COUNT ( encon(:,:) < 0.0_wp ) / SIZE ( encon(:,:) ), wp ) * 100.0, '%'

  !IF ( lwp ) WRITE(numout,*) REAL ( COUNT ( encon(:,:) > 0.0_wp ) / SIZE ( encon(:,:) ), wp ) * 100.0, '%'


    PRINT *, ' The whole depth energy flux GW             ' ,  SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 350.0 ) / 1.e9_wp


    IF ( lwp ) WRITE(numout,*) ' The whole depth energy flux GW             ' ,  SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 0.0 ) / 1.e9_wp

    IF ( lwp ) WRITE(numout,*) ' The energy flux depth greater than  500 m GW' , SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 500.0 ) / 1.e9_wp

    IF ( lwp ) WRITE(numout,*) ' The energy flux depth greater than 1000 m GW',  SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 1000.0 ) / 1.e9_wp
    
    IF ( lwp ) WRITE(numout,*) ' The energy flux depth greater than 1500 m GW',  SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 1500.0 ) / 1.e9_wp
    
    IF ( lwp ) WRITE(numout,*) ' The energy flux depth greater than 2000 m GW',  SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 2000.0 ) / 1.e9_wp
    
    IF ( lwp ) WRITE(numout,*) ' The energy flux depth greater than 2500 m GW',  SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 2500.0 ) / 1.e9_wp
    !
    IF ( lwp ) WRITE(numout,*) 

    IF ( lwp ) WRITE(numout,*) ' the contribution of the supercritical topography to the energy conversion rate are as follows:'

    IF ( lwp ) WRITE(numout,*) 

    IF ( lwp ) WRITE(numout,*) ' The whole depth energy flux GW             ',    SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 0.0 .AND. vsteep(:,:) >= 1.0 ) / 1.e9_wp

    IF ( lwp ) WRITE(numout,*) ' The energy flux depth greater  than 500 m GW' ,  SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 500.0 .AND. vsteep(:,:) >= 1.0 ) / 1.e9_wp

    IF ( lwp ) WRITE(numout,*) ' The energy flux depth greater than 1000 m GW',   SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 1000.0 .AND. vsteep(:,:) >= 1.0 ) / 1.e9_wp
    
    IF ( lwp ) WRITE(numout,*) ' The energy flux depth greater than 1500 m GW',   SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 1500.0 .AND. vsteep(:,:) >= 1.0 ) / 1.e9_wp
    
    IF ( lwp ) WRITE(numout,*) ' The energy flux depth greater than 2000 m GW',   SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 2000.0 .AND. vsteep(:,:) >= 1.0 ) / 1.e9_wp
    
    IF ( lwp ) WRITE(numout,*) ' The energy flux depth greater than 2500 m GW',   SUM( encon(:,:)* DA(:,:), mask = h(:,:) >= 2500.0 .AND. vsteep(:,:) >= 1.0 ) / 1.e9_wp
    !
    
  END SUBROUTINE sol_energy


  SUBROUTINE sol_velocity
    !!----------------------------------------------------------------------------------------
    !!                       *** ROUTINE sol_velocity  ***
    !!
    !! ** Purpose : Calculation of barotropic plus baroclinic tidal velocity on
    !!              topography points
    !!
    !! ** Method  : The calculations are carried out by computing convolution integrals in
    !!              the main routine. The outer loop is over all field points, i.e. all
    !!              points on the topographic grid where velocity of the internal tide is 
    !!              wanted. For each field point there is an inner loop over source points
    !!              this inner loop is the convolution integral.
    !! 
    !!  References :
    !!-----------------------------------------------------------------------------------------
    INTEGER  :: ji, jj, jk, jl         ! dummy integers used for loops
    INTEGER  :: is, iu, iw, ie         ! local integer used for loops
    INTEGER  :: iis, iiu, iiw, iie     ! local integer used for loops
    REAL(wp) :: zdist
    REAL(wp) :: zda
    REAL(wp) :: zxdif, zydif
    REAL(wp) :: zfcrm, zprefac
    REAL(wp) :: zureal, zuimag
    REAL(wp) :: zvreal, zvimag
    !!--------------------------------------------------------------------
    !
    CALL name_outvel
    !

    IF ( ln_bt ) CALL output_bt

    iw = nn_bfz + 2
    ie = SIZE(h,1) - nn_bfz
    is = nn_bfz + 2
    iu = SIZE(h,2) - nn_bfz
    !

    DO jj = is, iu

       !!PRINT *, jj, narea         !monitoring computation for each lattitude

       IF ( lwp ) WRITE(numout,*)  jj, narea       !monitoring computation for each lattitude

       DO ji = iw, ie
          IF ( vcut(ji,jj) > 0.0 ) THEN
             !vcut(ji,jj) = ABS(vcut(ji,jj))
             !
             iiw = MAX( 1, ji - nxinterval(ji,jj) )
             !                                                                                                                                                                             
             iie = MIN( nx - 1, ji + nxinterval(ji,jj) )
             !                                                                                                                                                                             
             iis = MAX( 1, jj - nyinterval(ji,jj) )
             !                                                                                                                                                                              
             iiu = MIN( ny - 1, jj + nyinterval(ji,jj) )
             !                                                                                                                                                                              
             DO jl = iis, iiu
                DO jk = iiw, iie

                   zda   = dx(jl) * dy(jl)

                   zxdif = dx(INT(( jj + jl) / 2.0 )) * REAL( ji -jk - 0.5, wp )

                   zydif = dy(INT((jj + jl) / 2.0 )) * REAL(jj - jl, wp)

                   zdist = SQRT( zxdif**2.0 + zydif**2.0 )

                   zdist = MIN( zdist, rn_trunc * vcut(ji,jj))

                   d2qdx2(ji,jj) = d2qdx2(ji,jj) + zda * dhdxs(jk,jl) * zxdif * ( 1.0 / vcut(ji,jj))**2.0 &
                        &   * grendr(zdist / vcut(ji,jj)) / zdist
                   !
                   zxdif = dx(INT((jj + jl) / 2.0)) * REAL(ji - jk, wp)

                   zydif = dy(INT((jj + jl) / 2.0)) * REAL( jj - jl - 0.50, wp)

                   zdist = SQRT( zxdif**2.0 + zydif**2.0)

                   zdist = MIN( zdist, rn_trunc * vcut(ji,jj))

                   d2qdy2(ji,jj) = d2qdy2(ji,jj) + zda * dhdys(jk,jl) * zydif * (1.0 / vcut(ji,jj))**2.0 &
                        &   * grendr(zdist / vcut(ji,jj)) / zdist
                   !
                   zxdif = dx(INT((jj + jl) / 2.0 )) * REAL( ji - jk - 0.5, wp)

                   zydif = dy(INT((jj + jl) / 2.0 )) * REAL( jj - jl - 0.5, wp)

                   zdist = SQRT(zxdif**2.0 + zydif**2.0)

                   zdist = MIN( zdist, rn_trunc * vcut(ji,jj))

                   d2qdxdy(ji,jj) = d2qdxdy(ji,jj) + zda * ( 1.0 / vcut(ji,jj))**2.0 * grendr( zdist / vcut(ji,jj)) * &
                        &   0.5_wp * ( dhdy(jk,jl) * zxdif + dhdx(jk,jl) * zydif ) / zdist
                   !
                ENDDO
             ENDDO

             d2qdx2(ji,jj) = d2qdx2(ji,jj) + ( h(ji-1,jj) - 2.0_wp * h(ji,jj) + h(ji+1,jj) )* vcorx(jj) / rpi

             d2qdy2(ji,jj) = d2qdy2(ji,jj) + ( h(ji,jj-1)-2.0_wp * h(ji,jj) + h(ji,jj+1)) * vcory(jj) / rpi

             d2qdxdy(ji,jj) = d2qdxdy(ji,jj)+ &
                  &  ( h(ji-1,jj-1) + h(ji+1,jj+1) - h(ji-1,jj+1) - h(ji+1,jj-1)) * vcorxy(jj) / (2.0_wp * rpi)

             zfcrm = fcor(jj) / tfrq(ncu)

             zprefac = bfbt(ji,jj) / tfrq(ncu) / SQRT(1.0 - zfcrm**2.0)

             zureal = zprefac * ( d2qdx2(ji,jj) * uimag(ji,jj) + d2qdy2(ji,jj) * vreal(ji,jj) * zfcrm &
                  &   + d2qdxdy(ji,jj) * ( vimag(ji,jj) + ureal(ji,jj) * zfcrm) )

             zuimag = zprefac * ( -d2qdx2(ji,jj) * ureal(ji,jj) + d2qdy2(ji,jj) * vimag(ji,jj) * zfcrm &
                  &   + d2qdxdy(ji,jj) * ( -vreal(ji,jj) + uimag(ji,jj) * zfcrm))

             zvreal = zprefac * ( -d2qdx2(ji,jj) * ureal(ji,jj) * zfcrm + d2qdy2(ji,jj) * vimag(ji,jj) &
                  &   + d2qdxdy(ji,jj) * ( uimag(ji,jj) - vreal(ji,jj) * zfcrm))

             zvimag = zprefac * ( -d2qdx2(ji,jj) * uimag(ji,jj) * zfcrm - d2qdy2(ji,jj) * vreal(ji,jj) &
                  &   + d2qdxdy(ji,jj) * ( -ureal(ji,jj) - vimag(ji,jj) * zfcrm))

             ureal(ji,jj) = ureal(ji,jj) + zureal !* tmask(ji,jj) 
             uimag(ji,jj) = uimag(ji,jj) + zuimag !* tmask(ji,jj)
             vreal(ji,jj) = vreal(ji,jj) + zvreal !* tmask(ji,jj)  
             vimag(ji,jj) = vimag(ji,jj) + zvimag !* tmask(ji,jj)
             !
          END IF
       ENDDO
    ENDDO
    !
    CALL output_velocity
    !
  END SUBROUTINE sol_velocity


  FUNCTION dstnce( pla1, phi1, pla2, phi2 )
    !!--------------------------------------------------------------------------
    !!                   *** ROUTINE  dstnce ***  
    !!
    !!--------------------------------------------------------------------------
    REAL(sp), INTENT(in) ::   pla1, phi1, pla2, phi2
    REAL(wp)             ::   zla1, zhi1, zla2, zhi2
    REAL(wp)             ::   dly1, dly2, dlx1, dlx2, dlx, dls, dld, dpi
    REAL(wp)             ::   dstnce
    !!--------------------------------------------------------------------------
    !
    dpi  = 2.* ASIN(1.)
    dls  = dpi / 180.

    zla1  = REAL( pla1, wp) 
    zhi1  = REAL( phi1, wp)
    zla2  = REAL( pla2, wp)
    zhi2  = REAL( phi2, wp)

    dly1 = zhi1 * dls
    dly2 = zhi2 * dls
    dlx1 = zla1 * dls
    dlx2 = zla2 * dls

    dlx = SIN(dly1) * SIN(dly2) + COS(dly1) * COS(dly2) * COS(dlx2-dlx1)

    IF( ABS(dlx) > 1.0 ) dlx = 1.0
    dld = ATAN(DSQRT( ( 1-dlx )/( 1+dlx ) )) * 222.24 / dls
    dstnce = dld * 1000.
    !
  END FUNCTION dstnce

  !!======================================================================================
END MODULE soltmx
