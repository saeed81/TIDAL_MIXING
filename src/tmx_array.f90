MODULE tmx_array
  !!===================================================================================
  !!                           *** MODULE tmx_array ***
  !! Array definition    : definition of array for tmx system 
  !!===================================================================================

  !!-----------------------------------------------------------------------------------
  !!  tmx_alloc    :  allocation of all the arrays in tmx system
  !!  tmx_dealloc  :  free up the memory using deallocation
  !!-----------------------------------------------------------------------------------  
  USE par_tide     ! TMX system paramaters and variables
  USE arraydef     ! definition of TMX array   
 
  IMPLICIT NONE  
  PRIVATE

  PUBLIC :: tmx_alloc           ! called by domini.f90 
  PUBLIC :: tmx_dealloc         ! called by tmx.f90

CONTAINS 

  SUBROUTINE tmx_alloc( kx, ky )
    !!------------------------------------------------------------
    !!                    *** ROUTINE tmx_alloc ***
    !! ** Purpose : allocation of arrays used in TMX system 
    !!-------------------------------------------------------------
    INTEGER, INTENT(in)    :: kx
    INTEGER, INTENT(in)    :: ky
    INTEGER                :: istat
    INTEGER, DIMENSION(80) :: istalloc
    !----------------------------------------------------------
    !
    istat = 0
    istalloc(:) = 0

    ALLOCATE ( ureal(kx,ky), uimag(kx,ky), vreal(kx,ky), vimag(kx,ky), stat = istat ); istalloc(1) = istat
    
    ALLOCATE ( h(kx,ky), tmask(kx,ky), stat = istat ); istalloc(2) = istat
    
    ALLOCATE ( nxinterval(kx,ky), nyinterval(kx,ky), vcut(kx,ky), stat = istat ); istalloc(3) = istat

    ALLOCATE ( dhdx(kx,ky), dhdy(kx,ky), stat = istat ); istalloc(4) = istat

    ALLOCATE ( dhdxs(kx,ky), dhdys(kx,ky), stat = istat ); istalloc(5) = istat

    ALLOCATE ( d2qdx2(kx,ky), d2qdy2(kx,ky), d2qdxdy(kx,ky), stat = istat ); istalloc(6) = istat

    ALLOCATE ( dqdx(kx,ky), dqdy(kx,ky), encon(kx,ky), stat = istat ); istalloc(7) = istat

    ALLOCATE ( dx(ky), dy(ky), stat = istat ); istalloc(8) = istat  

    ALLOCATE ( vlon(kx), vlat(ky), stat = istat ); istalloc(9) = istat

    ALLOCATE ( vlons(kx), vlats(ky), stat = istat ); istalloc(10) = istat

    ALLOCATE ( fcor(ky), fcors(ky), stat = istat ); istalloc(11) = istat  

    ALLOCATE ( asx(ky), asy(ky), stat = istat ); istalloc(12) = istat
    
    ALLOCATE ( vcorx(ky), vcory(ky), vcorxy(ky), stat = istat ); istalloc(13) = istat 
    
    IF ( .NOT. ln_modal ) ALLOCATE ( bfint(kx,ky), bfbt(kx,ky), stat = istat )
    IF ( .NOT. ln_modal ) istalloc(14) = istat
    
    IF (    ln_modal    ) ALLOCATE ( speed(kx,ky), sn(kx,ky), eps(kx,ky), stat = istat )
    IF (    ln_modal    ) istalloc(14) = istat

    ALLOCATE ( DA(kx,ky), stat = istat ); istalloc(15) = istat 
    
    ALLOCATE ( vsteep(kx,ky), stat = istat ); istalloc(16) = istat 
    
    ALLOCATE ( vhslope(kx,ky), stat = istat ); istalloc(17) = istat 
    
    ALLOCATE ( ftxx(kx,ky), stat = istat ); istalloc(18) = istat 
    
    ALLOCATE ( ftxy(kx,ky), stat = istat ); istalloc(19) = istat 
    
    ALLOCATE ( ftyy(kx,ky), stat = istat ); istalloc(20) = istat 
    
    
    ALLOCATE (vtime(1)); vtime = 0.0

    IF( ANY( istalloc(:) /= 0) ) THEN
       IF(lwp) WRITE(numout,*) 'memory storage problem'
       STOP  'memory problem'
    ELSE
       IF(lwp)  WRITE(numout,*)
       IF(lwp)  WRITE(numout,*) 'array allocation was successful'
       IF(lwp)  WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       IF(lwp)  WRITE(numout,*)
    END IF
    !
  END SUBROUTINE tmx_alloc


  SUBROUTINE tmx_dealloc
    !!-----------------------------------------------------------------------
    !!                     *** ROUTINE tmx_dealloc ***
    !!
    !! ** Purpose :  free up the memory
    !!-----------------------------------------------------------------------
    !
    IF( ALLOCATED(  ureal   ))  DEALLOCATE(    ureal    ) 
    IF( ALLOCATED(  uimag   ))  DEALLOCATE(    uimag    ) 
    IF( ALLOCATED(  vreal   ))  DEALLOCATE(    vreal    ) 
    IF( ALLOCATED(  uimag   ))  DEALLOCATE(    vimag    ) 
    IF( ALLOCATED(    h     ))  DEALLOCATE(     h       ) 
    IF( ALLOCATED(  bfint   ))  DEALLOCATE(    bfint    ) 
    IF( ALLOCATED(  bfbt    ))  DEALLOCATE(    bfbt     ) 
    IF( ALLOCATED(  tmask   ))  DEALLOCATE(    tmask    ) 
    IF( ALLOCATED(nxinterval))  DEALLOCATE( nxinterval  ) 
    IF( ALLOCATED(nyinterval))  DEALLOCATE( nyinterval  ) 
    IF( ALLOCATED(  dhdx    ))  DEALLOCATE(    dhdx     ) 
    IF( ALLOCATED(  dhdy    ))  DEALLOCATE(    dhdy     ) 
    IF( ALLOCATED(  dhdxs   ))  DEALLOCATE(    dhdxs    ) 
    IF( ALLOCATED(  dhdys   ))  DEALLOCATE(    dhdys    ) 
    IF( ALLOCATED(  vcut    ))  DEALLOCATE(    vcut     ) 
    IF( ALLOCATED(  d2qdx2  ))  DEALLOCATE(    d2qdx2   ) 
    IF( ALLOCATED(  d2qdy2  ))  DEALLOCATE(    d2qdy2   ) 
    IF( ALLOCATED(  d2qdxdy ))  DEALLOCATE(    d2qdxdy  )
    IF( ALLOCATED(  dqdx    ))  DEALLOCATE(    dqdx     ) 
    IF( ALLOCATED(  dqdy    ))  DEALLOCATE(    dqdy     ) 
    IF( ALLOCATED(  encon   ))  DEALLOCATE(    encon    ) 
    IF( ALLOCATED(  vlon    ))  DEALLOCATE(    vlon     ) 
    IF( ALLOCATED(  vlat    ))  DEALLOCATE(    vlat     ) 
    IF( ALLOCATED(  vlons   ))  DEALLOCATE(    vlons    ) 
    IF( ALLOCATED(  vlats   ))  DEALLOCATE(    vlats    ) 
    IF( ALLOCATED(  asx     ))  DEALLOCATE(    asx      ) 
    IF( ALLOCATED(  asy     ))  DEALLOCATE(    asy      ) 
    IF( ALLOCATED(  fcor    ))  DEALLOCATE(    fcor     ) 
    IF( ALLOCATED(  fcors   ))  DEALLOCATE(    fcors    ) 
    IF( ALLOCATED(  vcorx   ))  DEALLOCATE(    vcorx    )  
    IF( ALLOCATED(  vcory   ))  DEALLOCATE(    vcory    ) 
    IF( ALLOCATED(  vcorxy  ))  DEALLOCATE(    vcorxy   ) 
    IF( ALLOCATED(  speed   ))  DEALLOCATE(    speed    ) 
    IF( ALLOCATED(  eps     ))  DEALLOCATE(    eps      ) 
    IF( ALLOCATED(  sn      ))  DEALLOCATE(    sn       ) 
    
    !
    IF(lwp) WRITE(numout,*) 
    IF(lwp) WRITE(numout,*) 'freeing up the memory was successfuly done'
    IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    IF(lwp) WRITE(numout,*)
    !
  END SUBROUTINE tmx_dealloc

  !!======================================================================================
END MODULE tmx_array
