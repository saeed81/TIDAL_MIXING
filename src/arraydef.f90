MODULE arraydef
  !!===============================================================================
  !!                   *** MODULE arraydef***  
  !! Array Declaration : Defining of arrays for tmx system
  !!===============================================================================
  USE par_tide   ! type of variables
  
  IMPLICIT NONE
  PRIVATE

  REAL(sp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: ureal        !:                
  REAL(sp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: uimag        !:
  REAL(sp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: vreal        !:
  REAL(sp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: vimag        !:
  REAL(sp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: h            !:
  REAL(sp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: bfbt         !: 
  REAL(sp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: bfint        !:
  REAL(sp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: speed        !: 
  REAL(sp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: sn           !: 
  REAL(sp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: eps          !:
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: tmask        !:
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: vcut         !:
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: vsteep       !:
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: vhslope      !:
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: dhdx         !: 
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: dhdy         !: 
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: dhdxs        !:  
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: dhdys        !:  
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: d2qdx2       !:   
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: d2qdy2       !:   
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: d2qdxdy      !:   
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: dqdx         !:
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: dqdy         !:
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: encon        !: 
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: ftxx         !: 
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: ftxy         !: 
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: ftyy         !: 
  REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: DA           !:
  INTEGER , PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: nxinterval   !:    
  INTEGER , PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: nyinterval   !: 
  !!
  REAL(sp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: vlon         !:
  REAL(sp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: vlat         !:
  REAL(sp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: vlats        !:
  REAL(sp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: vlons        !:
  REAL(wp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: dx           !:
  REAL(wp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: dy           !:
  REAL(wp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: fcor         !:
  REAL(wp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: fcors        !:
  REAL(wp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: asx          !:
  REAL(wp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: asy          !:
  REAL(wp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: vcorx        !:
  REAL(wp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: vcory        !:
  REAL(wp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: vcorxy       !:
  REAL(wp), PUBLIC, DIMENSION(:),   ALLOCATABLE, SAVE :: vtime        !:

  !!==================================================================================
END MODULE arraydef
