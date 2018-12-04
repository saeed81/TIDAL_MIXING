MODULE par_tide
  !!====================================================================================
  !!                       *** MODULE par_tide ***
  !! TMX system parameter : definition of parameters and variables
  !!====================================================================================

  IMPLICIT NONE
  PRIVATE


  !!------------------------------------------------------------------------------------
  !!                       define the kind of real for the whole model
  !!------------------------------------------------------------------------------------
    INTEGER, PUBLIC, PARAMETER ::  &  !: Floating point section
       &   sp = 4,  &                 !: single precision (real 4)
       &   dp = 8,  &                 !: double precision (real 8)
       &   wp = dp                    !: working precision


  !!-----------------------------------------------------------------------------------
  !!                             domain size
  !!-----------------------------------------------------------------------------------
  INTEGER, PUBLIC             :: nx = 1000     !: size of domain in x direction 
  INTEGER, PUBLIC             :: ny = 1000     !: size of domain in y direction 
  INTEGER, PUBLIC             :: nz = 1        !: size of domain in z direction 
  INTEGER, PUBLIC             :: nt = 1        !: size of domain in t direction 

   
  !!-------------------------------------------------------------------------------------
  !!                              names of inputfile
  !!-------------------------------------------------------------------------------------
  CHARACTER(len=300) , PUBLIC  :: cfile_h     = ''   !:    
  CHARACTER(len=20 ) , PUBLIC  :: cvar_h      = ''   !:
  CHARACTER(len=300) , PUBLIC  :: cfile_bf    = ''   !:
  CHARACTER(len=20 ) , PUBLIC  :: cvar_bf     = ''   !:
  CHARACTER(len=300) , PUBLIC  :: cfile_intbf = ''   !:
  CHARACTER(len=20 ) , PUBLIC  :: cvar_intbf  = ''   !:
  CHARACTER(len=300) , PUBLIC  :: cfile_speed = ''   !:    
  CHARACTER(len=20 ) , PUBLIC  :: cvar_speed  = ''   !:
  CHARACTER(len=300) , PUBLIC  :: cfile_sn    = ''   !:    
  CHARACTER(len=20 ) , PUBLIC  :: cvar_sn     = ''   !:
  CHARACTER(len=300) , PUBLIC  :: cfile_eps   = ''   !:    
  CHARACTER(len=20 ) , PUBLIC  :: cvar_eps    = ''   !:
  CHARACTER(len=20 ) , PUBLIC  :: cvarlon     = ''   !:
  CHARACTER(len=20 ) , PUBLIC  :: cvarlat     = ''   !:
  
   
  !!-------------------------------------------------------------------------------------
  !!                                 directory of input and output
  !!--------------------------------------------------------------------------------------
  CHARACTER(len=300), PUBLIC  :: cn_dir = './' !: directory of input or output  files    


  !!-------------------------------------------------------------------------------------
  !!                          namgrn namelist parameter
  !!-------------------------------------------------------------------------------------
  INTEGER,  PUBLIC       :: nn_gr   =  0  !: namelist variable for choosing method
  REAL(wp), PUBLIC       :: rn_beta = 1.0 !:   

  
  !!-------------------------------------------------------------------------------------
  !!                          namtide namelist parameter
  !!-------------------------------------------------------------------------------------
  CHARACTER(len=2), PUBLIC   :: cn_cu  = 's2'  !:
  INTEGER,          PUBLIC   :: ncu    =  2    !:


  !!-------------------------------------------------------------------------------------
  !!                          namsol namelist parameter
  !!-------------------------------------------------------------------------------------
  LOGICAL, PUBLIC :: ln_modal    = .FALSE.    !:
  LOGICAL, PUBLIC :: ln_energy   = .FALSE.    !:
  LOGICAL, PUBLIC :: ln_velocity = .FALSE.    !:
  LOGICAL, PUBLIC :: ln_bt       = .FALSE.    !:
  LOGICAL, PUBLIC :: ln_cutoff   = .FALSE.    !:
  LOGICAL, PUBLIC :: ln_steep    = .FALSE.    !:
  LOGICAL, PUBLIC :: ln_hslope   = .FALSE.    !:
  INTEGER, PUBLIC :: nn_bfz      =  30        !:
  INTEGER, PUBLIC :: nn_corr     =  0         !:
  INTEGER, PUBLIC :: nn_mode     =  0         !:


  !!------------------------------------------------------------------------------------
  !!                      truncation variable for cutoff length
  !!------------------------------------------------------------------------------------
  REAL(wp), PUBLIC :: rn_trunc = 5.0_wp          !:
  

  !!------------------------------------------------------------------------------------
  !!                          Run control 
  !!------------------------------------------------------------------------------------
  LOGICAL, PUBLIC :: lwp = .FALSE.            !:


  !!------------------------------------------------------------------------------------
  !!                         area parameters
  !!------------------------------------------------------------------------------------
  CHARACTER(len=8), PUBLIC  :: carea  =  ''  !:
  CHARACTER(len=4), PUBLIC  :: cxarea =  ''  !:
  CHARACTER(len=4), PUBLIC  :: cyarea =  ''  !:
  INTEGER         , PUBLIC  :: narea  =  0   !:


  !!------------------------------------------------------------------------------------
  !!                           Logical units 
  !!------------------------------------------------------------------------------------
  INTEGER, PUBLIC           :: numout = -1  !: unit for tide.output file opened by processor one
  INTEGER, PUBLIC           :: numnam = -1  !: unit for namelist file
  INTEGER, PUBLIC           :: numgr  = -1  !: unit for green function file          
  INTEGER, PUBLIC           :: numdgr = -1  !: unit for dgreen function file

  !!====================================================================================
END MODULE par_tide

  

