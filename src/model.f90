PROGRAM model
  !!--------------------------------------------------------
  !!                    *** PROGRAM model ***
  !!
  !! ** Purpose : calcuate energy and velocity field for 
  !!              internal tidal mixing in the ocean
  !!
  !! ** Method  : -use tmx module and call tmx_model routine
  !!
  !!--------------------------------------------------------   
  USE tmx          !TMX system (tmx_model routine)
  !!--------------------------------------------------------
  !
  CALL tmx_model

END PROGRAM model
