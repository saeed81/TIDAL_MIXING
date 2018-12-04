MODULE lib_mpp

  IMPLICIT NONE

INCLUDE 'mpif.h'

CONTAINS


  SUBROUTINE mpp_init

    CALL mpi_comm_dup( mpi_comm_world, mpi_comm_tmx, code)
    IF( code /= MPI_SUCCESS ) THEN
       WRITE(*, *)'error'
       WRITE(*, *) ' lib_mpp: Error in routine mpi_comm_dup'
       CALL mpi_abort( mpi_comm_world, code, ierr )
    ENDIF
    
    CALL mpi_comm_rank( mpi_comm_tmx, mpprank, ierr )
    CALL mpi_comm_size( mpi_comm_tmx, mppsize, ierr )
    




  SUBROUTINE mppsync()
    !!----------------------------------------------------------------------
    !!                  ***  routine mppsync  ***
    !!
    !! ** Purpose :   Massively parallel processors, synchroneous
    !!
    !!-----------------------------------------------------------------------
    INTEGER :: ierror
    !!-----------------------------------------------------------------------
   !
    CALL mpi_barrier( mpi_comm_tmx, ierror )
    !
  END SUBROUTINE mppsync


  SUBROUTINE mppstop
    !!----------------------------------------------------------------------
    !!                  ***  routine mppstop  ***
    !!
    !! ** purpose :   Stop massilively parallel processors method
    !!
    !!----------------------------------------------------------------------
    INTEGER ::   info
    !!----------------------------------------------------------------------
    !
    CALL mppsync
    CALL mpi_finalize( info )
    !
  END SUBROUTINE mppstop


END MODULE lib_mpp
