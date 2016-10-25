!=======================================================================
!                        Version 1: 02/99, MRZ
!-----------------------------------------------------------------------
! MPI_param
!   This class wraps the MPI parameters obtained from the system-
! dependent include file so that it can be used by both free- and
! fixed-format Fortran90.
!-----------------------------------------------------------------------
! v1.0: Original implementation
!=======================================================================

module mpi_param_mod
                                                                                                     
  use kind_mod

  public

! MPI is enabled: include the system-dependent MPI include file

#include 'mpif.h

! MPI is enabled: include the system-dependent MPI include file

end module mpi_param_mod
