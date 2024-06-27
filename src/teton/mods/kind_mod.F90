! Kind Module:  define the kind required for the requested precision;
!               "adqt" is the default "adequate" precision.
!
!               Use the iso_c_binding C types in order to maximize compatibility
!               with C, as we are sharing array pointers in conduit nodes.
!
!    adqt:  default ("adequate") precision

module kind_mod
use iso_c_binding

private

  integer, parameter, public ::  adqt = C_DOUBLE

  integer, parameter, public ::  float = C_FLOAT
  integer, parameter, public ::  double = C_DOUBLE

end module kind_mod
