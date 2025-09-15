! Default iteration control values.

module default_iter_controls_mod

use iso_c_binding

private

! Iteration control defaults.
  real(C_DOUBLE), parameter, public :: &
     outer_temp_reltol = 1.0e-5_C_DOUBLE, &
     outer_intensity_reltol = 1.0e-5_C_DOUBLE, &
     grey_reltol = 1.0e-6_C_DOUBLE, &
     incident_flux_reltol = 5.0e-6_C_DOUBLE, &
     inner_nl_reltol = 1.0e-7_C_DOUBLE

  integer(C_INT), parameter, public :: &
     outer_max_it = 50_C_INT, &
     grey_max_sweeps = 21_C_INT, &
     incident_flux_max_it = 2_C_INT, &
     inner_nl_max_it = 100_C_INT

end module default_iter_controls_mod
