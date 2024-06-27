! Constant Module:  various commonly used constants.

module constant_mod

use kind_mod

private

  real(adqt), parameter, public ::                   &
                      zero     = 0.0_adqt,           &
                      eleventh = 1.0_adqt/11.0_adqt, &
                      ninth    = 1.0_adqt/9.0_adqt,  &
                      seventh  = 1.0_adqt/7.0_adqt,  &
                      sixth    = 1.0_adqt/6.0_adqt,  &
                      fifth    = 1.0_adqt/5.0_adqt,  &
                      fourth   = 0.25_adqt,          &
                      third    = 1.0_adqt/3.0_adqt,  &
                      half     = 0.5_adqt,   &
                      one      = 1.0_adqt,   &
                      two      = 2.0_adqt,   &
                      e        = 2.718281828459045235360287471352662_adqt, &
                      three    = 3.0_adqt,   &
                      pi       = 3.14159265358979323846264338327950_adqt,  &
                      four     = 4.0_adqt,   &
                      five     = 5.0_adqt,   &
                      six      = 6.0_adqt,   &
                      seven    = 7.0_adqt,   &
                      eight    = 8.0_adqt,   &
                      nine     = 9.0_adqt,   & 
                      ten      = 10.0_adqt,  &
                      eleven   = 11.0_adqt,  &
                      twenty   = 20.0_adqt,  &
                      sixty    = 60.0_adqt,  &
                      hundred  = 100.0_adqt, &
! The following are used in PWLD transport methods
                      one_3    = 1.0_adqt/3.0_adqt,  &
                      one_4    = 0.25_adqt,          &
                      one_6    = 1.0_adqt/6.0_adqt,  &
                      one_12   = 1.0_adqt/12.0_adqt, &
                      one_24   = 1.0_adqt/24.0_adqt, &
                      one_30   = 1.0_adqt/30.0_adqt, &
                      one_60   = 1.0_adqt/60.0_adqt, &
                      one_120  = 1.0_adqt/120.0_adqt


  real(adqt), parameter, public ::                     &
                      adqtTiny    = tiny(0.0_adqt),    &
                      adqtEpsilon = epsilon(0.0_adqt), &
                      adqtSmall   = 1.e-150_adqt

end module constant_mod
