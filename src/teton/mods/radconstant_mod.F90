! Radiation Constant Module:  commonly used radiation constants.
! Updated 1/18/18 by TAB to use PhysicsUtils constants.
!
module radconstant_mod

use kind_mod
use constant_mod, only: one, three, eight

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
use PhysicsUtilsConstants_ADiv, only: speed_of_light, &
                                      black_body_constant, &
                                      Thomson_cross_section, &
                                      energy_per_keV, &
                                      classical_electron_radius, &
                                      Planck_constant
use PhysicsUtilsConstants, only: electron_mass_energy_equivalent_in_keV, &
                                 pi
#else
use constant_mod, only: pi
#endif


private

  real(adqt), parameter, public ::                                       &
#if defined(TETON_ENABLE_MINIAPP_BUILD)
                  speed_light         = 299.792458_adqt,                 &
                  rad_constant        = 0.013720169264801055_adqt,       &
                  elecRestEnergy      = 510.9989499985809_adqt,          &
                  recipElecRestEnergy = one/510.9989499985809_adqt,      &
                  electronSigTh       = 6.652458732110024e-25_adqt,      &
                  hPlanck             = 1.0545718176461566e-35_adqt,     &
                  KeVToEnergy         = 1.602176634e-25_adqt,            &
                  classical_electron_radius = 2.8179403261914806e-13
#else
                  speed_light         = speed_of_light,                             &
                  rad_constant        = black_body_constant,                        &
                  elecRestEnergy      = electron_mass_energy_equivalent_in_keV,     &
                  recipElecRestEnergy = one/electron_mass_energy_equivalent_in_keV, &
                  electronSigTh       = Thomson_cross_section,                      &
                  hPlanck             = Planck_constant,                            &
                  KeVToEnergy         = energy_per_keV
#endif

  real(adqt), parameter, public ::                                 &
                  hC             = (hPlanck/KeVToEnergy)*          &
                                   speed_light,                    &
                  PlanckFactor   = eight*pi/(hC*hC*hC),            &
                  facForE        = PlanckFactor*KeVToEnergy*       &
                                   elecRestEnergy*elecRestEnergy*  &
                                   elecRestEnergy*elecRestEnergy,  &
                  electronRadius = classical_electron_radius

! TODO: Make these parameters tunable by the host codes? See TETON-131.
  real(adqt), parameter, public ::               &
                      sigmaCeiling = 1.e50_adqt, &
                      sigmaWarn    = 1.e40_adqt

end module radconstant_mod
