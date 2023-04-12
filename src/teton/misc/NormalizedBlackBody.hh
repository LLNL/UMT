#ifndef NORMALIZEDBLACKBODY_HH__
#define NORMALIZEDBLACKBODY_HH__

extern "C"
{
/// \file NormalizedBlackBody.hh
/// \brief Integrate functions of the normalized black body function over photon
/// energies.

/// \defgroup NormalizedBlackBody Normalized Black Body
/// \ingroup BlackBody
/// \brief Stateless functions to compute normalized black body emission as a
/// function of group bounds and material temperature
///
/// Compute the angle and energy integrated black body function using
/// the methods in this paper:
///
/// - Clark, Bradley A., "Computing Multigroup Radiation Integrals Using
///        Polylogarithm-Based Methods," *Journal of Compuational Physics*,
///        70, p. 311-329, 1987.
///
///
/// @{

/// Compute normalized emission for all groups at once with raw pointer
/// interface.
void NBB_integrateBlackBodyGroups(const double T,
                                  const double k,
                                  const double Bnorm,
                                  const int numGroups,
                                  const double *const __restrict__ groupBounds,
                                  double *__restrict__ const B);

/// Compute normalized emission and derivatives for all groups at once with raw
/// pointer interface
void NBB_integrateBlackBodyAndDerivGroups(const double T,
                                          const double k,
                                          const double Bnorm,
                                          const int numGroups,
                                          const double *__restrict__ const groupBounds,
                                          double *__restrict__ const B,
                                          double *__restrict__ const dBdT);

/// @}
}
#endif // NORMALIZEDBLACKBODY_H__
