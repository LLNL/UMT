#include <cmath>

////////////////////////////////////////////////////////////////////////////////

/// \addtogroup NormalizedBlackBody
/// @{

/// Returns the constant \f$ \pi^4/15 \f$
constexpr double NBB_pi4OverFifteen()
{
   return 6.493939402266829149096022L;
}
/// Returns the constant \f$ 15/\pi^4 \f$
constexpr double NBB_fifteenOverPi4()
{
   return 0.1539897338202650278372917L;
}

extern "C"
{
////////////////////////////////////////////////////////////////////////////////

/// This function evaluates the black body function
/// using a rational polynomial defined as Eq. 48
/// in Clark, JCP vol 70. no 2, p 311 for epsilon > 0.01 and Eq. 32
/// for small epsilon.
///
double NBB_cumulativeEmission(const double epsilon)
{
   if (epsilon < 0.01)
   {
      //    use small x approximation: Clark, JCP 70, no. 2 1987, p 316, eq 32
      const double a1 = 0.05132991127342032; // 1/3/(pi**4/15)
      const double a2 = 0.01924871672753262; // 1/8/(pi**4/15)

      return epsilon * epsilon * epsilon * (a1 - a2 * epsilon);
   }
   else
   {
      //    These coefs are 15/pi**4 times the ones in the pli function.
      //    Clark, JCP 70, no. 2 1987, p 316, eq 48
      const double a1 = 1.2807339766120354;
      const double a2 = 0.8578722513311724;
      const double a3 = 0.33288908614428098;
      const double a4 = 0.079984931563508915;
      const double a5 = 0.011878558806454416;

      const double b1 = 0.2807339758744;
      const double b2 = 0.07713864107538;

      double Int_epsilon_Inf = std::exp(-epsilon)
                               * (1.0
                                  + epsilon * (a1 + epsilon * (a2 + epsilon * (a3 + epsilon * (a4 + epsilon * a5)))))
                               / (1. + epsilon * (b1 + epsilon * b2));

      return 1.0 - Int_epsilon_Inf;
   }
}

////////////////////////////////////////////////////////////////////////////////

/// When computing the derivative of the black body function with respect
/// to temperature, we need to also take derivatives of the group bounds
/// that have been normalized by temperature.
///
/// We multiply the standard form by \f$ e^{-\epsilon}/e^{-\epsilon} \f$
/// to get a more stable numerical form, since \f$\epsilon\f$ can be large,
/// we'd rather work with small denormalized numbers than floating point
/// infinity.
///
/// \param[in] epsilon photon energy normalized by temperature
/// \return extra term in the derivative of the normalized black body function
/// due to the bounds
double NBB_cumulativeEmissionDerivPart(const double epsilon)
{
   const double eps2 = epsilon * epsilon;
   if (epsilon <= 0.0)
   {
      return 0.0;
   }
   else if (epsilon < 1.0e-3)
   {
      // For small values of epsilon, the formula below is inaccurate (and singular)
      return 0.0384974334550662569593229372517 * epsilon * eps2 - 0.0192487167275331284796614686258 * eps2 * eps2
             + 0.00320811945458885474661024477097 * eps2 * eps2 * epsilon
             - 0.0000534686575764809124435040795162 * eps2 * eps2 * eps2 * epsilon;
   }
   else // if (epsilon < 46.435)
   {
      const double z = std::exp(-epsilon);
      const double numerator = 0.25 * z * eps2 * eps2;
      const double denominator = 1.0 - z;
      return NBB_fifteenOverPi4() * numerator / denominator;
   }
   /*
  else
  {
    // Epsilon in NBB_cumulativeEmission where it is always 1.
    return 0.0;
  }
  */
}

////////////////////////////////////////////////////////////////////////////////

/// This function computes the black body emission in photon energy
/// groups for a given temperature.  Photon energy space is broken
/// into numGroups groups given by boundaries groupBounds.  The lowest
/// group bound (groupBounds[0]) is ignored and set assumed to be zero.  The
/// highest group bound (groupBounds[numGroups+1]) is ignored and assumed to be
/// infinity. For each group, we compute
/// \f[
///     B_g(T) = B_{\mathrm{norm}} \int_{E_g}^{E_{g+1}} \frac{E^3}{e^{E/k T} -1 } d E
/// \f]
/// where \f$ E_g \f$ are the group bounds, \f$ E = h \nu\f$ is the photon
/// energy,
/// \f$ T \f$ is the material temperature, and \f$k\f$ is Boltzmann's constant.
/// \f$ B_{\mathrm{norm}} \f$ is a normalization of the black body function. For
/// for methods that compute directional intensities, such as Sn, we could have
/// \f[
///       B_{\mathrm{norm}} = \frac{2\pi^4}{15 h^3 c^2}
/// \f]
/// For angle-integrated methods in terms of energy density, such as diffusion,
/// we have
/// \f[
///       B_{\mathrm{norm,diffusion}} = \frac{4\pi}{c} B_{\mathrm{norm}} =
///       \frac{8\pi^5}{15 h^3 c^3}
/// \f]
/// The integral at the upper and lower bounds is computed using
/// NBB_cumulativeEmission() using
/// \f[
///     B_g(T) = B_{\mathrm{norm}} \int_{E_g}^{E_{g+1}} \frac{E^3}{e^{E/k T} -1 } d E
///            = B_{\mathrm{norm}} (k T)^4 \left[ \frac{15}{\pi^4}
///            \int_{\epsilon_g}^{\epsilon_{g+1}} \frac{\epsilon^3}{e^{\epsilon}
///            -1 } d \epsilon \right],
/// \f]
/// where \f$ \epsilon = E / k T \f$.  (We could eliminate the multiply and
/// divide by
/// \f$\pi^4/15\f$, but then NBB_cumulativeEmission() would not be between zero and
/// one.)
///
/// \param[in] T Temperature of the material.
/// \param[in] k Boltzmann constant such that k*T is in the same units as the group boundaries.
/// \param[in] Bnorm is the normalization to get the intensity or energy density
///            in your units such that Bnorm*T^4 = sum_g B[g] = aT^4 (or multiplied by c for intensit)
/// \param[in] numGroups Number of groups
/// \param[in] groupBounds Group boundaries, must be numGroups+1 long.  First
///            and last ignored and assumed to be zero or infinity
/// \param[out] B Pre-allocated array, numGroups long, where group-wise emission
///            is stored.
void NBB_integrateBlackBodyGroups(const double T,
                                  const double k,
                                  const double Bnorm,
                                  const int numGroups,
                                  const double *__restrict__ const groupBounds,
                                  double *__restrict__ const B)
{
   const double kT = k * T;
   // Should be rare, if not, never, that we have a zero temperature,
   // but we're being safe.
   if (kT > 0.0)
   {
      const double Tsq = T * T;
      const double fullBNorm = Bnorm * Tsq * Tsq;

      double lowB = 0.0;
      int g;
      for (g = 0; g < numGroups - 1; ++g)
      {
         double highBound = groupBounds[g + 1] / kT;
         double highB = NBB_cumulativeEmission(highBound);
         B[g] = fullBNorm * (highB - lowB);
         lowB = highB;
      }
      B[numGroups - 1] = fullBNorm * (1.0 - lowB);
   }
   else
   {
      int g;
      for (g = 0; g < numGroups; ++g)
      {
         B[g] = 0.0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

/// In addition to computing everything that NBB_integrateBlackBodyGroups()
/// computes, this function also returns the derivative with respect to
/// temperature the emission in each group.
///
/// This uses the chain-rule-like formula
/// \f[
///       \frac{\partial}{\partial x} \int_{g(x)}^{h(x)} f(x) dx
///      = \int_{g(x)}^{h(x)} \frac{\partial}{\partial x} f(x) dx
///      + \frac{\partial h(x)}{\partial x} f(h(x))
///      - \frac{\partial g(x)}{\partial x} f(g(x))
/// \f]
/// to evaluate derivatives of the temperature normalized form of the integral.
/// Using the normalized emission form, we have
/// \f[
///     \frac{\partial B_g(T)}{\partial T}
///            = B_{\mathrm{norm}} (k T)^3 \left[4 \frac{15}{\pi^4}\right]
///          \left[ \int_{\epsilon_g}^{\epsilon_{g+1}}
///          \frac{\epsilon^3}{e^{\epsilon} -1 } d \epsilon
///             -\frac{1}{4}\frac{\epsilon_{g+1}^3}{e^{\epsilon_{g+1}} -1 }
///             +\frac{1}{4}\frac{\epsilon_{g}^3}{e^{\epsilon_{g}} -1 } \right]
/// \f]
/// where NBB_cumulativeEmissionDerivPart() is called to compute the derivatives
/// of the bounds.  By miracle of miracles, all the extra \f$k\f$'s that you
/// might expect here magically vanish.  And since you need the same integral as
/// NBB_integrateBlackBodyGroups() here too, it makes sense to evaluate these all
/// at the same time.
///
/// \todo Next time I need to make sure this is right, I should just write out
/// all the steps!
///
/// \param[in] T Temperature of the material.
/// \param[in] k Boltzmann constant such that k*T is in the same units as the group boundaries.
/// \param[in] Bnorm is the normalization to get the intensity or energy density
///            in your units such that Bnorm*T^4 = sum_g B[g] = aT^4 (or multiplied by c for intensit)
/// \param[in] numGroups Number of groups
/// \param[in] groupBounds Group boundaries, must be numGroups+1 long.  First
///            and last ignored and assumed to be zero or infinity
/// \param[out] B Pre-allocated array, numGroups long, where group-wise emission
///            is stored.
/// \param[out] dBdT Pre-allocated array, numGroups long, where group-wise
///            emission derivative is stored.
void NBB_integrateBlackBodyAndDerivGroups(const double T,
                                          const double k,
                                          const double Bnorm,
                                          const int numGroups,
                                          const double *__restrict__ const groupBounds,
                                          double *__restrict__ const B,
                                          double *__restrict__ const dBdT)
{
   const double kT = k * T;
   // Should be rare, if not, never, that we have a zero temperature,
   // but we're being safe.
   if (kT > 0.0)
   {
      const double Tsq = T * T;
      const double fullBNorm = Bnorm * Tsq * Tsq;
      const double fulldBdTNorm = 4.0 * Bnorm * Tsq * T;

      double lowB = 0.0;
      double lowdBdT = 0.0;
      int g;
      for (g = 0; g < numGroups - 1; ++g)
      {
         double highBound = groupBounds[g + 1] / kT;
         double highB = NBB_cumulativeEmission(highBound);
         double highdBdT = NBB_cumulativeEmissionDerivPart(highBound);

         // The order of operations is important here.  The terms that form delta (B)
         // tend to one at large values of the bounds.  But the dBdT terms are small.
         // It is possible that the B-dBdT ≈ B for epsilon > 49.  For group bounds
         // in this range deltaB≈0 too.  But we can have deltadBdT>0.
         //
         // So we first subtract the B terms.  And then the dBdT terms, and finally
         // combine them.  You might be tempted to think this doesn't matter, but
         // for some codes, it did.  Perhaps they should be more robust against
         // small numbers becoming zero, but if we can be more accurate, we should be.
         double deltaB = highB - lowB;
         double deltadBdT = highdBdT - lowdBdT;
         B[g] = fullBNorm * (deltaB);
         // The minus sign here is from Clark, Eq. 27.
         dBdT[g] = fulldBdTNorm * (deltaB - deltadBdT);
         lowB = highB;
         lowdBdT = highdBdT;
      }
      double deltaB = 1.0 - lowB;
      B[numGroups - 1] = fullBNorm * deltaB;
      dBdT[numGroups - 1] = fulldBdTNorm * (deltaB /*- 0.0 for highdBdT  */ + lowdBdT);
   }
   else
   {
      int g;
      for (g = 0; g < numGroups; ++g)
      {
         B[g] = 0.0;
         dBdT[g] = 0.0;
      }
   }
}
}

////////////////////////////////////////////////////////////////////////////////
/// @}
