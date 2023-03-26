#include "PorousFlowBrineSaltCO2Test.h"
#include "FluidPropertiesTestUtils.h"

/**
 * Verify that the correct values for the fluid phase and component indices are supplied
*/
TEST_F(PorousFlowBrineSaltCO2Test, indices)
{
  EXPECT_EQ((unsigned int)3, _fs->numPhases());
  EXPECT_EQ((unsigned int)3, _fs->numComponents());
  EXPECT_EQ((unsigned int)0, _fs->aqueousPhaseIndex());
  EXPECT_EQ((unsigned int)1, _fs->gasPhaseIndex());
//EXPECT_EQ((unsigned int)2, _fs->solidPhaseIndex()); doesn't work. probably because there's no solid phase enum
  EXPECT_EQ((unsigned int)0, _fs->aqueousComponentIndex());
  EXPECT_EQ((unsigned int)1, _fs->gasComponentIndex());
  EXPECT_EQ((unsigned int)2, _fs->saltComponentIndex());
}

/*
 * Verify the calculations of equilibrium mass fraction and its derivatives wrt
 * to the variables P,T,X,Z.

 * Note: All we're doing is re-computing the (equilibrium) mass fractions using some
 * specific 'test' variables based on the following 2 different approach:
 * 1) obtain the values directly using the functions we implemented in the source file
 * 2) manually compute the values w/o those functions.
 * After that, we compare the different results and see whether they match!
*/
TEST_F(PorousFlowBrineSaltCO2Test, equilibriumMassFraction)
{
  // Low temperature regime (T=350)
  DualReal p = 1.0e6;
  Moose::derivInsert(p.derivatives(), _pidx, 1.0);

  DualReal T = 350.0;
  Moose::derivInsert(T.derivatives(), _Tidx, 1.0);

  DualReal Xnacl = 0.1;
  Moose::derivInsert(Xnacl.derivatives(), _Xidx, 1.0);

  const Real dp = 1.0e-2;
  const Real dT = 1.0e-6;
  const Real dx = 1.0e-8;

  DualReal X, Y;
  _fs->equilibriumMassFractions(p, T, Xnacl, X, Y);
  //compare results:
  ABS_TEST(X.value(), 0.0035573020148, 1.0e-10);
  ABS_TEST(Y.value(), 0.0171977397214, 1.0e-10);

/// Verify the salt component in the solid phase (i.e., solid mass fraction, Snacl).
  // First, test the solubilities by comparing the implemented methods in the source file
  // and the results obtained manually:
  DualReal XEQ = _brine_fp->haliteSolubilityWater(350, 1e6);
  DualReal XEQG = _brine_fp->haliteSolubilityGas(350, 1e6);

  REL_TEST(XEQ, 0.104368, 2.0e-2);
  REL_TEST(XEQG, 3.8602e-11, 2.0e-2);

  // Next, test the salt component in the liquid phase:
  DualReal Ynacl = Xnacl * (XEQG/XEQ);
  REL_TEST(Ynacl, 3.6986e-11, 2.0e-2);

  // Finally, Verify the salt component in the solid phase:
  DualReal Snacl = 1 - Xnacl - Ynacl;
  REL_TEST(Snacl, 0.9, 2.0e-2);

/// Verify the derivatives of the mass fractions wrt to the variables
/// Note: the derivatives of Snacl wrt to the variables is not tested!

  DualReal X1, Y1, X2, Y2;
  // two different X,Y values from a small increment in pressure (dp) to compute the derivatives
  _fs->equilibriumMassFractions(p - dp, T, Xnacl, X1, Y1);
  _fs->equilibriumMassFractions(p + dp, T, Xnacl, X2, Y2);

  // Derivatives wrt pressure
  Real dX_dp_fd = (X2.value() - X1.value()) / (2.0 * dp);
  Real dY_dp_fd = (Y2.value() - Y1.value()) / (2.0 * dp);
  //compare results
  REL_TEST(X.derivatives()[_pidx], dX_dp_fd, 1.0e-6);
  REL_TEST(Y.derivatives()[_pidx], dY_dp_fd, 1.0e-6);

  // Derivatives wrt temperature (similar to pressure!)
  _fs->equilibriumMassFractions(p, T - dT, Xnacl, X1, Y1);
  _fs->equilibriumMassFractions(p, T + dT, Xnacl, X2, Y2);

  Real dX_dT_fd = (X2.value() - X1.value()) / (2.0 * dT);
  Real dY_dT_fd = (Y2.value() - Y1.value()) / (2.0 * dT);

  REL_TEST(X.derivatives()[_Tidx], dX_dT_fd, 1.0e-6);
  REL_TEST(Y.derivatives()[_Tidx], dY_dT_fd, 1.0e-6);

  // Derivative wrt salt mass fraction
  _fs->equilibriumMassFractions(p, T, Xnacl - dx, X1, Y1);
  _fs->equilibriumMassFractions(p, T, Xnacl + dx, X2, Y2);

  Real dX_dX_fd = (X2.value() - X1.value()) / (2.0 * dx);
  Real dY_dX_fd = (Y2.value() - Y1).value() / (2.0 * dx);

  REL_TEST(X.derivatives()[_Xidx], dX_dX_fd, 1.0e-6);
  REL_TEST(Y.derivatives()[_Xidx], dY_dX_fd, 1.0e-6);

/// High temperature and pressure regimes (T = 525.15)
  p = 10.0e6;
  Moose::derivInsert(p.derivatives(), _pidx, 1.0);

  T = 525.15;
  Moose::derivInsert(T.derivatives(), _Tidx, 1.0);

  _fs->equilibriumMassFractions(p, T, Xnacl, X, Y);

  ABS_TEST(X.value(), 0.016299479086, 1e-10);
  ABS_TEST(Y.value(), 0.249471400766, 1e-10);

/// Test the salt component in the solid phase (i.e., solid mass fraction, Snacl).
    // First, test the solubilities by comparing the implemented methods in the source file
    // and the results obtained manually:
  XEQ = _brine_fp->haliteSolubilityWater(525.15, 10e6);
  XEQG = _brine_fp->haliteSolubilityGas(350, 10e6);

    REL_TEST(XEQ, 0.139766, 2.0e-2);
    REL_TEST(XEQG, 0.00012207, 2.0e-2);

    // Next, test the salt component in the liquid phase:
    Ynacl = Xnacl * (XEQG/XEQ);
    REL_TEST(Ynacl, 8.73388e-5, 2.0e-2);

    // Finally, Verify the salt component in the solid phase:
    Snacl = 1 - Xnacl - Ynacl;
    REL_TEST(Snacl, 0.8999, 2.0e-2);
}


/*
* Verify the the massfraction method and its derivatives wrt to the variables P,T,X,Z.
*
* Note: This time, we're  re-computing the various mass fractions in each phase using
* specific 'test' variables. Again, we obtain the values using the methods implemented
* in the source file and then compare them with manually calculated results.
*/

TEST_F(PorousFlowBrineSaltCO2Test, MassFraction)
{
  DualReal p = 1.0e6;
  Moose::derivInsert(p.derivatives(), _pidx, 1.0);

  DualReal T = 350.0;
  Moose::derivInsert(T.derivatives(), _Tidx, 1.0);

  DualReal Xnacl = 0.1;
  Moose::derivInsert(Xnacl.derivatives(), _Xidx, 1.0);

  FluidStatePhaseEnum phase_state;
  const unsigned int np = _fs->numPhases();
  const unsigned int nc = _fs->numComponents();
  std::vector<FluidStateProperties> fsp(np, FluidStateProperties(nc));

  /// Liquid region (Z<_Zmin)
  DualReal Z = 0.0001;
  Moose::derivInsert(Z.derivatives(), _Zidx, 1.0);

  _fs->massFractions(p, T, Xnacl, Z, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::LIQUID);

  // Obtain the mass fractions using the methods from source file
  DualReal Xh2o =   fsp[0].mass_fraction[0];
  DualReal Xco2 =   fsp[0].mass_fraction[1];
  DualReal Xnacl2 = fsp[0].mass_fraction[2];  // halite (component 2) in liquid (phase 0)
  DualReal Yh2o =   fsp[1].mass_fraction[0];
  DualReal Yco2 =   fsp[1].mass_fraction[1];
  DualReal Ynacl  = fsp[1].mass_fraction[2];   // halite (component 2) in gas (phase 1)
  DualReal Snacl  = fsp[2].mass_fraction[2];;  // halite (component 2) in solid (phase 2)

  // Now, compare and verfify the values
  ABS_TEST(Xco2.value(), Z.value(), 1.0e-8);
  ABS_TEST(Yco2.value(), 0.0, 1.0e-8);
  ABS_TEST(Xh2o.value(), (1.0 - Z.value())/2.0, 1.0e-8);
  ABS_TEST(Yh2o.value(), 0.0, 1.0e-8);
  ABS_TEST(Xnacl2.value(), (1.0 - Z.value())/2.0, 1.0e-8);
  ABS_TEST(Ynacl.value(), 0.0, 1.0e-8);
  ABS_TEST(Snacl.value(), 0.0, 1.0e-8);

  /// Verify derivatives by comparing with manually computed results
  // note: derivatives of co2 in the gas phase (Yco2) is replaced by
  // water in the liquid phase (Xh2o) to ensure all tests are done in the liquid phase.
  ABS_TEST(Xco2.derivatives()[_pidx], 0.0, 1.0e-8);
  ABS_TEST(Xco2.derivatives()[_Tidx], 0.0, 1.0e-8);
  ABS_TEST(Xco2.derivatives()[_Xidx], 0.0, 1.0e-8);
  ABS_TEST(Xco2.derivatives()[_Zidx], 1.0, 1.0e-8);
  ABS_TEST(Xh2o.derivatives()[_pidx], 0.0, 1.0e-8);
  ABS_TEST(Xh2o.derivatives()[_Tidx], 0.0, 1.0e-8);
  ABS_TEST(Xh2o.derivatives()[_Xidx], 0.0, 1.0e-8);
  ABS_TEST(Xh2o.derivatives()[_Zidx], 0.5, 1.0e-8);
  ABS_TEST(Xnacl2.derivatives()[_pidx], 0.0, 1.0e-8);
  ABS_TEST(Xnacl2.derivatives()[_Tidx], 0.0, 1.0e-8);
  ABS_TEST(Xnacl2.derivatives()[_Xidx], 0.0, 1.0e-8);
  ABS_TEST(Xnacl2.derivatives()[_Zidx], 0.5, 1.0e-8);
//  ABS_TEST(Yco2.derivatives()[_pidx], 0.0, 1.0e-8);
//  ABS_TEST(Yco2.derivatives()[_Tidx], 0.0, 1.0e-8);
//  ABS_TEST(Yco2.derivatives()[_Xidx], 0.0, 1.0e-8);
//  ABS_TEST(Yco2.derivatives()[_Zidx], 0.0, 1.0e-8);

  // Gas region
  Z = 0.995;
  Moose::derivInsert(Z.derivatives(), _Zidx, 1.0);

  _fs->massFractions(p, T, Xnacl, Z, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::GAS);

  // Obtain mass fraction values directly
  DualReal Ynacl1;
  Xco2 = fsp[0].mass_fraction[1];
  Yco2 = fsp[1].mass_fraction[1];
  Xh2o = fsp[0].mass_fraction[0];
  Yh2o = fsp[1].mass_fraction[0];
  Xnacl = fsp[0].mass_fraction[2];  // halite (component 2) in liquid (phase 0)
  Ynacl1 = fsp[1].mass_fraction[2];    // halite (component 2) in liquid (phase 0)
  Snacl  = fsp[2].mass_fraction[2];;  // halite (component 2) in solid (phase 2)

  // Now, compare and verfify the values
  ABS_TEST(Xco2.value(), 0.0, 1.0e-8);
  ABS_TEST(Yco2.value(), Z.value(), 1.0e-8);
  ABS_TEST(Xh2o.value(), 0.0, 1.0e-8);
  ABS_TEST(Yh2o.value(), (1.0 - Z.value())/2.0, 1.0e-8);
  ABS_TEST(Xnacl.value(), 0.0, 1.0e-8);
  ABS_TEST(Ynacl1.value(), (1.0 - Z.value())/2.0, 1.0e-8);
  ABS_TEST(Snacl.value(), 0.0, 1.0e-8);

  /// Verify derivatives by comparing with manually computed results
  // note: derivatives of co2 in the liquid phase (Xco2) is replaced by
  // water in the gas phase (Yh2o) to ensure all tests are done in the gas phase.
  ABS_TEST(Yco2.derivatives()[_pidx], 0.0, 1.0e-8);
  ABS_TEST(Yco2.derivatives()[_Tidx], 0.0, 1.0e-8);
  ABS_TEST(Yco2.derivatives()[_Xidx], 0.0, 1.0e-8);
  ABS_TEST(Yco2.derivatives()[_Zidx], 1.0, 1.0e-8);
  ABS_TEST(Yco2.derivatives()[_pidx], 0.0, 1.0e-8);
  ABS_TEST(Yco2.derivatives()[_Tidx], 0.0, 1.0e-8);
  ABS_TEST(Yco2.derivatives()[_Xidx], 0.0, 1.0e-8);
  ABS_TEST(Yco2.derivatives()[_Zidx], 0.5, 1.0e-8);
  ABS_TEST(Ynacl1.derivatives()[_pidx], 0.0, 1.0e-8);
  ABS_TEST(Ynacl1.derivatives()[_Tidx], 0.0, 1.0e-8);
  ABS_TEST(Ynacl1.derivatives()[_Xidx], 0.0, 1.0e-8);
  ABS_TEST(Ynacl1.derivatives()[_Zidx], 0.5, 1.0e-8);
 // ABS_TEST(Xco2.derivatives()[_pidx], 0.0, 1.0e-8);
 // ABS_TEST(Xco2.derivatives()[_Tidx], 0.0, 1.0e-8);
 // ABS_TEST(Xco2.derivatives()[_Xidx], 0.0, 1.0e-8);
 // ABS_TEST(Xco2.derivatives()[_Zidx], 0.0, 1.0e-8);

/// Two phase region. In this region, the mass fractions and derivatives can
/// be verified using the equilibrium mass fraction derivatives that have
/// already been verified above
  Z = 0.45;
  Moose::derivInsert(Z.derivatives(), _Zidx, 1.0);

  _fs->massFractions(p, T, Xnacl, Z, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

  // compute new Eq. mass fractions using the already verified
  // 'equilibriumMassFractions' method.
  DualReal Xco2_eq, Yh2o_eq;
  _fs->equilibriumMassFractions(p, T, Xnacl, Xco2_eq, Yh2o_eq);

  //obtain the eq.mass fractions directly from source file.
  Xco2 = fsp[0].mass_fraction[1];
  Yco2 = fsp[1].mass_fraction[1];
  Xh2o = fsp[0].mass_fraction[0];
  Yh2o = fsp[1].mass_fraction[0];
  //compare them with newly computed ones Xco2_eq and Yh2o_eq.
  ABS_TEST(Xco2, Xco2_eq, 1.0e-8);
  ABS_TEST(Yco2, 1.0 - Yh2o_eq, 1.0e-8);
  ABS_TEST(Xh2o, 1.0 - Xco2_eq, 1.0e-8);
  ABS_TEST(Yh2o, Yh2o_eq, 1.0e-8);

/// Use finite differences to verify that the derivatives of the mass fraction
// wrt Z in the two-phase region is unaffected by Z.
// note: the derivatives of salt mass fraction is not tested.
  const Real dZ = 1.0e-8;
  _fs->massFractions(p, T, Xnacl, Z + dZ, phase_state, fsp);
  DualReal Xco21 = fsp[0].mass_fraction[1];
  DualReal Yco21 = fsp[1].mass_fraction[1];
  _fs->massFractions(p, T, Xnacl, Z - dZ, phase_state, fsp);
  DualReal Xco22 = fsp[0].mass_fraction[1];
  DualReal Yco22 = fsp[1].mass_fraction[1];

  ABS_TEST(Xco2.derivatives()[_Zidx], (Xco21.value() - Xco22.value()) / (2.0 * dZ), 1.0e-8);
  ABS_TEST(Yco2.derivatives()[_Zidx], (Yco21.value() - Yco22.value()) / (2.0 * dZ), 1.0e-8);
}



/*
 * Verify calculation of gas density, viscosity enthalpy, and derivatives. Note that as
 * these properties don't depend on (equilibrium) mass fraction because they are computed for
 * single phase regions, only the gas region needs to be tested (the calculations are identical
 * in the two phase region).
 *
 * Note: All we are doing is re-computing the property values using 2 different approaches:
 * 1).the functions from our source file and 2) directly w/o those functions. Then we
 * can compare the results and see whether they match!
 */
TEST_F(PorousFlowBrineSaltCO2Test, gasProperties)
{
  DualReal p = 1.0e6;
  Moose::derivInsert(p.derivatives(), _pidx, 1.0);

  DualReal T = 350.0;
  Moose::derivInsert(T.derivatives(), _Tidx, 1.0);

  DualReal Xnacl = 0.1;
  Moose::derivInsert(Xnacl.derivatives(), _Xidx, 1.0);

  FluidStatePhaseEnum phase_state;
  const unsigned int np = _fs->numPhases();
  const unsigned int nc = _fs->numComponents();
  std::vector<FluidStateProperties> fsp(np, FluidStateProperties(nc));

  // Gas region
  DualReal Z = 0.995;
  Moose::derivInsert(Z.derivatives(), _Zidx, 1.0);

  _fs->massFractions(p, T, Xnacl, Z, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::GAS);

  // Verify fluid density, viscosity and enthalpy
  // note: here, we use the methods/functions (i.e., gasProperties) from the source file.
  _fs->gasProperties(p, T, fsp);
  DualReal gas_density = fsp[1].density;
  DualReal gas_viscosity = fsp[1].viscosity;
  DualReal gas_enthalpy = fsp[1].enthalpy;

  // now, we re-compute these properties again (without the methods)!
  Real psat = _water_fp->vaporPressure(T.value());

  Real co2_density = _co2_fp->rho_from_p_T(p.value(), T.value());
  Real co2_viscosity = _co2_fp->mu_from_p_T(p.value(), T.value());
  Real co2_enthalpy = _co2_fp->h_from_p_T(p.value(), T.value());

  Real vapor_density =   _water_fp->rho_from_p_T((1.0 -Z.value()) * psat, T.value());
  Real vapor_viscosity = _water_fp->mu_from_p_T((1.0 - Z.value()) * psat, T.value());
  Real vapor_enthalpy =   _water_fp->h_from_p_T((1.0 - Z.value()) * psat, T.value());

  Real density = co2_density + vapor_density;
  Real viscosity = Z.value() * co2_viscosity + (1.0 - Z.value()) * vapor_viscosity;
  Real enthalpy = Z.value() * co2_enthalpy + (1.0 - Z.value()) * vapor_enthalpy;

 // now we can compare the two results:
  ABS_TEST(gas_density.value(), density, 1.0e-8);
  ABS_TEST(gas_viscosity.value(), viscosity, 1.0e-8);
  ABS_TEST(gas_enthalpy.value(), enthalpy, 1.0e-8);


  // Verify derivatives
  // note: the derivatives are w.r.t all variables (i.e., P,T,Z)

  // start with P.
  const Real dp = 1.0e-1; //small increment of P (change in P)
//  obtain parameters needed for the derivatives.
  _fs->gasProperties(p + dp, T, fsp);
  Real rho1 = fsp[1].density.value();
  Real mu1 = fsp[1].viscosity.value();
  Real h1 = fsp[1].enthalpy.value();

  _fs->gasProperties(p - dp, T, fsp);
  Real rho2 = fsp[1].density.value();
  Real mu2 = fsp[1].viscosity.value();
  Real h2 = fsp[1].enthalpy.value();
// now we compare the derivatives. One from the value obtained from the source file
// and the other (i.e., (prop1 - prop2) / (2.0 * dp)) computed here directly.
  REL_TEST(gas_density.derivatives()[_pidx], (rho1 - rho2) / (2.0 * dp), 1.0e-6);
  REL_TEST(gas_viscosity.derivatives()[_pidx], (mu1 - mu2) / (2.0 * dp), 1.0e-6);
  REL_TEST(gas_enthalpy.derivatives()[_pidx], (h1 - h2) / (2.0 * dp), 1.0e-6);

// do same for derivative w.r.t T
  const Real dT = 1.0e-3;

  _fs->gasProperties(p, T + dT, fsp);
  rho1 = fsp[1].density.value();
  mu1 = fsp[1].viscosity.value();
  h1 = fsp[1].enthalpy.value();

  _fs->gasProperties(p, T - dT, fsp);
  rho2 = fsp[1].density.value();
  mu2 = fsp[1].viscosity.value();
  h2 = fsp[1].enthalpy.value();

  REL_TEST(gas_density.derivatives()[_Tidx], (rho1 - rho2) / (2.0 * dT), 1.0e-6);
  REL_TEST(gas_viscosity.derivatives()[_Tidx], (mu1 - mu2) / (2.0 * dT), 1.0e-6);
  REL_TEST(gas_enthalpy.derivatives()[_Tidx], (h1 - h2) / (2.0 * dT), 1.0e-6);

 // do same for Z
 // Note: (total) mass fraction changes with Z
  const Real dZ = 1.0e-8;

  _fs->massFractions(p, T, Xnacl, Z + dZ, phase_state, fsp);
  _fs->gasProperties(p, T, fsp);
  rho1 = fsp[1].density.value();
  mu1 = fsp[1].viscosity.value();
  h1 = fsp[1].enthalpy.value();

  _fs->massFractions(p, T, Xnacl, Z - dZ, phase_state, fsp);
  _fs->gasProperties(p, T, fsp);
  rho2 = fsp[1].density.value();
  mu2 = fsp[1].viscosity.value();
  h2 = fsp[1].enthalpy.value();

  ABS_TEST(gas_density.derivatives()[_Zidx], (rho1 - rho2) / (2.0 * dZ), 1.0e-8);
  ABS_TEST(gas_viscosity.derivatives()[_Zidx], (mu1 - mu2) / (2.0 * dZ), 1.0e-8);
  ABS_TEST(gas_enthalpy.derivatives()[_Zidx], (h1 - h2) / (2.0 * dZ), 1.0e-6);
}



/*
 * Verify calculation of liquid properties (i.e., density, viscosity, enthalpy) and
 * their derivatives wrt the variables P,T,X and Z.
 *
 * Note: Again, all we are doing is re-computing these properties using 2 different
 * approach: 1).the functions from our source file and 2) directly w/o those functions.
 * Then we compare the results and see whether they match!
*/
TEST_F(PorousFlowBrineSaltCO2Test, liquidProperties)
{
  DualReal p = 1.0e6;
  Moose::derivInsert(p.derivatives(), _pidx, 1.0);

  DualReal T = 350.0;
  Moose::derivInsert(T.derivatives(), _Tidx, 1.0);

  DualReal Xnacl = 0.1;
  Moose::derivInsert(Xnacl.derivatives(), _Xidx, 1.0);

  FluidStatePhaseEnum phase_state;
  const unsigned int np = _fs->numPhases();
  const unsigned int nc = _fs->numComponents();
  std::vector<FluidStateProperties> fsp(np, FluidStateProperties(nc));

  // Liquid region
  DualReal Z = 0.0001;
  Moose::derivInsert(Z.derivatives(), _Zidx, 1.0);

// note: since Z < _Zmin, total mass fraction will be used.
// Also only liquid region will be considered (no multiphase!).
  _fs->massFractions(p, T, Xnacl, Z, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::LIQUID);

  // Verify fluid density and viscosity
  // note: here, we use the methods/functions (i.e.,liquidProperties) implemented in the source file.
  _fs->liquidProperties(p, T, Xnacl, fsp);
  DualReal liquid_density = fsp[0].density;
  DualReal liquid_viscosity = fsp[0].viscosity;
  DualReal liquid_enthalpy = fsp[0].enthalpy;

// now, we re-compute these properties again (without the methods)!
  Real co2_partial_density = _fs->partialDensityCO2(T).value();
  Real brine_density = _brine_fp->rho_from_p_T_X(p.value(), T.value(), Xnacl.value());

  Real density = 1.0 / (Z.value() / co2_partial_density + (1.0 - Z.value()) / brine_density);

  Real viscosity = _brine_fp->mu_from_p_T_X(p.value(), T.value(), Xnacl.value());

  Real brine_enthalpy = _brine_fp->h_from_p_T_X(p.value(), T.value(), Xnacl.value());
  Real hdis = _fs->enthalpyOfDissolution(T).value();
  Real co2_enthalpy = _co2_fp->h_from_p_T(p.value(), T.value());
  Real enthalpy = (1.0 - Z.value()) * brine_enthalpy + Z.value() * (co2_enthalpy + hdis);

  // now we can compare the two results:
  ABS_TEST(liquid_density.value(), density, 1.0e-12);
  ABS_TEST(liquid_viscosity.value(), viscosity, 1.0e-12);
  ABS_TEST(liquid_enthalpy.value(), enthalpy, 1.0e-12);


  // Verify derivatives
  // Derivatives wrt pressure
  const Real dp = 1.0;
  _fs->liquidProperties(p + dp, T, Xnacl, fsp);
  Real rho1 = fsp[0].density.value();
  Real mu1 = fsp[0].viscosity.value();
  Real h1 = fsp[0].enthalpy.value();

  _fs->liquidProperties(p - dp, T, Xnacl, fsp);
  Real rho2 = fsp[0].density.value();
  Real mu2 = fsp[0].viscosity.value();
  Real h2 = fsp[0].enthalpy.value();

  //compare results
  REL_TEST(liquid_density.derivatives()[_pidx], (rho1 - rho2) / (2.0 * dp), 1.0e-6);
  REL_TEST(liquid_viscosity.derivatives()[_pidx], (mu1 - mu2) / (2.0 * dp), 2.0e-6);
  REL_TEST(liquid_enthalpy.derivatives()[_pidx], (h1 - h2) / (2.0 * dp), 1.0e-6);

  // Derivatives wrt temperature
  const Real dT = 1.0e-4;
  _fs->liquidProperties(p, T + dT, Xnacl, fsp);
  rho1 = fsp[0].density.value();
  mu1 = fsp[0].viscosity.value();
  h1 = fsp[0].enthalpy.value();

  _fs->liquidProperties(p, T - dT, Xnacl, fsp);
  rho2 = fsp[0].density.value();
  mu2 = fsp[0].viscosity.value();
  h2 = fsp[0].enthalpy.value();

  //compare results
  REL_TEST(liquid_density.derivatives()[_Tidx], (rho1 - rho2) / (2.0 * dT), 1.0e-6);
  REL_TEST(liquid_viscosity.derivatives()[_Tidx], (mu1 - mu2) / (2.0 * dT), 1.0e-6);
  REL_TEST(liquid_enthalpy.derivatives()[_Tidx], (h1 - h2) / (2.0 * dT), 1.0e-6);

  // Derivatives wrt Xnacl
  const Real dx = 1.0e-8;
  _fs->liquidProperties(p, T, Xnacl + dx, fsp);
  rho1 = fsp[0].density.value();
  mu1 = fsp[0].viscosity.value();
  h1 = fsp[0].enthalpy.value();

  _fs->liquidProperties(p, T, Xnacl - dx, fsp);
  rho2 = fsp[0].density.value();
  mu2 = fsp[0].viscosity.value();
  h2 = fsp[0].enthalpy.value();

 //compare results
  REL_TEST(liquid_density.derivatives()[_Xidx], (rho1 - rho2) / (2.0 * dx), 1.0e-6);
  REL_TEST(liquid_viscosity.derivatives()[_Xidx], (mu1 - mu2) / (2.0 * dx), 1.0e-6);
  REL_TEST(liquid_enthalpy.derivatives()[_Xidx], (h1 - h2) / (2.0 * dx), 1.0e-6);

  // Derivatives wrt Z
  const Real dZ = 1.0e-8;
  _fs->massFractions(p, T, Xnacl, Z + dZ, phase_state, fsp);
  _fs->liquidProperties(p, T, Xnacl, fsp);
  rho1 = fsp[0].density.value();
  mu1 = fsp[0].viscosity.value();
  h1 = fsp[0].enthalpy.value();

  _fs->massFractions(p, T, Xnacl, Z - dZ, phase_state, fsp);
  _fs->liquidProperties(p, T, Xnacl, fsp);
  rho2 = fsp[0].density.value();
  mu2 = fsp[0].viscosity.value();
  h2 = fsp[0].enthalpy.value();

  // compare results
  REL_TEST(liquid_density.derivatives()[_Zidx], (rho1 - rho2) / (2.0 * dZ), 1.0e-6);
  ABS_TEST(liquid_viscosity.derivatives()[_Zidx], (mu1 - mu2) / (2.0 * dZ), 1.0e-6);
  REL_TEST(liquid_enthalpy.derivatives()[_Zidx], (h1 - h2) / (2.0 * dZ), 1.0e-6);


// Verify (the derivatives of) the Liquid properties in the multiphase region w.r.t
// the variables P,T,X and Z.
  Z = 0.045;
  Moose::derivInsert(Z.derivatives(), _Zidx, 1.0);

// Note: Z > _Zmin, so (equilibrium) mass fraction will be used and defined on the
// multiphase region. Here, we verify only the derivatives since the properties
// in the single phase regions are identical to the multiphase regions.
  _fs->massFractions(p, T, Xnacl, Z, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

// Initialize the properties:
  _fs->liquidProperties(p, T, Xnacl, fsp);

  liquid_density = fsp[0].density;
  liquid_viscosity = fsp[0].viscosity;
  liquid_enthalpy = fsp[0].enthalpy;

  // Derivatives wrt pressure in the multiphase region
  _fs->massFractions(p + dp, T, Xnacl, Z, phase_state, fsp);
  _fs->liquidProperties(p + dp, T, Xnacl, fsp);
  rho1 = fsp[0].density.value();
  mu1 = fsp[0].viscosity.value();
  h1 = fsp[0].enthalpy.value();

  _fs->massFractions(p - dp, T, Xnacl, Z, phase_state, fsp);
  _fs->liquidProperties(p - dp, T, Xnacl, fsp);
  rho2 = fsp[0].density.value();
  mu2 = fsp[0].viscosity.value();
  h2 = fsp[0].enthalpy.value();

  // compare results
  REL_TEST(liquid_density.derivatives()[_pidx], (rho1 - rho2) / (2.0 * dp), 1.0e-6);
  REL_TEST(liquid_viscosity.derivatives()[_pidx], (mu1 - mu2) / (2.0 * dp), 2.0e-6);
  REL_TEST(liquid_enthalpy.derivatives()[_pidx], (h1 - h2) / (2.0 * dp), 1.0e-6);

  // Derivatives wrt temperature in the multiphase region
  _fs->massFractions(p, T + dT, Xnacl, Z, phase_state, fsp);
  _fs->liquidProperties(p, T + dT, Xnacl, fsp);
  rho1 = fsp[0].density.value();
  mu1 = fsp[0].viscosity.value();
  h1 = fsp[0].enthalpy.value();

  _fs->massFractions(p, T - dT, Xnacl, Z, phase_state, fsp);
  _fs->liquidProperties(p, T - dT, Xnacl, fsp);
  rho2 = fsp[0].density.value();
  mu2 = fsp[0].viscosity.value();
  h2 = fsp[0].enthalpy.value();

// compare results
  REL_TEST(liquid_density.derivatives()[_Tidx], (rho1 - rho2) / (2.0 * dT), 1.0e-6);
  REL_TEST(liquid_viscosity.derivatives()[_Tidx], (mu1 - mu2) / (2.0 * dT), 1.0e-6);
  REL_TEST(liquid_enthalpy.derivatives()[_Tidx], (h1 - h2) / (2.0 * dT), 1.0e-6);

  // Derivatives wrt Xnacl in the multiphase region
  _fs->massFractions(p, T, Xnacl + dx, Z, phase_state, fsp);
  _fs->liquidProperties(p, T, Xnacl + dx, fsp);
  rho1 = fsp[0].density.value();
  mu1 = fsp[0].viscosity.value();
  h1 = fsp[0].enthalpy.value();

  _fs->massFractions(p, T, Xnacl - dx, Z, phase_state, fsp);
  _fs->liquidProperties(p, T, Xnacl - dx, fsp);
  rho2 = fsp[0].density.value();
  mu2 = fsp[0].viscosity.value();
  h2 = fsp[0].enthalpy.value();

// compare results
  REL_TEST(liquid_density.derivatives()[_Xidx], (rho1 - rho2) / (2.0 * dx), 1.0e-6);
  REL_TEST(liquid_viscosity.derivatives()[_Xidx], (mu1 - mu2) / (2.0 * dx), 1.0e-6);
  REL_TEST(liquid_enthalpy.derivatives()[_Xidx], (h1 - h2) / (2.0 * dx), 1.0e-6);

  // Derivatives wrt Z in the multiphase region
  _fs->massFractions(p, T, Xnacl, Z + dZ, phase_state, fsp);
  _fs->liquidProperties(p, T, Xnacl, fsp);
  rho1 = fsp[0].density.value();
  mu1 = fsp[0].viscosity.value();
  h1 = fsp[0].enthalpy.value();

  _fs->massFractions(p, T, Xnacl, Z - dZ, phase_state, fsp);
  _fs->liquidProperties(p, T, Xnacl, fsp);
  rho2 = fsp[0].density.value();
  mu2 = fsp[0].viscosity.value();
  h2 = fsp[0].enthalpy.value();

// compare results
  ABS_TEST(liquid_density.derivatives()[_Zidx], (rho1 - rho2) / (2.0 * dZ), 1.0e-6);
  ABS_TEST(liquid_viscosity.derivatives()[_Zidx], (mu1 - mu2) / (2.0 * dZ), 1.0e-6);
  ABS_TEST(liquid_enthalpy.derivatives()[_Zidx], (h1 - h2) / (2.0 * dZ), 1.0e-6);
}


/*
 * Verify calculation of the single-phase solid properties (density, viscosity, enthalpy)
 * and their derivatives wrt to the variables P,T,X and Z. Note that these properties
 * depend on (equilibrium) mass fractions!
*/
TEST_F(PorousFlowBrineSaltCO2Test, solidProperties)
{
  DualReal p = 1.0e6;
  Moose::derivInsert(p.derivatives(), _pidx, 1.0);

  DualReal T = 350.0;
  Moose::derivInsert(T.derivatives(), _Tidx, 1.0);

  DualReal Xnacl = 0.1;
  Moose::derivInsert(Xnacl.derivatives(), _Xidx, 1.0);

  FluidStatePhaseEnum phase_state;
  const unsigned int np = _fs->numPhases();
  const unsigned int nc = _fs->numComponents();
  std::vector<FluidStateProperties> fsp(np, FluidStateProperties(nc));

  // Single-phse Solid region
  DualReal Z = 0.995;
  Moose::derivInsert(Z.derivatives(), _Zidx, 1.0);

// note: Z > _Zmin so we should use (equilibrium) mass fraction. Here, we compute
// the solid phase properties in the multiphase (i.e., Two-Phase) region since
// solid phase enum is N/A. These properties are, however, identical in both regions!.
  _fs->massFractions(p, T, Xnacl, Z, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

  // Verify solid properties (density, viscosity and enthalpy) using the method
  // (i.e., solidProperties) implemented in the source file.
  _fs->solidProperties(p, T, fsp);
  DualReal solid_density = fsp[2].density;
  DualReal solid_viscosity = fsp[2].viscosity;
  DualReal solid_enthalpy = fsp[2].enthalpy;

  // now, we re-compute these properties again (without the methods)!
  Real halite_density = _brine_fp->halite_rho_from_p_T(p.value(), T.value());
  Real halite_enthalpy = _brine_fp->halite_h_from_p_T(p.value(), T.value());

  Real density = halite_density;
  Real enthalpy = halite_enthalpy;

 // now we can compare the two results:
  ABS_TEST(solid_density.value(), density, 1.0e-8);
  ABS_TEST(solid_enthalpy.value(), enthalpy, 1.0e-8);


// Verify derivatives
// note: the derivatives are w.r.t all variables (i.e., P,T,X,Z)
// start with P.
  const Real dp = 1.0e-1; //small increment of P (change in P)
  _fs->solidProperties(p + dp, T, fsp);
  Real rho1 = fsp[2].density.value();
  Real h1 = fsp[2].enthalpy.value();

  _fs->solidProperties(p - dp, T, fsp);
  Real rho2 = fsp[2].density.value();
  Real h2 = fsp[2].enthalpy.value();

// compare derivatives. One from the value obtained from the source file
// and the other (i.e., (prop1 - prop2) / (2.0 * dp)) computed directly here.
  REL_TEST(solid_density.derivatives()[_pidx], (rho1 - rho2) / (2.0 * dp), 1.0e-6);
  REL_TEST(solid_enthalpy.derivatives()[_pidx], (h1 - h2) / (2.0 * dp), 1.0e-6);

// do same for derivative w.r.t T
  const Real dT = 1.0e-3;
  _fs->solidProperties(p, T + dT, fsp);
  rho1 = fsp[2].density.value();
  h1 = fsp[2].enthalpy.value();

  _fs->solidProperties(p, T - dT, fsp);
  rho2 = fsp[2].density.value();
  h2 = fsp[2].enthalpy.value();

// compare derivatives.
  REL_TEST(solid_density.derivatives()[_Tidx], (rho1 - rho2) / (2.0 * dT), 1.0e-6);
  REL_TEST(solid_enthalpy.derivatives()[_Tidx], (h1 - h2) / (2.0 * dT), 1.0e-6);

 // do same for derivatives wrt Z
 // Note: (equilibrium) mass fraction changes with Z
  const Real dZ = 1.0e-8;

  _fs->massFractions(p, T, Xnacl, Z + dZ, phase_state, fsp);
  _fs->solidProperties(p, T, fsp);
  rho1 = fsp[2].density.value();
  h1 = fsp[2].enthalpy.value();

  _fs->massFractions(p, T, Xnacl, Z - dZ, phase_state, fsp);
  _fs->solidProperties(p, T, fsp);
  rho2 = fsp[2].density.value();
  h2 = fsp[2].enthalpy.value();

  ABS_TEST(solid_density.derivatives()[_Zidx], (rho1 - rho2) / (2.0 * dZ), 1.0e-8);
  ABS_TEST(solid_enthalpy.derivatives()[_Zidx], (h1 - h2) / (2.0 * dZ), 1.0e-6);
}



/*
 * Verify calculation of gas saturation and derivatives in the two-phase region
*/
TEST_F(PorousFlowBrineSaltCO2Test, saturationGAS)
{
  DualReal p = 1.0e6;
  Moose::derivInsert(p.derivatives(), _pidx, 1.0);

  DualReal T = 350.0;
  Moose::derivInsert(T.derivatives(), _Tidx, 1.0);

  DualReal Xnacl = 0.1;
  Moose::derivInsert(Xnacl.derivatives(), _Xidx, 1.0);

  FluidStatePhaseEnum phase_state;
  const unsigned int np = _fs->numPhases();
  const unsigned int nc = _fs->numComponents();
  std::vector<FluidStateProperties> fsp(np, FluidStateProperties(nc));

  // In the two-phase region, the mass fractions are the equilibrium values, so
  // a temporary value of Z can be used (as long as it corresponds to the two-phase
  // region)
  DualReal Z = 0.45; // Used in the massFraction method to determine the type of enum to be used.
  Moose::derivInsert(Z.derivatives(), _Zidx, 1.0);

  _fs->massFractions(p, T, Xnacl, Z, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

  // Calculate Z that gives a gas saturation of 0.25
  // note: we'll be comparing this value with a saturation value
  // computed using the saturationGAS method implemented in the source file !
  DualReal gas_saturation = 0.25;
  DualReal liquid_pressure = p - _pc->capillaryPressure(1.0 - gas_saturation);

  // Calculate gas density and liquid density
  _fs->gasProperties(p, T, fsp);
  _fs->liquidProperties(liquid_pressure, T, Xnacl, fsp);

  // The mass fraction that corresponds to a gas_saturation = 0.25
  DualReal Zc = (gas_saturation * fsp[1].density * fsp[1].mass_fraction[1] +
                 (1.0 - gas_saturation) * fsp[0].density * fsp[0].mass_fraction[1]) /
                (gas_saturation * fsp[1].density + (1.0 - gas_saturation) * fsp[0].density);

  // Calculate the gas saturation based on Zc!
  DualReal saturationGAS = _fs->saturationGAS(p, T, Xnacl, Zc, fsp);
  // compare the gas saturations!
  ABS_TEST(saturationGAS.value(), gas_saturation.value(), 1.0e-6);


  /// Now, verify the derivatives of the gas saturation
  gas_saturation = _fs->saturation(p, T, Xnacl, Z, fsp);

  // Derivative wrt pressure
  const Real dp = 1.0e-1;
  _fs->massFractions(p + dp, T, Xnacl, Z, phase_state, fsp);
  Real gsat1 = _fs->saturation(p + dp, T, Xnacl, Z, fsp).value();

  _fs->massFractions(p - dp, T, Xnacl, Z, phase_state, fsp);
  Real gsat2 = _fs->saturation(p - dp, T, Xnacl, Z, fsp).value();

  //compare the derivatives
  REL_TEST(gas_saturation.derivatives()[_pidx], (gsat1 - gsat2) / (2.0 * dp), 1.0e-6);

  // Derivative wrt temperature
  const Real dT = 1.0e-4;
  _fs->massFractions(p, T + dT, Xnacl, Z, phase_state, fsp);
  gsat1 = _fs->saturation(p, T + dT, Xnacl, Z, fsp).value();

  _fs->massFractions(p, T - dT, Xnacl, Z, phase_state, fsp);
  gsat2 = _fs->saturation(p, T - dT, Xnacl, Z, fsp).value();

  //compare the derivatives
  REL_TEST(gas_saturation.derivatives()[_Tidx], (gsat1 - gsat2) / (2.0 * dT), 1.0e-6);

  // Derivative wrt Xnacl
  const Real dx = 1.0e-8;
  _fs->massFractions(p, T, Xnacl + dx, Z, phase_state, fsp);
  gsat1 = _fs->saturation(p, T, Xnacl + dx, Z, fsp).value();

  _fs->massFractions(p, T, Xnacl - dx, Z, phase_state, fsp);
  gsat2 = _fs->saturation(p, T, Xnacl - dx, Z, fsp).value();

 //compare the derivatives
 REL_TEST(gas_saturation.derivatives()[_Xidx], (gsat1 - gsat2) / (2.0 * dx), 1.0e-6);

  // Derivative wrt Z
  const Real dZ = 1.0e-8;

  _fs->massFractions(p, T, Xnacl, Z, phase_state, fsp);
  gsat1 = _fs->saturation(p, T, Xnacl, Z + dZ, fsp).value();

  gsat2 = _fs->saturation(p, T, Xnacl, Z - dZ, fsp).value();

  //compare the derivatives
  REL_TEST(gas_saturation.derivatives()[_Zidx], (gsat1 - gsat2) / (2.0 * dZ), 1.0e-6);
}
