#include "BrineFluidPropertiesBetaTest.h"
#include "FluidPropertiesTestUtils.h"

#include "Water97FluidProperties.h"
#include "NaClFluidProperties.h"

/**
 * Verify the calculations made by the halitesolubility method implemented
 * in the source file using data from Driesner (2017):The system H2O-NaCl.
 * Part II: Correlations for molar volume, enthalpy, and isobaric heat capacity from 0
 * to 1000 C, 1 to 500 bar, and 0 to 1 Xnacl, Geochimica et Cosmochimica Acta 71, 4902-4919.
 *
* Note: All that is being done here is re-computing the solubility values using some
* specific 'test' variables (P and T) based on the following 2 different approaches:
* 1) obtain the values directly using the functions we implemented in the source file
* 2) manually compute the values w/o those functions.
* After that, we compare the different results and see whether they match!
*/
TEST_F(BrineFluidPropertiesBetaTest, solubilityWater)
{
  REL_TEST(_fp->haliteSolubilityWater(357.15,5e7), 0.106211, 2.0e-2);
  REL_TEST(_fp->haliteSolubilityWater(457.15,5e7), 0.12252, 2.0e-2);
  REL_TEST(_fp->haliteSolubilityWater(557.15,5e7), 0.151447, 2.0e-2);
}

/**
 * Verify calculation of halite solubility in (water) vapour using data from
 * Palliser and McKibbin (1998):A Model for Deep Geothermal Brines, I:
 *  T-p-X State-Space Description pg 78
 */
TEST_F(BrineFluidPropertiesBetaTest, solubilityGas)
{
  REL_TEST(_fp->haliteSolubilityGas(357.15,1e7), 0.00012207, 2.0e-2);
  REL_TEST(_fp->haliteSolubilityGas(357.15,8e6), 2.86217e-05, 2.0e-2);
  REL_TEST(_fp->haliteSolubilityGas(357.15,5e6), 1.3487e-06, 2.0e-2);
}

/**
 * Verify calculation of halite properties: density and enthalpy.
 * Density data from Brown, "The NaCl pressure standard", J. Appl. Phys., 86 (1999).
 *
 * Values for enthalpy are difficult to compare against. Instead, the values
 * provided by the BrineFluidProperties UserObject were compared against
 * simple correlations, eg. from NIST sodium chloride data
 */
TEST_F(BrineFluidPropertiesBetaTest, haliteProperties)
{

Real p0, p1, p2, T0, T1, T2; //test variables.
const Real tol = REL_TOL_EXTERNAL_VALUE;
p0 = 30.0e6;
p1 = 60.0e6;
p2 = 80.0e6;
T0 = 300.0;
T1 = 500.0;
T2 = 700.0;

// Density
REL_TEST(_fp->halite_rho_from_p_T(p0, T0), 2167.88, tol);
REL_TEST(_fp->halite_rho_from_p_T(p1, T1), 2116.0, tol);
REL_TEST(_fp->halite_rho_from_p_T(p2, T2), 2056.8, tol);

// Test enthalpy at the triple point pressure of water
Real pt = 611.657; //test variable (tripple pressure point)
ABS_TEST(_fp->halite_h_from_p_T(pt, 273.16), 0.0, tol);
REL_TEST(_fp->halite_h_from_p_T(pt, 573.15), 271.13e3, tol);
REL_TEST(_fp->halite_h_from_p_T(pt, 673.15), 366.55e3, tol);
 }


/**
 * Verify the computed derivatives of the halite properties by comparing them with finite
 * differences

 TEST_F(BrineFluidPropertiesBetaTest, derivatives)
 {
   const Real tol = REL_TOL_DERIVATIVE;

   const Real p = 30e6;
   const Real T = 300;

   DERIV_TEST(_fp->halite_rho_from_p_T, p, T, tol);
   DERIV_TEST(_fp->halite_h_from_p_T, p, T, tol);
 }
 */
