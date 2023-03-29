#pragma once

#include "BrineFluidProperties.h"
#include "Water97FluidProperties.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

/**
 * Updated the Previous BrineFluidProperties to include halite
 * solubility in water according to Driesner and Heinrich (2007) and
 *the vapor phaseaccording to Palliser and McKibbin (1998).
 **/
class BrineFluidPropertiesBeta : public BrineFluidProperties
{
public:
  static InputParameters validParams();

  BrineFluidPropertiesBeta(const InputParameters & parameters);

  /**
   * Updated the Solubility of halite (solid NaCl) in water
   * from Potter et al., (i.e., haliteSolubility above) to Driesner and Heinrich (2007):
   * The system H2O�NaCl. Part I pg 4891. This new solubility is referred to as
   * haliteSolubilityWater
   *
   * @param temperature temperature (K)
   * @return halite solubility (kg/kg)
   */
  Real haliteSolubilityWater(Real temperature, Real pressure) const;

  DualReal haliteSolubilityWater(DualReal temperature, DualReal pressure) const;

  /**
   *  Halite (solid NaCl) solubility in the gas phase from Palliser and McKibbin (1998):
   * A Model for Deep Geothermal Brines, I: T-p-X State-Space Description pg 78
   *
   * @param temperature temperature (K)
   * @return halite solubility (kg/kg)
   */
  Real haliteSolubilityGas(/*Real temperature,*/ Real pressure) const;

  DualReal haliteSolubilityGas(/*DualReal temperature,*/ DualReal pressure) const;
  /**
   *  Halite (solid NaCl) property, density (rho) and its derivative from
   *  Driesner and Heinrich (2007)
   **/
  virtual Real halite_rho_from_p_T(Real pressure, Real temperature) const;

  virtual DualReal halite_rho_from_p_T(DualReal pressure, DualReal temperature) const;

  virtual void halite_rho_from_p_T(
      Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const;

  /**
   *  Halite (solid NaCl) property, enthalpy (h) and its derivative from
   *  Driesner and Heinrich (2007).
   **/
  virtual Real halite_h_from_p_T(Real pressure, Real temperature) const;

  virtual DualReal halite_h_from_p_T(DualReal pressure, DualReal temperature) const;

  virtual void
  halite_h_from_p_T(Real pressure, Real temperature, Real & h, Real & dh_dp, Real & dh_dT) const;

protected:
  // Parameter for Driesner
  Real _alpha;

  // Parameters for Palliser
  Real _y;
  Real _XgSat;
  Real _Psat;

  // Params for Nacl
  const Real _T_triple;
  const Real _p_triple;
};

#pragma GCC diagnostic pop
