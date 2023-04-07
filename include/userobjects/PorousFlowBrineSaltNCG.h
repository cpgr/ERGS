#pragma once

#include "PorousFlowBrineSaltCO2.h"

class BrineFluidPropertiesBeta;

class PorousFlowBrineSaltNCG : public PorousFlowBrineSaltCO2
{
public:
  static InputParameters validParams();

  PorousFlowBrineSaltNCG(const InputParameters & parameters);

  virtual void thermophysicalProperties(Real pressure,
                                        Real temperature,
                                        Real Xnacl,
                                        Real Z,
                                        unsigned int qp,
                                        std::vector<FluidStateProperties> & fsp) const override;

  void equilibriumMassFractions(const DualReal & pressure,
                                const DualReal & temperature,
                                const DualReal & Xnacl,
                                const DualReal& Z,
                                DualReal & Xncg,
                                DualReal & Yh2o,
                                DualReal& Ynacl,
                                DualReal& Snacl,
      std::vector<FluidStateProperties>& fsp) const;

  void massFractions(const DualReal & pressure,
                     const DualReal & temperature,
                     const DualReal & Xnacl,
                     const DualReal & Z,
                     FluidStatePhaseEnum & phase_state,
                     std::vector<FluidStateProperties> & fsp) const;

  void gasProperties(const DualReal & pressure,
                     const DualReal & temperature,
                     const DualReal& Xnacl,
                    const DualReal& Z,
                     std::vector<FluidStateProperties> & fsp) const;

  void liquidProperties(const DualReal & pressure,
                        const DualReal & temperature,
                        const DualReal & Xnacl,
                        std::vector<FluidStateProperties> & fsp) const;

  void solidProperties(const DualReal& pressure,
      const DualReal& temperature,
      std::vector<FluidStateProperties>& fsp) const;

  DualReal saturationGAS(const DualReal & pressure,
                      const DualReal & temperature,
                      const DualReal & Xnacl,
                      const DualReal & Z,
                      std::vector<FluidStateProperties> & fsp) const;

  DualReal saturationSOLID(const DualReal& pressure,
      const DualReal& temperature,
      const DualReal& Xnacl,
      std::vector<FluidStateProperties>& fsp) const;

  void MultiPhaseProperties(const DualReal & pressure,
                          const DualReal & temperature,
                          const DualReal & Xnacl,
                          const DualReal & Z,
                          unsigned int qp,
                          std::vector<FluidStateProperties> & fsp) const;

  unsigned int solidPhaseIndex() const { return _solid_phase_number; };


  /**
   * Henry's constant of dissolution of gas phase NCG in brine. From
   * Battistelli et al, A fluid property module for the TOUGH2 simulator for saline brines
   * with non-condensible gas, Proc. Eighteenth Workshop on Geothermal Reservoir Engineering (1993)
   *
   * @param temperature fluid temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @return Henry's constant (Pa)
   */

  DualReal henryConstant(const DualReal& temperature, const DualReal& Xnacl) const;

  /**
 * Enthalpy of dissolution of gas phase CO2 in brine calculated using Henry's constant
 * From Himmelblau, Partial molal heats and entropies of solution for gases dissolved
 * in water from the freezing to the near critical point, J. Phys. Chem. 63 (1959).
 * Correction due to salinity from Battistelli et al, A fluid property module for the
 * TOUGH2 simulator for saline brines with non-condensible gas, Proc. Eighteenth Workshop
 * on Geothermal Reservoir Engineering (1993).
 *
 * @param temperature fluid temperature (K)
 * @param Xnacl NaCl mass fraction (kg/kg)
 * @return enthalpy of dissolution (J/kg)
 */

  DualReal enthalpyOfDissolutionNCG(const DualReal& temperature, 
                                     const DualReal& Xnacl) const;
  /**
 * Partial density of dissolved CO2
 * From Garcia, Density of aqueous solutions of CO2, LBNL-49023 (2001)
 *
 * @param temperature fluid temperature (K)
 * @return partial molar density (kg/m^3)
 */
  DualReal partialDensityNCG(const DualReal& temperature) const;

  /**
* Vapor Pressure Lowering (VPL) effect
* From Garcia, Density of aqueous solutions of CO2, LBNL-49023 (2001)
*
* @param temperature fluid temperature (K)
* @return VPL factor (-)
*/
  DualReal f_vpl(const DualReal& pressure,
      const DualReal& temperature,
      const DualReal& Xnacl,
      const DualReal& Z,
      std::vector<FluidStateProperties>& fsp) const;

  protected:
/* 
  /// Fluid properties UserObject for brine 
  const BrineFluidPropertiesBeta& _brine_fp;

  /// Fluid properties UserObject for water
  const SinglePhaseFluidProperties & _water_fp;

  /// Fluid properties UserObject for water (used to access Henry's law)
  const Water97FluidProperties& _water97_fp;

  /// Solid phase index
  const unsigned int _solid_phase_number;
 */

  /// Fluid properties UserObject for the NCG
  const SinglePhaseFluidProperties& _ncg_fp;
  /// Molar mass of water (kg/mol)
  const Real _Mh2o;
  /// Molar mass of non-condensable gas (kg/mol)
  const Real _Mncg;
  /// Triple point temperature of water (K)
  const Real _water_triple_temperature;
  /// Critical temperature of water (K)
  const Real _water_critical_temperature;
  /// Henry's coefficients for the NCG
  const std::vector<Real> _ncg_henry;
  /// Universal Gas Constant
  const Real _R;
};
