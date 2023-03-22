#pragma once

#include "PorousFlowBrineCO2.h"

//class BrineFluidProperties;
//class SinglePhaseFluidProperties;
//class Water97FluidProperties;

/**
 * Specialized class for brine and CO2 including calculation of mutual
 * solubility of the two fluids using a high-accuracy fugacity-based formulation.
 *
 * For temperatures 12C <= T <= 99C, the formulation is based on
 * Spycher, Pruess and Ennis-King, CO2-H2O mixtures in the geological
 * sequestration of CO2. I. Assessment and calculation of mutual solubilities
 * from 12 to 100C and up to 600 bar, Geochimica et Cosmochimica Acta, 67, 3015-3031 (2003),
 * and Spycher and Pruess, CO2-H2O mixtures in the geological sequestration of CO2. II.
 * Partitioning in chloride brine at 12-100C and up to 600 bar, Geochimica et
 * Cosmochimica Acta, 69, 3309-3320 (2005).
 *
 * For temperatures 109C <= T <= 300C, the formulation is based on
 * Spycher and Pruess, A Phase-Partitioning Model for CO2-Brine Mixtures at Elevated
 * Temperatures and Pressures: Application to CO2-Enhanced Geothermal Systems,
 * Transport in Porous Media, 82, 173-196 (2010)
 *
 * As the two formulations do not coincide at temperatures near 100C, a cubic
 * polynomial is used in the intermediate temperature range 99C < T < 109C to
 * provide a smooth transition from the two formulations in this region.
 *
 * Notation convention
 * Throughout this class, both mole fractions and mass fractions will be used.
 * The following notation will be used:
 * yk: mole fraction of component k in the gas phase
 * xk: mole fraction of component k in the liquid phase
 * Yk: mass fraction of component k in the gas phase
 * Xk: mass fraction of component k in the liquid phase
 */
class PorousFlowBrineSaltCO2 : public PorousFlowBrineCO2
{
public:
  static InputParameters validParams();

  PorousFlowBrineSaltCO2(const InputParameters & parameters);

  virtual void thermophysicalProperties(Real pressure,
                                        Real temperature,
                                        Real Xnacl,
                                        Real Z,
                                        unsigned int qp,
                                        std::vector<FluidStateProperties> & fsp) const override;

  /**
   * Mass fractions of CO2 in brine and water vapor in CO2 at equilibrium
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param[out] Xco2 mass fraction of CO2 in liquid (kg/kg)
   * @param[out] Yh2o mass fraction of H2O in gas (kg/kg)
   */
  void equilibriumMassFractions(const DualReal & pressure,
                                const DualReal & temperature,
                                const DualReal & Xnacl,
                                DualReal & Xco2,
                                DualReal & Yh2o) const;

  /**
   * Mass fractions of CO2 and H2O in both phases, as well as derivatives wrt
   * PorousFlow variables. Values depend on the phase state (liquid, gas or two phase)
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param Z total mass fraction of CO2 component
   * @param[out] PhaseStateEnum current phase state
   * @param[out] FluidStateProperties data structure
   */
  void massFractions(const DualReal & pressure,
                     const DualReal & temperature,
                     const DualReal & Xnacl,
                     const DualReal & Z,
                     FluidStatePhaseEnum & phase_state,
                     std::vector<FluidStateProperties> & fsp) const;

  /**
   * Thermophysical properties of the gaseous state
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (K)
   * @param[out] FluidStateProperties data structure
   */
  void gasProperties(const DualReal & pressure,
                     const DualReal & temperature,
                     std::vector<FluidStateProperties> & fsp) const;

  /**
   * Thermophysical properties of the liquid state
   *
   * @param pressure liquid pressure (Pa)
   * @param temperature temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param[out] FluidStateProperties data structure
   */
  void liquidProperties(const DualReal & pressure,
                        const DualReal & temperature,
                        const DualReal & Xnacl,
                        std::vector<FluidStateProperties> & fsp) const;

  /**
   * Thermophysical properties of the solid state (pure Halite)
   *
   * @param pressure solid pressure (Pa)
   * @param temperature temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param[out] FluidStateProperties data structure
   */
  void solidProperties(const DualReal& pressure,
      const DualReal& temperature,
      std::vector<FluidStateProperties>& fsp) const;

  /**
   * Gas saturation in the two-phase region
   *
   * @param pressure gas pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param Z total mass fraction of CO2 component
   * @param FluidStateProperties data structure
   * @return gas saturation (-)
   */
  DualReal saturationGAS(const DualReal & pressure,
                      const DualReal & temperature,
                      const DualReal & Xnacl,
                      const DualReal & Z,
                      std::vector<FluidStateProperties> & fsp) const;

  /**
   * Gas and liquid properties in the two-phase region
   *
   * @param pressure gas pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param Z total mass fraction of NCG component
   * @param qp quadpoint for capillary presssure
   * @param[out] FluidStateProperties data structure
   */
  void MultiPhaseProperties(const DualReal & pressure,
                          const DualReal & temperature,
                          const DualReal & Xnacl,
                          const DualReal & Z,
                          unsigned int qp,
                          std::vector<FluidStateProperties> & fsp) const;

  /// Fluid properties UserObject for brine
  const BrineFluidPropertiesBeta & _brine_fp;
  /// Fluid properties UserObject for water
  const SinglePhaseFluidProperties & _water_fp;

  /// Solid phase index
  const unsigned int _solid_phase_number;

  /// initial non-zero salt saturation
  const Real _saturationSOLID;
};
