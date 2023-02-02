#pragma once

#include "PorousFlowWaterNCG.h"

class PorousFlowWaterAir : public PorousFlowWaterNCG
{
public:
  static InputParameters validParams();

  PorousFlowWaterAir(const InputParameters & parameters);

 
    virtual std::string fluidStateName() const override;

    void thermophysicalProperties(Real pressure,
        Real temperature,
        Real Xnacl,
        Real Z,
        unsigned int qp,
        std::vector<FluidStateProperties>& fsp) const override;

    void equilibriumMassFractions(const DualReal& pressure,
        const DualReal& temperature,
        DualReal& Xair,
        DualReal& Yh2o) const;

    void massFractions(const DualReal& pressure,
        const DualReal& temperature,
        const DualReal& Z,
        FluidStatePhaseEnum& phase_state,
        std::vector<FluidStateProperties>& fsp) const;

    void gasProperties(const DualReal& pressure,
        const DualReal& temperature,
        std::vector<FluidStateProperties>& fsp) const;

    void liquidProperties(const DualReal& pressure,
        const DualReal& temperature,
        std::vector<FluidStateProperties>& fsp) const;

    DualReal liquidDensity(const DualReal& pressure, const DualReal& temperature) const;

    DualReal gasDensity(const DualReal& pressure,
        const DualReal& temperature,
        std::vector<FluidStateProperties>& fsp) const;

    DualReal saturation(const DualReal& pressure,
        const DualReal& temperature,
        const DualReal& Z,
        std::vector<FluidStateProperties>& fsp) const;

    void twoPhaseProperties(const DualReal& pressure,
        const DualReal& temperature,
        const DualReal& Z,
        unsigned int qp,
        std::vector<FluidStateProperties>& fsp) const;

    DualReal enthalpyOfDissolution(const DualReal& temperature) const;

    virtual Real totalMassFraction(
        Real pressure, Real temperature, Real Xnacl, Real saturation, unsigned int qp) const override;

protected:
    DualReal moleFractionToMassFraction(const DualReal& xmol) const;

    void checkVariables(Real temperature) const;
  
    
    /// Fluid properties UserObject for the Air (In this case, air is Ideal gas)
    const SinglePhaseFluidProperties& _air_fp;
    /// Molar mass of air (i.e., Ideal gas (kg/mol))
    const Real _Mair;
    /// Henry's coefficients for Air (Ideal gas)
    const std::vector<Real> _air_henry;
    /// Universal Gas Constant
    const Real _R;



    using PorousFlowWaterNCG::_Mh2o;
    using PorousFlowWaterNCG::_water_triple_temperature;
    using PorousFlowWaterNCG::_water_critical_temperature;
    using PorousFlowWaterNCG::_water_fp;
    using PorousFlowWaterNCG::_water97_fp;
  };

