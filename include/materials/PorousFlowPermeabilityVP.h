#pragma once

#include "PorousFlowPermeabilityBase.h"

/**
 * Material designed to provide a constant permeability tensor based on Verma and Pruess (1988)
 */
class PorousFlowPermeabilityVP : public PorousFlowPermeabilityBase
{
public:
  static InputParameters validParams();

  PorousFlowPermeabilityVP(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// Constant value of permeability tensor
  const RealTensorValue _input_initial_permeability;
  
  /// solid saturation (Ss)
  const VariableValue& _solid_sat;

  // independent geometric parameters in Verma and Pruess
  const Real _Gamma;
  const Real _theta_r;
  //Const
};
