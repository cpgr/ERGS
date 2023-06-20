#pragma once

#include "PorousFlowPorosityConst.h"

/**
 * Material that provides an 'effective' constant porosity value. It is a modification of the
 * PorousFlowPorosityConst material. The main modification is the addition of the solid-phase
 * saturation from a mixture containing liquid, gas and solid phases. The formulation is as 
 * follows: phi_{eff} = phi * (1 - Ss). Where Ss is the solid-phase saturation, phi and phi_{eff}
 * are the total and effective porosity values, respectively. 
 * 
 * NOTE: The total porosity can be specified either by a constant value in the input file,  
 * or taken from an aux variable. Ss is an auxvariable and affected by the dissolution/precipitation
 * of the solid phase. This material assumes that the total porosity remains constant throughout the
 *  simulation, so the coupled aux variable porosity must also remain constant.
 */
class PorousFlowEffectivePorosityConst : public PorousFlowPorosityConst
{
public:
  static InputParameters validParams();

  PorousFlowEffectivePorosityConst(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// solid saturation (Ss)
  const VariableValue & _solid_sat;
};
