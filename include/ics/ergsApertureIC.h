#pragma once

#include "InitialCondition.h"
#include "Material.h"

class SinglePhaseFluidProperties;

/**
 * ergsApertureIC calculates an initial value for fracture aperture (b) based on 
 * a lagged approach consisting of the old state of the initial fracture aperture (b_old),
 * mineral precipitation rate coefficient (rm), salt mass fraction (Xnacl), salt solubility limit (XEQ) 
 * and liquid saturation (satLIQUID).
 */

class ergsApertureIC : public InitialCondition
{
public:
  static InputParameters validParams();

  ergsApertureIC(const InputParameters & parameters);

  virtual Real value(const Point & p) override;

protected:

  /// initial fracture aperture
  const MaterialProperty<Real> & _b0;
  


};
