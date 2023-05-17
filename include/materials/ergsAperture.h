#pragma once

#include "PorousFlowMaterial.h"

/**
 * The ergsAperture material computes the initial fracture aperture (b0) based on
 * a lagged approach consisting of the old state of the initial fracture aperture (b0_old),
 * mineral precipitation rate coefficient (rm), salt mass fraction (Xnacl), salt solubility 
 * limit (XEQ) and liquid saturation (satLIQUID).
 */

class ergsAperture : public PorousFlowMaterial
{
public:
  static InputParameters validParams();

  ergsAperture(const InputParameters & parameters);

  virtual void initQpStatefulProperties();
 // virtual void propagateQpStatefulProperties();


protected:
	virtual void computeQpProperties() override;

	/// initial fracture aperture
	MaterialProperty<Real>& _b;
	const MaterialProperty<Real>& _b_old;

	/// liquid saturation
	const VariableValue& _satLIQUID;

	/// salt mass fraction
	const VariableValue& _Xnacl;

	/// initial aperture
	const VariableValue& _aperture;

	/// salt solubility limit
	const Real _XEQ;

	/// mineral precipitation rate coefficient
	const VariableValue& _rm;
//	const Real _rm;

	const VariableValue& _Dt;
};
