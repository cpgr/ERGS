#pragma once

#include "PorousFlowPermeabilityBase.h"
#include "RankTwoTensor.h"
#include <Eigen/Geometry>

/**
 * Material designed to provide the permeability based on cubic law
 */
class PorousFlowCubicLaw : public PorousFlowPermeabilityBase
{
public:
	static InputParameters validParams();

	PorousFlowCubicLaw(const InputParameters& parameters);

protected:
	void computeQpProperties() override;

	/// Computed strain in the fracture normal vector direction as a material property
	const MaterialProperty<Real>& _aperture;

//	const VariableValue& _aperture;
};