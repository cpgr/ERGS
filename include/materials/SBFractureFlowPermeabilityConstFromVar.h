#pragma once

#include "PorousFlowPermeabilityBase.h"

/**
 * Material to provide permeability calculated from a variable representing a fracture aperture.
 * This material is primarily designed for use with heterogeneous fracture aperture models
 * where the aperture is provided by an
 * elemental aux variables that does not change.
 * The three diagonal entries corresponding to the x, y, and z directions
 * are assumed to be equal and calculated using the cubic law.
 */

class SBFractureFlowPermeabilityConstFromVar : public PorousFlowPermeabilityBase
{
public:
	static InputParameters validParams();

	SBFractureFlowPermeabilityConstFromVar(const InputParameters& parameters);

protected:
	void computeQpProperties() override;

	/// Fracture aperture
	const VariableValue& _aperture;
};