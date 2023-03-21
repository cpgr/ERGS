#pragma once

#include "AuxKernel.h"
#include "PorousFlowDictator.h"

/**
 * Computes a component of the Permeability Tensor:
 */
class EFPComponents : public AuxKernel
{
public:
	static InputParameters validParams();

	EFPComponents(const InputParameters& parameters);

protected:
	virtual Real computeValue();

	/// Permeability of porous material
	const MaterialProperty<RealTensorValue>& _EmbeddedFracturePermeability;

	/// PorousFlowDicatator UserObject
	const PorousFlowDictator& _dictator;

	/// Index of the fluid phase
	const unsigned int _ph;

	/// Desired spatial component
	unsigned int _i;
	unsigned int _j;
};
