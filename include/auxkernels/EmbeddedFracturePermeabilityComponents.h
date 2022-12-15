#pragma once

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "PorousFlowDictator.h"

class EmbeddedFracturePermeabilityComponents : public AuxKernel
{
public:
	static InputParameters validParams();
	EmbeddedFracturePermeabilityComponents(const InputParameters& parameters);

protected:
	virtual Real computeValue();

	const MaterialProperty<RankTwoTensor>& _EmbeddedFracturePermeability;
	const unsigned int _i;
	const unsigned int _j;


    /// PorousFlowDicatator UserObject
    const PorousFlowDictator& _dictator;

	/// Index of the fluid phase
	const unsigned int _ph;
};

