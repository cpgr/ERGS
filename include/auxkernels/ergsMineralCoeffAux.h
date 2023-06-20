#pragma once

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class ergsMineralCoeffAux : public AuxKernel
{
public:
	static InputParameters validParams();
	ergsMineralCoeffAux(const InputParameters& parameters);

protected:
	virtual Real computeValue();

	const VariableValue& _satLIQUID;

	const VariableValue& _rm;

	const MaterialProperty<Real>& _aperture;
};

