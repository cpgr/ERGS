#pragma once

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class RandomAnglesAux : public AuxKernel
{
public:
	static InputParameters validParams();
	RandomAnglesAux(const InputParameters& parameters);

protected:
	virtual Real computeValue();
	const MaterialProperty<Real>& _random_radXY;
};

