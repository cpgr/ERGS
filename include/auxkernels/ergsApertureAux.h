#pragma once

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class ergsApertureAux : public AuxKernel
{
public:
	static InputParameters validParams();
	ergsApertureAux(const InputParameters& parameters);

protected:
	virtual Real computeValue();
	const MaterialProperty<Real>& _aperture;
};

