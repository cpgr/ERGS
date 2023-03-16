#pragma once

#include "IdealGasFluidProperties.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

/**
 * Air fluid properties. This class is essentially Ideal gas with the addition of
 * Henry's Coefficient. Default parameters are for air at atmospheric pressure and temperature
 */
class IdealGasAirFluidProperties : public IdealGasFluidProperties
{
public:
	static InputParameters validParams();

	IdealGasAirFluidProperties(const InputParameters& parameters);

	virtual std::vector<Real> henryCoefficients() const override;
};

#pragma GCC diagnostic pop