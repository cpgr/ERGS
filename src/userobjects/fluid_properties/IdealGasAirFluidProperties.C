#include "IdealGasAirFluidProperties.h"
#include "Conversion.h"

#include "metaphysicl/raw_type.h"

registerMooseObject("FluidPropertiesApp", IdealGasAirFluidProperties);

InputParameters
IdealGasAirFluidProperties::validParams()
{
  InputParameters params = IdealGasFluidProperties::validParams();

  return params;
}

IdealGasAirFluidProperties::IdealGasAirFluidProperties(const InputParameters & parameters)
  : IdealGasFluidProperties(parameters){}

std::vector<Real>
IdealGasAirFluidProperties::henryCoefficients() const
{
  return {-8.55445, 4.01195, 9.52345};
}
