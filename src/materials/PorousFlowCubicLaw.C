#include "PorousFlowCubicLaw.h"


registerMooseObject("PorousFlowApp", PorousFlowCubicLaw);

InputParameters
PorousFlowCubicLaw::validParams()
{
  InputParameters params = PorousFlowPermeabilityBase::validParams();
  params.addClassDescription(
    " This permeability material calculates the permeability based on the cubic Law .");
//    params.addRequiredCoupledVar("Aperture", "Aperture as an auxvariable.");
  return params;
}

PorousFlowCubicLaw::PorousFlowCubicLaw(const InputParameters & parameters)
  : PorousFlowPermeabilityBase(parameters),
    _aperture(getMaterialProperty<Real>("initial_fracture_aperture_qp"))
//  _aperture(coupledValue("Aperture"))
{
}


void
PorousFlowCubicLaw::computeQpProperties()
{
     _permeability_qp[_qp] = std::pow(_aperture[_qp],3.0)/12.0;
}
