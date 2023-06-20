#include "ergsApertureIC.h"
/*
#include "SinglePhaseFluidProperties.h"
#include <MaterialProperty.h>
*/
registerMooseObject("PorousFlowApp", ergsApertureIC);

InputParameters
ergsApertureIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addClassDescription("An IC to compute the initial fracture aperture (b0)"
                              "which will serve as an auxvariable");
  return params;
}

ergsApertureIC::ergsApertureIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _b0(getMaterialProperty<Real>("initial_fracture_aperture_qp"))
{
}


Real
ergsApertureIC::value(const Point & /*p*/)
{

  return _b0[_qp];
}
