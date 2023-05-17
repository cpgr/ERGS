#include "ergsApertureAux.h"

registerMooseObject("PorousFlowApp", ergsApertureAux);

InputParameters
ergsApertureAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  return params;
}

ergsApertureAux::ergsApertureAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _aperture(getMaterialProperty<Real>("initial_fracture_aperture_qp"))
{
}

Real
ergsApertureAux::computeValue()
{
  return _aperture[_qp];
}
