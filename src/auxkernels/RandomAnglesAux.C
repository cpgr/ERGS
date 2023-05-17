#include "RandomAnglesAux.h"

registerMooseObject("PorousFlowApp", RandomAnglesAux);

InputParameters
RandomAnglesAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  return params;
}

RandomAnglesAux::RandomAnglesAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _random_radXY(getMaterialProperty<Real>("random_xy_rotation_angle_for_each_element_qp"))
{
}

Real
RandomAnglesAux::computeValue()
{
  return _random_radXY[_qp];
}
