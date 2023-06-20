#include "ergsMineralCoeffAux.h"

registerMooseObject("PorousFlowApp", ergsMineralCoeffAux);

InputParameters
ergsMineralCoeffAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("satLIQUID", "liquid saturation");
  params.addRequiredCoupledVar("rm", "mineral precipitation coefficient");
  return params;
}

ergsMineralCoeffAux::ergsMineralCoeffAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _satLIQUID(coupledValue("satLIQUID")),
    _rm(coupledValue("rm")),
    _aperture(getMaterialProperty<Real>("initial_fracture_aperture_qp"))
{
}

Real
ergsMineralCoeffAux::computeValue()
{
  if (_satLIQUID[_qp] > 0.001)
  {
    return  _aperture[_qp] * _satLIQUID[_qp] * 1251 * _rm[_qp];
  }
    else
  return _aperture[_qp] *  1251 * _rm[_qp];
}
