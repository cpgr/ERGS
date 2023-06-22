#include "PorousFlowPermeabilityVP.h"

registerMooseObject("PorousFlowApp", PorousFlowPermeabilityVP);

InputParameters
PorousFlowPermeabilityVP::validParams()
{
  InputParameters params = PorousFlowPermeabilityBase::validParams();
  params.addRequiredCoupledVar(
      "solid_sat",
      "The saturation of the solid-phase in a three (multi) phase mixture of liquid-gas-solid."
      "This value is a constant monomial variable (not a linear lagrange or any variable).");
  params.addParam<Real>("Gamma", 0.8, "The fractional length of the pore bodies");
  params.addParam<Real>("theta_r", 0.8, "The fraction of original porosity at which"
                        "permeability is reduced to zero");
  params.addRequiredParam<RealTensorValue>("permeability0",
      "The initial permeability tensor (usually in m^2)");
  params.addClassDescription(
      "This Material calculates the permeability tensor based on Verma and Pruess (1988)");
  return params;
}

PorousFlowPermeabilityVP::PorousFlowPermeabilityVP(const InputParameters & parameters)
  : PorousFlowPermeabilityBase(parameters),
   _solid_sat(coupledValue("solid_sat")),
   _Gamma(getParam<Real>("Gamma")),
   _theta_r(getParam<Real>("theta_r")),
   _input_initial_permeability(getParam<RealTensorValue>("permeability0"))
{
}

void
PorousFlowPermeabilityVP::computeQpProperties()
{
   Real theta = (1.0 - _solid_sat[0] - _theta_r)/(1.0 -_theta_r);
   Real omega = 1.0 + ((1/_Gamma)/((1.0/_theta_r) - 1.0));
   Real numer = std::pow(theta,2.0) * 1.0 - _Gamma + (_Gamma/std::pow(omega, 2.0));
   Real deno = 1.0 - _Gamma + _Gamma * (std::pow((theta/(theta + omega - 1.0)), 2.0));
  _permeability_qp[_qp] = _input_initial_permeability * (numer/deno);
  (*_dpermeability_qp_dvar)[_qp].assign(_num_var, RealTensorValue());
  (*_dpermeability_qp_dgradvar)[_qp].resize(LIBMESH_DIM);
  for (unsigned i = 0; i < LIBMESH_DIM; ++i)
    (*_dpermeability_qp_dgradvar)[_qp][i].assign(_num_var, RealTensorValue());
}
