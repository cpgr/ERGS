#include "PorousFlowEffectivePorosityConst.h"

registerMooseObject("PorousFlowApp", PorousFlowEffectivePorosityConst);

InputParameters
PorousFlowEffectivePorosityConst::validParams()
{
  InputParameters params = PorousFlowPorosityConst::validParams();
  params.addRequiredCoupledVar(
      "solid_sat",
      "The saturation of the solid-phase in a three (multi) phase mixture of liquid-gas-solid."
      "This value is a constant monomial variable (not a linear lagrange or any variable).");
  params.addClassDescription("This Material computes the 'effective' porosity as follows: "
                             "phi_{eff} = phi * (1 - Ss). Where Ss is the solid-phase saturation,"
                             "phi and phi_{eff} are the total and effective porosity values,"
                             "respectively. The total porosity is assumed to be constant");
  return params;
}

PorousFlowEffectivePorosityConst::PorousFlowEffectivePorosityConst(const InputParameters & parameters)
  : PorousFlowPorosityConst(parameters),
   _solid_sat(coupledValue("solid_sat"))
{
}

void
PorousFlowEffectivePorosityConst::initQpStatefulProperties()
{
  // note the [0] below: _phi0 is a constant monomial and we use [0] regardless of _nodal_material
  _porosity[_qp] = _input_porosity[0] * (1.0 - _solid_sat[0]);
//   _console << "solid_sat = " << _solid_sat[0] << std::endl;
}

void
PorousFlowEffectivePorosityConst::computeQpProperties()
{
  initQpStatefulProperties();

  // The derivatives are zero for all time
  _dporosity_dvar[_qp].assign(_num_var, 0.0);
  _dporosity_dgradvar[_qp].assign(_num_var, RealGradient());
}
