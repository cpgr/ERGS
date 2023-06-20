#include "ergsAperture.h"

registerMooseObject("PorousFlowApp", ergsAperture);

InputParameters
ergsAperture::validParams()
{
  InputParameters params = PorousFlowMaterial::validParams();
  params.addRequiredCoupledVar("satLIQUID", "liquid saturation");
  params.addRequiredCoupledVar("Xnacl", "Fluid temperature");
  params.addRequiredCoupledVar("rm", "mineral precipitation coefficient");
  params.addRequiredCoupledVar("Dt", "time step");
  params.addParam<Real>("XEQ", 0.277, "solubility limit");
  params.addCoupledVar("initial_fracture_aperture_input", "aperture field or aperture variable");         //
  params.addClassDescription("An IC to compute the initial fracture aperture (b0)"
                              "based on a lagged approach");
  return params;
}

ergsAperture::ergsAperture(const InputParameters & parameters)
  : PorousFlowMaterial(parameters),

    _b(_nodal_material
                     ? declareProperty<Real>("initial_fracture_aperture_nodal")
                     : declareProperty<Real>("initial_fracture_aperture_qp")),
    _b_old(_nodal_material
                     ? getMaterialPropertyOld<Real>("initial_fracture_aperture_nodal")
                     : getMaterialPropertyOld<Real>("initial_fracture_aperture_qp")),

    _satLIQUID(coupledValue("satLIQUID")),
    _Xnacl(coupledValue("Xnacl")),
    _rm(coupledValue("rm")),
    _Dt(coupledValue("Dt")),
    _XEQ(getParam<Real>("XEQ")),
  // _aperture_old(coupledValueOld("initial_fracture_aperture_input")),
    _aperture(coupledValue("initial_fracture_aperture_input"))
  //  _XEQ(getParam<Real>("temperature_unit") == 0 ? 0.0 : 273.15)
{
}


 void
 ergsAperture::initQpStatefulProperties()
 {
  _b[_qp] = _aperture[_qp];
 }


void
 ergsAperture::computeQpProperties()
 {
   if (_rm[_qp] == 0.0)
   {
     _b[_qp] = _b_old[_qp];  //_aperture_old[_qp];  //
   }
   else
  _b[_qp] = _b_old[_qp] - ( _b_old[_qp] * _satLIQUID[_qp] * 0.5765 * /*_dt*/ _Dt[_qp] * _rm[_qp] * (_Xnacl[_qp] -_XEQ));
//  _b[_qp] = _aperture_old[_qp] - ( _aperture_old[_qp] * _satLIQUID[_qp] * 0.5765 * /*_dt*/ _Dt[_qp] * _rm[_qp] * (_Xnacl[_qp] -_XEQ));
}
