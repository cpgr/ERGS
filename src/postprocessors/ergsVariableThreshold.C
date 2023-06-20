#include "ergsVariableThreshold.h"

registerMooseObject("PorousFlowApp", ergsVariableThreshold);

InputParameters
ergsVariableThreshold::validParams()
{
  InputParameters params = NodalPostprocessor::validParams();
  params.addRequiredCoupledVar("variable", "The name of the variable that shall be examined.");
  params.addRequiredParam<Real>("threshold", "A threshold value to compare variable value against.");
  params.addClassDescription(" A PostProcessor to obtain either 1 or 0 at any node when the value "
  " at that node exceeds a threshold value.");
  return params;
}

ergsVariableThreshold::ergsVariableThreshold(const InputParameters & parameters)
  : NodalPostprocessor(parameters),
  _variable(coupledValue("variable")),
  _threshold(getParam<Real>("threshold")),
  _flag(0)
{
}

 void
 ergsVariableThreshold::initialize()
 {
  _flag = 0.0;
 }

 void
 ergsVariableThreshold::finalize()
 {
  gatherMax(_flag);
 }

 void
ergsVariableThreshold::execute()
 {
  if (_variable[_qp] > _threshold)
    _flag = 1.0;
  else
    _flag = 0.0;
 }

 Real
ergsVariableThreshold::getValue()
 {
  return _flag;
 }

 void
ergsVariableThreshold::threadJoin(const UserObject & y)
 {
  const ergsVariableThreshold & pps = static_cast<const ergsVariableThreshold &>(y);
  if (pps._flag == 1.0)
   _flag = 1.0;
 }
