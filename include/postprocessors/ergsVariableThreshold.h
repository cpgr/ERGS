#pragma once

#include "NodalPostprocessor.h"
#include "MooseVariableInterface.h"

/**
 * A PostProcessor to obtain either 1 or 0 at any node when the value at that 
 * node exceeds a threshold value.
 */

class ergsVariableThreshold : public NodalPostprocessor
{
public:
  static InputParameters validParams();

  ergsVariableThreshold(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void finalize() override;
  virtual void execute() override;
  virtual Real getValue() override;
  virtual void threadJoin(const UserObject& y) override;

protected:
	const VariableValue& _variable;
	const Real& _threshold;
	Real _flag;

};
