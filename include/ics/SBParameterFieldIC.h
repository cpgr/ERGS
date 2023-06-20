#pragma once

#include "InitialCondition.h"
#include "Material.h"
#include "VariableField.h"

/* 
class InputParameters;
class SBParameterFieldIC;
namespace libMesh
{
	class Point;
}

template <typename T>
InputParameters validParams();

template <>
InputParameters validParams<SBParameterFieldIC>();
*/

class SBParameterFieldIC : public InitialCondition
{
public:
	static InputParameters validParams();

	SBParameterFieldIC(const InputParameters& parameters);
	
	virtual Real value(const Point& p) override;

private:
	VariableField _var_field;
};
