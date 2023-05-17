#pragma once

#include "AuxKernel.h"
#include "VariableField.h"

class SBRotatedParameterFieldAux : public AuxKernel
{
public:
	static InputParameters validParams();
	SBRotatedParameterFieldAux(const InputParameters& parameters);
	
	virtual ~SBRotatedParameterFieldAux() {};

protected:
	virtual Real computeValue();

private:
	Real _x_min = std::numeric_limits<Real>::max();
	Real _x_max = -std::numeric_limits<Real>::max();
	Real _y_min = std::numeric_limits<Real>::max();
	Real _y_max = -std::numeric_limits<Real>::max();

	VariableField _var_field;
	Real _block_id;
	Real _alphaX;
	Real _alphaY;
	Point _orig;
	Real _x_scale;
	Real _y_scale;
	Real _scale;
	Point _mesh_origin;

};