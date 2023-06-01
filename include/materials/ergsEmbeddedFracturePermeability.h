/******************************************************************************/
/*         PERGS - Permeability for Enhanced RockSalt Geothermal Systems      */
/*                                                                            */
/*          Copyright (C) 2022 by Ishmael Dominic Yevugah                     */
/*        University of Manitoba, Price Faculty of Engineering                */
/*                                                                            */
/*        Special Thanks to Guillaume Giudicelli, Chris Green                 */
/*        and the rest of the Moose Team for helping on the model             */
/*                                                                            */
/*       This program is free software: you can redistribute it and/or modify */
/*    it under the terms of the GNU General Public License as published by    */
/*      the Free Software Foundation, either version 3 of the License, or     */
/*                     (at your option) any later version.                    */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program.  If not, see <http://www.gnu.org/licenses/>    */
/******************************************************************************/

#pragma once

#include "PorousFlowEmbeddedFracturePermeability.h"
#include "RankTwoTensor.h"
#include <Eigen/Geometry>

/**
 * Derived material class from PorousFlowEmbeddedFracturePermeability that obtains the initial 
 * fracture aperture as a coupled variable instead of an ordinary parameter. This initial
 * fracture aperture affects the permeability.
 */
class ergsEmbeddedFracturePermeability : public PorousFlowEmbeddedFracturePermeability
{
public:
	static InputParameters validParams();

	ergsEmbeddedFracturePermeability(const InputParameters& parameters);

	virtual void initQpStatefulProperties();

	void computeQpProperties() override;

protected:

	/// initial fracture aperture
	MaterialProperty<Real>& _b;
	const MaterialProperty<Real>& _b_old;

	/// liquid saturation
	const VariableValue& _satLIQUID;

	/// salt mass fraction
	const VariableValue& _Xnacl;

	/// salt solubility limit
	const Real _XEQ;

	/// mineral precipitation rate coefficient
	const VariableValue& _rm;
	//	const Real _rm;

	const VariableValue& _Dt;
};