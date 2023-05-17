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

#include "PorousFlowPermeabilityBase.h"
#include "RankTwoTensor.h"
#include <Eigen/Geometry>

/**
 * Material designed to provide the permeability tensor which is a function of
 * the computed strain. This permeability material is based on Embedded Fractures.
 * See Zill et. al.(2021): Hydro-mechanical continuum modelling of fluid percolation
 * through rock. The permeability is given as follows:
 *
 * k = (k_m * I_{ij}) + (b/a*[(b^2/12 - k_m)]*(I-M_{ij}))
 *
 * where b is the fracture aperture given by: b_0 + /Delta{b}
 * /Delta{b} depends on the strain (/epsilon) as follows:
 * /Delta{b} = a * 〈/epsilon_n -/epsilon_0〉. Here, /epsilon_0 is a threshold strain
 * and /epsilon_n is the computed strain normal to the embedded fracture plane/surface.
 * a is mean fracture distance, b is fracture aperture, K_m is matrix/intrinsic permeability.
 * I_{ij} is the identity tensor and M_{ij} is a structure tensor given as n⊗n. n is a vector normal
 * to the fracture.
 */
class PorousFlowEmbeddedFracturePermeability : public PorousFlowPermeabilityBase
{
public:
	static InputParameters validParams();

	PorousFlowEmbeddedFracturePermeability(const InputParameters& parameters);

protected:

	void computeQpProperties() override;

	/// optional parameter that allows multiple mechanics materials to be defined
	const std::string _base_name;

	/// mean fracture distance
	const Real _a;

	/// Initial fracture aperture
	const Real _b0;

	/// Threshold strain
	const Real _e0;

	/// matrix/intrinsic permeability
	const Real _km;

	/// whether normal vector to fracture is constant or not
	bool _n_const;

	/// whether to use random rotation angle or fix
	bool _Random_field;

	/// normal vector to fracture surface
	RealVectorValue _nVec;

	/// string to hold the name of the fracture rotation angle around xy
	std::string _phi_xy;

	/// string to hold the name of the fracture rotation angle around yz
	std::string _phi_yz;

	/// Structure tensor M_{ij} = n⊗n {not used}
	const RankTwoTensor _M;

	///  Identity RankTwotensor I_{ij}
	const RankTwoTensor _identity_two;

	/// get the stress tensor
	const MaterialProperty<RankTwoTensor>& _stress;

	/// get the strain tensor (actually, this is the creep strain)
	const MaterialProperty<RankTwoTensor>& _strain;

	/// fracture rotation angles about xy and yz (in radians)
	const Real& _fix_rad_xy;
	const Real& _fix_rad_yz;

	/// random rotation_angle_for_each_element
	MaterialProperty<Real>& _randm_rad_xy;
	MaterialProperty<Real>& _randm_rad_yz;

	/// Computed strain in the fracture normal vector direction as a material property
	MaterialProperty<Real>& _en;

	const VariableValue& _rotXY;
	const VariableValue& _rotYZ;

};