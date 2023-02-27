/******************************************************************************/
/*         PERGS - Permeability for Enhanced RockSalt Geothermal Systems      */
/*                                                                            */
/*          Copyright (C) 2022 by Ishmael Dominic Yevugah                     */
/*      University of Manitoba, Price Faculty of Engineering                  */
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

#include "PFEMBase.h"

/**
 * Material designed to provide the permeability tensor which is a function of
 * the computed strain. This permeability material is based on Embedded Orthotropic Fractures.
 * See Zill et. al.(2021): Hydro-mechanical continuum modelling of fluid percolation
 * through rock. The permeability is given as follows:
 *
 * k = k_m*I +  [b_i/a_i * (\frac{b_{i}^2}{12} - k_m)*(I-M_i)].
 *
 * where i=summation Over number of fracture planes/surfaces. b is the fracture
 *  aperture given by: b_{i0} + /Delta{b_i} /Delta{b_i} depends on the strain (/epsilon) as follows:
 * /Delta{b_i} = a_i * 〈/epsilon :M_i - /epsilon_{0i}〉. Here, /epsilon_{0i} is a threshold strain of
 * the material in each of the fracture normal vector direction. /epsilon is total strain
 * a_i and b_i are fracture distance and fracture aperture respectively in each fracture
 * normal vector direction. K_m is the matrix/intrinsic permeability. I_{ij} is the identity tensor
 * and M_{ij} is a structure tensor given as n_i⊗n_i. n_i is a vector normal to each fracture plane.
 */

class PFOrthoEM : public PFEMBase
{
public:
	static InputParameters validParams();

	PFOrthoEM(const InputParameters& parameters);

protected:
	void computeQpProperties() ;

	/// mean fracture distance in all 3 directions
    std::vector<double> _alpha;

	/// Threshold strain in all 3 directions
     std::vector<double> _eps;

	/// Tensor that holds the normal vector to the fracture surface in each direction
	RealTensorValue _NVec;
};
