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

#include "PorousFlowEmbeddedOrthotropicFracturePermeability.h"

/**
 * Alternative version of PorousFlowEmbeddedOrthotropicFracturePermeability. This material obtains the
 * Automatic-Differentiation of the total stress and total creep strain.
 */

class ADPorousFlowEmbeddedOrthotropicFracturePermeability : public PorousFlowEmbeddedOrthotropicFracturePermeability
{
public:
	static InputParameters validParams();

	ADPorousFlowEmbeddedOrthotropicFracturePermeability(const InputParameters& parameters);

protected:
  void computeQpProperties() override;

  /// get the AD of total stress tensor
	const ADMaterialProperty<RankTwoTensor>& _stress;

	/// get the AD of total creep strain tensor (actually, this is the creep strain)
	const ADMaterialProperty<RankTwoTensor>& _strain;
};
