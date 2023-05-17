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
/*         Time was heavily invested in these codes. Please, show              */
/*               support by properly referencing it. Thanks.                  */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program.  If not, see <http://www.gnu.org/licenses/>    */
/******************************************************************************/

#include "ADergsEmbeddedFracturePermeability.h"


registerMooseObject("ergsApp", ADergsEmbeddedFracturePermeability);

InputParameters
ADergsEmbeddedFracturePermeability::validParams()
{
  InputParameters params = ergsEmbeddedFracturePermeability::validParams();
  params.addClassDescription(
    " Alternative version of ergsEmbeddedFracturePermeability. This material obtains the"
    " Automatic Differentiation of the total stress and total creep strain.");
  return params;
}


ADergsEmbeddedFracturePermeability::ADergsEmbeddedFracturePermeability(const InputParameters & parameters)
  : ergsEmbeddedFracturePermeability(parameters),
   _stress(getADMaterialProperty<RankTwoTensor>(_base_name + "stress")),
   _strain(getADMaterialProperty<RankTwoTensor>("creep_strain"))
{
}


void
ADergsEmbeddedFracturePermeability::computeQpProperties()
{
// This code block describes how the 'normal vector' (n) wrt the fracture face should
// be computed. if the components of n is known (e.g., sigma_xx, tau_xy and tau_zx),
// then it should be specify and obtain from the input file. Otherwise, n is computed as the
// direction (eigenvector) of the principal stress vector corresponding to the normal direction.

      RealVectorValue _n;
      if (_n_const)
     {
         _n =_nVec;
      }
        else
          // Eigenvectors were derived from the total stress obtained from the tensor mech. action.
          // Then, the eigenvector in the second column corresponding to the third principal stress
          // vector was computed.
      {
        RankTwoTensor eigvec;
        std::vector<Real> eigvals;
        raw_value(_stress[_qp]).symmetricEigenvaluesEigenvectors(eigvals, eigvec);
        _n = eigvec.column(0);
      }


  if (_Random_field)
  // Get the spatially random rotation angle (in radians) for each element at each timestep either
  // from a specific distribution or randomly.
      {
       _randm_rad_xy[_qp] = _rotXY[_qp];
      _randm_rad_yz[_qp] = _rotYZ[_qp];
      }

      else
  // Get the fixed (or user-specified) rotation angle for all elements in the domain at each timestep.
      {
        _randm_rad_xy[_qp] = _fix_rad_xy;
        _randm_rad_yz[_qp] = _fix_rad_yz;
      }

  // The fracture normal vector (n) is rotated around the Z-axis (i.e., X-Y plane) during random
  // rotation of the material using the correct rotation matrix (rotMat_xy and angle rad_xy).
  // (See Zill et al., 2022 for why material is randomly rotated).

    RankTwoTensor rotMat_xy;
    rotMat_xy(0, 0) = std::cos(_randm_rad_xy[_qp]);
    rotMat_xy(0, 1) = -std::sin(_randm_rad_xy[_qp]);
    rotMat_xy(0, 2) = 0;
    rotMat_xy(1, 0) = std::sin(_randm_rad_xy[_qp]);
    rotMat_xy(1, 1) = std::cos(_randm_rad_xy[_qp]);
    rotMat_xy(1, 2) = 0;
    rotMat_xy(2, 0) = 0;
    rotMat_xy(2, 1) = 0;
    rotMat_xy(2, 2) = 1;

// Similarly, the fracture normal is rotated around the X-axis (i.e., Y-Z plane) during random
// rotation of the material using the correct rotation matrix (rotMat_yz and angle rad_xy).

     RankTwoTensor rotMat_yz;
     rotMat_yz(0, 0) = 1;
     rotMat_yz(0, 1) = 0;
     rotMat_yz(0, 2) = 0;
     rotMat_yz(1, 0) = 0;
     rotMat_yz(1, 1) = std::cos(_randm_rad_yz[_qp]);
     rotMat_yz(1, 2) = -std::sin(_randm_rad_yz[_qp]);
     rotMat_yz(2, 0) = 0;
     rotMat_yz(2, 1) = std::sin(_randm_rad_yz[_qp]);
     rotMat_yz(2, 2) = std::cos(_randm_rad_yz[_qp]);

  // Rotation of the  fracture normal vector using the rotation matrices
     RealVectorValue n_r = rotMat_xy * rotMat_yz * _n;

  // strain in the normal fracture direction
    _en[_qp] = n_r *(raw_value(_strain[_qp]) * n_r);

  // H_de is the heaviside function that implements the macaulay-bracket in Zill et al.
  // since _e0 is the initial/threshold strain state of the material, and strain is always
  // increasing in n-direction, _en should always be bigger than _e0. otherwise, change in
  // fracture aperture = 0

     Real H_de = (_en[_qp] > _e0) ? 1.0 : 0.0;

  // initial fracture aperture is sqrt(12 * k_m) in the literature
  //   Real _b0 = std::sqrt(12. * _km);

  // change in fracture aperture
     Real b_f = _b0evol[_qp] + (H_de * _a * (_en[_qp] - _e0));

     Real coeff =  H_de * (b_f / _a) * ((b_f * b_f / 12.0) - _km);

     RankTwoTensor I = _identity_two;
     auto _M = RankTwoTensor::selfOuterProduct(n_r);

  // Finally, the permeability
     _permeability_qp[_qp] = (_km*I) + (coeff *(I -_M));
}
