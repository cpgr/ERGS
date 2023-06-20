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
/*           Time was heavily invested in these codes. Please, show           */
/*                 support by properly referencing it. Thanks.                */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program.  If not, see <http://www.gnu.org/licenses/>    */
/******************************************************************************/



#include "PorousFlowEmbeddedOrthotropicFracturePermeability.h"

registerMooseObject("PorousFlowApp", PorousFlowEmbeddedOrthotropicFracturePermeability);

InputParameters
PorousFlowEmbeddedOrthotropicFracturePermeability::validParams()
{
  InputParameters params = PorousFlowEmbeddedFracturePermeability::validParams();
  params.addClassDescription(
    " This permeability material calculates the permeability tensor based on attributes of the "
    " Orthotropic Embedded Fractures in Zill et. al.(2021): Hydro-mechanical continuum modelling of "
    " fluid percolation through rock. The permeability is given as follows: "
    " k = k_m*I +  [b_i/a_i * (\frac{b_{i}^2}{12} - k_m)*(I-M_i)]. where i=summation Over number of "
    " fracture planes/surfaces. b is the fracture aperture given by: b_{i0} + /Delta{b_i} "
    " /Delta{b_i} depends on the strain (/epsilon) as follows "
    " /Delta{b_i} = a_i * 〈/epsilon :M_i - /epsilon_{0i}〉. Here, /epsilon_{0i} is a threshold strain of "
    " the material in each of the fracture normal vector direction. /epsilon is total strain."
    " a_i and b_i are fracture distance and fracture aperture respectively in each fracture "
    " normal vector direction. K_m is the matrix/intrinsic permeability. I_{ij} is the identity tensor"
    " and M_{ij} is a structure tensor given as n_i⊗n_i. n_i is a vector normal to each fracture plane.");

    params.addRequiredParam<std::vector<double>>("alpha", "Mean fracture distance value in all 3 directions");
    params.addRequiredParam<std::vector<double>>("eps0", "threshold strain");
    params.addParam<RealTensorValue>("N",
                           "normal vector wrt to fracture surface");
  return params;
}

PorousFlowEmbeddedOrthotropicFracturePermeability::PorousFlowEmbeddedOrthotropicFracturePermeability(
    const InputParameters & parameters)
  : PorousFlowEmbeddedFracturePermeability(parameters),
    _alpha(getParam<std::vector<double>>("alpha")),
    _eps(getParam<std::vector<double>>("eps0")),
    _NVec(parameters.isParamValid("N")
                  ? getParam<RealTensorValue>("N")
                  : RealTensorValue(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
{
}


void
PorousFlowEmbeddedOrthotropicFracturePermeability::computeQpProperties()
{
// This code block describes how the 'normal vector' (n) wrt each (of the 3) fracture face
// should be computed. if the components of n is known (e.g., sigma_xx, tau_xy and tau_zx),
// then it should be specify in the input file. Otherwise, n is computed as the
// direction (eigenvector) of the all the three principal stresses. The assumption here is
// that the three fracture planes lie within the principal stresses.

    RankTwoTensor _n;

    if (_n_const)
      {
         _n =_NVec;
      }
    else
    // Eigenvectors were derived from the total stress obtained from the tensor mech. action.
    // Then, the eigenvector in each column, corresponding to the direction of each principal
    // stress vector was computed.
      {
        RankTwoTensor eigvec;
        std::vector<Real> eigvals;
        _stress[_qp].symmetricEigenvaluesEigenvectors(eigvals, eigvec);
       _n = eigvec;
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

  // The fracture normal vectors captured in the Tensor (_n) are rotated around both Z-axis
  // (i.e., X-Y plane) and X-axis (i.e., Y-Z plane) during random rotation of the material using
  // the  rotation matrices below. (See Zill et al. for why material is randomly rotated)

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

 // The permeability is computed by first initializing it:
    RankTwoTensor I = _identity_two;
    _permeability_qp[_qp] = _km*I;

 // The final/total permeability is the summation over the permeability due to each
 // individual strain corresponding to its unique rotated fracture normal vector.
 // Note that each column of the _n tensor corresponding to the fracture normal vector
 // refers to each principal stress direction.

 for (int i = 0; i < 3; i++)
  {
  // rotate each fracture normal vector captured in the _n tensor
   RealVectorValue n_r = rotMat_xy * rotMat_yz * _n.column(i);

  // strain due to each fracture normal vector direction
    _en[_qp] = n_r*(_strain[_qp] * n_r);

  // The heaviside function (H_de) that implements the macaulay-bracket in Zill et al.
   Real H_de = (_en[_qp] > _eps[i]) ? 1.0 : 0.0;

  // change in fracture aperture
   Real b_f = _b0 + (H_de * _alpha[i] * (_en[_qp] - _eps[i]));

   Real coeff = H_de * (b_f / _alpha[i]) * ((b_f * b_f / 12.0) - _km);

   RankTwoTensor I = _identity_two;

   auto _M = RankTwoTensor::selfOuterProduct(n_r);


   _permeability_qp[_qp] += coeff * (I - _M);
 }
}
