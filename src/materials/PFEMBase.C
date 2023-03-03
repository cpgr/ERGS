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
/*         Time was heavily invested in this codes. Please, show              */
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

#include "PFEMBase.h"


registerMooseObject("PorousFlowApp", PFEMBase);

InputParameters
PFEMBase::validParams()
{
  InputParameters params = PorousFlowPermeabilityBase::validParams();
  params.addClassDescription(
    " This permeability material calculates the permeability tensor based on attributes of "
    " Embedded Fractures. See Zill et. al.(2021): Hydro-mechanical continuum modelling of "
    " fluid percolation through rock. The permeability is given as follows: "
    " k = (k_m * I_{ij}) + (b/a*[(b^2/12 - k_m)]*(I-M_{ij})) "
    " where b is the fracture aperture given by: b_0 + /Delta{b} "
    " /Delta{b} depends on the strain (/epsilon) as follows "
    " /Delta{b} = a * 〈/epsilon_n -/epsilon_0〉. Here, /epsilon_0 is a threshold strain "
    " and /epsilon_n is the computed strain normal to the embedded fracture direction."
    " a is mean fracture distance, b is fracture aperture, K_m is matrix/intrinsic permeability."
    " I_{ij} is the identity tensor and M_{ij} is a structure tensor given as n⊗n. n is a vector normal"
    " to the fracture.");
    params.addRangeCheckedParam<Real>("a", 1,
                                      "a > 0",
                                      "Mean (scalar) fracture distance value");
    params.addParam<Real>("e0", 1,"threshold strain");
    params.addParam<Real>("b0", 0,"Initial fracture aperture");
    params.addParam<Real>("km", 1, "matrix/intrinsic permeability");
    params.addParam<Real>("fix_rad_xy", 0, "fix fracture rotation angle in radians");
    params.addParam<Real>("fix_rad_yz", 0, "fix fracture rotation angle in radians");
    params.addParam<RealVectorValue>("n",
                           "normal vector wrt to fracture surface");
    params.addParam<std::string>("base_name",
                             "Optional parameter that allows the user to define "
                             "multiple mechanics material systems on the same "
                             "block, i.e. for multiple phases");
    params.addParam<bool>("normal_vector_to_fracture_is_constant",
                      true,
                    "Whether the normal vector wrt to fracture surface is constant/known or not.");
    params.addParam<bool>("Random_field",
                      true,
                      "Whether to use spatially random angle of rotation at each timestep or a fixed"
                      "angle of rotation. Set true if you want random field, false otherwise");
    params.addRequiredCoupledVar("rotation_angleXY", "The random rotation field_XY.");
    params.addRequiredCoupledVar("rotation_angleYZ", "The random rotation field_YZ.");
  return params;
}

PFEMBase::PFEMBase(
    const InputParameters & parameters)
  : PorousFlowPermeabilityBase(parameters),
    _a(getParam<Real>("a")),
    _b0(getParam<Real>("b0")),
    _e0(getParam<Real>("e0")),
    _km(getParam<Real>("km")),
    _nVec(parameters.isParamValid("n")
                  ? getParam<RealVectorValue>("n")
                  : RealVectorValue(1.0, 0.0, 0.0)),
    _identity_two(RankTwoTensor::initIdentity),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _n_const(parameters.get<bool>("normal_vector_to_fracture_is_constant")),
    _stress(getMaterialProperty<RankTwoTensor>(_base_name + "stress")),
    _fix_rad_xy(getParam<Real>("fix_rad_xy")),
    _fix_rad_yz(getParam<Real>("fix_rad_yz")),
    _strain(getMaterialProperty<RankTwoTensor>("total_strain")),
    _en (_nodal_material ? declareProperty<Real>("fracture_normal_strain_nodal")
                              : declareProperty<Real>("fracture_normal_strain_qp")),
     _Random_field(parameters.get<bool>("Random_field")),
     _randm_rad_xy(_nodal_material ? declareProperty<Real>("random_xy_rotation_angle_for_each_element")
                              : declareProperty<Real>("random_xy_rotation_angle_for_each_element_qp")),
     _randm_rad_yz(_nodal_material ? declareProperty<Real>("random_yz_rotation_angle_for_each_element")
                              : declareProperty<Real>("random_yz_rotation_angle_for_each_element_qp")),
     _rotXY(coupledValue("rotation_angleXY")),
     _rotYZ(coupledValue("rotation_angleYZ"))
{
}


void
PFEMBase::computeQpProperties()
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
        _stress[_qp].symmetricEigenvaluesEigenvectors(eigvals, eigvec);
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
    _en[_qp] = n_r *(_strain[_qp] * n_r);

  // H_de is the heaviside function that implements the macaulay-bracket in Zill et al.
  // since _e0 is the initial/threshold strain state of the material, and strain is always
  // increasing in n-direction, _en should always be bigger than _e0. otherwise, change in
  // fracture aperture = 0

     Real H_de = (_en[_qp] > _e0) ? 1.0 : 0.0;

  // initial fracture aperture is sqrt(12 * k_m) in the literature
  //   Real _b0 = std::sqrt(12. * _km);

  // change in fracture aperture
     Real b_f = _b0 + (H_de * _a * (_en[_qp] - _e0));

     Real coeff =  H_de * (b_f / _a) * ((b_f * b_f / 12.0) - _km);

     RankTwoTensor I = _identity_two;
     auto _M = RankTwoTensor::selfOuterProduct(n_r);

  // Finally, the permeability
     _permeability_qp[_qp] = (_km*I) + (coeff *(I -_M));
}
