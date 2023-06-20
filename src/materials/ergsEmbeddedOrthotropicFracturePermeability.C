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


#include "ergsEmbeddedOrthotropicFracturePermeability.h"

registerMooseObject("PorousFlowApp", ergsEmbeddedOrthotropicFracturePermeability);

InputParameters
ergsEmbeddedOrthotropicFracturePermeability::validParams()
{
  InputParameters params = PorousFlowEmbeddedOrthotropicFracturePermeability::validParams();
  params.addClassDescription(
  " Derived material class from PorousFlowEmbeddedOrthotropicFracturePermeability that computes"
  " the permeability using extra aperture change due to halite dissolution. This extra aperture"
  " change is computed from various coupled variables.");
  params.addParam<Real>("rho_w", 1250, "water density");
  params.addParam<Real>("rho_m", 2170, "Density of halite mineral");
/*
  params.addParam<Real>("r_plus", 0.277, "Halite dissolution rate in an infinitely dilute solution");
  params.addParam<Real>("R", 0.277, "Gas constant");
  params.addParam<Real>("T", 0.277, "Absolute temperature");
  params.addParam<Real>("sig", 0.277, "Temkin's average stoichiometric number");
  params.addParam<Real>("K", 0.277, "Thermodynamic equilibrium constant");
  params.addParam<Real>("gamma_Na", 0.277, "Thermodynamic equilibrium constant");
  params.addParam<Real>("gamma_Cl", 0.277, "Thermodynamic equilibrium constant");
  params.addParam<Real>("c_Na", 0.277, "concentration of Sodium ions");
  params.addParam<Real>("c_Cl", 0.277, "concentration of Chloride");
  params.addRequiredCoupledVar("aperture_old", "initial aperture field");
*/
  params.addRequiredCoupledVar("sw", "water saturation");
  params.addRequiredCoupledVar("Xnacl", "Fluid temperature");
  params.addRequiredCoupledVar("rm", "mineral precipitation coefficient");
  params.addRequiredCoupledVar("Dt", "time step");
  params.addParam<Real>("XEQ", 0.277, "solubility limit");
  return params;
}

ergsEmbeddedOrthotropicFracturePermeability::ergsEmbeddedOrthotropicFracturePermeability(
    const InputParameters & parameters)
  : PorousFlowEmbeddedOrthotropicFracturePermeability(parameters),
      _b(_nodal_material
                       ? declareProperty<Real>("initial_fracture_aperture_nodal")
                       : declareProperty<Real>("initial_fracture_aperture_qp")),
      _b_old(_nodal_material
                       ? getMaterialPropertyOld<Real>("initial_fracture_aperture_nodal")
                       : getMaterialPropertyOld<Real>("initial_fracture_aperture_qp")),
      _rho_w(getParam<Real>("rho_w")),
      _rho_m(getParam<Real>("rho_m")),
/*
      _r_plus(getParam<Real>("r_plus")),
      _R(getParam<Real>("R")),
      _T(getParam<Real>("T")),
      _sig(getParam<Real>("sig")),
      _K(getParam<Real>("K")),
      _gamma_Na(getParam<Real>("gamma_Na")),
      _gamma_Cl(getParam<Real>("gamma_Cl")),
      _c_Na(getParam<Real>("c_Na")),
      _c_Cl(getParam<Real>("c_Cl")),
      _aperture_old(coupledValue("aperture_old")),
*/
      _sw(coupledValue("sw")),
      _Xnacl(coupledValue("Xnacl")),
      _rm(coupledValue("rm")),
      _Dt(coupledValue("Dt")),
      _XEQ(getParam<Real>("XEQ"))
{
}


void
ergsEmbeddedOrthotropicFracturePermeability::initQpStatefulProperties()
{
 _b[_qp] = 0.0;
}


void
ergsEmbeddedOrthotropicFracturePermeability::computeQpProperties()
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

/*
  /* compute the halite dissolution rate according to Seales et. al. (2016)*/
    // activity of the various ions (i.e., Na and Cl)
//    Real a_Na = _gamma_Na * c_Na;
//    Real a_Cl = _gamma_Cl * c_Cl;
    // compute the reaction ionic activity product
//    Real Q = a_Na + a_Cl;
    // compute the chemical affinity
//    Real A = -(_R *_T)*std::log(Q/_K);
    //finally, the reaction rate
//    Real r = _r_plus * (1-std::exp(-A/(_sigT*_R*_T)));


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

  // change in fracture aperture according to initial aperture and strain
   Real b_f = _b0 + (H_de * _alpha[i] * (_en[_qp] - _eps[i]));

  // final aperture evolution, accounting for the halite dissolution
   _b[_qp] = b_f + (_b_old[_qp] * (1-( 1 * _sw[_qp] * (_rho_w/_rho_m) * _rm[_qp] * (_Xnacl[_qp] - _XEQ) * /*_dt*/ _Dt[_qp])));
//   _b[_qp] = b_f + (_b_old[_qp] * (1-( 1 * _sw[_qp] * (_rho_w/_rho_m) * _r * (_XEQ-_Xnacl[_qp]) * _dt /*_Dt[_qp]*/ )));

   Real coeff = H_de * (_b[_qp] / _alpha[i]) * ((_b[_qp] * _b[_qp] / 12.0) - _km);

   RankTwoTensor I = _identity_two;

   auto _M = RankTwoTensor::selfOuterProduct(n_r);


   _permeability_qp[_qp] += coeff * (I - _M);
 }
}
