#include "PorousFLowPermeabilityEmbeddedFractures.h"

registerMooseObject("PorousFlowApp", PorousFLowPermeabilityEmbeddedFractures);

InputParameters
PorousFLowPermeabilityEmbeddedFractures::validParams()
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
    params.addRangeCheckedParam<Real>("a",
                                      "a > 0",
                                      "Mean (scalar) fracture distance value");
    params.addRequiredParam<Real>("e0", "threshold strain");
    params.addRequiredParam<Real>("km", "matrix/intrinsic permeability");
    params.addRequiredParam<Real>("b0", "initial fracture aperture");
    params.addRequiredParam<Real>("rad_xy", "fracture rotation angle in radians");
    params.addRequiredParam<Real>("rad_yz", "fracture rotation angle in radians");
    params.addParam<RealVectorValue>("n",
                           "normal vector wrt to fracture surface");
    params.addParam<std::string>("base_name",
                             "Optional parameter that allows the user to define "
                             "multiple mechanics material systems on the same "
                             "block, i.e. for multiple phases");
    params.addParam<bool>("normal_vector_to_fracture_is_constant",
                      false,
                      "Whether the normal vector wrt to fracture surface is constant/known or not.");
  return params;
}

PorousFLowPermeabilityEmbeddedFractures::PorousFLowPermeabilityEmbeddedFractures(
    const InputParameters & parameters)
  : PorousFlowPermeabilityBase(parameters),
    _a(getParam<Real>("a")),
    _e0(getParam<Real>("e0")),
    _b0(getParam<Real>("b0")),
    _km(getParam<Real>("km")),
    _nVec(getParam<Real>("n")),
    _identity_two(RankTwoTensor::initIdentity),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _n_const(parameters.get<bool>("normal_vector_to_fracture_is_constant")),
    _stress(getMaterialProperty<RankTwoTensor>(_base_name + "stress")),
    _rad_xy(getParam<Real>("rad_xy")),
    _rad_yz(getParam<Real>("rad_yz")),
    _strain(getMaterialProperty<RankTwoTensor>("creep_strain"))
  //  _dvol_strain_qp_dvar(_mechanical ? &getMaterialProperty<std::vector<RealGradient>>(
  //                              "dPorousFlow_total_volumetric_strain_qp_dvar")
  //                                                   : nullptr),
{

// should be included if the derivatives/jacobian of this material is computed as well.
 //  _dictator.usePermDerivs(true);
}

void
PorousFLowPermeabilityEmbeddedFractures::computeQpProperties()
{
// This code block describes how the 'normal vector' (n) wrt the fracture face should
// be computed. if the components of n is known (e.g., sigma_xx, tau_xy and tau_zx),
// then it should be obtain from the input file. Otherwise, n is computed as the
// direction (eigenvector) of the third principal stress vector.

   if (_n_const)
       {
          RealVectorValue _n =_nVec;
       }
 // Eigenvectors were derived from the total stress obtained from the tensor mech. action.
 // Then, the eigenvector in the second column corresponding to the third principal stress
 // vector was computed.

      RankTwoTensor eigvec;
      std::vector<Real> eigvals;
      _stress[_qp].symmetricEigenvaluesEigenvectors(eigvals, eigvec);
      RealVectorValue _n = eigvec.column(2);

// The Z axis of the material is rotated around the X-Y plane.
    RankTwoTensor _transformation_matrix;

    _transformation_matrix(0, 0) = std::cos(_rad_xy);
    _transformation_matrix(0, 1) = std::sin(_rad_xy);
    _transformation_matrix(0, 2) = 0;
    _transformation_matrix(1, 0) = -std::sin(_rad_xy);
    _transformation_matrix(1, 1) = std::cos(_rad_xy);
    _transformation_matrix(1, 2) = 0;
    _transformation_matrix(2, 0) = 0;
    _transformation_matrix(2, 1) = 0;
    _transformation_matrix(2, 2) = 1;

// Similarly, X-axis is rotated around the Y-Z plane.
    _transformation_matrix(0, 0) = 1;
    _transformation_matrix(0, 1) = 0;
    _transformation_matrix(0, 2) = 0;
    _transformation_matrix(1, 0) = 0;
    _transformation_matrix(1, 1) = std::cos(_rad_yz);
    _transformation_matrix(1, 2) = std::sin(_rad_yz);
    _transformation_matrix(2, 0) = 0;
    _transformation_matrix(2, 1) = -std::sin(_rad_yz);
    _transformation_matrix(2, 2) = std::cos(_rad_yz);

    RealVectorValue rotMat_xy = _transformation_matrix * _n;
    RealVectorValue rotMat_yz = _transformation_matrix * _n;

  // finally, the normal vector is re-computed to always remain normal even during rotation
  // of the material (i.e., wrt to the rotated fracture plane (n_r)) by multiplying n with the
  // rotated material about the fix x-y and y-z axis (rotMat_xy and rotMat_yz)

    RealVectorValue n_r = rotMat_xy * rotMat_yz * _n;

  // strain in the normal fracture direction (i.e., strain dependent on the
  // fracture normal vector

   Real e_n = (_strain[_qp] * n_r)*(n_r);

  // H_de implements the macaulay-bracket in Zill et al. since _e0 is the initial/threshold
  // strain state of the material, and strain is always increasing in n-direction, e_n should
  // always be bigger than e_0. otherwise, H_de = 0
     Real H_de = (e_n > _e0) ? 1.0 : 0.0;

  // initial fracture aperture is sqrt(12 * k_m) in the literature
     Real _b0 = std::sqrt(12. * _km);

  // change in fracture aperture
     Real b_f = _b0 + (H_de * _a * (e_n - _e0));

     Real coeff = H_de * (b_f / _a) * ((b_f * b_f / 12.0) - _km);

     RankTwoTensor I = _identity_two;
     auto _M = RankTwoTensor::selfOuterProduct(n_r);

  // Finally, the permeability
    _permeability_qp[_qp] = (_km*I) + (coeff * (I - _M));

 // The computation of this material does not include its derivatives (jacobian)
}
