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
    params.addRequiredParam<RealEigenMatrix>("n",
                           "normal vector wrt to fracture surface");
  // params.addRequiredParam<Eigen::Matrix<double, 3, 1>>("n","normal vector wrt to fracture surface");
    params.addParam<std::string>("base_name",
                             "Optional parameter that allows the user to define "
                             "multiple mechanics material systems on the same "
                             "block, i.e. for multiple phases");
    params.addParam<std::string>("fracture_rotation_xy",
                            "name for the input that hold the angle (in radius)"
                            "by which the fracture normal is rotated in the xy-plane");
    params.addParam<std::string>("fracture_rotation_yz",
                             "name for the input that hold the angle (in radius)"
                             "by which the fracture normal is rotated in the yz-plane");
    params.addParam<bool>("normal_vector_to_fracture_is_constant",
                      false,
                      "Whether the normal vector wrt to fracture surface is constant or not.");
  return params;
}

PorousFLowPermeabilityEmbeddedFractures::PorousFLowPermeabilityEmbeddedFractures(
    const InputParameters & parameters)
  : PorousFlowPermeabilityBase(parameters),
    _a(getParam<Real>("a")),
    _e0(getParam<Real>("e0")),
    _b0(getParam<Real>("b0")),
    _km(getParam<Real>("km")),
    _identity_two(RankTwoTensor::initIdentity),
    _n(getParam<RealEigenMatrix>("n")),
    //_n(getParam<Eigen::Matrix<double, 3, 1>>("n")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _n_const(parameters.get<bool>("normal_vector_to_fracture_is_constant")),
    _stress(getMaterialProperty<RankTwoTensor>(_base_name + "stress")),
    _phi_xy(getParam<std::string>("fracture_rotation_xy")),
    _phi_yz(getParam<std::string>("fracture_rotation_yz")),
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
  // This code block describes how the 'normal vector' (n) wrt the fracture face
  // should be computed. if the components of n is known (e.g., sigma_xx, tau_xy and tau_zx),
  // then it should be obtain from the input file. Otherwise, n is computed as the
  // direction (eigenvector) of the third principal stress vector.
  Eigen::Matrix<double, 3, 1> const n = [&]
   {
       if (_n_const)
       {
           return _n;
       }

       // Here, eigenvectors were derived from the total stress obtained from
       // the tensor mech. action. Then, the eigenvector in the second column
       // corresponding to the third principal stress vector was computed.
       // Similar to line 69 to 72 of the OGS code.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> e_s(_stress);
        return (Eigen::Matrix<double, 3, 1>)e_s.eigenvectors().col(2);
   }();

   // Here, the Z axis of the material is rotated around fixed X-Y axis based on
   // quad point (_qp) and time (_dt). Both _qp and _dt are passed as an argument to _phi_xy.
   // Same is true for X axis.
    Real const rotMat_xy =
    Eigen::AngleAxisd(_phi_xy(_dt, _qp)[0], Eigen::Vector3d::UnitZ());

    Real const rotMat_yz =
    Eigen::AngleAxisd(_phi_yz(_dt, _qp)[0], Eigen::Vector3d::UnitX());

    // finally, the normal vector is re-computed to always remain normal even during rotation
    // of the material (i.e., wrt to the rotated fracture plane (n_r)) by multiplying n with the
    // rotated material about the fix x-y and y-z axis (rotMat_xy and rotMat_yz)
     Eigen::Matrix<double, 3, 1> const n_r = rotMat_yz * (rotMat_xy * _n);

    // Here, the strain was obtained and multiplied by the normal vector. see line 85 of OGS code.
    // Not sure the way I obtained the strain is the best.
     Real e_n = (_strain * n_r).dot(n_r.transpose());

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
     RankTwoTensor _M = n_r * n_r.transpose();

     // Finally, the permeability
     _permeability_qp[_qp] = (_km*I) + (coeff * (I - M));

 // The computation of this material does not include its derivatives (jacobian)
}
