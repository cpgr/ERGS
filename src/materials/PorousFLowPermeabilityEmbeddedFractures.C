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
    // 'n' has to be obtained as a 3x1 vector normal to the fracture surface from the
    //  input file but i'm not sure if this next line of code can do that.
    params.addRequiredParam<Eigen::Matrix<double, 3, 1>>("n","fracture_normal");

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
    _km(getParam<Real>("km")),
    _identity_two(RankTwoTensor::initIdentity),

    // see the comment above. not sure of this next line of code
    _n(getParam<Eigen::Matrix<double, 3, 1>>("n")),

    _n_const(parameters.get<bool>("normal_vector_to_fracture_is_constant")),
    _vol_strain_qp(_mechanical ? &getMaterialProperty<Real>("PorousFlow_total_volumetric_strain_qp")
                                               : nullptr),
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
  // not entirely sure how to go about the implementation here. the code block below is supposed
  // to be implemented as a lambda function, see line 63 to 73 of the opengeosys (OGS) code. The code block describes
  // how the 'normal vector' (n) wrt the fracture face should be computed. if n is constant, then
  // it should be obtain from the input file, otherwise, it is computed by obtaining the applied stress as shown below.
  Eigen::Matrix<double, 3, 1> const n = [&]
   {
       if (_n_const)
       {
           return _n;
       }
       auto const sigma = // I need to obtain the applied stress and convert it to
                          // a symmetric tensor here, but don't know how to obtain the stress.
                          // After obtaining the symmetric stress tensor, the 'Eigen library' functions
                          // was used to obtain a column 2 of its eigen vectors as 'n'. See line 69 to
                          // 72 of the OGS code.

   }();

   // here, the Z axis vector is rotated together with the material around fixed X-Y axis using
   // (am assuming) a position vector (pos) and time (t) passed as an argument to _phi_xy on line 76.
   // of the OGS code. I think these arguments should be in the InputParameters list (I might be wrong).
    Real const rotMat_xy =
    Eigen::AngleAxisd(_phi_xy(t, pos)[0], Eigen::Vector3d::UnitZ());

    // Similarly, here, the X axis vector is rotated together with the material around fixed Y-Z axis using
    // pos and t as arguments.
    Real const rotMat_yz =
    Eigen::AngleAxisd(_phi_yz(t, pos)[0], Eigen::Vector3d::UnitX());

    // finally, the normal vector is re-computed to always remain normal even during rotation
    // of the material (i.e., wrt to the rotated fracture plane (n_r)) by multiplying n with the
    // rotated material about the fix x-y and y-z axis (rotMat_xy and rotMat_yz)
     Eigen::Matrix<double, 3, 1> const n_r = rotMat_yz * (rotMat_xy * n);

    // Here, the strain was obtained and multiplied by the normal vector. see line 85 of OGS code.
    // I'm not sure whether the way I obtained the strain is the best.
     Real e_n = (_vol_strain_qp * n_r).dot(n_r.transpose());

    // H_de implements the macaulay-bracket in Zill et al. since _e0 is the initial/threshold
    // strain state of the material, and strain is always increasing in n-direction, e_n should
    // always be bigger than e_0. otherwise, H_de = 0
     Real H_de = (e_n > _e0) ? 1.0 : 0.0;

     // initial fracture aperture is sqrt(12 * k_m) in the literature
     Real _b0 = std::sqrt(12. * _km);

     // change in fracture aperture
     Real b_f = _b0 + H_de * _a * (e_n - _e0);

     Real coeff = H_de * (b_f / _a) * ((b_f * b_f / 12.0) - _k);

     // not sure how to compute the tensor I here. either one of these two might suffix
     Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
     RankTwoTensor I = _identity_two;


     // finally, the permeability
     _permeability_qp[_qp] =
     (_km*I) + (coeff * (I - n_r * n_r.transpose())
      (_k_anisotropy * _km*_identity_two) + (b/_a)*(std::pow(b,2)/12.-_km)*(_identity_two);


 // The computation of this material does not include its derivatives (jacobian)
  }
}
