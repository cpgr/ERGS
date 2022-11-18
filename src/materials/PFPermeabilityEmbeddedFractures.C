#include "PFPermeabilityEmbeddedFractures.h"

registerMooseObject("PorousFlowApp", PFPermeabilityEmbeddedFractures);

InputParameters
PFPermeabilityEmbeddedFractures::validParams()
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
    params.addParam<std::string>("base_name",
                             "Optional parameter that allows the user to define "
                             "multiple mechanics material systems on the same "
                             "block, i.e. for multiple phases");
  return params;
}

PFPermeabilityEmbeddedFractures::PFPermeabilityEmbeddedFractures(
    const InputParameters & parameters)
  : PorousFlowPermeabilityBase(parameters),
    _a(getParam<Real>("a")),
    _e0(getParam<Real>("e0")),
    _b0(getParam<Real>("b0")),
    _km(getParam<Real>("km")),
    _identity_two(RankTwoTensor::initIdentity),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _strain(getMaterialProperty<RankTwoTensor>("creep_strain"))
{

 // should be included if the derivatives/jacobian of this material is computed as well. 
 //  _dictator.usePermDerivs(true);
}

void
PFPermeabilityEmbeddedFractures::computeQpProperties()
{
    VectorValue<Real> _n(0.0, 0.0, 1.0);
    RankTwoTensor _M;
    _M.vectorOuterProduct(_n,_n);

    Real e_n = (_strain[_qp] * _n).dot(_n.transpose());

    // H_de implements the macaulay-bracket. Since _e0 is the initial/threshold
    // strain state of the material, and strain is always increasing in n-direction,
    // e_n should always be bigger than e_0. otherwise, H_de = 0
     Real H_de = (e_n > _e0) ? 1.0 : 0.0;

     // initial fracture aperture is sqrt(12 * k_m) in the literature
     Real _b0 = std::sqrt(12. * _km);

     // change in fracture aperture
     Real b_f = _b0 + (H_de * _a * (e_n - _e0));

     Real coeff = H_de * (b_f / _a) * ((b_f * b_f / 12.0) - _km);

     RankTwoTensor I = _identity_two;

     // Finally, the permeability
     _permeability_qp[_qp] = (_km*I) + (coeff * (I - _M));

 // The computation of this material does not include its derivatives (jacobian)
}
