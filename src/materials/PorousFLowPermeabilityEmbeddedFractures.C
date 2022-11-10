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
    params.addRequiredParam<Real>("eps0", "threshold strain");
    params.addRequiredParam<Real>("km", "matrix/intrinsic permeability");
    params.addParam<RealTensorValue>("k_anisotropy",
                                     "A tensor to multiply the calculated scalar "
                                     "permeability, in order to obtain anisotropy if "
                                     "required. Defaults to isotropic permeability "
                                     "if not specified.");
  return params;
}

PorousFLowPermeabilityEmbeddedFractures::PorousFLowPermeabilityEmbeddedFractures(
    const InputParameters & parameters)
  : PorousFlowPermeabilityBase(parameters),
    _a(getParam<Real>("a")),
    _eps0(getParam<Real>("eps0")),
    _km(getParam<Real>("km")),
    _identity_two(RankTwoTensor::initIdentity),
    _k_anisotropy(parameters.isParamValid("k_anisotropy")
                      ? getParam<RealTensorValue>("k_anisotropy")
                      : RealTensorValue(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)),
    _vol_strain_qp(_mechanical ? &getMaterialProperty<Real>("PorousFlow_total_volumetric_strain_qp")
                                               : nullptr),
    _dvol_strain_qp_dvar(_mechanical ? &getMaterialProperty<std::vector<RealGradient>>(
                                "dPorousFlow_total_volumetric_strain_qp_dvar")
                                                     : nullptr),
{

 // use derivatives in the Jacobian calculations
  _dictator.usePermDerivs(true);
}

void
PorousFLowPermeabilityEmbeddedFractures::computeQpProperties()
{
  Real b = std::sqrt(12. *_km);
  _permeability_qp[_qp] =
      (_k_anisotropy * _km*_identity_two) + (b/_a)*(std::pow(b,2)/12.-_km)*(_identity_two);

  _dpermeability_qp_dvar[_qp].resize(_num_var, RealTensorValue());
  for (unsigned int v = 0; v < _num_var; ++v)
    _dpermeability_qp_dvar[_qp][v] = _dporosity_qp_dvar[_qp][v] * _permeability_qp[_qp] *
                                     (_n / _porosity_qp[_qp] + _m / (1.0 - _porosity_qp[_qp]));

  _dpermeability_qp_dgradvar[_qp].resize(LIBMESH_DIM);
  for (unsigned i = 0; i < LIBMESH_DIM; ++i)
  {
    _dpermeability_qp_dgradvar[_qp][i].resize(_num_var, RealTensorValue());
    for (unsigned int v = 0; v < _num_var; ++v)
      _dpermeability_qp_dgradvar[_qp][i][v] =
          _dporosity_qp_dgradvar[_qp][v](i) * _permeability_qp[_qp] *
          (_n / _porosity_qp[_qp] + _m / (1.0 - _porosity_qp[_qp]));
  }
}
