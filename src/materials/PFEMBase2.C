#include "PFEMBase2.h"
#include "MooseRandom.h"
#include "Distribution.h"
#include <RandomInterface.h>

registerMooseObject("PorousFlowApp", PFEMBase2);

#include "libmesh/point.h"

namespace
{
inline Real
valueHelper(dof_id_type id, MooseRandom & generator, std::map<dof_id_type, Real> & map)
{
  auto it_pair = map.lower_bound(id);

  // Do we need to generate a new number?
  if (it_pair == map.end() || it_pair->first != id)
    it_pair = map.emplace_hint(it_pair, id, generator.rand(id));

  return it_pair->second;
}
}

InputParameters
PFEMBase2 ::validParams()
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
      " I_{ij} is the identity tensor and M_{ij} is a structure tensor given as n⊗n. n is a vector "
      "normal"
      " to the fracture.");
  params.addRangeCheckedParam<Real>("a", 1, "a > 0", "Mean (scalar) fracture distance value");
  params.addParam<Real>("e0", 1, "threshold strain");
  params.addParam<Real>("b0", 0, "Initial fracture aperture");
  params.addParam<Real>("km", 1, "matrix/intrinsic permeability");
  params.addParam<Real>("fix_rad_xy", 0, "fix fracture rotation angle in radians");
  params.addParam<Real>("fix_rad_yz", 0, "fix fracture rotation angle in radians");
  params.addParam<RealVectorValue>("n", "normal vector wrt to fracture surface");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addParam<bool>(
      "normal_vector_to_fracture_is_constant",
      true,
      "Whether the normal vector wrt to fracture surface is constant/known or not.");
  params.addParam<unsigned int>("seed", 0, "Seed value for the random number generator");
  params.addParam<Real>(
      "min", 0, "Lower bound of randomly or uniformly distributed random generated values");
  params.addParam<Real>(
      "max", 1.57, "Upper bound of randomly or uniformly distributed random generated values");
  params.addParam<DistributionName>(
      "distribution", "Name of distribution defining distribution of randomly generated values");
  params.addParam<bool>(
      "Random_field",
      true,
      "Whether to use spatially random angle of rotation at each timestep or a fixed"
      "angle of rotation. Set true if you want random field, false otherwise");
  return params;
}

PFEMBase2 ::PFEMBase2(const InputParameters & parameters)
  : PorousFlowPermeabilityBase(parameters),
    _distribution(nullptr),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _a(getParam<Real>("a")),
    _b0(getParam<Real>("b0")),
    _e0(getParam<Real>("e0")),
    _km(getParam<Real>("km")),
    _n_const(parameters.get<bool>("normal_vector_to_fracture_is_constant")),
    _Random_field(parameters.get<bool>("Random_field")),
    _nVec(parameters.isParamValid("n") ? getParam<RealVectorValue>("n")
                                       : RealVectorValue(1.0, 0.0, 0.0)),

    _identity_two(RankTwoTensor::initIdentity),
    _stress(getMaterialProperty<RankTwoTensor>(_base_name + "stress")),
    _strain(getMaterialProperty<RankTwoTensor>("total_strain")),
    _fix_rad_xy(getParam<Real>("fix_rad_xy")),
    _fix_rad_yz(getParam<Real>("fix_rad_yz")),

    _randm_rad_xy(_nodal_material
                      ? declareProperty<Real>("random_xy_rotation_angle_for_each_element")
                      : declareProperty<Real>("random_xy_rotation_angle_for_each_element_qp")),
    _randm_rad_yz(_nodal_material
                      ? declareProperty<Real>("random_yz_rotation_angle_for_each_element")
                      : declareProperty<Real>("random_yz_rotation_angle_for_each_element_qp")),

    _en(_nodal_material ? declareProperty<Real>("fracture_normal_strain_nodal")
                        : declareProperty<Real>("fracture_normal_strain_qp")),
    _min(getParam<Real>("min")),
    _max(getParam<Real>("max")),
    _current_node(nullptr),
    _elem_random_generator(nullptr),
    _node_random_generator(nullptr)
{
  unsigned int processor_seed = getParam<unsigned int>("seed");

  MooseRandom::seed(processor_seed);
  _elem_random_data =
      std::make_unique<RandomData>(_fe_problem, false, EXEC_INITIAL, MooseRandom::randl());
  _node_random_data =
      std::make_unique<RandomData>(_fe_problem, true, EXEC_INITIAL, MooseRandom::randl());

  _elem_random_generator = &_elem_random_data->getGenerator();
  _node_random_generator = &_node_random_data->getGenerator();

  if (_min >= _max)
    paramError("min", "Min >= Max for PFEMBase2 !");

  if (parameters.isParamSetByUser("distribution"))
  {
    _distribution = &getDistributionByName(getParam<DistributionName>("distribution"));
    if (parameters.isParamSetByUser("min") || parameters.isParamSetByUser("max"))
      paramError("distribution", "Cannot use together with 'min' or 'max' parameter");
  }
}

Real
PFEMBase2 ::generateRandom()
{
  Real rand_num;

  if (_current_node)
  {
    rand_num = valueHelper(_current_node->id(), *_node_random_generator, _node_numbers);
  }
  else if (_current_elem)
  {
    rand_num = valueHelper(_current_elem->id(), *_elem_random_generator, _elem_numbers);
  }
  else
  {
    mooseError("fail to generate parallel consistent random numbers. Please, check your input "
               "parameters or contact the MOOSE team for help");
  }

  return rand_num;
}

void
PFEMBase2 ::computeQpProperties()
{
  // This code block describes how the 'normal vector' (n) wrt the fracture face should
  // be computed. if the components of n is known (e.g., sigma_xx, tau_xy and tau_zx),
  // then it should be specify and obtain from the input file. Otherwise, n is computed as the
  // direction (eigenvector) of the principal stress vector corresponding to the normal direction.

  RealVectorValue _n;
  if (_n_const)
  {
    _n = _nVec;
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

  Real _rad_xy;
  Real _rad_yz;
  if (_Random_field)
  // Get the spatially random rotation angle (in radians) for each element at each timestep either
  // from a specific distribution or randomly.
  {
    if (_distribution)
    // Get random field from the specified distribution
    {
      _randm_rad_xy[_qp] = _distribution->quantile(generateRandom());
      _randm_rad_yz[_qp] = _distribution->quantile(generateRandom());
    }
    else
    // Get random field between min and max.
    {
      _randm_rad_xy[_qp] = generateRandom() * (_max - _min) + _min;
      _randm_rad_yz[_qp] = generateRandom() * (_max - _min) + _min;
    }

    _rad_xy = _randm_rad_xy[_qp]; // generateRandom() * (_max - _min) + _min;
    _rad_yz = _randm_rad_yz[_qp]; // generateRandom() * (_max - _min) + _min;
  }
  else
  // Get the fixed (or user-specified) rotation angle for all elements in the domain at each
  // timestep.
  {
    _rad_xy = _fix_rad_xy;
    _rad_yz = _fix_rad_yz;
  }

  // The fracture normal vector (n) is rotated around the Z-axis (i.e., X-Y plane) during random
  // rotation of the material using the correct rotation matrix (rotMat_xy and angle rad_xy).
  // (See Zill et al., 2022 for why material is randomly rotated).

  RankTwoTensor rotMat_xy;
  rotMat_xy(0, 0) = std::cos(_rad_xy);
  rotMat_xy(0, 1) = -std::sin(_rad_xy);
  rotMat_xy(0, 2) = 0;
  rotMat_xy(1, 0) = std::sin(_rad_xy);
  rotMat_xy(1, 1) = std::cos(_rad_xy);
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
  rotMat_yz(1, 1) = std::cos(_rad_yz);
  rotMat_yz(1, 2) = -std::sin(_rad_yz);
  rotMat_yz(2, 0) = 0;
  rotMat_yz(2, 1) = std::sin(_rad_yz);
  rotMat_yz(2, 2) = std::cos(_rad_yz);

  // Rotation of the  fracture normal vector using the rotation matrices
  RealVectorValue n_r = rotMat_xy * rotMat_yz * _n;

  // strain in the normal fracture direction
  _en[_qp] = n_r * (_strain[_qp] * n_r);

  // H_de is the heaviside function that implements the macaulay-bracket in Zill et al.
  // since _e0 is the initial/threshold strain state of the material, and strain is always
  // increasing in n-direction, _en should always be bigger than _e0. otherwise, change in
  // fracture aperture = 0

  Real H_de = (_en[_qp] > _e0) ? 1.0 : 0.0;

  // initial fracture aperture is sqrt(12 * k_m) in the literature
  //   Real _b0 = std::sqrt(12. * _km);

  // change in fracture aperture
  Real b_f = _b0 + (H_de * _a * (_en[_qp] - _e0));

  Real coeff = H_de * (b_f / _a) * ((b_f * b_f / 12.0) - _km);

  RankTwoTensor I = _identity_two;
  auto _M = RankTwoTensor::selfOuterProduct(n_r);

  // Finally, the permeability
  _permeability_qp[_qp] = (_km * I) + (coeff * (I - _M));
}
