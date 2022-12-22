#include "PFEMBase1.h"

#include "Moose.h"
#include "RandomInterface.h"
#include "RandomData.h"
#include "MooseRandom.h"
#include <FEProblemBase.h>
#include "Assembly.h"

#include "libmesh/fe_interface.h"
#include "libmesh/quadrature.h"


#include "MooseTypes.h"

registerMooseObject("PorousFlowApp", PFEMBase1);

InputParameters
PFEMBase1::validParams()
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

    params.addParam<unsigned int>("seed", 0, "The seed for the master random number generator");
    params.addParamNamesToGroup("seed", "Advanced");
/*
    params.addParam<Real>(
        "min", 0, "Lower bound of randomly or uniformly distributed random generated values");
    params.addParam<Real>(
        "max", 1.57, "Upper bound of randomly or uniformly distributed random generated values");
    params.addParam<DistributionName>(
    "distribution", "Name of distribution defining distribution of randomly generated values");
*/

    params.addParam<bool>("Random_field",
                      true,
                      "Whether to use spatially random angle of rotation at each timestep or a fixed"
                      "angle of rotation. Set true if you want random field, false otherwise");
  return params;
}

PFEMBase1::PFEMBase1(const InputParameters & parameters)
  : PorousFlowPermeabilityBase(parameters),
/*
    _min(getParam<Real>("min")),
    _max(getParam<Real>("max")),
    _distribution(nullptr),
*/

    _random_data(nullptr),
    _generator(nullptr),
    _ri_problem(problem),
    _ri_name(parameters.get<std::string>("_object_name")),
    _master_seed(parameters.get<unsigned int>("seed")),
    _is_nodal(is_nodal),
    _reset_on(EXEC_LINEAR),

    problem(*getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")),
    tid(getParam<THREAD_ID>("_tid")),

    _curr_node(problem.assembly(tid).node()),
    _curr_element(problem.assembly(tid).elem()),



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
                              : declareProperty<Real>("random_yz_rotation_angle_for_each_element_qp"))

{
  /*
   unsigned int processor_seed = getParam<unsigned int>("seed");
   MooseRandom::seed(processor_seed);
     _elem_random_data =
          std::make_unique<RandomData>(_fe_problem, false, EXEC_INITIAL, MooseRandom::randl());
     _node_random_data =
          std::make_unique<RandomData>(_fe_problem, true, EXEC_INITIAL, MooseRandom::randl());

//     _elem_random_generator = &_elem_random_data->getGenerator();
//     _node_random_generator = &_node_random_data->getGenerator();

     _elem_random_generator = &_elem_random_data->getRandmFieldLong();
     _node_random_generator = &_node_random_data->getRandmFieldLong();

  if (_min >= _max)
  paramError("min", "Min >= Max for PFEMBase1!");

  if (parameters.isParamSetByUser("distribution"))
  {
    _distribution = &getDistributionByName(getParam<DistributionName>("distribution"));
    if (parameters.isParamSetByUser("min") || parameters.isParamSetByUser("max"))
    paramError("distribution", "Cannot use together with 'min' or 'max' parameter");
  }
*/
}



PFEMBase1::~PFEMBase1() {}

void
PFEMBase1::setRandomResetFrequency(ExecFlagType exec_flag)
{
  _reset_on = exec_flag;
  _ri_problem.registerRandomInterface(*this, _ri_name);
}


void
PFEMBase1::setRandomDataPointer(RandomData * random_data)
{
  _random_data = random_data;
  _generator = &_random_data->getGenerator();
}

unsigned int
PFEMBase1::getSeed(std::size_t id)
{
  mooseAssert(_random_data, "RandomData object is NULL!");

  return _random_data->getSeed(id);
}


unsigned long
PFEMBase1::getRandmFieldLong() const
{
  mooseAssert(_generator, "Random Generator is NULL, did you call setRandomResetFrequency()?");

  dof_id_type id;
 if (_is_nodal)
    id = _curr_node->id();
 else
    id = _curr_element->id();

  return _generator->randl(id);
}


Real
PFEMBase1::getRandmFieldReal() const
{
  mooseAssert(_generator, "Random Generator is NULL, did you call setRandomResetFrequency()?");

  dof_id_type id;
  if (_is_nodal)
    id = _curr_node->id();
  else
    id = _curr_element->id();

  return _generator->rand(id);
}



void
PFEMBase1::computeQpProperties()
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


      Real  _rad_xy;
      Real  _rad_yz;
    if (_Random_field)
  // Get the spatially random rotation angle (in radians) for each element at each timestep either
  // from a specific distribution or randomly.
      {
//          if (_distribution)
//        // Get random field from the specified distribution
//          {
//            _randm_rad_xy[_qp] = _distribution->quantile(generateRandom());
//            _randm_rad_yz[_qp] = _distribution->quantile(generateRandom());
//          }
//          else
//         // Get random field between min and max.
//          {
           _randm_rad_xy[_qp] = getRandmFieldReal();
           _randm_rad_yz[_qp] = getRandmFieldReal();
//           }
        _rad_xy = _randm_rad_xy[_qp];
        _rad_yz = _randm_rad_yz[_qp];
        }
      else
  // Get the fixed (or user-specified) rotation angle for all elements in the domain at each timestep.
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
