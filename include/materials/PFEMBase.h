#pragma once

#include "PorousFlowPermeabilityBase.h"
#include "RankTwoTensor.h"
#include <Eigen/Geometry>
// #include <RandomInterface.h>
// #include "RandomData.h"
// #include <InitialConditionBase.h>
#include "MooseTypes.h"
#include "MooseEnumItem.h"

/**
 * Material designed to provide the permeability tensor which is a function of
 * the computed strain. This permeability material is based on Embedded Fractures.
 * See Zill et. al.(2021): Hydro-mechanical continuum modelling of fluid percolation
 * through rock. The permeability is given as follows:
 *
 * k = (k_m * I_{ij}) + (b/a*[(b^2/12 - k_m)]*(I-M_{ij}))
 *
 * where b is the fracture aperture given by: b_0 + /Delta{b}
 * /Delta{b} depends on the strain (/epsilon) as follows:
 * /Delta{b} = a * 〈/epsilon_n -/epsilon_0〉. Here, /epsilon_0 is a threshold strain
 * and /epsilon_n is the computed strain normal to the embedded fracture plane/surface.
 * a is mean fracture distance, b is fracture aperture, K_m is matrix/intrinsic permeability.
 * I_{ij} is the identity tensor and M_{ij} is a structure tensor given as n⊗n. n is a vector normal
 * to the fracture.
 */
class PFEMBase : public PorousFlowPermeabilityBase
{
public:
  static InputParameters validParams();

  PFEMBase(const InputParameters & parameters);

protected:
  void computeQpProperties() override;
  Distribution const * _distribution;

  /// optional parameter that allows multiple mechanics materials to be defined
  const std::string _base_name;

  /// mean fracture distance
  const Real _a;

  /// Initial fracture aperture
  const Real _b0;

  /// Threshold strain
  const Real _e0;

  /// matrix/intrinsic permeability
  const Real _km;

  /// whether normal vector to fracture is constant or not
  bool _n_const;

  /// whether to use random rotation angle or fix
  bool _Random_field;

  /// normal vector to fracture surface
  RealVectorValue _nVec;

  /// string to hold the name of the fracture rotation angle around xy
  std::string _phi_xy;

  /// string to hold the name of the fracture rotation angle around yz
  std::string _phi_yz;

  /// Structure tensor M_{ij} = n⊗n {not used}
  const RankTwoTensor _M;

  ///  Identity RankTwotensor I_{ij}
  const RankTwoTensor _identity_two;

  /// get the stress tensor
  const MaterialProperty<RankTwoTensor> & _stress;

  /// get the strain tensor (actually, this is the creep strain)
  const MaterialProperty<RankTwoTensor> & _strain;

  /// fracture rotation angles about xy and yz (in radians)
  const Real _fix_rad_xy;
  const Real _fix_rad_yz;

  /// random rotation_angle_for_each_element
  MaterialProperty<Real> & _randm_rad_xy;
  MaterialProperty<Real> & _randm_rad_yz;

  /// Computed strain in the fracture normal vector direction as a material property
  MaterialProperty<Real> & _en;

  /// Lower and Upper bound of the 'randomly' or the 'uniformly distributed random' generated values
  const Real _min;
  const Real _max;

  std::size_t _seed;

  const Real _rotXY;
  const Real _rotYZ;
};
