#pragma once

#include "PorousFlowPermeabilityBase.h"
#include "RankTwoTensor.h"
#include "libmesh/utility.h"

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
class PFPermeabilityEmbeddedFractures : public PorousFlowPermeabilityBase
{
public:
	static InputParameters validParams();

	PFPermeabilityEmbeddedFractures(const InputParameters& parameters);

	using Material::_qp;
	using Material::_dt;
  using Material::_q_point;

protected:
	void computeQpProperties() override;

	/// optional parameter that allows multiple mechanics materials to be defined
  const std::string _base_name;

	/// mean fracture distance
	const Real _a;

	/// Threshold strain
	const Real _e0;

	/// matrix/intrinsic permeability
	const Real _km;

	/// Initial fracture aperture b_0 = sqrt(12 * k_m)
	const Real _b0;

	///  Identity RankTwotensor I_{ij}
	const RankTwoTensor _identity_two;

  /// get the strain tensor (actually, this is the creep strain)
  const MaterialProperty<RankTwoTensor> & _strain;

};
