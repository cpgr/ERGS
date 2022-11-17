#pragma once

#include "PorousFlowPermeabilityBase.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include "RankTwoTensor.h"

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
class PorousFLowPermeabilityEmbeddedFractures : public PorousFlowPermeabilityBase
{
public:
	static InputParameters validParams();

	PorousFLowPermeabilityEmbeddedFractures(const InputParameters& parameters);

	using Material::_qp;
	using Material::_dt;
  using Material::_q_point;

protected:
	void computeQpProperties() override;

	/// mean fracture distance
	const Real _a;

	/// Threshold strain
	const Real _e0;

	/// matrix/intrinsic permeability
	const Real _km;

	/// Initial fracture aperture b_0 = sqrt(12 * k_m)
	const Real _b0;

	/// whether normal vector to fracture is constant or not
	bool _n_const;

	/// normal vector to fracture surface
	RealEigenMatrix _n;

	/// string to hold the name of the fracture rotation angle around xy
	std::string _phi_xy;

	/// string to hold the name of the fracture rotation angle around yz
	std::string _phi_yz;

	/// Structure tensor M_{ij} = n⊗n {not used}
  const RankTwoTensor _M;

	///  Identity RankTwotensor I_{ij}
	const RankTwoTensor _identity_two;

	/// stress
  const MaterialProperty<RankTwoTensor> & _stress;

	/// Strain (first const means we never want to dereference and change the value, second means we'll always be pointing to the same address after initialization (like a reference))
  const MaterialProperty<Real> * const _vol_strain_qp;


};
