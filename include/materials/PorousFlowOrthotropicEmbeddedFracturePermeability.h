#pragma once

#include "PorousFlowEmbeddedFracturePermeabilityBase.h"
#include "RankTwoTensor.h"
#include <Eigen/Geometry>
#include "libmesh/utility.h"

/**
 * Material designed to provide the permeability tensor which is a function of
 * the computed strain. This permeability material is based on Embedded Orthotropic Fractures.
 * See Zill et. al.(2021): Hydro-mechanical continuum modelling of fluid percolation
 * through rock. The permeability is given as follows:
 *
 * k = k_m*I +  [b_i/a_i * (\frac{b_{i}^2}{12} - k_m)*(I-M_i)].
 *
 * where i=summation Over number of fracture planes/surfaces. b is the fracture
 *  aperture given by: b_{i0} + /Delta{b_i} /Delta{b_i} depends on the strain (/epsilon) as follows:
 * /Delta{b_i} = a_i * 〈/epsilon :M_i - /epsilon_{0i}〉. Here, /epsilon_{0i} is a threshold strain of
 * the material in each of the fracture normal vector direction. /epsilon is total strain
 * a_i and b_i are fracture distance and fracture aperture respectively in each fracture
 * normal vector direction. K_m is the matrix/intrinsic permeability. I_{ij} is the identity tensor
 * and M_{ij} is a structure tensor given as n_i⊗n_i. n_i is a vector normal to each fracture plane.
 */
class PorousFlowOrthotropicEmbeddedFracturePermeability : public PorousFlowEmbeddedFracturePermeabilityBase
{
public:
	static InputParameters validParams();

	PorousFlowOrthotropicEmbeddedFracturePermeability(const InputParameters& parameters);

	using Material::_qp;
	using Material::_dt;
    using Material::_q_point;

protected:
	virtual void initQpStatefulProperties() override;
	void computeQpProperties() override;

	/// optional parameter that allows multiple mechanics materials to be defined
    const std::string _base_name;

	/// mean fracture distance in all 3 directions
    std::vector<double> _alpha;

	/// Threshold strain in all 3 directions
     std::vector<double> _eps;

	/// matrix/intrinsic permeability
	const Real _km;

	/// whether normal vector to fracture is constant or not
	bool _n_const;

	/// Tensor that holds the normal vector to the fracture surface in each direction
	RealTensorValue _NVec;

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
  	 const Real _rad_xy;
  	 const Real _rad_yz;

    /// Jacobian factor
	 const Real _jf;

    /// Computed strain in the fracture normal vector direction as a material property
	 MaterialProperty<Real>& _en;

    /// Old value of strain in the fracture normal vector direction
	 const MaterialProperty<Real>& _en_old;
};
