#pragma once

#include "RadialReturnCreepStressUpdateBase.h"
#include<algorithm>
#include <complex>
#include<cmath>

/**
 * This is a basic creep model for rocksalt based on the Lubby2CreepBase model. See Zhang
 * and Nagel (2020): Error-controlled implicit time integration of elasto-visco-plastic
 * constitutive models for rock salt. It uses the StressUpdate class in the radialreturn
 * isotropic creep model.
 */

template <bool is_ad>
class Lubby2CreepBaseTempl : public RadialReturnCreepStressUpdateBaseTempl<is_ad>
{
public:
  static InputParameters validParams();
  Lubby2CreepBaseTempl(const InputParameters & parameters);


  virtual bool substeppingCapabilityEnabled() override;
  virtual void resetIncrementalMaterialProperties() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;

  virtual void
  computeStressInitialize(const GenericReal<is_ad> & effective_trial_stress,
                          const GenericRankFourTensor<is_ad> & elasticity_tensor) override;
  virtual void
  computeStressFinalize(const GenericRankTwoTensor<is_ad> & plastic_strain_increment) override;


  //Declare residuals when automatic_differentiation = false

  virtual GenericReal<is_ad> computeResidual(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar) override
  {
    return computeResidualInternal<GenericReal<is_ad>>(effective_trial_stress, scalar);
  }


  ///Declare Derivatives when automatic_differentiation = false
  virtual GenericReal<is_ad> computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                               const GenericReal<is_ad> & scalar) override;


  // Declare computeResidual And Derivatives for automatic_differentiation:
  virtual GenericChainedReal<is_ad>
  computeResidualAndDerivative(const GenericReal<is_ad> & effective_trial_stress,
                               const GenericChainedReal<is_ad> & scalar) override
  {
    return computeResidualInternal<GenericChainedReal<is_ad>>(effective_trial_stress, scalar);
  }


  /// whether the creep test is uniaxial or triaxial (von-mises stress is computed differently)
  bool _uniaxial_test;


  /// whether the applied minimum/confining stress depends on time or fix
  bool _smin_function;
  const VariableValue& _smin_ramp;
  Real _smin_fix;

  /// whether the maximum applied stress depends on time or fix
  bool _smax_function;
  const VariableValue& _smax_ramp;
  Real _smax_fix;

  /// Maxwell initial viscosity
  const Real _etaM0;
  /// Maxwell viscosity parameter
  const Real _mvM;
  /// Kelvin ViscoParameter
  const Real _mvK;
  /// Kelvin Elastic Parameter
  const Real _mk;
  /// Initial Kelvin Viscosity
  const Real _etaK0;
  /// Initial Kelvin Shear Modulus
  const Real _GK0;

  GenericMaterialProperty<Real, is_ad> & _kelvin_creep_strain;
  const MaterialProperty<Real> & _kelvin_creep_strain_old;

  // The (applied) stress responsible for the creep (or strain-rate) as a single value
  const RankTwoTensor _Stensor;

  // bool to determine whether the (applied) stress responsible for the creep (or strain-rate) is a function or not
 // const Function* const _StensorFn;

  // Whether the applied stress is a single value or a function
 // const bool _function_Stensor;


  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_qp;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_dt;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_t;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_three_shear_modulus;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_creep_strain;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_creep_strain_old;

private:

template <typename ScalarType>
  ScalarType computeResidualInternal(const GenericReal<is_ad> & effective_trial_stress,
                                     const ScalarType & scalar);

};

typedef Lubby2CreepBaseTempl<false> Lubby2CreepBase;
typedef Lubby2CreepBaseTempl<true> ADLubby2CreepBase;
