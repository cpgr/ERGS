/******************************************************************************/
/*               The Lubby2-Creep Model (Part of DRACO-ERGS)                  */
/*                                                                            */
/*   DRACO-ERGS is an abreviation for Damage,Recovery And Creep Of Enhanced   */
/*   RockSalt Geothermal Systems. It is a compilation of the different creep  */
/*   models available in the literature that describe the mechanical          */
/*   behaviour mechanical behaviour of rock salt. The aim is to use these     */
/*   models to describe how rock salt behave when utilized as a geothermal    */
/*   repo                                                                     */
/*                                                                            */
/*          Copyright (C) 2022 by Ishmael Dominic Yevugah                     */
/*      University of Manitoba, Price Faculty of Engineering                  */
/*                                                                            */
/*        Special Thanks to Guillaume Giudicelli and the Moose Team           */
/*                                                                            */
/*       This program is free software: you can redistribute it and/or modify */
/*    it under the terms of the GNU General Public License as published by    */
/*      the Free Software Foundation, either version 3 of the License, or     */
/*                     (at your option) any later version.                    */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program.  If not, see <http://www.gnu.org/licenses/>    */
/******************************************************************************/

#include "Lubby2CreepBase.h"
#include <RankTwoTensor.h>

registerMooseObject("PorousFlowApp", Lubby2CreepBase);
registerMooseObject("PorousFlowApp", ADLubby2CreepBase);

template <bool is_ad>
InputParameters
Lubby2CreepBaseTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription(
    "This is a basic creep model for rocksalt based on the Lubby2CreepBase model. See Zhang"
    "and Nagel (2020): Error-controlled implicit time integration of elasto-visco-plastic"
    "constitutive models for rock salt. It uses the StressUpdate class in the radialreturn"
    "isotropic creep model. The creep results were compared with the analytical solution from"
    "Heusermann et. al. (2003): Nonlinear finite-element analysis of solution mined storage caverns"
    "in rock salt using the Lubby2CreepBase constitutive model. The results correspond very well");

  // Maxwell parameters
  params.addParam<Real>("mvM", "Maxwell viscosity parameter");
  params.addRequiredParam<Real>("etaM0", "Initial Maxwell Viscosity");
  // Kelvin parameters
  params.addParam<Real>("mvK", "Kelvin ViscoParameter");
  params.addParam<Real>("mk", "Kelvin Elastic Parameter");
  params.addRequiredParam<Real>("etaK0", "Initial Kelvin Viscosity");
  params.addParam<Real>("GK0", "Initial Kelvin Shear Modulus");
  // Applied stress
  params.addParam<RealTensorValue>("Stensor", "Applied Equivalent Stress");
params.addRequiredParam<Real>("smin_fix", "model parameter smin");
params.addRequiredCoupledVar("smin_ramp", "applied minimum stress as a function");
params.addParam<bool>("smin_function",
    true,
    "Whether to the applied minimum/confining stress is a function or not");
params.addParam<Real>("smax_fix", "fixed (scalar) applied maximum stress");
params.addRequiredCoupledVar("smax_ramp", "applied maximum stress as a function");
params.addParam<bool>("smax_function",
    false,
    "Whether to apply the maximum stress as a function or not");
params.addParam<bool>("uniaxial_test",
    true,
"Whether the test is uniaxial or triaxial");
 return params;
}

template <bool is_ad>
Lubby2CreepBaseTempl<is_ad>::Lubby2CreepBaseTempl(const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
    _etaM0(this->template getParam<Real>("etaM0")),
    _mvM(this->template getParam<Real>("mvM")),
    _mvK(this->template getParam<Real>("mvK")),
    _mk(this->template getParam<Real>("mk")),
    _etaK0(this->template getParam<Real>("etaK0")),
    _GK0(this->template getParam<Real>("GK0")),
    _kelvin_creep_strain(this->template declareGenericProperty<Real, is_ad>("kelvin_creep_strain")),
    _kelvin_creep_strain_old(this->template getMaterialPropertyOld<Real>("kelvin_creep_strain")),
    _uniaxial_test(parameters.get<bool>("uniaxial_test")),
    _smin_function(parameters.get<bool>("smin_function")),
    _smin_ramp(coupledValue("smin_ramp")),
    _smin_fix(this->template getParam<Real>("smin_fix")),
    _smax_function(parameters.get<bool>("smax_function")),
    _smax_ramp(coupledValue("smax_ramp")),
    _smax_fix(this->template getParam<Real>("smax_fix")),
    _Stensor(parameters.isParamValid("Stensor")
                 ? this->template getParam<RealTensorValue>("Stensor")
                 : RealTensorValue(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)) // RankTwoTensor
{
  if (_etaM0 == 0.0 && _etaK0 == 0.0)
    mooseError("Lubby2CreepBase: at least one of the creep should be active.");
}


template <bool is_ad>
void
Lubby2CreepBaseTempl<is_ad>::initQpStatefulProperties()
{
  _kelvin_creep_strain[_qp] = 0.0;
}


template <bool is_ad>
void
Lubby2CreepBaseTempl<is_ad>::propagateQpStatefulProperties()
{
  _kelvin_creep_strain[_qp] = _kelvin_creep_strain_old[_qp];
  RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}

template <bool is_ad>
void
Lubby2CreepBaseTempl<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & /* effective_trial_stress */,
    const GenericRankFourTensor<is_ad> & /* elasticity_tensor */)
{
  _kelvin_creep_strain[_qp] = _kelvin_creep_strain_old[_qp];
}


/// Compute the Residuals when automatic_differentiation = false
template <bool is_ad>
template <typename ScalarType>
ScalarType
Lubby2CreepBaseTempl<is_ad>::computeResidualInternal(const GenericReal<is_ad> & effective_trial_stress,
                                                 const ScalarType & scalar)
{
  // minimum stress
  ScalarType _smin;
  if (_smin_function)
    {
      _smin = _smin_ramp[_qp];
    }
    else
   {
    _smin = _smin_fix;
   }

 // maximum stress
   ScalarType _smax;
   if (_smax_function)
     {
       _smax = _smax_ramp[_qp];
     }
     else
    {
     _smax = _smax_fix;
   }

  ScalarType eqv_sigma =  _smax - _smin;

  RankTwoTensor _sigma = _Stensor;
  RankTwoTensor deviator_stress = _sigma.deviatoric();
  // compute the J2 stress
  ScalarType dev_stress_squared =
      deviator_stress.doubleContraction(deviator_stress);
  ScalarType J2 = dev_stress_squared/2.0;
  ScalarType eff = std::sqrt(J2); //Use von mises stress instead of the equivalent/effective von-mises stress (sqrt (3*J2))

  // final von mises equivalent stress based on whether test is uniaxial or triaxial
  ScalarType sigma_V;
  if (_uniaxial_test)
    {
      sigma_V = eff;
    }
    else
   {
    sigma_V = eqv_sigma;
   }

  const ScalarType stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);

  const ScalarType etaM = _etaM0 * std::exp(_mvM * sigma_V);
  const ScalarType etaK = _etaK0 * std::exp(_mvK * sigma_V);
  const ScalarType GK = _GK0 * std::exp(_mk * sigma_V);

  _kelvin_creep_strain[_qp] =  _kelvin_creep_strain_old[_qp] + MetaPhysicL::raw_value(scalar);

if (_etaM0 != 0.0 && _etaK0 != 0.0)
  {
    // Maxwell and Kelvin
  const ScalarType M_creep_rate = stress_delta / (3.0 * etaM);
 const ScalarType K_creep_rate = (stress_delta / (3.0 * etaK)) - ((GK *_kelvin_creep_strain[_qp])/(etaK));

return ((M_creep_rate + K_creep_rate) * _dt) - scalar;
   }
    else if (_etaM0 != 0.0 && _etaK0 == 0.0)
  // Maxwell
  {
   const ScalarType creep_rate = stress_delta / (3.0 * etaM);

return (creep_rate * _dt) - scalar;
  }
   else
  // Kelvin
  {
  const ScalarType creep_rate =
 (stress_delta / (3.0 * etaK)) - ((GK*_kelvin_creep_strain[_qp])/(etaK));

 return (creep_rate * _dt) - scalar;
  }
}


/// Activate the Derivatives when automatic_differentiation = false
template <bool is_ad>
GenericReal<is_ad>
Lubby2CreepBaseTempl<is_ad>::computeDerivative(
    const GenericReal<is_ad> & /* effective_trial_stress */,
    const GenericReal<is_ad> & /* scalar */)
{
  return 0.0;
}

// Note: The plastic_strain_increment is the inelastic_strain_increment coming from the newton iteration.
// Here, it is corrected by a factor to ensure a good match with the analytical solution of EnhancedLubby2.
// This correction factor is necessary because the initial strain_increment (i.e.,total strain increment)
// required for the newton iteration is extremely large, resulting in large creep strains.
template <bool is_ad>
void
Lubby2CreepBaseTempl<is_ad>::computeStressFinalize(
const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment * 1e-3;
}

template <bool is_ad>
void
Lubby2CreepBaseTempl<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
Lubby2CreepBaseTempl<is_ad>::substeppingCapabilityEnabled()
{
  return this->template getParam<bool>("use_substep");
}

template class Lubby2CreepBaseTempl<false>;
template class Lubby2CreepBaseTempl<true>;

template Real Lubby2CreepBaseTempl<false>::computeResidualInternal<Real>(const Real &, const Real &);
template ADReal Lubby2CreepBaseTempl<true>::computeResidualInternal<ADReal>(const ADReal &,
                                                                     const ADReal &);
template ChainedReal
Lubby2CreepBaseTempl<false>::computeResidualInternal<ChainedReal>(const Real &, const ChainedReal &);

template ChainedADReal
Lubby2CreepBaseTempl<true>::computeResidualInternal<ChainedADReal>(const ADReal &, const ChainedADReal &);
