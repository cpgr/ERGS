#include "PorousFlowBrineSaltCO2.h"
#include "BrineFluidPropertiesBeta.h"
#include "SinglePhaseFluidProperties.h"
#include "MathUtils.h"
#include "Conversion.h"
#include "Water97FluidProperties.h"
#include "NaClFluidProperties.h"

registerMooseObject("PorousFlowApp", PorousFlowBrineSaltCO2);

InputParameters
PorousFlowBrineSaltCO2::validParams()
{
  InputParameters params = PorousFlowBrineCO2::validParams();
  params.addParam<unsigned int>("salt_component", 2, "The component number of salt");
  params.addParam<unsigned int>("solid_phase_number", 2, "The phase number for salt");
  params.addRequiredParam<UserObjectName>("water_fp", "The name of the user object for water");
  params.addParam<Real>("saturationSOLID", 0.01, "A small non-zero initial saturation"
                             "of halite in the multiphase system");
  params.addClassDescription("Fluid state class for brine, salt and Co2. Includes the"
                              "dissolution/precipitation of the solid-salt/halite");
  return params;
}

PorousFlowBrineSaltCO2::PorousFlowBrineSaltCO2(const InputParameters & parameters)
  : PorousFlowBrineCO2(parameters),
    _brine_fp(getUserObject<BrineFluidPropertiesBeta>("brine_fp")),
    _water_fp(getUserObject<SinglePhaseFluidProperties>("water_fp")),
    _solid_phase_number(getParam<unsigned int>("solid_phase_number")),
    _saturationSOLID(getParam<Real>("saturationSOLID"))
{
  // Set the number of phases and components, and their indexes
  _num_phases = 3;
  _num_components = 3;
  _gas_phase_number = 3 - _aqueous_phase_number - _solid_phase_number;
  _gas_fluid_component = 3 - _aqueous_fluid_component - _salt_component;

  // Check that _aqueous_phase_number is <= total number of phases
  if (_aqueous_phase_number >= _num_phases)
    paramError("liquid_phase_number",
               "This value is larger than the possible number of phases ",
               _num_phases);

  // Check that the solid phase number/index is not identical to the aqueous phase number/index
  if (_solid_phase_number == _aqueous_phase_number)
    paramError( "solid_phase_number",
          "The value provided must be different from the value entered in aqueous_phase_number");

  // Check that solid_phase number is <= total number of phases
  if (_solid_phase_number >= _num_phases)
    paramError("solid_phase_number",
      "The value provided is larger than the possible number of phases",
                _num_phases);

  // Check that _aqueous_fluid_component is <= total number of fluid components
  if (_aqueous_fluid_component >= _num_components)
    paramError("liquid_fluid_component",
               "This value is larger than the possible number of fluid components",
               _num_components);

  // Check that the salt component index is not identical to the liquid fluid component
  if (_salt_component == _aqueous_fluid_component)
    paramError(
        "salt_component",
        "The value provided must be different from the value entered in liquid_fluid_component");

  // Check that _salt_component is <= total number of fluid components
  if (_salt_component >= _num_components)
    paramError("salt_component",
               "The value provided is larger than the possible number of fluid components",
               _num_components);

  _empty_fsp = FluidStateProperties(_num_components);
}

void
PorousFlowBrineSaltCO2::thermophysicalProperties(Real pressure,
                                                 Real temperature,
                                                 Real Xnacl,
                                                 Real Z,
                                                 unsigned int qp,
                                                 std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];
  FluidStateProperties & solid = fsp[_solid_phase_number];

  // Check whether the input temperature is within the region of validity
  checkVariables(pressure, temperature);

  // AD versions of primary variables
  DualReal p = pressure;
  Moose::derivInsert(p.derivatives(), _pidx, 1.0);
  DualReal T = temperature;
  Moose::derivInsert(T.derivatives(), _Tidx, 1.0);
  DualReal Zco2 = Z;
  Moose::derivInsert(Zco2.derivatives(), _Zidx, 1.0);
  DualReal X = Xnacl;
  Moose::derivInsert(X.derivatives(), _Xidx, 1.0);

  // Clear all of the FluidStateProperties data
  clearFluidStateProperties(fsp);

  FluidStatePhaseEnum phase_state;
  massFractions(p, T, X, Zco2, phase_state, fsp);

  switch (phase_state)
  {
    case FluidStatePhaseEnum::GAS:
    {
      _console << "phase state gas " << std::endl;

      // Set the gas saturations
      // note: since it is single gas phase, gas saturation is 1. It
      // is not used in the gas properties, instead the total mass
      // fraction is used.
      gas.saturation = 1.0;

      // Calculate gas properties
      gasProperties(p, T, fsp);

      break;
    }

    case FluidStatePhaseEnum::LIQUID:
    {
      _console << "phase state liquid " << std::endl;

      // Set the gas saturations
      // note: since it is single liquid phase, gas and solid saturation are 0,
      // although they are not used in the liquid properties.
      gas.saturation = 0.0;

      // Calculate the liquid properties
      // note: here, the liquid properties depend on liquid pressure!
      const DualReal liquid_pressure = p - _pc.capillaryPressure(1.0, qp);
      liquidProperties(liquid_pressure, T, X, fsp);

      break;
    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      _console << "phase state multiphase " << std::endl;

      // Calculate the gas and liquid properties in the two phase region
      // note: Here, we have multiphase. Therefore, we need the saturation
      // (computed internally) together with the various equilibrium mass fractions
      MultiPhaseProperties(p, T, X, Zco2, qp, fsp);

      break;
    }
  }

  // set the solid/halite saturation (initialized to a small non-zero number)
  solid.saturation = _saturationSOLID;
  // Liquid saturations can now be set
  liquid.saturation = 1.0 - gas.saturation - solid.saturation;

  // Save pressures to FluidStateProperties object
  gas.pressure = p;
  liquid.pressure = p - _pc.capillaryPressure(liquid.saturation, qp);
  solid.pressure = p; // pressure in solid phase is the same as gas phase. no Pc
}

void
PorousFlowBrineSaltCO2::massFractions(const DualReal & pressure,
                                      const DualReal & temperature,
                                      const DualReal & Xnacl,
                                      const DualReal & Z,
                                      FluidStatePhaseEnum & phase_state,
                                      std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];
  FluidStateProperties solid = fsp[_solid_phase_number];

  DualReal Xco2 = 0.0;
  DualReal Yh2o = 0.0;
  DualReal Yco2 = 0.0;
  DualReal Ynacl;  //halite in gas phase
  DualReal Snacl;  //halite in solid phase

 // halite solubility in liquid phase
  const DualReal XEQ = _brine_fp.haliteSolubilityWater(temperature,pressure);
  _console << "XEQ " << XEQ << std::endl;

  // If the amount of CO2 is less than the smallest solubility,then all CO2 will
  // be dissolved/dissapear, and the equilibrium mass fractions do not need to be computed
  if (Z < _Zmin)
  // note: Zmin is minimum amount of co2 that could exist in the gas phase. Hence,
  // the above info means that if CO2 is less than the smallest amount set, then there
  // is no Co2 in the gas phase and all CO2 is concentrated in the liquid phase.
  // Therefore, there is no multiphase.
  {
    phase_state = FluidStatePhaseEnum::LIQUID;
  }
  else if (Xnacl > XEQ)
     {
  // the amount of halite in the liquid phase is greater than or above its solubility/concentration
  // (limit) due to high amount of water evaporating from the liquid phase into the gas phase .
  // Gas phase is fully developed while liquid phase dissapears
       phase_state = FluidStatePhaseEnum::GAS;
     }
  else
  {
 // Equilibrium mass fraction of CO2 in liquid (Yco2) and H2O in gas phases (Xh20) are computed here.
 //note: Also, the contribution of salt in both liquid and gas is now included, i.e., Xnacl and Ynacl
    equilibriumMassFractions(pressure, temperature, Xnacl, Xco2, Yh2o);

    Yco2 = 1.0 - Yh2o - Ynacl ; // CO2 in liquid corrected with halite in the liquid phase

    // Determine which phases are present based on the value of z
    phaseState(Z.value(), Xco2.value(), Yco2.value(), phase_state); // LATER
  }

  // The equilibrium mass fractions calculated above are only correct in the two phase
  // state. If only liquid or gas phases are present, the mass fractions are given by
  // the total mass fraction z which includes the salt/halite mass fraction too

  DualReal Xh2o = 0.0;

  switch (phase_state)
  {
    case FluidStatePhaseEnum::LIQUID:
    {
    DualReal Xnacl;
    DualReal Ynacl;
    DualReal Snacl;
        Xco2 = Z;
        Yco2 = 0.0;
// note: the total mass fraction of 'gas', Z, assigned is taken by co2 in the liquid phase
// since solid phase is included, the mass fraction of h20 and halite that is
// left after the entire gas phase is taken by co2 should be shared among both h2o and halite:
        Xh2o = (1.0 - Z)/2.0;
        Yh2o = 0.0;
        Xnacl = (1.0 - Z)/2.0;
        Ynacl = 0.0;
        Snacl = 0.0;
      Moose::derivInsert(Xco2.derivatives(), _pidx, 0.0);
      Moose::derivInsert(Xco2.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(Xco2.derivatives(), _Xidx, 0.0);
      Moose::derivInsert(Xco2.derivatives(), _Zidx, 1.0);
      Moose::derivInsert(Yco2.derivatives(), _pidx, 0.0);
      Moose::derivInsert(Yco2.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(Yco2.derivatives(), _Xidx, 0.0);
      break;
    }

    case FluidStatePhaseEnum::GAS:
    {
    DualReal Xnacl;
    DualReal Ynacl;
    DualReal Snacl;
      Xco2 = 0.0;
      Yco2 = Z;
 // same here, total mass fraction of 'gas' assign is taken by co2 in the gas phase.
 // the rest should be shared.
 // note: Xh2o is not included here because it is already initialize to zero above.
      Yh2o = 1.0 - Z;
      Xnacl = 0.0;
      Ynacl = (1.0 - Z)/2;
      Snacl = 0.0;
      Moose::derivInsert(Xco2.derivatives(), _pidx, 0.0);
      Moose::derivInsert(Xco2.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(Xco2.derivatives(), _Xidx, 0.0);
      Moose::derivInsert(Yco2.derivatives(), _pidx, 0.0);
      Moose::derivInsert(Yco2.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(Yco2.derivatives(), _Xidx, 0.0);
      Moose::derivInsert(Yco2.derivatives(), _Zidx, 1.0);
      break;

    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      // Keep equilibrium mass fractions
      Xh2o = 1.0 - Xco2 - Xnacl; //  H2o in gas phase corrected with halite in the gas phase
      break;
    }
  }

/// Save the mass fractions in the FluidStateProperties object
// Note: There is no mass fraction of h20 and co2 in the solid phase (i.e., Sh2o & SCo2)
// because adsorption of gas and liquid onto the solid mass is negligible.
  liquid.mass_fraction[_aqueous_fluid_component] = Xh2o;
  liquid.mass_fraction[_gas_fluid_component] = Xco2;
  liquid.mass_fraction[_aqueous_fluid_component] = Xnacl;
  gas.mass_fraction[_aqueous_fluid_component] = Yh2o;
  gas.mass_fraction[_gas_fluid_component] = Yco2;
  gas.mass_fraction[_gas_fluid_component] = Ynacl;
  solid.mass_fraction[_salt_component] = Snacl;
}

void
PorousFlowBrineSaltCO2::gasProperties(const DualReal & pressure,
                                      const DualReal & temperature,
                                     std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

 /// included correction factor to account for water vapor in the gas phase
 /// note: since halite is purely solid, it's density and enthalpy does not
 /// affect the gas properties/behaviour.
  const DualReal psat = _water_fp.vaporPressure(temperature);

  const DualReal Yco2 = gas.mass_fraction[_gas_fluid_component];
  const DualReal Xco2 = liquid.mass_fraction[_gas_fluid_component];

  // Gas density, viscosity and enthalpy are calculated using partial pressure
  // Yco2 * gas_poreressure (Dalton's law)
  DualReal co2_density, co2_viscosity;
  _co2_fp.rho_mu_from_p_T(Yco2 * pressure, temperature, co2_density, co2_viscosity);

  DualReal co2_enthalpy = _co2_fp.h_from_p_T(Yco2 * pressure, temperature);

  // Vapor density, viscosity and enthalpy calculated using partial pressure
  // X1 * psat (Raoult's law)
  DualReal vapor_density, vapor_viscosity;

  _water_fp.rho_mu_from_p_T((1.0 - Xco2) * psat, temperature, vapor_density, vapor_viscosity);
  DualReal vapor_enthalpy = _water_fp.h_from_p_T((1.0 - Xco2) * psat, temperature);

  /// Save the values to the FluidStateProperties object. Note that derivatives wrt z are 0
  // Density is just the sum of individual component densities
  gas.density = co2_density + vapor_density;
  // Viscosity of the gas phase is a weighted sum of the individual viscosities
  gas.viscosity = Yco2 * co2_viscosity + (1.0 - Yco2) * vapor_viscosity;
  // Enthalpy of the gas phase is a weighted sum of the individual enthalpies
  gas.enthalpy = Yco2 * co2_enthalpy + (1.0 - Yco2) * vapor_enthalpy;

  //  Internal energy of the gas phase (e = h - pv)
  mooseAssert(gas.density.value() > 0.0, "Gas density must be greater than zero");
  gas.internal_energy = gas.enthalpy - pressure / gas.density;
}

void
PorousFlowBrineSaltCO2::liquidProperties(const DualReal & pressure,
                                         const DualReal & temperature,
                                         const DualReal & Xnacl,
                                         std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];

  // The liquid density includes the density increase due to dissolved CO2 (liquid)
  /// note: since halite is purely solid, it's density and enthalpy does not
  /// affect the liquid properties/behaviour.
  const DualReal brine_density = _brine_fp.rho_from_p_T_X(pressure, temperature, Xnacl);

  // Mass fraction of CO2 in liquid phase
  const DualReal Xco2 = liquid.mass_fraction[_gas_fluid_component];

  // The liquid density
  const DualReal co2_partial_density = partialDensityCO2(temperature);

  const DualReal liquid_density = 1.0 / (Xco2 / co2_partial_density + (1.0 - Xco2) / brine_density);

  // Assume that liquid viscosity is just the brine viscosity
  const DualReal liquid_viscosity = _brine_fp.mu_from_p_T_X(pressure, temperature, Xnacl);

  // Liquid enthalpy (including contribution due to the enthalpy of dissolution)
  const DualReal brine_enthalpy = _brine_fp.h_from_p_T_X(pressure, temperature, Xnacl);

  // Enthalpy of CO2
  const DualReal co2_enthalpy = _co2_fp.h_from_p_T(pressure, temperature);

  // Enthalpy of dissolution
  const DualReal hdis = enthalpyOfDissolution(temperature);

  const DualReal liquid_enthalpy = (1.0 - Xco2) * brine_enthalpy + Xco2 * (co2_enthalpy + hdis);

  // Save the values to the FluidStateProperties object
  liquid.density = liquid_density;
  liquid.viscosity = liquid_viscosity;
  liquid.enthalpy = liquid_enthalpy;

  mooseAssert(liquid.density.value() > 0.0, "Liquid density must be greater than zero");
  liquid.internal_energy = liquid.enthalpy - pressure / liquid.density;
}


void
PorousFlowBrineSaltCO2::solidProperties(const DualReal & pressure,
                                        const DualReal & temperature,
                                        std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & solid = fsp[_solid_phase_number];

  /// note: Pure solid properties. No modifications/corrections needed to account
  /// for Co2 and Brine since they do not exist in solid state.

  // halite density
  const DualReal halite_density = _brine_fp.halite_rho_from_p_T(pressure, temperature);

  // halite enthalpy
  const DualReal halite_enthalpy = _brine_fp.halite_h_from_p_T(pressure, temperature);

  // Save the solid thermodynamic properties to the solid phase in the FluidStateProperties object
  solid.density  = halite_density;
  solid.enthalpy = halite_enthalpy;

  mooseAssert(solid.density.value() > 0.0, "solid density must be greater than zero");
  solid.internal_energy = solid.enthalpy - pressure / solid.density;
}


DualReal
PorousFlowBrineSaltCO2::saturationGAS(const DualReal & pressure,
                               const DualReal & temperature,
                               const DualReal & Xnacl,
                               const DualReal & Z,
                               std::vector<FluidStateProperties> & fsp) const
{
  auto & gas = fsp[_gas_phase_number];
  auto & liquid = fsp[_aqueous_fluid_component];

  // Approximate liquid density as saturation isn't known yet, by using the gas
  // pressure rather than the liquid pressure. This does result in a small error
  // in the calculated saturation, but this is below the error associated with
  // the correlations. A more accurate saturation could be found iteraviely,
  // at the cost of increased computational expense
  // note: Gas saturation is computed from vaporMassfraction, gas_density and liquid_density

  // compute the Gas density
  const DualReal gas_density = _co2_fp.rho_from_p_T(pressure, temperature);

  // compute an approximate liquid density as saturation isn't known yet
  const DualReal brine_density = _brine_fp.rho_from_p_T_X(pressure, temperature, Xnacl);

  // Mass fraction of CO2 in liquid phase
  const DualReal Xco2 = liquid.mass_fraction[_gas_fluid_component];

  const DualReal co2_partial_density = partialDensityCO2(temperature);

    // finally the liquid density
  const DualReal liquid_density = 1.0 / (Xco2 / co2_partial_density + (1.0 - Xco2) / brine_density);

 /// compute the saturation using vapor_mass_fraction, liquid and gas density
  // mass fraction of co2 in the gas phase
  const DualReal Yco2 = gas.mass_fraction[_gas_fluid_component];

  // Set mass equilibrium constants used in the calculation of vapor mass fraction
  const DualReal K0 = Yco2 / Xco2;
  const DualReal K1 = (1.0 - Yco2) / (1.0 - Xco2);
  //vapor mass fraction
  const DualReal vapor_mass_fraction = vaporMassFraction(Z, K0, K1);

  // The gas saturation in the two phase case
  const DualReal saturation = vapor_mass_fraction * liquid_density /
                              (gas_density + vapor_mass_fraction * (liquid_density - gas_density));

  return saturation;
}

void
PorousFlowBrineSaltCO2::MultiPhaseProperties(const DualReal & pressure,
                                       const DualReal & temperature,
                                       const DualReal & Xnacl,
                                       const DualReal & Z,
                                       unsigned int qp,
                                       std::vector<FluidStateProperties> & fsp) const
{
  auto & gas = fsp[_gas_phase_number];
//  auto & solid = fsp[_solid_phase_number];
//  auto & liquid = fsp[_aqueous_fluid_component];

  // The gas saturation in the multiphase system
  gas.saturation = saturationGAS(pressure, temperature, Xnacl, Z, fsp);
  _console << "gas.saturation " << gas.saturation << std::endl;

  /// Calculate  the gas phase properties
  // note: since the gas properties are corrected with the dissolved water vapor,
  // the pressure is also affected by the gas saturation:
  const DualReal gas_pressure = pressure - _pc.capillaryPressure(gas.saturation, qp);
  gasProperties(gas_pressure, temperature, fsp);

  // The liquid pressure and properties can now be calculated
  // const DualReal liquid_pressure = pressure - _pc.capillaryPressure(1.0 - gas.saturation - solid.saturation, qp);
  const DualReal liquid_pressure = pressure - _pc.capillaryPressure(1 - gas.saturation - _saturationSOLID, qp);
  liquidProperties(liquid_pressure, temperature, Xnacl, fsp);

  // The solid properties (note that the effect of liquid and gas adsorptions is neglected)
  // Pc occurs only in fluids
  solidProperties(pressure, temperature, fsp);
}

void
PorousFlowBrineSaltCO2::equilibriumMassFractions(const DualReal & pressure,
                                                 const DualReal & temperature,
                                                 const DualReal & Xnacl,
                                                 DualReal & Xco2,
                                                 DualReal & Yh2o) const
{
  // Mole fractions at equilibrium
  DualReal xco2, yh2o;
  equilibriumMoleFractions(pressure, temperature, Xnacl, xco2, yh2o);

  // The mass fraction of H2O in gas (assume no salt in gas phase) and derivatives
  // wrt p, T, and X
  Yh2o = yh2o * _Mh2o / (yh2o * _Mh2o + (1.0 - yh2o) * _Mco2);

  // NaCl molality (mol/kg)
  const DualReal mnacl = Xnacl / (1.0 - Xnacl) / _Mnacl;

  // The molality of CO2 in 1kg of H2O
  const DualReal mco2 = xco2 * (2.0 * mnacl + _invMh2o) / (1.0 - xco2);
  // The mass fraction of CO2 in brine is then
  const DualReal denominator = (1.0 + mnacl * _Mnacl + mco2 * _Mco2);
  Xco2 = mco2 * _Mco2 / denominator;

  // Compute the mass fraction of NaCl/halite in liquid and vapor phase to account
  // for precipitation and dissolution of halite due to persistent boiling of the
  // brine solution

  // halite solubility in liquid phase
  const DualReal XEQ = _brine_fp.haliteSolubilityWater(temperature,pressure);
  // halite solubility in water vapour (gas) phase
  const DualReal XEQG = _brine_fp.haliteSolubilityGas(temperature, pressure);

  // The mass fraction of halite in the (water) vapor phase(Ynacl), computed from
  // the mass fraction in liquid phase (Xnacl).
  // (note: derivatives wrt to p, T, and X is not accounted for):
  const DualReal Ynacl  =  Xnacl * (XEQG/XEQ);

  // halite in the solid phase (Snacl)
  // note: since adsorption is not allowed, Sh2o and Sco2 equal 0 or neglected.
  const DualReal Snacl  =  1 - Xnacl - Ynacl;
}
