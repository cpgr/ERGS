#include "PorousFlowWaterAir.h"
#include "SinglePhaseFluidProperties.h"
#include "Water97FluidProperties.h"
#include "Conversion.h"

registerMooseObject("PorousFlowApp", PorousFlowWaterAir);

InputParameters
PorousFlowWaterAir::validParams()
{
  InputParameters params = PorousFlowWaterNCG::validParams();
  params.addParam<Real>("R", 8314.56,"Universal Gas Constant");
  params.addClassDescription("Fluid state class for water and air (Ideal gas),"
                              "with vapor pressure lowering capability");
  return params;
}


PorousFlowWaterAir::PorousFlowWaterAir(const InputParameters & parameters)
  : PorousFlowWaterNCG(parameters),
    _air_fp(getUserObject<SinglePhaseFluidProperties>("gas_fp")),
    _Mair(_air_fp.molarMass()),
    _air_henry(_air_fp.henryCoefficients()),
    _R(getParam<Real>("R"))

{
  // Check that the correct FluidProperties UserObjects have been provided
  if (_water_fp.fluidName() != "water")
    paramError("water_fp", "A valid water FluidProperties UserObject must be provided in water_fp");

  // Set the number of phases and components, and their indexes
  _num_phases = 2;
  _num_components = 2;
  _gas_phase_number = 1 - _aqueous_phase_number;
  _gas_fluid_component = 1 - _aqueous_fluid_component;

  // Check that _aqueous_phase_number is <= total number of phases
  if (_aqueous_phase_number >= _num_phases)
    paramError("liquid_phase_number",
               "This value is larger than the possible number of phases ",
               _num_phases);

  // Check that _aqueous_fluid_component is <= total number of fluid components
  if (_aqueous_fluid_component >= _num_components)
    paramError("liquid_fluid_component",
               "This value is larger than the possible number of fluid components",
               _num_components);

  _empty_fsp = FluidStateProperties(_num_components);
}

std::string
PorousFlowWaterAir::fluidStateName() const
{
  return "water-air";
}


void
PorousFlowWaterAir::thermophysicalProperties(Real pressure,
                                             Real temperature,
                                             Real /* Xnacl */,
                                             Real Z,
                                             unsigned int qp,
                                             std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Check whether the input temperature is within the region of validity
  checkVariables(temperature);

  // AD versions of primary variables
  DualReal p = pressure;
  Moose::derivInsert(p.derivatives(), _pidx, 1.0);
  DualReal T = temperature;
  Moose::derivInsert(T.derivatives(), _Tidx, 1.0);
  DualReal Zair = Z;
  Moose::derivInsert(Zair.derivatives(), _Zidx, 1.0);
  // saturated vapor pressure
  DualReal psat;

  // Clear all of the FluidStateProperties data
  clearFluidStateProperties(fsp);

  FluidStatePhaseEnum phase_state;
  massFractions(p, T, Zair, phase_state, fsp);

  switch (phase_state)
  {
    case FluidStatePhaseEnum::GAS:
    {
      // Set the gas saturations
      gas.saturation = 1.0;

      // Calculate gas properties
      gasProperties(p, T, fsp);

      break;
    }

    case FluidStatePhaseEnum::LIQUID:
    {
      // Calculate the liquid properties
      const DualReal liquid_pressure = p - _pc.capillaryPressure(1.0, qp);
      liquidProperties(liquid_pressure, T, fsp);

      break;
    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      // Calculate the gas and liquid properties in the two phase region
      twoPhaseProperties(p, T, Zair, qp, fsp);

      break;
    }
  }

  // Liquid saturations can now be set
  liquid.saturation = 1.0 - gas.saturation;

  // Vapor pressure (with lowering capability)
  DualReal numer = _Mh2o *_pc.capillaryPressure(liquid.saturation, qp) * liquid.saturation;
  DualReal deno  = liquid.density * _R * T;
  DualReal f_vpl = std::exp(numer/deno);
  DualReal pv = f_vpl * psat;

  // Save pressures to FluidStateProperties object
  gas.pressure = p + pv;
//  liquid.pressure = p - _pc.capillaryPressure(liquid.saturation, qp);
  liquid.pressure = gas.pressure - _pc.capillaryPressure(liquid.saturation, qp);
}

void
PorousFlowWaterAir::massFractions(const DualReal & pressure,
                                  const DualReal & temperature,
                                  const DualReal & Z,
                                  FluidStatePhaseEnum & phase_state,
                                  std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Equilibrium mass fraction of Air in liquid and H2O in gas phases
  DualReal Xair, Yh2o;
  equilibriumMassFractions(pressure, temperature, Xair, Yh2o);

  DualReal Yair = 1.0 - Yh2o;

  // Determine which phases are present based on the value of Z
  phaseState(Z.value(), Xair.value(), Yair.value(), phase_state);

  // The equilibrium mass fractions calculated above are only correct in the two phase
  // state. If only liquid or gas phases are present, the mass fractions are given by
  // the total mass fraction Z.
  DualReal Xh2o = 0.0;

  switch (phase_state)
  {
    case FluidStatePhaseEnum::LIQUID:
    {
      Xair = Z;
      Yair = 0.0;
      Xh2o = 1.0 - Z;
      Yh2o = 0.0;
      Moose::derivInsert(Xair.derivatives(), _pidx, 0.0);
      Moose::derivInsert(Xair.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(Xair.derivatives(), _Zidx, 1.0);
      break;
    }

    case FluidStatePhaseEnum::GAS:
    {
      Xair = 0.0;
      Yair = Z;
      Yh2o = 1.0 - Z;
      Moose::derivInsert(Yair.derivatives(), _pidx, 0.0);
      Moose::derivInsert(Yair.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(Yair.derivatives(), _Zidx, 1.0);
      break;
    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      // Keep equilibrium mass fractions
      Xh2o = 1.0 - Xair;
      break;
    }
  }

  // Save the mass fractions in the FluidStateMassFractions object
  liquid.mass_fraction[_aqueous_fluid_component] = Xh2o;
  liquid.mass_fraction[_gas_fluid_component] = Xair;
  gas.mass_fraction[_aqueous_fluid_component] = Yh2o;
  gas.mass_fraction[_gas_fluid_component] = Yair;
}

void
PorousFlowWaterAir::gasProperties(const DualReal & pressure,
                                  const DualReal & temperature,
                                  std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  const DualReal psat = _water_fp.vaporPressure(temperature);

  const DualReal Yair = gas.mass_fraction[_gas_fluid_component];
  const DualReal Xair = liquid.mass_fraction[_gas_fluid_component];

  // Air density, viscosity and enthalpy calculated using partial pressure
  // Yair * gas_poreressure (Dalton's law)
  DualReal air_density, air_viscosity;
  _air_fp.rho_mu_from_p_T(Yair * pressure, temperature, air_density, air_viscosity);
  DualReal air_enthalpy = _air_fp.h_from_p_T(Yair * pressure, temperature);

  // Vapor density, viscosity and enthalpy calculated using partial pressure
  // X1 * psat (Raoult's law)
  DualReal vapor_density, vapor_viscosity;

  _water_fp.rho_mu_from_p_T((1.0 - Xair) * psat, temperature, vapor_density, vapor_viscosity);
  DualReal vapor_enthalpy = _water_fp.h_from_p_T((1.0 - Xair) * psat, temperature);

  // Density is just the sum of individual component densities
  gas.density = air_density + vapor_density;

  // Viscosity of the gas phase is a weighted sum of the individual viscosities
  gas.viscosity = Yair * air_viscosity + (1.0 - Yair) * vapor_viscosity;

  // Enthalpy of the gas phase is a weighted sum of the individual enthalpies
  gas.enthalpy = Yair * air_enthalpy + (1.0 - Yair) * vapor_enthalpy;

  //  Internal energy of the gas phase (e = h - pv)
  mooseAssert(gas.density.value() > 0.0, "Gas density must be greater than zero");
  gas.internal_energy = gas.enthalpy - pressure / gas.density;
}

void
PorousFlowWaterAir::liquidProperties(const DualReal & pressure,
                                     const DualReal & temperature,
                                     std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];

  // Calculate liquid density and viscosity if in the two phase or single phase
  // liquid region, assuming they are not affected by the presence of dissolved
  // Air (why? effect of dissolved air should be considered!!! Need to go over this!).
  // Note: the (small) contribution due to derivative of capillary pressure
  // wrt pressure (using the chain rule) is not implemented.
  DualReal liquid_density, liquid_viscosity;
  _water_fp.rho_mu_from_p_T(pressure, temperature, liquid_density, liquid_viscosity);

  liquid.density = liquid_density;
  liquid.viscosity = liquid_viscosity;

  // Enthalpy does include a contribution due to the enthalpy of dissolution
  const DualReal hdis = enthalpyOfDissolution(temperature);

  const DualReal water_enthalpy = _water_fp.h_from_p_T(pressure, temperature);
  const DualReal air_enthalpy = _air_fp.h_from_p_T(pressure, temperature);

  const DualReal Xair = liquid.mass_fraction[_gas_fluid_component];
  liquid.enthalpy = (1.0 - Xair) * water_enthalpy + Xair * (air_enthalpy + hdis);

  //  Internal energy of the liquid phase (e = h - pv)
  mooseAssert(liquid.density.value() > 0.0, "Liquid density must be greater than zero");
  liquid.internal_energy = liquid.enthalpy - pressure / liquid.density;
}

DualReal
PorousFlowWaterAir::liquidDensity(const DualReal & pressure, const DualReal & temperature) const
{
  return _water_fp.rho_from_p_T(pressure, temperature);
}

DualReal
PorousFlowWaterAir::gasDensity(const DualReal & pressure,
                               const DualReal & temperature,
                               std::vector<FluidStateProperties> & fsp) const
{
  auto & liquid = fsp[_aqueous_phase_number];
  auto & gas = fsp[_gas_phase_number];

  DualReal psat = _water_fp.vaporPressure(temperature);

  const DualReal Yair = gas.mass_fraction[_gas_fluid_component];
  const DualReal Xair = liquid.mass_fraction[_gas_fluid_component];

  DualReal air_density = _air_fp.rho_from_p_T(Yair * pressure, temperature);
  DualReal vapor_density = _water_fp.rho_from_p_T((1.0 - Xair) * psat, temperature);

  // Density is just the sum of individual component densities
  return air_density + vapor_density;
}

DualReal
PorousFlowWaterAir::saturation(const DualReal & pressure,
                               const DualReal & temperature,
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

  // The gas and liquid densities
  const DualReal gas_density = gasDensity(pressure, temperature, fsp);
  const DualReal liquid_density = liquidDensity(pressure, temperature);

  // Set mass equilibrium constants used in the calculation of vapor mass fraction
  const DualReal Xair = liquid.mass_fraction[_gas_fluid_component];
  const DualReal Yair = gas.mass_fraction[_gas_fluid_component];

  const DualReal K0 = Yair / Xair;
  const DualReal K1 = (1.0 - Yair) / (1.0 - Xair);
  const DualReal vapor_mass_fraction = vaporMassFraction(Z, K0, K1);

  // The gas saturation in the two phase case
  const DualReal saturation = vapor_mass_fraction * liquid_density /
                              (gas_density + vapor_mass_fraction * (liquid_density - gas_density));

  return saturation;
}

void
PorousFlowWaterAir::twoPhaseProperties(const DualReal & pressure,
                                       const DualReal & temperature,
                                       const DualReal & Z,
                                       unsigned int qp,
                                       std::vector<FluidStateProperties> & fsp) const
{
  auto & gas = fsp[_gas_phase_number];

  // Calculate all of the gas phase properties, as these don't depend on saturation
  gasProperties(pressure, temperature, fsp);

  // The gas saturation in the two phase case
  gas.saturation = saturation(pressure, temperature, Z, fsp);

  // The liquid pressure and properties can now be calculated
  const DualReal liquid_pressure = pressure - _pc.capillaryPressure(1.0 - gas.saturation, qp);
  liquidProperties(liquid_pressure, temperature, fsp);
}

void
PorousFlowWaterAir::equilibriumMassFractions(const DualReal & pressure,
                                             const DualReal & temperature,
                                             DualReal & Xair,
                                             DualReal & Yh2o) const
{
  // Equilibrium constants for each component (Henry's law for the Air
  // component, and Raoult's law for water).
  const DualReal Kh = _water97_fp.henryConstant(temperature, _air_henry);
  const DualReal psat = _water_fp.vaporPressure(temperature);

  const DualReal Kair = Kh / pressure;
  const DualReal Kh2o = psat / pressure;

  // The mole fractions for the Air (Ideal gas) component in the two component
  // case can be expressed in terms of the equilibrium constants only
  const DualReal xair = (1.0 - Kh2o) / (Kair - Kh2o);
  const DualReal yair = Kair * xair;

  // Convert mole fractions to mass fractions
  Xair = moleFractionToMassFraction(xair);
  Yh2o = 1.0 - moleFractionToMassFraction(yair);
}

DualReal
PorousFlowWaterAir::moleFractionToMassFraction(const DualReal & xmol) const
{
  return xmol * _Mair / (xmol * _Mair + (1.0 - xmol) * _Mh2o);
}

void
PorousFlowWaterAir::checkVariables(Real temperature) const
{
  // Check whether the input temperature is within the region of validity of this equation
  // of state (T_triple <= T <= T_critical)
  if (temperature < _water_triple_temperature || temperature > _water_critical_temperature)
    mooseException(name() + ": temperature " + Moose::stringify(temperature) +
                   " is outside range 273.16 K <= T <= 647.096 K");
}

DualReal
PorousFlowWaterAir::enthalpyOfDissolution(const DualReal & temperature) const
{
  // Henry's constant
  const DualReal Kh = _water97_fp.henryConstant(temperature, _air_henry);

  DualReal hdis = -_R * temperature * temperature * Kh.derivatives()[_Tidx] / Kh / _Mair;

  // Derivative of enthalpy of dissolution wrt temperature requires the second derivative of
  // Henry's constant wrt temperature. For simplicity, approximate this numerically
  const Real dT = temperature.value() * 1.0e-8;
  const DualReal t2 = temperature + dT;
  const DualReal Kh2 = _water97_fp.henryConstant(t2, _air_henry);

  const Real dhdis_dT =
      (-_R * t2 * t2 * Kh2.derivatives()[_Tidx] / Kh2 / _Mair - hdis).value() / dT;

  hdis.derivatives() = temperature.derivatives() * dhdis_dT;

  return hdis;
}

Real
PorousFlowWaterAir::totalMassFraction(
    Real pressure, Real temperature, Real /* Xnacl */, Real saturation, unsigned int qp) const
{
  // Check whether the input temperature is within the region of validity
  checkVariables(temperature);

  // As we do not require derivatives, we can simply ignore their initialisation
  const DualReal p = pressure;
  const DualReal T = temperature;

  // FluidStateProperties data structure
  std::vector<FluidStateProperties> fsp(_num_phases, FluidStateProperties(_num_components));
  auto & liquid = fsp[_aqueous_phase_number];
  auto & gas = fsp[_gas_phase_number];

  // Calculate equilibrium mass fractions in the two-phase state
  DualReal Xair, Yh2o;
  equilibriumMassFractions(p, T, Xair, Yh2o);

  // Save the mass fractions in the FluidStateMassFractions object to calculate gas density
  const DualReal Yair = 1.0 - Yh2o;
  liquid.mass_fraction[_aqueous_fluid_component] = 1.0 - Xair;
  liquid.mass_fraction[_gas_fluid_component] = Xair;
  gas.mass_fraction[_aqueous_fluid_component] = Yh2o;
  gas.mass_fraction[_gas_fluid_component] = Yair;

  // Gas density
  const Real gas_density = gasDensity(p, T, fsp).value();

  // Liquid density
  const DualReal liquid_pressure = p - _pc.capillaryPressure(1.0 - saturation, qp);
  const Real liquid_density = liquidDensity(liquid_pressure, T).value();

  // The total mass fraction of air (Z) can now be calculated
  const Real Z = (saturation * gas_density * Yair.value() +
                  (1.0 - saturation) * liquid_density * Xair.value()) /
                 (saturation * gas_density + (1.0 - saturation) * liquid_density);

  return Z;
}
