#include "PorousFlowBrineSaltNCG.h"
#include "BrineFluidPropertiesBeta.h"
#include "SinglePhaseFluidProperties.h"
#include "MathUtils.h"
#include "Conversion.h"
#include "Water97FluidProperties.h"
#include "NaClFluidProperties.h"

registerMooseObject("ERGSApp", PorousFlowBrineSaltNCG);

InputParameters
PorousFlowBrineSaltNCG::validParams()
{
  InputParameters params = PorousFlowBrineSaltCO2::validParams();
 params.addRequiredParam<UserObjectName>(
    "gas_fp", "The name of the user object for the non-condensable gas");
    params.addParam<Real>("R", 8314.56,"Universal Gas Constant");
  params.addClassDescription("Fluid state class for brine, salt and non-condensable gas."
                            "Includes the dissolution/precipitation of the solid-salt/halite");
  return params;
}

PorousFlowBrineSaltNCG::PorousFlowBrineSaltNCG(const InputParameters & parameters)
  : PorousFlowBrineSaltCO2(parameters),
   _ncg_fp(getUserObject<SinglePhaseFluidProperties>("gas_fp")),
   _Mh2o(_water_fp.molarMass()),
   _Mncg(_ncg_fp.molarMass()),
   _water_triple_temperature(_water_fp.triplePointTemperature()),
   _water_critical_temperature(_water_fp.criticalTemperature()),
   _ncg_henry(_ncg_fp.henryCoefficients()),
   _R(getParam<Real>("R"))
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
PorousFlowBrineSaltNCG::thermophysicalProperties(Real pressure,
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
  DualReal Zncg = Z;
  Moose::derivInsert(Zncg.derivatives(), _Zidx, 1.0);
  DualReal X = Xnacl;
  Moose::derivInsert(X.derivatives(), _Xidx, 1.0);

  // Clear all of the FluidStateProperties data
  clearFluidStateProperties(fsp);

  FluidStatePhaseEnum phase_state;
  massFractions(p, T, X, Zncg, phase_state, fsp);

  switch (phase_state)
  {
    case FluidStatePhaseEnum::GAS:
    {
//      _console << "phase state gas " << std::endl;

      // Set the gas saturations
      // note: since it is single gas phase, gas saturation is 1. It
      // is not used in the gas properties, instead the total mass
      // fraction is used.
      gas.saturation = 1.0;

      // Calculate gas properties
      gasProperties(p, T, Xnacl,Z, fsp);

      break;
    }

    case FluidStatePhaseEnum::LIQUID:
    {
//      _console << "phase state liquid " << std::endl;

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
//      _console << "phase state multiphase " << std::endl;

      // Calculate the gas and liquid properties in the two phase region
      // note: Here, we have multiphase. Therefore, we need the saturation
      // (computed internally) together with the various equilibrium mass fractions
      MultiPhaseProperties(p, T, X, Zncg, qp, fsp);

      break;
    }
  }

  // set the solid/halite saturation
  solid.saturation = saturationSOLID(pressure, temperature, Xnacl, fsp);
  // Liquid saturations can now be set
  liquid.saturation = 1.0 - gas.saturation - solid.saturation;

  // Save pressures to FluidStateProperties object
  gas.pressure = p;
  liquid.pressure = p - _pc.capillaryPressure(liquid.saturation, qp);
  solid.pressure = p; // pressure in solid phase is the same as gas phase. no Pc
}

void
PorousFlowBrineSaltNCG::massFractions(const DualReal & pressure,
                                      const DualReal & temperature,
                                      const DualReal & Xnacl,
                                      const DualReal & Z,
                                      FluidStatePhaseEnum & phase_state,
                                      std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];
  FluidStateProperties solid = fsp[_solid_phase_number];

  DualReal Xncg = 0.0;
  DualReal Yh2o = 0.0;
  DualReal Yncg = 0.0;
  DualReal Ynacl = 0.0;  //halite in gas phase
  DualReal Snacl = 0.0;  //halite in solid phase
  DualReal Xh2o = 0.0;

 // halite solubility in liquid phase
 DualReal XEQ = _brine_fp.haliteSolubilityWater(temperature,pressure);
 //  _console << "XEQ " << XEQ << std::endl;
 const DualReal R = 1.0 - Z;
 /*
 /Check the salt composition in the TWOPHASE state. That is, whether the amount
 / of salt (or salt mass fraction) in the TWOPHASE is less than its concentration
 / or solubility limit even after evaporation and precipitation of brine.
 */
 if (Xnacl < XEQ)
   {
  /*
  / If true, the amount of salt present is not enough to change the TWOPHASE state.
  / TWOPHASE state is still liquid & gas (i.e., Brine-ncg). Proceed to check the
  / whether only liquid brine or gaseous ncg exists.
  /
  / check if the total pressure in the system is greater than the boiling
  / pressure. If true, boiling does not start and there is no gas phase.
  / ONlY liquid phase exist. Note: boiling pressure is the pressure exerted by
  / the gas phase due when boiling starts.
  */
     if (Z < _Zmin)
   // only liquid phase exist.
    {
      phase_state = FluidStatePhaseEnum::LIQUID;
    }
  // if false, probably there exist either gas or both liquid and gas. check the
  // brine vapor pressure against the steam (i.e., water vapour) pressure. if the
  // former is greater, there is no liquid. only gas phase:
     else if (Z > R)
   {
     // only gas phase exist.
       phase_state = FluidStatePhaseEnum::GAS;
   }
     else
   {
  // There is TWOPHASE. The equilibrium mass fraction of ncg in liquid (Yncg) and H2O
  // in gas phases (Xh20) can now be computed. Note: There is NO contribution of salt
  // in both liquid and gas phases (i.e., Xnacl and Ynacl is not accounted for)!
  equilibriumMassFractions(pressure, temperature, Xnacl, Z, Xncg, Yh2o,Ynacl,Snacl,fsp);
    Yncg = 1.0 - Yh2o;
    Xh2o = 1.0 - Xncg;
  // Determine which phases are present based on the value of Z
    phaseState(Z.value(), Xncg.value(), Yncg.value(), phase_state);
   }
 }
else
// Xnacl > or = XEQ. Hence, the amount of salt present after evaporation and
// precipitation is enough to include the solid phase. The system could be
// either TWO-PHASE (i.e., S+G, S+L) or THREE-PHASE (i.e., S+L+G)
{
  // Set Salt mass fraction equal to the solid phase saturation to
  // maintain the presence of salt.
  //solid.saturation = saturationSOLID(pressure, temperature, Xnacl, fsp);
  //const DualReal Xnacl = solid.saturation;
  //_console << "solid.saturation = " << solid.saturation << std::endl;
  //_console << "Xnacl = " << Xnacl << std::endl;

  /** Check whether there is S+L using the persistent variable (Z)!**/
   if (Z < _Zmin)
   {
  // The existing two-phases are solid and liquid. No Gas. Compute h2o in the
  // liquid phase. Note: salt component (i.e. Xnacl) is accounted for and Snacl =1.
  equilibriumMassFractions(pressure, temperature, Xnacl, Z, Xncg, Yh2o,Ynacl,Snacl,fsp);
   Xh2o = 1.0 - Xncg -  Xnacl;
   // Determine which phases are present based on the value of z
   // note: Snacl = 1 and since there is no ncg in the liquid phase, Sncg = 0.
   phaseState(Z.value(), Xncg.value(), 0.0, phase_state);
   }
   /** No S+l? Check for S+G **/
   else if (Z > R)
  {
   // Only two-phase (i.e., solid and gas) exist. No liquid-phase. Use the
   // equilibrium mass fraction to compute the ncg in the gas phase. Note: salt
   // component (i.e., Ynacl) is account for  and Snacl is already initialize to 1.
  equilibriumMassFractions(pressure, temperature, Xnacl, Z, Xncg, Yh2o,Ynacl,Snacl,fsp);
   Yncg = 1.0 - Yh2o - Ynacl;
   // Determine which phases are present based on the value of z
   // note: Snacl = 1 and since there is no ncg in the solid phase, Sncg = 0.
   phaseState(Z.value(), 0.0, Yncg.value(), phase_state);
   }
    else
  {
  // Equilibrium mass fraction of ncg in liquid (Yncg) and h20 in gas phases (Xh20)
  // can now be computed. Note: Here,the contribution of salt in both liquid and gas
  // (i.e., Xnacl and Ynacl) is now included. Snacl is already initialize to 1.
  equilibriumMassFractions(pressure, temperature, Xnacl, Z, Xncg, Yh2o,Ynacl,Snacl,fsp);
  Yncg = 1.0 - Yh2o - Ynacl; // ncg in gas phase corrected with halite in the gas phase
  Xh2o = 1.0 - Xncg - Xnacl; // h2o in liquid corrected with halite in the liquid phase
  // Determine which phases are present based on the value of z
  phaseState(Z.value(), Xncg.value(), Yncg.value(), phase_state);
  }
}
  // The equilibrium mass fractions calculated above are only correct in the two phase
  // state. If only liquid or gas phases are present, the mass fractions are given by
  // the total mass fraction z which includes the salt/halite mass fraction too

  switch (phase_state)
  {
    case FluidStatePhaseEnum::LIQUID:
    {
  //  DualReal Xnacl;
// note: the total mass fraction of 'gas', Z, assigned to the liquid phase is taken
// by ncg (i.e., Xncg). The mass fraction of h20 left in the liquid phase (i.e., Xh20)
// is 1-Z. Salt in the liquid phase (Xnacl) is already initialized. The rest of the
// phases do not exist in this pure liquid phase case.
        Xncg = Z;
        Yncg = 0.0;
        Xh2o = 1.0 - Z;
        Yh2o = 0.0;
      Moose::derivInsert(Xncg.derivatives(), _pidx, 0.0);
      Moose::derivInsert(Xncg.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(Xncg.derivatives(), _Xidx, 0.0);
      Moose::derivInsert(Xncg.derivatives(), _Zidx, 1.0);
      Moose::derivInsert(Yncg.derivatives(), _pidx, 0.0);
      Moose::derivInsert(Yncg.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(Yncg.derivatives(), _Xidx, 0.0);
      break;
    }

    case FluidStatePhaseEnum::GAS:
    {
    DualReal Xnacl;
    DualReal Ynacl;
    DualReal Snacl;
      Xncg = 0.0;
      Yncg = Z;
 // same here, total mass fraction of 'gas' assign is taken by ncg in the gas phase.
 // the rest is taken by water in gas phase.
      Yh2o = 1.0 - Z;
      Xh2o = 0.0;
      Ynacl = 0.0;     //note: no salt in the (pure) gas phase.
      Snacl = 0.0;     //salt in liquid phase (Xnacl) is already initialized to zero.
      Moose::derivInsert(Xncg.derivatives(), _pidx, 0.0);
      Moose::derivInsert(Xncg.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(Xncg.derivatives(), _Xidx, 0.0);
      Moose::derivInsert(Yncg.derivatives(), _pidx, 0.0);
      Moose::derivInsert(Yncg.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(Yncg.derivatives(), _Xidx, 0.0);
      Moose::derivInsert(Yncg.derivatives(), _Zidx, 1.0);
      break;

    }
    case FluidStatePhaseEnum::TWOPHASE:
    {
  // Keep equilibrium mass fraction
      break;
    }
  }
 /// Save the mass fractions in the FluidStateProperties object
 // Note: There is no mass fraction of h20 and ncg in the solid phase (i.e., Sh2o & Sncg)
 // because adsorption of gas and liquid onto the solid mass is negligible.
  liquid.mass_fraction[_aqueous_fluid_component] = Xh2o;
  liquid.mass_fraction[_gas_fluid_component] = Xncg;
  liquid.mass_fraction[_salt_component] = Xnacl;
  gas.mass_fraction[_aqueous_fluid_component] = Yh2o;
  gas.mass_fraction[_gas_fluid_component] = Yncg;
  gas.mass_fraction[_salt_component] = Ynacl;
  solid.mass_fraction[_salt_component] = Snacl;
}

void
PorousFlowBrineSaltNCG::gasProperties(const DualReal & pressure,
                                      const DualReal & temperature,
                                      const DualReal & Xnacl,
                                      const DualReal & Z,
                                     std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

 /// included correction factor to account for water vapor in the gas phase
 /// note: since halite is purely solid, it's density and enthalpy does not
 /// affect the gas properties/behaviour.

  const DualReal psat = _brine_fp.vaporPressure(temperature,Xnacl);
// note: vapour pressure lowering effect (i.e.,  f_vpl * psat) not included yet.

  const DualReal Yncg = gas.mass_fraction[_gas_fluid_component];
  const DualReal Xncg = liquid.mass_fraction[_gas_fluid_component];

  // Gas density, viscosity and enthalpy are calculated using partial pressure
  // Yncg * gas_poreressure (Dalton's law)
  DualReal ncg_density, ncg_viscosity;
  _ncg_fp.rho_mu_from_p_T(Yncg * pressure, temperature, ncg_density, ncg_viscosity);

  DualReal ncg_enthalpy = _ncg_fp.h_from_p_T(Yncg * pressure, temperature);

  // Vapor density, viscosity and enthalpy calculated using partial pressure
  // X1 * psat (Raoult's law)

  DualReal vapor_density, vapor_viscosity;

  _water_fp.rho_mu_from_p_T((1.0 - Xncg) * psat, temperature, vapor_density, vapor_viscosity);
  DualReal vapor_enthalpy = _water_fp.h_from_p_T((1.0 - Xncg) * psat, temperature);

  /// Save the values to the FluidStateProperties object. Note that derivatives wrt z are 0
  // Density is just the sum of individual component densities
  gas.density = ncg_density + vapor_density;
  // Viscosity of the gas phase is a weighted sum of the individual viscosities
  gas.viscosity = Yncg * ncg_viscosity + (1.0 - Yncg) * vapor_viscosity;
  // Enthalpy of the gas phase is a weighted sum of the individual enthalpies
  gas.enthalpy = Yncg * ncg_enthalpy + (1.0 - Yncg) * vapor_enthalpy;

  //  Internal energy of the gas phase (e = h - pv)
  mooseAssert(gas.density.value() > 0.0, "Gas density must be greater than zero");
  gas.internal_energy = gas.enthalpy - pressure / gas.density;
}

void
PorousFlowBrineSaltNCG::liquidProperties(const DualReal & pressure,
                                         const DualReal & temperature,
                                         const DualReal & Xnacl,
                                         std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];

  // The liquid density includes the density increase due to dissolved ncg (liquid)
  /// note: since halite is purely solid, it's density and enthalpy does not
  /// affect the liquid property/behaviour.
  const DualReal brine_density = _brine_fp.rho_from_p_T_X(pressure, temperature, Xnacl);

  // Mass fraction of ncg in liquid phase
  const DualReal Xncg = liquid.mass_fraction[_gas_fluid_component];

  // The liquid density
  const DualReal ncg_partial_density = partialDensityNCG(temperature);

  const DualReal liquid_density = 1.0 / (Xncg / ncg_partial_density + (1.0 - Xncg) / brine_density);

  // Assume that liquid viscosity is just the brine viscosity
  const DualReal liquid_viscosity = _brine_fp.mu_from_p_T_X(pressure, temperature, Xnacl);

  // Liquid enthalpy (including contribution due to the enthalpy of dissolution)
  const DualReal brine_enthalpy = _brine_fp.h_from_p_T_X(pressure, temperature, Xnacl);

  // Enthalpy of ncg
  const DualReal ncg_enthalpy = _ncg_fp.h_from_p_T(pressure, temperature);

  // Enthalpy of dissolution
/// note: Unlike PorousFlowBrineSaltCO2, the enthalpy of dissolution depends
/// on henry's constant!
  const DualReal hdis = enthalpyOfDissolutionNCG(temperature,Xnacl);

  const DualReal liquid_enthalpy = (1.0 - Xncg) * brine_enthalpy + Xncg * (ncg_enthalpy + hdis);

  // Save the values to the FluidStateProperties object
  liquid.density = liquid_density;
  liquid.viscosity = liquid_viscosity;
  liquid.enthalpy = liquid_enthalpy;

  mooseAssert(liquid.density.value() > 0.0, "Liquid density must be greater than zero");
  liquid.internal_energy = liquid.enthalpy - pressure / liquid.density;
}


void
PorousFlowBrineSaltNCG::solidProperties(const DualReal & pressure,
                                        const DualReal & temperature,
                                        std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & solid = fsp[_solid_phase_number];

  /// note: Pure solid properties. No modifications/corrections needed to account
  /// for ncg and Brine since they do not exist in solid state.

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
PorousFlowBrineSaltNCG::saturationGAS(const DualReal & pressure,
                                      const DualReal & temperature,
                                      const DualReal & Xnacl,
                                      const DualReal & Z,
                                      std::vector<FluidStateProperties> & fsp) const
{
  auto & gas = fsp[_gas_phase_number];
  auto & liquid = fsp[_aqueous_fluid_component];
  FluidStateProperties & solid = fsp[_solid_phase_number];

  // Approximate liquid density as saturation isn't known yet, by using the gas
  // pressure rather than the liquid pressure. This does result in a small error
  // in the calculated saturation, but this is below the error associated with
  // the correlations. A more accurate saturation could be found iteraviely,
  // at the cost of increased computational expense
// note: Gas saturation is computed from vaporMassfraction, gas_density and liquid_density

  // compute the Gas density
// note: Unlike PorousFlowBrineSaltCO2, here, gas density includes vapor_density
// Also: the vapour pressure here is not affected by vpl effect because vpl
// depends gas saturation.

/* const DualReal gas_density = _ncg_fp.rho_from_p_T(pressure, temperature);*/
  const DualReal psat = _brine_fp.vaporPressure(temperature,Xnacl);

  const DualReal Yncg = gas.mass_fraction[_gas_fluid_component];
  const DualReal Xncg = liquid.mass_fraction[_gas_fluid_component];

  DualReal ncg_density = _ncg_fp.rho_from_p_T(Yncg * pressure, temperature);
  DualReal vapor_density = _water_fp.rho_from_p_T((1.0 - Xncg) * psat, temperature);
  const DualReal gas_density = ncg_density + vapor_density;

  // compute an approximate liquid density as saturation isn't known yet
  const DualReal brine_density = _brine_fp.rho_from_p_T_X(pressure, temperature, Xnacl);

  const DualReal ncg_partial_density = partialDensityNCG(temperature);

  // finally the liquid density
  const DualReal liquid_density = 1.0 / (Xncg / ncg_partial_density + (1.0 - Xncg) / brine_density);

/// compute the saturation using vapor_mass_fraction, liquid and gas density.
  // first, set mass equilibrium constants used in the calculation of vapor mass fraction
  const DualReal K0 = Yncg / Xncg;
  const DualReal K1 = (1.0 - Yncg) / (1.0 - Xncg);
  //vapor mass fraction
  const DualReal vapor_mass_fraction = vaporMassFraction(Z, K0, K1);

  // The gas saturation in the two phase case
  const DualReal saturationGAS = vapor_mass_fraction * liquid_density /
                              (gas_density + vapor_mass_fraction * (liquid_density - gas_density));

  return saturationGAS;
}


DualReal
PorousFlowBrineSaltNCG::saturationSOLID(const DualReal & pressure,
                                        const DualReal & temperature,
                                        const DualReal & Xnacl,
                                        std::vector<FluidStateProperties> & fsp) const
{
  auto & liquid = fsp[_aqueous_fluid_component];
  FluidStateProperties & solid = fsp[_solid_phase_number];

  // Solid phase saturation is computed from Xnacl, XEQ, brine and halite densities:
  // compute the halite density
  const DualReal halite_density = _brine_fp.halite_rho_from_p_T(pressure, temperature);

  // compute brine/liquid density
  const DualReal brine_density = _brine_fp.rho_from_p_T_X(pressure, temperature, Xnacl);

  // Mass fraction of ncg in liquid phase
  // note: Xncg is used to compute the saturationSOLID!
  const DualReal Xncg = liquid.mass_fraction[_gas_fluid_component];

  // compute the halite solubility in the liquid phase (XEQ)
  const DualReal XEQ = _brine_fp.haliteSolubilityWater(temperature,pressure);

  // The solid saturation:
  const DualReal saturationSOLID = ((Xnacl - XEQ) * brine_density * (1.0 - Xncg))/
                                    ((halite_density)*(1.0 - Xnacl));

  return saturationSOLID;
}


void
PorousFlowBrineSaltNCG::MultiPhaseProperties(const DualReal & pressure,
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
//  _console << "gas.saturation " << gas.saturation << std::endl;

  /// Calculate  the gas phase properties
  // note: since the gas properties are corrected with the dissolved water vapor,
  // the pressure is also affected by the gas saturation:
  const DualReal gas_pressure = pressure - _pc.capillaryPressure(gas.saturation, qp);
  gasProperties(gas_pressure, temperature, Xnacl, Z, fsp);

  DualReal solid_saturation = saturationSOLID(pressure, temperature, Xnacl, fsp);
//  _console << "solid.saturation " << solid_saturation << std::endl;

  // The liquid pressure and properties can now be calculated
  const DualReal liquid_pressure = pressure - _pc.capillaryPressure(1 - gas.saturation - solid_saturation, qp);
  liquidProperties(liquid_pressure, temperature, Xnacl, fsp);

  // The solid properties (note that the effect of liquid and gas adsorptions is neglected)
  // Pc occurs only in fluids
  solidProperties(pressure, temperature, fsp);
}


void
PorousFlowBrineSaltNCG::equilibriumMassFractions(const DualReal & pressure,
                                                 const DualReal & temperature,
                                                 const DualReal & Xnacl,
                                                 const DualReal & Z,
                                                 DualReal & Xncg,
                                                 DualReal & Yh2o,
                                                 DualReal & Ynacl,
                                                 DualReal & Snacl,
                                                 std::vector<FluidStateProperties> & fsp) const
{
  // Equilibrium constants for each component (Henry's law for the NCG
  // component, and Raoult's law for water).
  const DualReal Kh = _brine_fp.henryConstant(temperature, _ncg_henry);
  const DualReal psat = _brine_fp.vaporPressure(temperature,Xnacl);

  const DualReal Kncg = Kh / pressure;
  const DualReal Kh2o = psat / pressure;

  // The mole fractions for the NCG component in the two component
  // case can be expressed in terms of the equilibrium constants only
/// note: Unlike PorousFlowBrineSaltCO2, the mole fractions are computed
///  differently using Henry's constants!
  const DualReal xncg = (1.0 - Kh2o) / (Kncg - Kh2o);
  const DualReal yncg = Kncg * xncg;
  DualReal yh2o = 1.0 - yncg;
  DualReal xh2o = 1.0 - xncg;
/*
  // Mole fractions at equilibrium
  DualReal xncg, yh2o;
  equilibriumMoleFractions(pressure, temperature, Xnacl, xncg, yh2o);
*/

  // The mass fraction of H2O in gas (assume no salt in gas phase) and
  // derivatives wrt p, T, and X
  Yh2o = yh2o * _Mh2o / (yh2o * _Mh2o + (1.0 - yh2o) * _Mncg);

  // NaCl molality (mol/kg)
  const DualReal mnacl = Xnacl / (1.0 - Xnacl) / _Mnacl;

  // The molality of ncg in 1kg of H2O
  const DualReal mncg = xncg * (2.0 * mnacl + _invMh2o) / (1.0 - xncg);
  // The mass fraction of ncg in brine is then
  const DualReal denominator = (1.0 + mnacl * _Mnacl + mncg * _Mncg);
  Xncg = mncg * _Mncg / denominator;

  // Compute the mass fraction of NaCl/halite in liquid and vapor phase to account
  // for precipitation and dissolution of halite due to persistent boiling of the
  // brine solution

  // halite solubility in liquid phase
  DualReal XEQ = _brine_fp.haliteSolubilityWater(temperature, pressure);
  // halite solubility in water vapour (gas) phase
  DualReal XEQG = _brine_fp.haliteSolubilityGas(pressure);

  // The mass fraction of halite in the (water) vapor phase(Ynacl), computed from
  // the mass fraction in liquid phase (Xnacl).
  // (note: derivatives wrt to p, T, and X is not accounted for):
  //const DualReal
  Ynacl  =  Xnacl * (XEQG/XEQ);

  // halite in the solid phase (Snacl)
  // note: since adsorption is not allowed, Sh2o and Sncg equal 0 or neglected.
  Snacl  =  1.0;
}


DualReal
PorousFlowBrineSaltNCG::henryConstant(const DualReal & temperature, const DualReal & Xnacl) const
{
  // Henry's constant (for dissolution in brine!)
  const DualReal Kh_h2o = _brine_fp.henryConstant(temperature, _ncg_henry);

  // The correction to salt is obtained through the salting out coefficient
  const std::vector<Real> b{1.19784e-1, -7.17823e-4, 4.93854e-6, -1.03826e-8, 1.08233e-11};

  // Need temperature in Celsius
  const DualReal Tc = temperature - _T_c2k;

  DualReal kb = 0.0;
  for (unsigned int i = 0; i < b.size(); ++i)
    kb += b[i] * std::pow(Tc, i);

  // Need salt mass fraction in molality
  const DualReal xmol = Xnacl / (1.0 - Xnacl) / _Mnacl;

  // Henry's constant and its derivative wrt temperature and salt mass fraction
  return Kh_h2o * std::pow(10.0, xmol * kb);
}

 DualReal
 PorousFlowBrineSaltNCG::enthalpyOfDissolutionNCG(const DualReal & temperature,
                                             const DualReal & Xnacl) const
{
  // Henry's constant. note: This includes the effect of Xnacl!
  const DualReal Kh = henryConstant(temperature, Xnacl);

  DualReal hdis = -_R * temperature * temperature * Kh.derivatives()[_Tidx] / Kh /_Mncg;

  // Derivative of enthalpy of dissolution wrt temperature and xnacl requires the second
  // derivatives of Henry's constant. For simplicity, approximate these numerically
  const Real dT = temperature.value() * 1.0e-8;
  const DualReal T2 = temperature + dT;
  const DualReal Kh2 = henryConstant(T2, Xnacl);

  const Real dhdis_dT =
      (-_R * T2 * T2 * Kh2.derivatives()[_Tidx] / Kh2 / _Mncg - hdis).value() / dT;

  const Real dX = Xnacl.value() * 1.0e-8;
  const DualReal X2 = Xnacl + dX;
  const DualReal Kh3 = henryConstant(temperature, X2);

  const Real dhdis_dX =
      (-_R * temperature * temperature * Kh3.derivatives()[_Tidx] / Kh3 / _Mncg - hdis).value() /
      dX;

  hdis.derivatives() = temperature.derivatives() * dhdis_dT + Xnacl.derivatives() * dhdis_dX;

  return hdis;
}

DualReal
PorousFlowBrineSaltNCG::partialDensityNCG(const DualReal & temperature) const
{
  // This correlation uses temperature in C
  const DualReal Tc = temperature - _T_c2k;
  // The parial molar volume
  const DualReal V = 37.51 - 9.585e-2 * Tc + 8.74e-4 * Tc * Tc - 5.044e-7 * Tc * Tc * Tc;

  return 1.0e6 * _Mncg / V;
}


DualReal
PorousFlowBrineSaltNCG::f_vpl(const DualReal & pressure,
                              const DualReal & temperature,
                              const DualReal & Xnacl,
                              const DualReal & Z,
                              std::vector<FluidStateProperties> & fsp) const
{
  auto & gas = fsp[_gas_phase_number];
  auto & liquid = fsp[_aqueous_fluid_component];
  FluidStateProperties & solid = fsp[_solid_phase_number];

  // estimate the active liquid saturation from the gas and solid saturations:
  gas.saturation = saturationGAS(pressure, temperature, Xnacl, Z, fsp);
  liquid.saturation = 1.0 - gas.saturation;
  solid.saturation = saturationSOLID(pressure, temperature, Xnacl, fsp);

  DualReal S_La =  liquid.saturation/(1.0 - solid.saturation);

  // estimate liquid density:
  DualReal brine_density = _brine_fp.rho_from_p_T_X(pressure, temperature, Xnacl);

  DualReal Xncg = liquid.mass_fraction[_gas_fluid_component];

  DualReal ncg_partial_density = partialDensityNCG(temperature);

  liquid.density = 1.0 / (Xncg / ncg_partial_density + (1.0 - Xncg) / brine_density);

 // now, compute the vpl factor using the 'active' liquid saturation:
  DualReal numer = _Mh2o *_pc.capillaryPressure(S_La) ; //note: no qp!
  DualReal deno  = liquid.density * _R ; // * (temperature)

  return std::exp(numer/deno);
}
