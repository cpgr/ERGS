#include "BrineFluidPropertiesBeta.h"

registerMooseObject("FluidPropertiesApp", BrineFluidPropertiesBeta);

InputParameters
BrineFluidPropertiesBeta::validParams()
{
  InputParameters params = BrineFluidProperties::validParams();
  params.addParam<Real>("alpha", 2.47260e-2," Constant Driesner");
  params.addParam<Real>("y", 6.5," correction factor for temperature"
                         "values less than 350 degrees");
  params.addParam<Real>("XgSat",1," Mass fraction of halite in the"
                          "halite-saturated gas");
  params.addParam<Real>("Psat",400,"Saturated Pressure (in bars)"
                         "at the highest temperature");
  params.addClassDescription("Updated the Previous BrineFluidPropertiesBeta to include halite"
                              "solubility in water (Driesner and Heinrich (2007)) and halite "
                              "solubility in saturated watervapour (Palliser and McKibbin (1998))");
  return params;
}

BrineFluidPropertiesBeta::BrineFluidPropertiesBeta(const InputParameters & parameters)
  : BrineFluidProperties(parameters),
    _alpha(Real(getParam<Real>("alpha"))),
    _y(Real(getParam<Real>("y"))),
    _XgSat(Real(getParam<Real>("XgSat"))),
    _Psat(Real(getParam<Real>("Psat"))),
    _T_triple(1073.85), // note, temp is in K. need to convert to dCelsius in function
    _p_triple(50.0)    // note, pressure is in Pa. need to convert to bar in the function

{
//note: Below, the fluid property objects are assigned new names because a fluid
// property object already exist (i.e., BrineFluidProperties) with these combination
// of names (i.e., brine-water, brine-nacl). To differentiate it with the new one  being
// created (i.e., BrineFluidPropertiesBeta), new names are given to water and nacl to
// obtain a different combination of names (i.e., brine-water_new and brine_nacl_new).
  const std::string water_name = name() + ":water_new";
  {
    const std::string class_name = "Water97FluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    if (_tid == 0)
      _fe_problem.addUserObject(class_name, water_name, params);
  }
  _water97_fp = &_fe_problem.getUserObject<Water97FluidProperties>(water_name);

  if (parameters.isParamSetByUser("water_fp"))
  {
    // SinglePhaseFluidPropertiesPT UserObject for water
    _water_fp = &getUserObject<SinglePhaseFluidProperties>("water_fp");

    // Check that a water userobject has actually been supplied
    if (_water_fp->fluidName() != "water_new")
      paramError("water_fp", "A water FluidProperties UserObject must be supplied");
  }
  else
  {
    // Construct a SinglePhaseFluidProperties UserObject for water
    _water_fp = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(water_name);
  }

  // SinglePhaseFluidProperties UserObject for NaCl
  const std::string nacl_name = name() + ":nacl_new";
  {
//    const std::string class_name = "BrineFluidPropertiesBeta";
    const std::string class_name = "NaClFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    if (_tid == 0)
      _fe_problem.addUserObject(class_name, nacl_name, params);
  }
  _nacl_fp = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(nacl_name);

  // Molar mass of NaCl and H20
  _Mnacl = _nacl_fp->molarMass();
  _Mh2o = _water_fp->molarMass();
}

BrineFluidPropertiesBeta::~BrineFluidPropertiesBeta() {}

const SinglePhaseFluidProperties &
BrineFluidPropertiesBeta::getComponent(unsigned int component) const
{
  switch (component)
  {
    case WATER:
      return *_water_fp;

    case NACL:
      return *_nacl_fp;

    default:
      mooseError("BrineFluidPropertiesBeta::getComponent has been provided an incorrect component");
  }
}

std::string
BrineFluidPropertiesBeta::fluidName() const
{
  //note: Actual name of this new userobject is still brine.
  return "brine";
}

FPDualReal
BrineFluidPropertiesBeta::molarMass(const FPDualReal & xnacl) const
{
  return 1.0 / (xnacl / _Mnacl + (1.0 - xnacl) / _Mh2o);
}

Real
BrineFluidPropertiesBeta::molarMass(Real xnacl) const
{
  return 1.0 / (xnacl / _Mnacl + (1.0 - xnacl) / _Mh2o);
}

Real
BrineFluidPropertiesBeta::molarMassNaCl() const
{
  return _Mnacl;
}

Real
BrineFluidPropertiesBeta::molarMassH2O() const
{
  return _Mh2o;
}

FPDualReal
BrineFluidPropertiesBeta::rho_from_p_T_X(const FPDualReal & pressure,
                                     const FPDualReal & temperature,
                                     const FPDualReal & xnacl) const
{
  // The correlation requires the pressure in bar, not Pa.
  FPDualReal pbar = pressure * 1.0e-5;
  FPDualReal pbar2 = pbar * pbar;
  FPDualReal pbar3 = pbar2 * pbar;

  // The correlation requires mole fraction
  const FPDualReal Xnacl = massFractionToMoleFraction(xnacl);

  const FPDualReal n11 = -54.2958 - 45.7623 * std::exp(-9.44785e-4 * pbar);
  const FPDualReal n21 = -2.6142 - 2.39092e-4 * pbar;
  const FPDualReal n22 = 0.0356828 + 4.37235e-6 * pbar + 2.0566e-9 * pbar2;
  const FPDualReal n1x1 = 330.47 + 0.942876 * std::sqrt(pbar) + 0.0817193 * pbar -
                          2.47556e-8 * pbar2 + 3.45052e-10 * pbar3;
  const FPDualReal n2x1 = -0.0370751 + 0.00237723 * std::sqrt(pbar) + 5.42049e-5 * pbar +
                          5.84709e-9 * pbar2 - 5.99373e-13 * pbar3;
  const FPDualReal n12 = -n1x1 - n11;
  const FPDualReal n20 = 1.0 - n21 * std::sqrt(n22);
  const FPDualReal n23 = n2x1 - n20 - n21 * std::sqrt(1.0 + n22);

  // The temperature Tv where the brine has the same molar volume as pure water
  // Note: correlation uses temperature in Celcius
  const FPDualReal n1 = n1x1 + n11 * (1.0 - Xnacl) + n12 * (1.0 - Xnacl) * (1.0 - Xnacl);
  const FPDualReal n2 = n20 + n21 * std::sqrt(Xnacl + n22) + n23 * Xnacl;
  const FPDualReal Tv = n1 + n2 * (temperature - _T_c2k);

  // The density of water at temperature Tv
  // Note: convert Tv to Kelvin to calculate water density
  FPDualReal water_density;
  if (_water_fp_derivs)
  {
    Real rho, drho_dp, drho_dT;
    _water_fp->rho_from_p_T(pressure.value(), Tv.value() + _T_c2k, rho, drho_dp, drho_dT);
    water_density = rho;

    water_density.derivatives() = pressure.derivatives() * drho_dp + Tv.derivatives() * drho_dT;
  }
  else
    water_density = _water_fp->rho_from_p_T(pressure.value(), Tv.value() + _T_c2k);

  // The brine density is given by the water density scaled by the ratio of
  // brine molar mass to pure water molar mass
  return water_density * molarMass(xnacl) / _Mh2o;
}

Real
BrineFluidPropertiesBeta::rho_from_p_T_X(Real pressure, Real temperature, Real xnacl) const
{
  // Initialise the AD value (no derivatives required)
  FPDualReal p = pressure;
  FPDualReal T = temperature;
  FPDualReal x = xnacl;

  _water_fp_derivs = false;
  FPDualReal ad_rho = this->rho_from_p_T_X(p, T, x);

  return ad_rho.value();
}

void
BrineFluidPropertiesBeta::rho_from_p_T_X(Real pressure,
                                     Real temperature,
                                     Real xnacl,
                                     Real & rho,
                                     Real & drho_dp,
                                     Real & drho_dT,
                                     Real & drho_dx) const
{
  // Initialise the AD value and derivatives
  FPDualReal p = pressure;
  Moose::derivInsert(p.derivatives(), 0, 1.0);
  FPDualReal T = temperature;
  Moose::derivInsert(T.derivatives(), 1, 1.0);
  FPDualReal x = xnacl;
  Moose::derivInsert(x.derivatives(), 2, 1.0);

  _water_fp_derivs = true;
  FPDualReal ad_rho = this->rho_from_p_T_X(p, T, x);

  rho = ad_rho.value();
  drho_dp = ad_rho.derivatives()[0];
  drho_dT = ad_rho.derivatives()[1];
  drho_dx = ad_rho.derivatives()[2];
}

Real
BrineFluidPropertiesBeta::mu_from_p_T_X(Real pressure, Real temperature, Real xnacl) const
{
  // Correlation requires molal concentration (mol/kg)
  const Real mol = massFractionToMolalConc(xnacl);
  const Real mol2 = mol * mol;
  const Real mol3 = mol2 * mol;

  // Correlation requires temperature in C
  const Real Tc = temperature - _T_c2k;

  const Real a = 1.0 + 0.0816 * mol + 0.0122 * mol2 + 0.128e-3 * mol3 +
                 0.629e-3 * Tc * (1.0 - std::exp(-0.7 * mol));

  const Real water_viscosity = _water_fp->mu_from_p_T(pressure, temperature);

  return a * water_viscosity;
}

void
BrineFluidPropertiesBeta::mu_from_p_T_X(Real pressure,
                                    Real temperature,
                                    Real xnacl,
                                    Real & mu,
                                    Real & dmu_dp,
                                    Real & dmu_dT,
                                    Real & dmu_dx) const
{
  // Viscosity of water and derivatives wrt pressure and temperature
  Real muw, dmuw_dp, dmuw_dT;
  _water_fp->mu_from_p_T(pressure, temperature, muw, dmuw_dp, dmuw_dT);

  // Correlation requires molal concentration (mol/kg)
  Real mol = massFractionToMolalConc(xnacl);
  Real dmol_dx = 1.0 / ((1.0 - xnacl) * (1.0 - xnacl) * _Mnacl);
  Real mol2 = mol * mol;
  Real mol3 = mol2 * mol;

  // Correlation requires temperature in C
  Real Tc = temperature - _T_c2k;

  Real a = 1.0 + 0.0816 * mol + 0.0122 * mol2 + 0.128e-3 * mol3 +
           0.629e-3 * Tc * (1.0 - std::exp(-0.7 * mol));
  Real da_dx =
      (0.0816 + 0.0244 * mol + 3.84e-4 * mol2 + 4.403e-4 * Tc * std::exp(-0.7 * mol)) * dmol_dx;
  Real da_dT = 0.629e-3 * (1.0 - std::exp(-0.7 * mol));

  mu = a * muw;
  dmu_dp = a * dmuw_dp;
  dmu_dx = da_dx * muw;
  dmu_dT = da_dT * muw + a * dmuw_dT;
}

FPDualReal
BrineFluidPropertiesBeta::h_from_p_T_X(const FPDualReal & pressure,
                                   const FPDualReal & temperature,
                                   const FPDualReal & xnacl) const
{
  FPDualReal q1, q2, q10, q11, q12, q20, q21, q22, q23, q1x1, q2x1, Th;

  // The correlation requires the pressure in bar, not Pa.
  const FPDualReal pbar = pressure * 1.0e-5;
  const FPDualReal pbar2 = pbar * pbar;

  // The correlation requires mole fraction
  const FPDualReal Xnacl = massFractionToMoleFraction(xnacl);

  q11 = -32.1724 + 0.0621255 * pbar;
  q21 = -1.69513 - 4.52781e-4 * pbar - 6.04279e-8 * pbar2;
  q22 = 0.0612567 + 1.88082e-5 * pbar;

  q1x1 = 47.9048 - 9.36994e-3 * pbar + 6.51059e-6 * pbar2;
  q2x1 = 0.241022 + 3.45087e-5 * pbar - 4.28356e-9 * pbar2;

  q12 = -q11 - q1x1;
  q10 = q1x1;

  q20 = 1.0 - q21 * std::sqrt(q22);
  q23 = q2x1 - q20 - q21 * std::sqrt(1.0 + q22);

  q1 = q10 + q11 * (1.0 - Xnacl) + q12 * (1.0 - Xnacl) * (1.0 - Xnacl);
  q2 = q20 + q21 * std::sqrt(Xnacl + q22) + q23 * Xnacl;
  // The temperature Th where the brine has the same enthalpy as pure water
  // Note: correlation uses temperature in Celcius
  Th = q1 + q2 * (temperature - _T_c2k);

  // The brine enthalpy is then given by the enthalpy of water at temperature Th
  // Note: water enthalpy requires temperature in Kelvin
  FPDualReal enthalpy;
  if (_water_fp_derivs)
  {
    Real h, dh_dp, dh_dT;
    _water_fp->h_from_p_T(pressure.value(), Th.value() + _T_c2k, h, dh_dp, dh_dT);
    enthalpy = h;

    enthalpy.derivatives() = pressure.derivatives() * dh_dp + Th.derivatives() * dh_dT;
  }
  else
    enthalpy = _water_fp->h_from_p_T(pressure.value(), Th.value() + _T_c2k);

  return enthalpy;
}

Real
BrineFluidPropertiesBeta::h_from_p_T_X(Real pressure, Real temperature, Real xnacl) const
{
  // Initialise the AD value (no derivatives required)
  FPDualReal p = pressure;
  FPDualReal T = temperature;
  FPDualReal x = xnacl;

  _water_fp_derivs = false;
  return h_from_p_T_X(p, T, x).value();
}

void
BrineFluidPropertiesBeta::h_from_p_T_X(Real pressure,
                                   Real temperature,
                                   Real xnacl,
                                   Real & h,
                                   Real & dh_dp,
                                   Real & dh_dT,
                                   Real & dh_dx) const
{
  // Initialise the AD value and derivatives
  FPDualReal p = pressure;
  Moose::derivInsert(p.derivatives(), 0, 1.0);
  FPDualReal T = temperature;
  Moose::derivInsert(T.derivatives(), 1, 1.0);
  FPDualReal x = xnacl;
  Moose::derivInsert(x.derivatives(), 2, 1.0);

  _water_fp_derivs = true;
  FPDualReal ad_h = h_from_p_T_X(p, T, x);

  h = ad_h.value();
  dh_dp = ad_h.derivatives()[0];
  dh_dT = ad_h.derivatives()[1];
  dh_dx = ad_h.derivatives()[2];
}

Real
BrineFluidPropertiesBeta::cp_from_p_T_X(Real pressure, Real temperature, Real xnacl) const
{
  Real q1, q2, q10, q11, q12, q20, q21, q22, q23, q1x1, q2x1, Th;

  // The correlation requires the pressure in bar, not Pa.
  Real pbar = pressure * 1.0e-5;
  Real pbar2 = pbar * pbar;

  // The correlation requires mole fraction
  Real Xnacl = massFractionToMoleFraction(xnacl);

  q11 = -32.1724 + 0.0621255 * pbar;
  q21 = -1.69513 - 4.52781e-4 * pbar - 6.04279e-8 * pbar2;
  q22 = 0.0612567 + 1.88082e-5 * pbar;

  q1x1 = 47.9048 - 9.36994e-3 * pbar + 6.51059e-6 * pbar2;
  q2x1 = 0.241022 + 3.45087e-5 * pbar - 4.28356e-9 * pbar2;

  q12 = -q11 - q1x1;
  q10 = q1x1;

  q20 = 1.0 - q21 * std::sqrt(q22);
  q23 = q2x1 - q20 - q21 * std::sqrt(1.0 + q22);

  q1 = q10 + q11 * (1.0 - Xnacl) + q12 * (1.0 - Xnacl) * (1.0 - Xnacl);
  q2 = q20 + q21 * std::sqrt(Xnacl + q22) + q23 * Xnacl;
  // The temperature Th where the brine has the same isobaric heat capacity
  // as pure water. Note: correlation uses temperature in Celcius
  Th = q1 + q2 * (temperature - _T_c2k);

  // The brine isobaric heat capacity is then given by the isobaric heat
  // capacity of water at temperature Th multiplied by q2
  // Note: water isobaric heat capacity requires temperature in Kelvin
  return q2 * _water_fp->cp_from_p_T(pressure, Th + _T_c2k);
}

FPDualReal
BrineFluidPropertiesBeta::e_from_p_T_X(const FPDualReal & pressure,
                                   const FPDualReal & temperature,
                                   const FPDualReal & xnacl) const
{
  FPDualReal enthalpy = h_from_p_T_X(pressure, temperature, xnacl);
  FPDualReal density = rho_from_p_T_X(pressure, temperature, xnacl);

  return enthalpy - pressure / density;
}

Real
BrineFluidPropertiesBeta::e_from_p_T_X(Real pressure, Real temperature, Real xnacl) const
{
  Real enthalpy = h_from_p_T_X(pressure, temperature, xnacl);
  Real density = rho_from_p_T_X(pressure, temperature, xnacl);

  return enthalpy - pressure / density;
}

void
BrineFluidPropertiesBeta::e_from_p_T_X(Real pressure,
                                   Real temperature,
                                   Real xnacl,
                                   Real & e,
                                   Real & de_dp,
                                   Real & de_dT,
                                   Real & de_dx) const
{
  // Initialise the AD value and derivatives
  FPDualReal p = pressure;
  Moose::derivInsert(p.derivatives(), 0, 1.0);
  FPDualReal T = temperature;
  Moose::derivInsert(T.derivatives(), 1, 1.0);
  FPDualReal x = xnacl;
  Moose::derivInsert(x.derivatives(), 2, 1.0);

  _water_fp_derivs = true;
  FPDualReal ad_e = e_from_p_T_X(p, T, x);

  e = ad_e.value();
  de_dp = ad_e.derivatives()[0];
  de_dT = ad_e.derivatives()[1];
  de_dx = ad_e.derivatives()[2];
}

Real
BrineFluidPropertiesBeta::k_from_p_T_X(Real pressure, Real temperature, Real xnacl) const
{
  // Correlation requires molal concentration (mol/kg)
  Real mol = massFractionToMolalConc(xnacl);
  // Correlation requires temperature in C
  Real Tc = temperature - _T_c2k;

  Real S = 100.0 * _Mnacl * mol / (1.0 + _Mnacl * mol);
  Real lambdaw = _water_fp->k_from_p_T(pressure, temperature);
  Real lambda = 1.0 - (2.3434e-3 - 7.924e-6 * Tc + 3.924e-8 * Tc * Tc) * S +
                (1.06e-5 - 2.0e-8 * Tc - 1.2e-10 * Tc * Tc) * S * S;

  return lambda * lambdaw;
}

Real
BrineFluidPropertiesBeta::vaporPressure(Real temperature, Real xnacl) const
{
  // Correlation requires molal concentration (mol/kg)
  Real mol = massFractionToMolalConc(xnacl);
  Real mol2 = mol * mol;
  Real mol3 = mol2 * mol;

  Real a = 1.0 + 5.93582e-6 * mol - 5.19386e-5 * mol2 + 1.23156e-5 * mol3;
  Real b = 1.1542e-6 * mol + 1.41254e-7 * mol2 - 1.92476e-8 * mol3 - 1.70717e-9 * mol * mol3 +
           1.0539e-10 * mol2 * mol3;

  // The temperature of pure water at the same pressure as the brine is given by
  Real th20 = std::exp(std::log(temperature) / (a + b * temperature));

  // The brine vapour pressure is then found by evaluating the saturation pressure for pure water
  // using this effective temperature
  return _water_fp->vaporPressure(th20);
}

Real
BrineFluidPropertiesBeta::haliteSolubility(Real temperature) const
{
  // This correlation requires temperature in Celcius
  Real Tc = temperature - _T_c2k;

  return (26.18 + 7.2e-3 * Tc + 1.06e-4 * Tc * Tc) / 100.0;
}


Real
BrineFluidPropertiesBeta::haliteSolubilityWater(Real temperature, Real pressure) const
{
  // This correlation requires temperature in Celcius
  Real Tc = temperature - _T_c2k;
  Real tripleTc = _T_triple - _T_c2k;

 // This correlation requires pressure in bar
  Real pbar = pressure * 1.0e-5;
  Real triplePbar = _p_triple * 1.0e-5;

  Real T_hm = tripleTc + _alpha * (pressure -triplePbar);

 // Parameters from Driesner (2017):
  Real  e0,e1, e2, e3, e4, e5;

  e0 = 0.0989944 + 3.30796 * 1e-6 * pbar - 4.71759 * 1e-10* pbar * pbar;
  e1 = 0.00947257 - 8.66460 * 1e-6 * pbar + 1.69417 * 1e-9 * pbar * pbar;
  e2 = 0.610863 - 1.51716 * 1e-5 * pbar + 1.19290 * 1e-8 * pbar * pbar;
  e3 = -1.64994 + 2.03441 * 1e-4 * pbar - 6.46015 * 1e-8 * pbar * pbar;
  e4 = 3.36474 - 1.54023 * 1e-4 * pbar + 8.17048 * 1e-8 * pbar * pbar;
  e5 = 1.0 - e0 - e1 - e2 - e3 - e4;

  std::vector<Real> E {e0, e1, e2, e3, e4, e5};

Real X_LSol = 0.0;
  for (int i = 0; i < 7; i++)
     {
       X_LSol += E[i] * std::pow((Tc/T_hm),i);
     }

   return X_LSol;
}


DualReal
BrineFluidPropertiesBeta::haliteSolubilityWater(DualReal temperature, DualReal pressure) const
{
  // This correlation requires temperature in Celcius
  DualReal Tc = temperature - _T_c2k;
  DualReal tripleTc = _T_triple - _T_c2k;

 // This correlation requires pressure in bar
  DualReal pbar = pressure * 1.0e-5;
  DualReal triplePbar = _p_triple * 1.0e-5;

  DualReal T_hm = tripleTc + _alpha * (pressure -triplePbar);

 // Parameters from Driesner (2017):
  DualReal  e0,e1, e2, e3, e4, e5;

  e0 = 0.0989944 + 3.30796 * 1e-6 * pbar - 4.71759 * 1e-10* pbar * pbar;
  e1 = 0.00947257 - 8.66460 * 1e-6 * pbar + 1.69417 * 1e-9 * pbar * pbar;
  e2 = 0.610863 - 1.51716 * 1e-5 * pbar + 1.19290 * 1e-8 * pbar * pbar;
  e3 = -1.64994 + 2.03441 * 1e-4 * pbar - 6.46015 * 1e-8 * pbar * pbar;
  e4 = 3.36474 - 1.54023 * 1e-4 * pbar + 8.17048 * 1e-8 * pbar * pbar;
  e5 = 1.0 - e0 - e1 - e2 - e3 - e4;

  std::vector<DualReal> E {e0, e1, e2, e3, e4, e5};

DualReal X_LSol = 0.0;
  for (int i = 0; i < 7; i++)
     {
       X_LSol += E[i] * std::pow((Tc/T_hm),i);
     }

   return X_LSol;
}


Real
BrineFluidPropertiesBeta::haliteSolubilityGas(Real temperature, Real pressure) const
{
  // This correlation requires temperature in Celcius
  Real Tc = temperature - _T_c2k;
 // This correlation requires pressure in bar
  Real pbar = pressure * 1.0e-5;

  return _XgSat * std::pow((pbar/_Psat),_y);
}


DualReal
BrineFluidPropertiesBeta::haliteSolubilityGas(DualReal temperature, DualReal pressure) const
{
  // This correlation requires temperature in Celcius
  DualReal Tc = temperature - _T_c2k;
 // This correlation requires pressure in bar
  DualReal pbar = pressure * 1.0e-5;

  return _XgSat * std::pow((pbar/_Psat),_y);
}



Real
BrineFluidPropertiesBeta::halite_rho_from_p_T(Real pressure, Real temperature) const
{
  // Correlation needs pressure in bar
  Real pbar = pressure * 1.0e-5;
  // Correlation requires temperature in Celcius
  Real Tc = temperature - _T_c2k;

  // Halite density at 0 Pa
  Real density_P0 = 2.17043e3 - 2.4599e-1 * Tc - 9.5797e-5 * Tc * Tc;

  // Halite density as a function of pressure
  Real l = 5.727e-3 + 2.715e-3 * std::exp(Tc / 733.4);

  return density_P0 + l * pbar;
}


DualReal
BrineFluidPropertiesBeta::halite_rho_from_p_T(DualReal pressure, DualReal temperature) const
{
  // Correlation needs pressure in bar
  DualReal pbar = pressure * 1.0e-5;
  // Correlation requires temperature in Celcius
  DualReal Tc = temperature - _T_c2k;

  // Halite density at 0 Pa
  DualReal density_P0 = 2.17043e3 - 2.4599e-1 * Tc - 9.5797e-5 * Tc * Tc;

  // Halite density as a function of pressure
  DualReal l = 5.727e-3 + 2.715e-3 * std::exp(Tc / 733.4);

  return density_P0 + l * pbar;
}


void
BrineFluidPropertiesBeta::halite_rho_from_p_T(
    Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho = this->halite_rho_from_p_T(pressure, temperature);

  // Correlation needs pressure in bar
  Real pbar = pressure * 1.0e-5;
  // Correlation requires temperature in Celcius
  Real Tc = temperature - _T_c2k;

  // Halite density at 0 Pa
  Real ddensity_P0_dT = -2.4599e-1 - 1.91594e-4 * Tc;

  Real l = 5.727e-3 + 2.715e-3 * std::exp(Tc / 733.4);
  Real dl_dT = 2.715e-3 * std::exp(Tc / 733.4) / 733.4;

  drho_dp = l * 1.0e-5;
  drho_dT = ddensity_P0_dT + dl_dT * pbar;
}


Real
BrineFluidPropertiesBeta::halite_h_from_p_T(Real pressure, Real temperature) const
{
  // Correlation needs pressure in bar
  Real pbar = pressure * 1.0e-5;
  // Correlation requires temperature in Celcius
  Real Tc = temperature - _T_c2k;
  // Triple point temperature of water (in C)
  Real Tt = 273.16 - _T_c2k;
  // Triple point presure of water (in bar)
  Real pt = 611.657 * 1.0e-5;

  // Note: the enthalpy of halite is 0 at the triple point of water
  return 8.7664e2 * (Tc - Tt) + 6.4139e-2 * (Tc * Tc - Tt * Tt) +
         8.8101e-5 * (Tc * Tc * Tc - Tt * Tt * Tt) + 44.14 * (pbar - pt);
}


DualReal
BrineFluidPropertiesBeta::halite_h_from_p_T(DualReal pressure, DualReal temperature) const
{
  // Correlation needs pressure in bar
  DualReal pbar = pressure * 1.0e-5;
  // Correlation requires temperature in Celcius
  DualReal Tc = temperature - _T_c2k;
  // Triple point temperature of water (in C)
  DualReal Tt = 273.16 - _T_c2k;
  // Triple point presure of water (in bar)
  DualReal pt = 611.657 * 1.0e-5;

  // Note: the enthalpy of halite is 0 at the triple point of water
  return 8.7664e2 * (Tc - Tt) + 6.4139e-2 * (Tc * Tc - Tt * Tt) +
         8.8101e-5 * (Tc * Tc * Tc - Tt * Tt * Tt) + 44.14 * (pbar - pt);
}


void
BrineFluidPropertiesBeta::halite_h_from_p_T(
    Real pressure, Real temperature, Real & h, Real & dh_dp, Real & dh_dT) const
{
  // Correlation needs pressure in bar
  Real pbar = pressure * 1.0e-5;
  // Correlation requires temperature in Celcius
  Real Tc = temperature - _T_c2k;
  // Triple point temperature of water (in C)
  Real Tt = 273.16 - _T_c2k;
  // Triple point presure of water (in bar)
  Real pt = 611.657 * 1.0e-5;

  // Note: the enthalpy of halite is 0 at the triple point of water
  h = 8.7664e2 * (Tc - Tt) + 6.4139e-2 * (Tc * Tc - Tt * Tt) +
      8.8101e-5 * (Tc * Tc * Tc - Tt * Tt * Tt) + 44.14 * (pbar - pt);

  dh_dp = 44.14 * 1.0e-5;
  dh_dT = 8.7664e2 + 2.0 * 6.4139e-2 * Tc + 3.0 * 8.8101e-5 * Tc * Tc;
}

Real
BrineFluidPropertiesBeta::massFractionToMolalConc(Real xnacl) const
{
  return xnacl / ((1.0 - xnacl) * _Mnacl);
}

Real
BrineFluidPropertiesBeta::massFractionToMoleFraction(Real xnacl) const
{
  // The average molar mass of brine from the mass fraction
  Real Mbrine = molarMass(xnacl);
  // The mole fraction is then
  return xnacl * Mbrine / _Mnacl;
}

FPDualReal
BrineFluidPropertiesBeta::massFractionToMoleFraction(const FPDualReal & xnacl) const
{
  // The average molar mass of brine from the mass fraction
  FPDualReal Mbrine = molarMass(xnacl);
  // The mole fraction is then
  return xnacl * Mbrine / _Mnacl;
}

Real
BrineFluidPropertiesBeta::henryConstant(Real temperature, const std::vector<Real> & coeffs) const
{
  return _water97_fp->henryConstant(temperature, coeffs);
}

void
BrineFluidPropertiesBeta::henryConstant(Real temperature,
                                    const std::vector<Real> & coeffs,
                                    Real & Kh,
                                    Real & dKh_dT) const
{
  _water97_fp->henryConstant(temperature, coeffs, Kh, dKh_dT);
}

DualReal
BrineFluidPropertiesBeta::henryConstant(const DualReal & temperature,
                                    const std::vector<Real> & coeffs) const
{
  Real Kh, dKh_dT;
  henryConstant(temperature.value(), coeffs, Kh, dKh_dT);

  DualReal henry = Kh;
  henry.derivatives() = temperature.derivatives() * dKh_dT;

  return henry;
}
