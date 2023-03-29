#include "BrineFluidPropertiesBeta.h"

registerMooseObject("ERGSApp", BrineFluidPropertiesBeta);

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
    _p_triple(50.0)     // note, pressure is in Pa. need to convert to bar in the function
{
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

  Real T_hm = tripleTc + _alpha * (pbar - triplePbar);

 // Parameters from Driesner (2017):
  Real  e0,e1, e2, e3, e4, e5;

  e0 = 0.0989944 + (3.30796 * 1e-6 * pbar) - (4.71759 * 1e-10* pbar * pbar);
  e1 = 0.00947257 - (8.66460 * 1e-6 * pbar) + (1.69417 * 1e-9 * pbar * pbar);
  e2 = 0.610863 - (1.51716 * 1e-5 * pbar) + (1.19290 * 1e-8 * pbar * pbar);
  e3 = -1.64994 + (2.03441 * 1e-4 * pbar) - (6.46015 * 1e-8 * pbar * pbar);
  e4 = 3.36474 - (1.54023 * 1e-4 * pbar) + (8.17048 * 1e-8 * pbar * pbar);
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
BrineFluidPropertiesBeta::haliteSolubilityGas(/* Real temperature,*/ Real pressure) const
{
 // This correlation requires pressure in bar
  Real pbar = pressure * 1.0e-5;

  return _XgSat * std::pow((pbar/_Psat),_y);
}


DualReal
BrineFluidPropertiesBeta::haliteSolubilityGas(/*DualReal temperature,*/ DualReal pressure) const
{
  // This correlation requires temperature in Celcius
//  DualReal Tc = temperature - _T_c2k;
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
