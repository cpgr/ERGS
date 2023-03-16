#pragma once

#include "BrineFluidProperties.h"
#include "Water97FluidProperties.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Woverloaded-virtual"

/**
 * Updated the Previous BrineFluidProperties to include halite
 * solubility in water according to Driesner and Heinrich (2007) and 
 *the vapor phaseaccording to Palliser and McKibbin (1998).
 **/
 class BrineFluidPropertiesBeta : public BrineFluidProperties
{
public:
  static InputParameters validParams();

  BrineFluidPropertiesBeta(const InputParameters & parameters);
  virtual ~BrineFluidPropertiesBeta();

 
  virtual std::string fluidName() const override;

  Real molarMass(Real xnacl) const;
  FPDualReal molarMass(const FPDualReal & xnacl) const;

  Real molarMassNaCl() const;

  Real molarMassH2O() const;
  
  virtual Real rho_from_p_T_X(Real pressure, Real temperature, Real xnacl) const override;
  using BrineFluidProperties::rho_from_p_T_X;

  FPDualReal rho_from_p_T_X(const FPDualReal & pressure,
                            const FPDualReal & temperature,
                            const FPDualReal & xnacl) const;

  virtual void rho_from_p_T_X(Real pressure,
                              Real temperature,
                              Real xnacl,
                              Real & rho,
                              Real & drho_dp,
                              Real & drho_dT,
                              Real & drho_dx) const override;

  virtual Real mu_from_p_T_X(Real pressure, Real temperature, Real xnacl) const override;
  using BrineFluidProperties::mu_from_p_T_X;

  virtual void mu_from_p_T_X(Real pressure,
                             Real temperature,
                             Real xnacl,
                             Real & mu,
                             Real & dmu_dp,
                             Real & dmu_dT,
                             Real & dmu_dx) const override;

  FPDualReal h_from_p_T_X(const FPDualReal & pressure,
                          const FPDualReal & temperature,
                          const FPDualReal & xnacl) const;

  virtual Real h_from_p_T_X(Real pressure, Real temperature, Real xnacl) const override;
  using BrineFluidProperties::h_from_p_T_X;

  virtual void h_from_p_T_X(Real pressure,
                            Real temperature,
                            Real xnacl,
                            Real & h,
                            Real & dh_dp,
                            Real & dh_dT,
                            Real & dh_dx) const override;

  virtual Real cp_from_p_T_X(Real pressure, Real temperature, Real xnacl) const override;

  FPDualReal e_from_p_T_X(const FPDualReal & pressure,
                          const FPDualReal & temperature,
                          const FPDualReal & xnacl) const;

  virtual Real e_from_p_T_X(Real pressure, Real temperature, Real xnacl) const override;

  virtual void e_from_p_T_X(Real pressure,
                            Real temperature,
                            Real xnacl,
                            Real & e,
                            Real & de_dp,
                            Real & de_dT,
                            Real & de_dx) const override;

  virtual Real k_from_p_T_X(Real pressure, Real temperature, Real xnacl) const override;

     Real vaporPressure(Real temperature, Real xnacl) const;


     Real henryConstant(Real temperature, const std::vector<Real>& coeffs) const;
     void
         henryConstant(Real temperature, const std::vector<Real>& coeffs, Real& Kh, Real& dKh_dT) const;
     DualReal henryConstant(const DualReal& temperature, const std::vector<Real>& coeffs) const;

     /// Fluid component numbers for water and NaCl
     static const unsigned int WATER = 0;
     static const unsigned int NACL = 1;

     virtual const SinglePhaseFluidProperties& getComponent(unsigned int component) const override;


  /** 
      Old halite solubility from Potter et al. (1977): A new method for determining the solubility of salts 
      in aqueous solutions at elevated temperatures. J. Res. US Geol. Surv., 5:389–395, 1977
  **/
      Real haliteSolubility(Real temperature) const;

  /**
   * Updated the Solubility of halite (solid NaCl) in water
   * from Potter et al., (i.e., haliteSolubility above) to Driesner and Heinrich (2007):
   * The system H2O–NaCl. Part I pg 4891. This new solubility is referred to as haliteSolubilityWater
   *
   * @param temperature temperature (K)
   * @return halite solubility (kg/kg)
   */
      Real haliteSolubilityWater (Real temperature, Real pressure) const;
 
  /**
  *  Halite (solid NaCl) solubility in the gas phase from Palliser and McKibbin (1998):
  * A Model for Deep Geothermal Brines, I: T-p-X State-Space Description pg 78
  *
  * @param temperature temperature (K)
  * @return halite solubility (kg/kg)
  */
     Real haliteSolubilityGas(Real temperature, Real pressure) const;

protected:
 
  Real massFractionToMolalConc(Real xnacl) const;

  Real massFractionToMoleFraction(Real xnacl) const;
 
  FPDualReal massFractionToMoleFraction(const FPDualReal & xnacl) const;
  
/// Flag to indicate whether to calculate derivatives in water_fp
  mutable bool _water_fp_derivs;

  using BrineFluidProperties::_Mh2o;
  using BrineFluidProperties::_Mnacl;
  using BrineFluidProperties::_nacl_fp;
  using BrineFluidProperties::_water_fp;
  using BrineFluidProperties::_water97_fp;
  


  //Parameter for Driesner
  Real _alpha;

  //Parameters for Palliser
  Real _y;
  Real _XgSat;
  Real _Psat;

  //Params for Nacl
  const Real _T_triple;
  const Real _p_triple;

};

#pragma GCC diagnostic pop
