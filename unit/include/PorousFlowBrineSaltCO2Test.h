#pragma once

#include "MooseObjectUnitTest.h"
#include "PorousFlowCapillaryPressureVG.h"
#include "PorousFlowBrineSaltCO2.h"
#include "BrineFluidPropertiesBeta.h"
#include "Water97FluidProperties.h"
#include "NaClFluidProperties.h"
#include "CO2FluidProperties.h"


class PorousFlowBrineSaltCO2Test : public MooseObjectUnitTest
{
public:
  PorousFlowBrineSaltCO2Test() :MooseObjectUnitTest("ERGSApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters pc_params = _factory.getValidParams("PorousFlowCapillaryPressureVG");
    pc_params.set<Real>("m") = 0.5;
    pc_params.set<Real>("alpha") = 1.0e-4;
    _fe_problem->addUserObject("PorousFlowCapillaryPressureVG", "pc", pc_params);
    _pc = &_fe_problem->getUserObject<PorousFlowCapillaryPressureVG>("pc");

    InputParameters brine_params = _factory.getValidParams("BrineFluidPropertiesBeta");
    _fe_problem->addUserObject("BrineFluidPropertiesBeta", "brine_fp", brine_params);
    _brine_fp = &_fe_problem->getUserObject<BrineFluidPropertiesBeta>("brine_fp");

    InputParameters water_params = _factory.getValidParams("Water97FluidProperties");
    _fe_problem->addUserObject("Water97FluidProperties", "water_fp", water_params);
    _water_fp = &_fe_problem->getUserObject<Water97FluidProperties>("water_fp");

    InputParameters co2_params = _factory.getValidParams("CO2FluidProperties");
    _fe_problem->addUserObject("CO2FluidProperties", "co2_fp", co2_params);
    _co2_fp = &_fe_problem->getUserObject<SinglePhaseFluidProperties>("co2_fp");

    InputParameters uo_params = _factory.getValidParams("PorousFlowBrineSaltCO2");
    uo_params.set<UserObjectName>("water_fp") = "water_fp";
    uo_params.set<UserObjectName>("brine_fp") = "brine_fp";
    uo_params.set<UserObjectName>("co2_fp") = "co2_fp";
    uo_params.set<UserObjectName>("capillary_pressure") = "pc";
    _fe_problem->addUserObject("PorousFlowBrineSaltCO2", "fs", uo_params);
    _fs = &_fe_problem->getUserObject<PorousFlowBrineSaltCO2>("fs");

    // Indices for derivatives
    _pidx = _fs->getPressureIndex();
    _Tidx = _fs->getTemperatureIndex();
    _Zidx = _fs->getZIndex();
    _Xidx = _fs->getXIndex();
  }
 
  const PorousFlowCapillaryPressureVG * _pc;
  const Water97FluidProperties * _water_fp;
  const SinglePhaseFluidProperties * _co2_fp;
  const PorousFlowBrineSaltCO2* _fs;
  const BrineFluidPropertiesBeta* _brine_fp;

  unsigned int _pidx;
  unsigned int _Tidx;
  unsigned int _Zidx;
  unsigned int _Xidx;
};
