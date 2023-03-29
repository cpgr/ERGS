#pragma once

#include "MooseObjectUnitTest.h"

#include "BrineFluidPropertiesBeta.h"


class BrineFluidPropertiesBeta;
class SinglePhaseFluidProperties;

class BrineFluidPropertiesBetaTest : public MooseObjectUnitTest
{
public:
    BrineFluidPropertiesBetaTest() : MooseObjectUnitTest("ERGSApp") { buildObjects(); }

protected:
    void buildObjects()
    {
        InputParameters uo_pars = _factory.getValidParams("BrineFluidPropertiesBeta");
        _fe_problem->addUserObject("BrineFluidPropertiesBeta", "fp", uo_pars);
        _fp = &_fe_problem->getUserObject<BrineFluidPropertiesBeta>("fp");
        _water_fp = &_fp->getComponent(BrineFluidPropertiesBeta::WATER);
    }

    const BrineFluidPropertiesBeta* _fp;
    const SinglePhaseFluidProperties* _water_fp;
};