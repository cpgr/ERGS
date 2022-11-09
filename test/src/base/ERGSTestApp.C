//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ERGSTestApp.h"
#include "ERGSApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
ERGSTestApp::validParams()
{
  InputParameters params = ERGSApp::validParams();
  return params;
}

ERGSTestApp::ERGSTestApp(InputParameters parameters) : MooseApp(parameters)
{
  ERGSTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

ERGSTestApp::~ERGSTestApp() {}

void
ERGSTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  ERGSApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"ERGSTestApp"});
    Registry::registerActionsTo(af, {"ERGSTestApp"});
  }
}

void
ERGSTestApp::registerApps()
{
  registerApp(ERGSApp);
  registerApp(ERGSTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
ERGSTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ERGSTestApp::registerAll(f, af, s);
}
extern "C" void
ERGSTestApp__registerApps()
{
  ERGSTestApp::registerApps();
}
