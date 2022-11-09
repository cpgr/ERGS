#include "ERGSApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
ERGSApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

ERGSApp::ERGSApp(InputParameters parameters) : MooseApp(parameters)
{
  ERGSApp::registerAll(_factory, _action_factory, _syntax);
}

ERGSApp::~ERGSApp() {}

void
ERGSApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"ERGSApp"});
  Registry::registerActionsTo(af, {"ERGSApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
ERGSApp::registerApps()
{
  registerApp(ERGSApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
ERGSApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ERGSApp::registerAll(f, af, s);
}
extern "C" void
ERGSApp__registerApps()
{
  ERGSApp::registerApps();
}
