#include "EmbeddedFracturePermeabilityComponents.h"

registerMooseObject("PorousFlowApp", EmbeddedFracturePermeabilityComponents);

InputParameters
EmbeddedFracturePermeabilityComponents::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Access components of the permeability tensor.");
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addParam<unsigned int>("fluid_phase", 0, "The index corresponding to the fluid phase");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_i",
      "index_i >= 0 & index_i <= 2",
      "The index i of ij for the permeability tensor (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_j",
      "index_j >= 0 & index_j <= 2",
      "The index j of ij for the permeability tensor (0, 1, 2)");
  return params;
}

EmbeddedFracturePermeabilityComponents::EmbeddedFracturePermeabilityComponents(
    const InputParameters & parameters)
  : AuxKernel(parameters),
    _EmbeddedFracturePermeability(getMaterialProperty<RankTwoTensor>("PorousFlow_permeability_qp")),
    _i(getParam<unsigned int>("index_i")),
    _j(getParam<unsigned int>("index_j")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _ph(getParam<unsigned int>("fluid_phase"))
{
  if (_ph >= _dictator.numPhases())
    paramError("fluid_phase",
               "The Dictator proclaims that the maximum phase index in this simulation is ",
               _dictator.numPhases() - 1,
               " whereas you have used ",
               _ph,
               ". Remember that indexing starts at 0. The Dictator is watching you, to "
               "ensure your wellbeing.");
}

Real
EmbeddedFracturePermeabilityComponents::computeValue()
{
  return _EmbeddedFracturePermeability[_qp](_i, _j);
}
