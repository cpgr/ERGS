#include "SBFractureFlowPermeabilityConstFromVar.h"

registerMooseObject("PorousFlowApp", SBFractureFlowPermeabilityConstFromVar);

InputParameters
SBFractureFlowPermeabilityConstFromVar::validParams()
{
  InputParameters params = PorousFlowPermeabilityBase::validParams();
  params.addRequiredCoupledVar("aperture", "The aperture of the fracture");
  params.addClassDescription(
      "This Material calculates the permeability tensor given by the input variable for aperture");
  return params;
}

SBFractureFlowPermeabilityConstFromVar::SBFractureFlowPermeabilityConstFromVar(
    const InputParameters & parameters)
  : PorousFlowPermeabilityBase(parameters),
   _aperture(coupledValue("aperture"))
{
}

void
SBFractureFlowPermeabilityConstFromVar::computeQpProperties()
{

  Real perm;
  perm = std::pow(_aperture[_qp], 3.0) / 12.0;

  RealTensorValue permeability(perm, 0.0, 0.0, 0.0, perm, 0.0, 0.0, 0.0, perm);

  _permeability_qp[_qp] = permeability;
  (*_dpermeability_qp_dvar)[_qp].resize(_num_var, RealTensorValue());
  (*_dpermeability_qp_dgradvar)[_qp].resize(LIBMESH_DIM);

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    (*_dpermeability_qp_dgradvar)[_qp][i].resize(_num_var, RealTensorValue());
}
