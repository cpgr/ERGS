# vars: pgas zi Xnacl

[Mesh]
[gen]
  type = GeneratedMeshGenerator
  dim = 2
  xmin = 5
  xmax = 1000
  nx = 100
  bias_x = 1.01
  ymin = -500 
  ymax = 0 
  ny = 1
  []
  coord_type = RZ
[]

[FluidProperties]
  [brine]
    type = BrineFluidPropertiesBeta
  []
  [co2]
    type = CO2FluidProperties
  []
  [water]
    type = Water97FluidProperties
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zi Xnacl'
    number_fluid_phases = 3
    number_fluid_components = 3
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  []
  [fs]
    type = PorousFlowBrineSaltCO2
    brine_fp = brine
    water_fp = water
    co2_fp = co2
    capillary_pressure = pc
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0.8  0.8  0.8'
  temperature_unit = Celsius
[]

[Variables]
  [pgas]
    initial_condition = 6e6
  []
  [zi]
  []
  [Xnacl]
    initial_condition = 0.3 
  []
[]

[ICs]
  [zi]
    type = PorousFlowFluidStateIC
    variable = zi
    saturation = 0.45
    gas_porepressure = pgas
    temperature = temperature
    fluid_state = fs
  []
[]

[AuxVariables]
  [temperature]
    initial_condition = 100
  []
  [sgas]
    order = CONSTANT
    family = MONOMIAL
  []
  [sbrine]
   order = CONSTANT
   family = MONOMIAL
  []
  [ssolid]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [sgas]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 1
    variable = sgas
    execute_on = 'TIMESTEP_END'
  []
  [sbrine]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = sbrine
    execute_on = 'TIMESTEP_END'
  []
  [ssolid]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 2
    variable = ssolid
    execute_on = 'TIMESTEP_END'
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    variable = pgas
    fluid_component = 0
  []
  [adv0]
    type = PorousFlowAdvectiveFlux
    variable = pgas
    fluid_component = 0
  []
  [mass1]
    type = PorousFlowMassTimeDerivative
    variable = zi
    fluid_component = 1
  []
  [adv1]
    type = PorousFlowAdvectiveFlux
    variable = zi
    fluid_component = 1
  []
  [mass2]
    type = PorousFlowMassTimeDerivative
    variable = Xnacl
    fluid_component = 2
  []
  [adv2]
    type = PorousFlowAdvectiveFlux
    variable = Xnacl
    fluid_component = 2
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = temperature
  []
  [brineCo2Properties]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = zi
    temperature = temperature
    xnacl = Xnacl
    capillary_pressure = pc
    fluid_state = fs
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.05
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '50e-15 0 0   0 50e-15 0   0 0 50e-15'
  []
  [relperm0]
    type = PorousFlowRelativePermeabilityCorey
    n = 4
    s_res = 0.30
    sum_s_res = 0.35
    phase = 0
  []
  [relperm1]
    type = PorousFlowRelativePermeabilityCorey
    n = 4
    s_res = 0.30
    sum_s_res = 0.35
    phase = 1
  []
  [relperm2]
    type = PorousFlowRelativePermeabilityConst
    kr = 0
    phase = 2
  []
[]


[Postprocessors]
  [Pgas]
    type = ElementAverageValue
    variable = pgas
    execute_on = 'initial TIMESTEP_END'
  []
  [sgas]
    type = ElementAverageValue
    variable = sgas
    execute_on = 'initial TIMESTEP_END'
  []
  [sbrine]
    type = ElementAverageValue
    variable = sbrine
    execute_on = 'initial TIMESTEP_END'
  []
  [ssolid]
    type = ElementAverageValue
    variable = ssolid
    execute_on = 'initial TIMESTEP_END'
  []
  [Xnacl]
    type = ElementAverageValue
    variable = Xnacl
    execute_on = 'initial TIMESTEP_END'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]         
  type = Transient
  solve_type = NEWTON
  dt = 1
  end_time = 1e6  
[]

[Outputs]
[]


