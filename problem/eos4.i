#TOUGH2 EOS3/4 example

[Mesh]
  [efmcube]
  type = FileMeshGenerator
  file = eos3.msh
  []
[]


[FluidProperties]
  [air]
    type = IdealGasAirFluidProperties
  []
  [water]
    type = Water97FluidProperties
  []
[]


[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas sgas temperature'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    m = 0.45
    alpha = 1e1 
  []
  [fs]
    type = PorousFlowWaterAir
    water_fp = water
    gas_fp = air
    capillary_pressure = pc
  []
[]


[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]


[Variables]
  [pgas]
    initial_condition = 1.0E5
  []
  [sgas]
    initial_condition = 0.2
  []
  [temperature]
    initial_condition = 291.15
  []
[]


[AuxVariables]
  [swater]
   order = CONSTANT
   family = MONOMIAL
  []
  [xAir]
   order = CONSTANT
    family = MONOMIAL
  []
  [Z]
    order = FIRST
    family = LAGRANGE
  []
[]


[AuxKernels]
  [swater]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = swater
  []
  [xAir]
    type = PorousFlowPropertyAux
    variable = xAir
    property = mass_fraction
    phase = 0
    fluid_component = 1
    execute_on = timestep_end
  []
 [Z]
    type = FunctionAux
    variable = Z
    function = sol_variable
 []
[]


[Functions]
  [sol_variable]
   type = ParsedFunction 
   value = (x/sqrt(t))
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
   [disp0]
    type = PorousFlowDispersiveFlux
    variable = pgas
    disp_trans = '0 0'
    disp_long = '0 0'
    fluid_component = 0
    use_displaced_mesh = false
  []
  [mass1]
    type = PorousFlowMassTimeDerivative
    variable = sgas
    fluid_component = 1
  []
  [adv1]
    type = PorousFlowAdvectiveFlux
    variable = sgas
    fluid_component = 1
  []
   [disp1]
    type = PorousFlowDispersiveFlux
    variable = sgas
    disp_trans = '0 0'
    disp_long = '0 0'
    fluid_component = 1
    use_displaced_mesh = false
  []
  [energy]
    type = PorousFlowEnergyTimeDerivative
    variable = temperature
  []
  [heat_adv]
    type = PorousFlowHeatAdvection
    variable = temperature
  []  
  [conduction]
    type = PorousFlowHeatConduction
    variable = temperature
  []  
  [heat_source]
    type = HeatSource
    variable = temperature
    value = 3e3
  []
[]


[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = temperature
  []
  [WaterAirProperties]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = sgas
    temperature = temperature
    capillary_pressure = pc
    fluid_state = fs
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '20.e-15 0 0   0 20.e-15 0   0 0 20.e-15'
  []
  [relperm0]
    type = PorousFlowRelativePermeabilityVG
    m = 0.45
    phase = 0
  []
  [relperm1]
    type = PorousFlowRelativePermeabilityVG
    m = 0.45
    phase = 1
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.01
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 800
    density = 2550
  [] 
  [diffusivity]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '2.13E-5 0 2.13E-5 0'
    tortuosity = '0.25 0.25'
  []
  [rock_thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2 0 0  0 2 0  0 0 2'
  []
[]


[VectorPostprocessors]
  [variables]
    type = LineValueSampler
    variable = 'xAir pgas temperature swater'
    sort_by = id
    start_point = '0 0 0'
    end_point = '1000 0 0'
    num_points = 100
    execute_on = 'timestep_end'
  []
[]


[Postprocessors]
  [Z]
    type = PointValue
    point = '1 0 0'
    variable = Z
  []
  [xAir]
    type = PointValue
    point =  '1 0 0'
    variable = xAir
  []
  [Pgas]
    type = PointValue
   point =  '1 0 0'
    variable = pgas
  []
  [T]
    type = PointValue
    point = '1 0 0'
    variable = temperature
  []
  [swater]
    type = PointValue
    point =  '1 0 0'
    variable = swater
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
  end_time = 10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  [csv]
  type = CSV
  execute_on = 'initial timestep_end'
  []
[]




