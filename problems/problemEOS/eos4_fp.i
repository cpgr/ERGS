#TOUGH2 EOS3/4 example

[Mesh]
 # [efmcube]
 # type = FileMeshGenerator
 # file = eos3.msh
 # []
  type = GeneratedMesh
  dim = 2
  nx = 500
  xmax = 100 #1000
  ny = 1
  ymax = 4.5
  bias_x = 1.01
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
    porous_flow_vars = 'pwater sgas temperature'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    alpha = 8e-5
    m = 0.45
    sat_lr = 1e-3
    pc_max = 5e-8
  []
[]


[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]


[Variables]
  [pwater]
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
  [pgas]
    family = MONOMIAL
    order = CONSTANT
  []
  [swater]
    family = MONOMIAL
    order = CONSTANT
  []
  [xAir]
    initial_condition = 1 
  []
   [yAir]
    initial_condition = 0 
  []
[]


[AuxKernels]
  [pgas]
    type = PorousFlowPropertyAux
    property = pressure
    phase = 1
    variable = pgas
  []
  [swater]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = swater
  []
#  [xAir]
#    type = PorousFlowPropertyAux
#    variable = xAir
#    property = mass_fraction
#    phase = 0
#    fluid_component = 1
#    execute_on = timestep_end
#  []
#  [yWater]
#    type = PorousFlowPropertyAux
#    variable = yWater
#    property = mass_fraction
#    phase = 1
#    fluid_component = 0
#    execute_on = timestep_end
#  []
[]


[Kernels]
  [mass1]
    type = PorousFlowMassTimeDerivative
    variable = pwater
    fluid_component = 0
  []
  [adv1]
    type = PorousFlowAdvectiveFlux
    variable = pwater
    fluid_component = 0
  []
   [disp1]
    type = PorousFlowDispersiveFlux
    variable = pwater
    disp_trans = '0 0'
    disp_long = '0 0'
    fluid_component = 0
    use_displaced_mesh = false
  []
  [mass0]
    type = PorousFlowMassTimeDerivative
    variable = sgas
    fluid_component = 1
  []
  [adv0]
    type = PorousFlowAdvectiveFlux
    variable = sgas
    fluid_component = 1
  []
   [disp0]
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
[]


[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = temperature
  []
  [ppss]
    type = PorousFlow2PhasePS
    phase0_porepressure = pwater
    phase1_saturation = sgas
    capillary_pressure = pc
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'xAir yAir'
  []
  [water]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
  []
  [gas]
    type = PorousFlowSingleComponentFluid
    fp = air
    phase = 1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '20.e-15 0 0   0 20.e-15 0   0 0 20.e-15'
  []
  [relperm0]
    type = PorousFlowRelativePermeabilityVG
    m = 0.45
    phase = 0
    s_res = 9.6e-4
    sum_s_res = 9.6e-4
  []
  [relperm1]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 1
    s_res = 0.05
    sum_s_res = 0.35
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.01
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 800
    density = 2550.0
  [] 
  [diffusivity]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '2.13E-5  2.13E-5   2.13E-5 2.13E-5'
    tortuosity = '0.25 0.25'
  []
  [rock_thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2 0 0  0 2 0  0 0 2'
  []
[]


[BCs]
  [flux]
    type = PorousFlowSink
    boundary = 'left'
    variable = temperature
    use_mobility = false
    use_relperm = false
    mass_fraction_component = 1
    fluid_phase = 0
    flux_function = -3e3
  []
[]


[VectorPostprocessors]
  [vars]
    type = NodalValueSampler
    sort_by = id
    variable = 'pwater sgas temperature'
    execute_on = 'timestep_end'
    outputs = spatial
  []
  [auxvars]
    type = ElementValueSampler
    sort_by = id
    variable = 'swater'
    execute_on = 'timestep_end'
    outputs = spatial
  []
[]


[Postprocessors]
#  [xAir]
#    type = PointValue
#    point =  '1 0 0'
#    variable = xAir
#  []
 # [yWater]
 #   type = PointValue
 #   point =  '1 0 0'
 #   variable = yWater
 # []
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
  [SWater]
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
  dt = 1e4
  num_steps = 1e6
  end_time = 3.56E8
  nl_abs_tol = 1e-20
  nl_max_its = 25
  l_max_its = 100
  dtmax = 5e6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
  []
[]


[Outputs]
  exodus = true
  sync_times = '1e5 2.592e6 3.64e7 3e8'
  [time]
    type = CSV
  []
  [spatial]
    type = CSV
    sync_only = true
  []
[]




