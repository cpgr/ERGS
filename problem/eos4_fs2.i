# TOUGH2 EOS3/4 example

# NOTE: The difference between eos_fs and eos_fs2 is the way pc and relative permeability are computed.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 500
  xmax = 10000    #100
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
    porous_flow_vars = 'pgas  zi temperature'
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
  [fs]
    type = PorousFlowWaterAir
    water_fp = water
    gas_fp =   air
    capillary_pressure = pc
  []
[]


[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]


[Variables]
  [pgas]
    initial_condition = 1E5
  []
  [zi]
    initial_condition = 0.01 #0.20 #0.14 
  []
  [temperature]
    initial_condition = 291.15 
  []
[]


[AuxVariables]
  [sgas]
   order = CONSTANT
   family = MONOMIAL
  []
  [swater]
   order = CONSTANT
   family = MONOMIAL
  []
  [xAir]
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
  []
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
    variable = zi
    fluid_component = 1
  []
  [adv1]
    type = PorousFlowAdvectiveFlux
    variable = zi
    fluid_component = 1
  []
   [disp1]
    type = PorousFlowDispersiveFlux
    variable = zi
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
  [WaterAirProperties]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = zi
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
    s_res = 9.6e-4
    sum_s_res = 9.6e-4
    phase = 0
  []
  [relperm1]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    s_res = 0.05
    sum_s_res = 0.35
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
    diffusion_coeff = ' 0 2.13E-5 0 2.13E-5 '
    tortuosity = '0.25 0.25'
  []
  [rock_thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2 0 0  0 2 0  0 0 2'
  []
[]


[BCs]
  #[newton]
  #  type = PorousFlowPiecewiseLinearSink
  #  variable = temperature
  #  boundary = left
  #  pt_vals = '0 1 2'
  #  multipliers = '-1 0 1'
  #  flux_function = -3
  #[] 
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
    variable = 'pgas temperature'
    execute_on = 'timestep_end'
    outputs = spatial
  []
  [auxvars]
    type = ElementValueSampler
    sort_by = id
    variable = 'xAir swater'
    execute_on = 'timestep_end'
    outputs = spatial
  []
[]

[Postprocessors]
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
  dt = 1e4
  num_steps = 1e6
  end_time = 3.56E8
  nl_abs_tol = 1e-20
  nl_max_its = 40 #25
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



