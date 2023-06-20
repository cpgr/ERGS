#TOUGH2 econ2N example

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 500
  xmax = 10000
  bias_x = 1.01
  coord_type = 'RZ'
  rz_coord_axis = Y
[]


[FluidProperties]
  [co2]
    type = CO2FluidProperties
  []
  [brine]
    type = BrineFluidProperties
  []
[]


[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zi' 
    number_fluid_phases = 2
    number_fluid_components = 3
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    alpha = 5.099e-5
    m = 0.457
    sat_lr = 0.0
    pc_max = 1e7
  []
  [fs]
   type = PorousFlowBrineCO2
   brine_fp = brine
   co2_fp = co2
  capillary_pressure = pc
  []
[]


[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]


[Variables]
  [pgas]
    initial_condition = 12E6
  []
  [zi]
    initial_condition =  0.0 
    scaling = 1e4
  []
[]


[AuxVariables]
  [Xnacl]
    initial_condition = 0.15
  []
 [sgas]
   order = CONSTANT
   family = MONOMIAL
 []
  [sbrine]
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
  [sbrine]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = sbrine
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
[]


[DiracKernels]
  [heat_source]
    type = PorousFlowSquarePulsePointSource
    point = '0 0 0'
    mass_flux = 1
    variable = zi
  []
[]


[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = 45          
  []
  [brineCo2Properties]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = zi
    temperature_unit = Celsius          
    xnacl = Xnacl
    capillary_pressure = pc
    fluid_state = fs
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.12
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-13 0 0   0 1E-13 0   0 0 1E-13'
  []
  [relperm0]
    type = PorousFlowRelativePermeabilityVG
    m = 0.457
    phase = 0
    s_res = 0.3
    sum_s_res = 0.35
  []
  [relperm1]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 1
    s_res = 0.05
    sum_s_res = 0.35
  []
[]


[BCs]
  [rightwater]
    type = PorousFlowPiecewiseLinearSink
    boundary = 'right'
    variable = pgas
    use_mobility = true
    PorousFlowDictator = dictator
    fluid_phase = 0
    multipliers = '0 1e9'
    PT_shift = '12e6'
    pt_vals = '0 1e9'
    mass_fraction_component = 0
    use_relperm = true
  []
  [rightco2]
    type = PorousFlowPiecewiseLinearSink
    variable = zi
    boundary = 'right'
    use_mobility = true
    PorousFlowDictator = dictator
    fluid_phase = 1
    multipliers = '0 1e9'
    PT_shift = '12e6'
    pt_vals = '0 1e9'
    mass_fraction_component = 1
    use_relperm = true
  []
[]


[VectorPostprocessors]
  [vars]
    type = NodalValueSampler
    sort_by = x
    variable = 'pgas zi Xnacl'
    execute_on = 'timestep_end'
  []
  [auxvars]
    type = ElementValueSampler
    sort_by = x
    variable = 'sgas'
    execute_on = 'timestep_end'
  []
[]


[Postprocessors]
  [Pgas]
    type = PointValue
   point =  '5 0 0'
    variable = pgas
  []
  [sgas]
    type = PointValue
    point =  '5 0 0'
    variable = sgas
  []
  [sbrine]
    type = PointValue
    point =  '5 0 0'
    variable = sbrine
  []
[]


[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres bjacobi lu NONZERO'
  []
[]


[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 1e4
  num_steps = 1e6
  end_time = 8.64E8
  nl_abs_tol = 1e-12
  nl_max_its = 25
  l_max_its = 100
  dtmax = 5e6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 100
  []
[]


[Outputs]
  exodus = true
  sync_times = '5e4 2.64e6 8.64e6 8.64e7 8.64e8'
  [time]
    type = CSV
    sync_only = true
  []
  [spatial]
    type = CSV
    sync_only = true
  []
[]



