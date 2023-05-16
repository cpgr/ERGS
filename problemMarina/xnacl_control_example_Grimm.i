[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  xmax = 0.28028
  ymax = 0.320276
  boundary_id = '0 1 2 3'
  boundary_name = 'bottom right top left'
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zi xnacl'
    number_fluid_phases = 2
    number_fluid_components = 3
  [../]

  [./pc]
   type = PorousFlowCapillaryPressureConst
   pc = 0
  [../]

  [./fs]
    type = PorousFlowBrineCO2
    brine_fp = brine
    co2_fp = co2
    capillary_pressure = pc
  [../]
[]


[Modules]
  [./FluidProperties]
    [./co2]
      type = CO2FluidProperties
    [../]
    [./brine]
      type = BrineFluidProperties
    [../]
  [../]
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]


[Variables]
  [./pgas]
   initial_condition = 15e6
  [../]
  [./zi]
   initial_condition = 0
  [../]
  [./xnacl]
   initial_condition = 0.2500
  [../]
[]

[Functions]
 [./conditional_function]
  type = ParsedFunction
  vars = 'xnacl_sol'
  vals = 'xnacl'
  value = 'xnacl_sol <= 0.27'
 [../]
[]

[AuxVariables]

 [./r_m]
  initial_condition = 5.0
 [../]

 [./aperture]
  order = CONSTANT
  family = MONOMIAL
 [../]

 [./salt_change_rate]
 [../]

 [./xnacl_control]
  initial_condition = 0.27000000
 [../]

[]

[ICs]
  [./apertureread]
   type = ParameterFieldIC
   variable = aperture
   file_name = gap000_N18Ss_pcUO.txt
  [../]
[]


[Kernels]
  [./mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pgas
  [../]
  [./mass1]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = zi
  [../]
   [./mass2]
    type = PorousFlowMassTimeDerivative
    fluid_component = 2
    variable = xnacl
  [../]


  [./diff0]
    type = PorousFlowDispersiveFlux
    fluid_component = 0
    variable = pgas
    disp_trans = '0 0'
    disp_long = '0 0'
  [../]
  [./diff1]
    type = PorousFlowDispersiveFlux
    fluid_component = 1
    variable = zi
    disp_trans = '0 0'
    disp_long = '0 0'
  [../]
  [./diff2]
    type = PorousFlowDispersiveFlux
    fluid_component = 2
    variable = xnacl
    disp_trans = '0 0'
    disp_long = '0 0'
  [../]


  [./flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pgas
  [../]
  [./flux1]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = zi
  [../]
  [./flux2]
    type = PorousFlowAdvectiveFlux
    fluid_component = 2
    variable = xnacl
  [../]


  [./react2]
    type = PorousFlowHeatMassTransfer
    variable = xnacl
    v = xnacl_control
    transfer_coefficient = r_m
    save_in = salt_change_rate
  [../]

[]

[Controls]
 [./xnacl_threshold]
  type = ConditionalFunctionEnableControl
  conditional_function = conditional_function
  disable_objects = 'Kernel::react2'
  execute_on = 'initial timestep_begin'
 [../]
[]


[Materials]
  [./temperature]
    type = PorousFlowTemperature
  temperature = 323.15
  [../]
  [./temperature_nodal]
  type = PorousFlowTemperature
  temperature = 323.15
  at_nodes = true
  [../]
  [./brineco2]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = zi
    temperature_unit = Kelvin
    xnacl = xnacl
    capillary_pressure = pc
    fluid_state = fs
    at_nodes = true
  [../]
  [./brineco2_qp]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = zi
    temperature_unit = Kelvin
    xnacl = xnacl
    capillary_pressure = pc
    fluid_state = fs
  [../]
  [./porosity]
    type = PorousFlowPorosityConst
    porosity = aperture
  [../]
  [./permeability]
    type = FractureFlowPermeabilityConstFromVar
    aperture = aperture
  [../]
  [./diff]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '1.e-9 1.e-9 1.e-9 1.e-9 1.e-9 1.e-9'
    tortuosity = '1.0 1.0'
  [../]
  [./relperm_water]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 0
    s_res = 0.1
    sum_s_res = 0.1
  [../]
  [./relperm_gas]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 1
  [../]
[]

[BCs]
  [./co2_injection]
    type = PorousFlowSink
    boundary = top
    variable = zi
    flux_function = -2.2E-4
  [../]

  [./right_liquid_water]
    type = PorousFlowPiecewiseLinearSink
    boundary = bottom
    variable = pgas
    fluid_phase = 0
    pt_vals = '0 1E9'
    multipliers = '0 1E9'
    PT_shift = 15E6
    mass_fraction_component = 0
    use_mobility = true
    use_relperm = true
    flux_function = 1
  [../]
  [./right_gas_water]
    type = PorousFlowPiecewiseLinearSink
    boundary = bottom
    variable = pgas
    fluid_phase = 1
    pt_vals = '0 1E9'
    multipliers = '0 1E9'
    PT_shift = 15E6
    mass_fraction_component = 0
    use_mobility = true
    use_relperm = true
    flux_function = 1
  [../]

  [./right_liquid_co2]
    type = PorousFlowPiecewiseLinearSink
    boundary = bottom
    variable = zi
    fluid_phase = 0
    pt_vals = '0 1E9'
    multipliers = '0 1E9'
    PT_shift = 15E6
    mass_fraction_component = 1
    use_mobility = true
    use_relperm = true
    flux_function = 1
  [../]
  [./right_gas_co2]
    type = PorousFlowPiecewiseLinearSink
    boundary = bottom
    variable = zi
    fluid_phase = 1
    pt_vals = '0 1E9'
    multipliers = '0 1E9'
    PT_shift = 15E6
    mass_fraction_component = 1
    use_mobility = true
    use_relperm = true
    flux_function = 1
  [../]
 [./right_liquid_nacl]
 type = PorousFlowPiecewiseLinearSink
 boundary = bottom
 variable = xnacl
 fluid_phase = 0
 pt_vals = '0 1E9'
 multipliers = '0 1E9'
 PT_shift = 15E6
 mass_fraction_component = 2
 use_mobility = true
 use_relperm = true
 flux_function = 1
 [../]
[]

[Preconditioning]
  active = 'preferred'
  [./basic]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = 'gmres asm lu NONZERO 2'
  [../]
  [./preferred]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = 'lu mumps'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 1e5
  dtmax = 1e4
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.5
    optimal_iterations = 6
  [../]
  l_max_its = 1000
  l_tol = 1e-10
  nl_max_its = 100
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-8
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  csv = true
  [./out]
  type = Exodus
  execute_on = 'initial timestep_end'
  [../]
[]
