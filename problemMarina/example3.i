[Mesh]
    uniform_refine = 1
 [file]
 type = FileMeshGenerator
 file = aperture_permcubicgap004_out.e
 use_for_exodus_restart = true
  []
 []


[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]


[UserObjects]

 [./dictator]
  type = PorousFlowDictator
  porous_flow_vars = 'pgas zi xnacl'
  number_fluid_phases = 2
  number_fluid_components = 3
 [../]

   [./pc]
    type = PorousFlowCapillaryPressureVG
    alpha = 0.001803
    m = 0.69235
    sat_lr = 0.005
    log_extension =true
    pc_max = 1e5
  [../]
  # [./pc]
  #   type = PorousFlowCapillaryPressureBC
  #   lambda = 1
  #   pe = 500e3
  # [../]
# Dummy Pc to pass to PorousFlowFluidState material to make the derivs 0
 [./pcdummy]
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
    [./co2sw]
      type = CO2FluidProperties
    [../]
    [./co2]
      type = TabulatedBicubicFluidProperties
      fp = co2sw
      save_file = false
    [../]
    [./water]
      type = Water97FluidProperties
    [../]
    [./watertab]
      type = TabulatedBicubicFluidProperties
      fp = water
      save_file = false
    [../]
    [./brine]
      type = BrineFluidProperties
      water_fp = watertab
    [../]
  [../]
[]

[Functions]
  [./time_func]
    type = ParsedFunction
    vals = 'get_dtime_var'
    vars = 'dtime_var'
    value = 'dtime_var'
  [../]
[]

[AuxVariables]
 [./saturation_gas]
  order = CONSTANT
  family = MONOMIAL
 [../]
 [./saturation_water]
  order = CONSTANT
  family = MONOMIAL
  initial_condition = 1.
 [../]

 [./frac_co2_in_liquid]
  order = CONSTANT
  family = MONOMIAL
 [../]
 [./frac_nacl_in_liquid]
  order = CONSTANT
  family = MONOMIAL
 [../]
 [./frac_water_in_liquid]
  order = CONSTANT
  family = MONOMIAL
 [../]
 [./frac_water_in_gas]
  order = CONSTANT
  family = MONOMIAL
 [../]

 [./co2_mass_ingas_flux_out]
 [../]
 [./co2_mass_inliquid_flux_out]
 [../]
 [./co2_mass_flux_in]
 [../]
 [./brine_mass_inliquid_flux_out]
 [../]
 [./brine_mass_ingas_flux_out]
 [../]
 [./nacl_mass_inliquid_flux_out]
 [../]
 [./nacl_mass_ingas_flux_out]
 [../]

 [./pressure_gas]
 order = CONSTANT
 family = MONOMIAL
 [../]

 [./pressure_liquid]
 order = CONSTANT
 family = MONOMIAL
 [../]

 [./capillary_pressure]
 order = CONSTANT
 family = MONOMIAL
 [../]

 [./unit]
  initial_condition = 1.
 [../]

 [./r_m]
  order = CONSTANT
  family = MONOMIAL
  initial_condition = 0.00
 [../]

 [./aperture]
  order = CONSTANT
  family = MONOMIAL
  initial_from_file_var = aperture
  initial_from_file_timestep = "LATEST"
 [../]
 [./permeability]
  order = CONSTANT
  family = MONOMIAL
  initial_from_file_var = permeability
  initial_from_file_timestep = "LATEST"
 [../]
 [./aperture_old]
  order = CONSTANT
  family = MONOMIAL
 [../]

 [./coeff_c]
  order = CONSTANT
  family = MONOMIAL
  initial_condition = 0.0
 [../]

 [./perm_view]
  order = CONSTANT
  family = MONOMIAL
 [../]
 [./poro_view]
  order = CONSTANT
  family = MONOMIAL
 [../]

 [./xnacl_control]
  initial_condition = 0.367
 [../]

 [./salt_change_reaction]
 [../]
 [./salt_change_diffusive]
 [../]
 [./salt_change_advective]
 [../]

 [./dtime_var]
  initial_condition = 0.0000000000001
 [../]

 [./brine_density]
  family = MONOMIAL
  order = CONSTANT
 [../]

[]

[AuxKernels]
  [./get_coeff_c]
    type = ParsedAux
    variable = coeff_c
    function = if(saturation_water>0.001,aperture*saturation_water*1251*r_m,aperture*1251*r_m)
    args = 'aperture saturation_water r_m'
    execute_on = 'timestep_begin'
  [../]

  [./change_aperture]
    type = ParsedAux
    variable = aperture
    function = if(r_m=0.0,aperture_old,(aperture_old-aperture_old*saturation_water*0.5765*dtime_var*r_m*(xnacl-0.367)))
    args = 'aperture_old aperture r_m dtime_var saturation_water xnacl'
    execute_on = 'timestep_begin linear'
  [../]

  [./change_permeability]
    type = ParsedAux
    variable = permeability
    function =  aperture*aperture*aperture/12
    args = 'aperture'
    execute_on = 'timestep_begin'
  [../]

  [./density]
    type = PorousFlowPropertyAux
    property = density
    phase = 0
    variable = brine_density
  [../]

  [./get_aperture_old]
    type = ParsedAux
    variable = aperture_old
    function = aperture
    args = 'aperture'
    execute_on = 'initial timestep_end'
  [../]

  [./get_r_m]
    type = ParsedAux
    variable = r_m
    function = if(xnacl<=0.3670,0.0,10)
    args = 'xnacl'
    execute_on = 'timestep_begin timestep_end'
  [../]

  [./get_time]
    type = FunctionAux
    variable = dtime_var
    function = time_func
    execute_on = 'initial linear timestep_end'
  [../]

  [./perm_view]
    type = PorousFlowPropertyAux
    property = permeability
    variable = perm_view
    execute_on = timestep_end
  [../]
  [./poro_view]
    type = PorousFlowPropertyAux
    property = porosity
    variable = poro_view
    execute_on = 'timestep_end'
  [../]

  [./saturation_gas]
    type = PorousFlowPropertyAux
    variable = saturation_gas
    property = saturation
    phase = 1
    execute_on = timestep_end
  [../]

  [./saturation_water]
    type = PorousFlowPropertyAux
    variable = saturation_water
    property = saturation
    phase = 0
    execute_on = 'timestep_begin linear'
  [../]

 [./pressure_gas]
 type = PorousFlowPropertyAux
 variable = pressure_gas
 property = pressure
 phase = 1
 [../]

 [./pressure_liquid]
 type = PorousFlowPropertyAux
 variable = pressure_liquid
 property = pressure
 phase = 0
 [../]

 [./capillary_pressure]
 type = ParsedAux
 args = 'pressure_gas pressure_liquid'
 function = 'pressure_gas - pressure_liquid'
 variable = capillary_pressure
 [../]

  [./frac_co2_in_liquid]
    type = PorousFlowPropertyAux
    variable = frac_co2_in_liquid
    property = mass_fraction
    phase = 0
    fluid_component = 1
    execute_on = timestep_end
  [../]
  [./frac_nacl_in_liquid]
    type = PorousFlowPropertyAux
    variable = frac_nacl_in_liquid
    property = mass_fraction
    phase = 0
    fluid_component = 2
    execute_on = timestep_end
  [../]
  [./frac_water_in_liquid]
    type = PorousFlowPropertyAux
    variable = frac_water_in_liquid
    property = mass_fraction
    phase = 0
    fluid_component = 0
    execute_on = timestep_end
  [../]
  [./frac_water_in_gas]
    type = PorousFlowPropertyAux
    variable = frac_water_in_gas
    property = mass_fraction
    phase = 1
    fluid_component = 0
    execute_on = timestep_end
  [../]
[]

[Variables]
  [./pgas]
    initial_condition = 15e6
  [../]
  [./zi]
    initial_condition = 0
  [../]
  [./xnacl]
    initial_condition = 0.340
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
    save_in = salt_change_diffusive
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
    save_in = salt_change_advective
  [../]

  [./react2]
    type = PorousFlowHeatMassTransfer
    variable = xnacl
    v = xnacl_control
    transfer_coefficient = coeff_c
    save_in = salt_change_reaction
  [../]
[]


[Materials]
  [./temperature]
    type = PorousFlowTemperature
    temperature = 323.15
  [../]
  [./brineco2]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = zi
    temperature_unit = Kelvin
    xnacl = xnacl
    capillary_pressure = pcdummy
    fluid_state = fs
  [../]
  [./porosity]
    type = PorousFlowPorosityConst
    porosity = aperture
  [../]
  [./permeability]
    type = PorousFlowPermeabilityTensorFromVar
    perm = permeability
  [../]
  [./diff]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '1.e-10 1.e-10 1.e-10 1.e-10 1.e-10 1.e-10'
    tortuosity = '1.0 1.0'
  [../]
  [./relperm_water]
    type = PorousFlowRelativePermeabilityCorey
    n = 4
    phase = 0
    s_res = 0.4
    sum_s_res = 0.4
  [../]
  [./relperm_gas]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 1
    sum_s_res = 0.4
  [../]
[]


[BCs]
  [./co2_injection]
    type = PorousFlowSink
    boundary = top
    variable = zi
    flux_function = -1.3649362222E-05
    save_in = co2_mass_flux_in
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
    # flux_function = 1000
    save_in = brine_mass_inliquid_flux_out
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
    # flux_function = 1000
    save_in = brine_mass_ingas_flux_out
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
    # flux_function = 1000
    save_in = co2_mass_inliquid_flux_out
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
    # flux_function = 1000
    save_in = co2_mass_ingas_flux_out
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
    # flux_function = 1000
    save_in = nacl_mass_inliquid_flux_out
  [../]
  [./right_gas_nacl]
    type = PorousFlowPiecewiseLinearSink
    boundary = bottom
    variable = xnacl
    fluid_phase = 1
    pt_vals = '0 1E9'
    multipliers = '0 1E9'
    PT_shift = 15E6
    mass_fraction_component = 2
    use_mobility = true
    use_relperm = true
    # flux_function = 1000
    save_in = nacl_mass_ingas_flux_out
  [../]
[]

[Preconditioning]
  active = 'basic'
  [./basic]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
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
  end_time = 9.1e4
  dtmax = 100
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.05
    growth_factor = 1.5
    # optimal_iterations = 10
  [../]
  dtmin = 1e-16
  l_max_its = 1000
  l_tol = 1e-10
  nl_max_its = 25
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-7
  automatic_scaling = true
  compute_scaling_once = false
[]

[Postprocessors]

 [./get_dtime_var]
  type = TimestepSize
  execute_on = 'initial linear timestep_end'
 [../]

 [./ave_zi_bottom]
  type = SideAverageValue
  variable = zi
  boundary = bottom
 [../]

 [./ave_pp_top]
  type = SideAverageValue
  variable = pgas
  boundary = top
 [../]

 [./ave_pp_bottom]
  type = SideAverageValue
  variable = pgas
  boundary = bottom
 [../]

 [./ave_sgas_bottom]
  type = SideAverageValue
  variable = saturation_gas
  boundary = bottom
 [../]

 [./ave_xnacl_top]
  type = SideAverageValue
  variable = xnacl
  boundary = top
 [../]

 [./ave_xnacl_bottom]
  type = SideAverageValue
  variable = xnacl
  boundary = bottom
 [../]

 [./total_mass_co2]
  type = PorousFlowFluidMass
  fluid_component = 1
  phase = '0 1'
  execute_on = 'initial timestep_end'
 [../]

 [./gas_mass_co2]
  type = PorousFlowFluidMass
  fluid_component = 1
  phase = 1
  execute_on = 'initial timestep_end'
 [../]

 [./sg_pos]
  type = ElementIntegralVariablePostprocessor
  variable = saturation_gas
  execute_on = 'initial timestep_end'
 [../]

 [./sw_pos]
  type = ElementIntegralVariablePostprocessor
  variable = saturation_water
  execute_on = 'initial timestep_end'
 [../]

 [./total_mass_nacl]
  type = PorousFlowFluidMass
  fluid_component = 2
  phase = '0 1'
  execute_on = 'initial timestep_end'
 [../]

 [./liq_mass_nacl]
  type = PorousFlowFluidMass
  fluid_component = 2
  phase = 0
  execute_on = 'initial timestep_end'
 [../]

 [./total_mass_water]
  type = PorousFlowFluidMass
  fluid_component = 0
  phase = '0 1'
  execute_on = 'initial timestep_end'
 [../]

 [./liquid_mass_water]
  type = PorousFlowFluidMass
  fluid_component = 0
  phase = 0
  execute_on = 'initial timestep_end'
 [../]

 [./injection_area]
  type = SideIntegralVariablePostprocessor
  boundary = top
  variable = aperture
  execute_on = 'initial'
 [../]

 [./total_fracture_volume]
  type = ElementIntegralVariablePostprocessor
  variable = aperture
  execute_on = 'initial timestep_end'
 [../]

 [./production_area]
  type = SideIntegralVariablePostprocessor
  boundary = bottom
  variable = aperture
  execute_on = 'initial'
 [../]

 [./production_area_unit]
  type = SideIntegralVariablePostprocessor
  boundary = bottom
  variable = unit
  execute_on = 'initial'
 [../]

 [./water_outlet_flux_inliquid]
  type = NodalSum
  boundary = bottom
  variable = brine_mass_inliquid_flux_out
  execute_on = 'timestep_end'
 [../]

 [./water_outlet_flux_ingas]
  type = NodalSum
  boundary = bottom
  variable = brine_mass_ingas_flux_out
  execute_on = 'timestep_end'
 [../]

 [./co2_outlet_flux_ingas]
  type = NodalSum
  boundary = bottom
  variable = co2_mass_ingas_flux_out
  execute_on = 'timestep_end'
 [../]

 [./co2_outlet_flux_inliquid]
  type = NodalSum
  boundary = bottom
  variable = co2_mass_inliquid_flux_out
  execute_on = 'timestep_end'
 [../]

 [./co2_inlet_flux]
  type = NodalSum
  boundary = top
  variable = co2_mass_flux_in
  execute_on = 'timestep_end'
 [../]

 [./nacl_outlet_flux_inliquid]
  type = NodalSum
  boundary = bottom
  variable = nacl_mass_inliquid_flux_out
  execute_on = 'timestep_end'
 [../]

 [./nacl_outlet_flux_ingas]
  type = NodalSum
  boundary = bottom
  variable = nacl_mass_ingas_flux_out
  execute_on = 'timestep_end'
 [../]

 [./brine_density_bottom]
  type = SideAverageValue
  variable = brine_density
  boundary = bottom
 [../]

[]

[Outputs]
  print_linear_residuals = false
  perf_graph = true
  csv = true
  checkpoint = true
  [./out]
    type = Exodus
    execute_on = 'initial timestep_end'
  [../]
[]

[Debug]
  show_var_residual = 'pgas zi xnacl'
[]
