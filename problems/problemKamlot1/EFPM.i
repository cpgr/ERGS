[Mesh]
  [annular]
    type = AnnularMeshGenerator
    nr = 10
    rmin = 1.0
    rmax = 10
    growth_r = 1.4
    nt = 4
    dmin = 0
    dmax = 90
  []

  [./extrude]
    type = MeshExtruderGenerator
    input = 'annular'
    num_layers = 3
    extrusion_vector = '0 0 12' 
   # existing_subdomains = '3'
   # layers = '3'
   # new_ids = 2
   # show_info = false
    bottom_sideset = 'bottom'
    top_sideset = 'top'
  []
[]

[Variables]
  [pwater]
    initial_condition = 1.01325e5       #atm
  []
  [T]
    initial_condition = 330
    scaling = 1E-5
  []
  [disp_x]
    scaling = 1E-5
  []
  [disp_y]
    scaling = 1E-5
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pwater disp_x disp_y'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    alpha = 1E-6
    m = 0.6
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  gravity = '0 0 0'
  biot_coefficient = 1.0
  PorousFlowDictator = dictator
[]

[Kernels]
  [mass_water_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pwater
  []
  [flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    use_displaced_mesh = false
    variable = pwater
  []
  [vol_strain_rate_water]
    type = PorousFlowMassVolumetricExpansion
    fluid_component = 0
    variable = pwater
  []
  [grad_stress_x] 
    type = StressDivergenceTensors 
    variable = disp_x
    use_displaced_mesh = false
    component = 0
  []
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_x
    use_displaced_mesh = false
    component = 0
  []
  [grad_stress_y]
    type = StressDivergenceTensors
    variable = disp_y
    use_displaced_mesh = false
    component = 1
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_y
    use_displaced_mesh = false
    component = 1
  []
[]

[AuxVariables]
  [disp_z]
  []
  [effective_fluid_pressure]
    family = MONOMIAL
    order = CONSTANT
  []
  [mass_frac_phase0_species0]
    initial_condition = 1 # all water in phase=0
  []
  [swater]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_rr]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_tt]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_zz]
    family = MONOMIAL
    order = CONSTANT
  []
  [porosity]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [effective_fluid_pressure]
    type = ParsedAux
    args = 'pwater swater'
    function = 'pwater * swater'
    variable = effective_fluid_pressure
  []
  [swater]
    type = PorousFlowPropertyAux
    variable = swater
    property = saturation
    phase = 0
    execute_on = timestep_end
  []
  [stress_rr]
    type = RankTwoScalarAux
    variable = stress_rr
    rank_two_tensor = stress
    scalar_type = RadialStress
    point1 = '0 0 0'
    point2 = '0 0 1'
    execute_on = timestep_end
  []
  [stress_tt]
    type = RankTwoScalarAux
    variable = stress_tt
    rank_two_tensor = stress
    scalar_type = HoopStress
    point1 = '0 0 0'
    point2 = '0 0 1'
    execute_on = timestep_end
  []
  [stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  []
  [porosity]
    type = PorousFlowPropertyAux
    variable = porosity
    property = porosity
    execute_on = timestep_end
  []
[]

[BCs]
  [roller_tmax]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = dmax        #sideset in max theta (rollers@ side/Y-dir.)
  []
  [roller_tmin]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = dmin         #sideset in min theta (rollers@ side/X-dir.)
  []
  [pinned_top_bottom_x]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'bottom'      #bottom-sideset (fixed displacement @bottom-x)
  []
  [pinned_top_bottom_y]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'bottom'        #bottom-sideset(fixed displacememnt @bottom-y)
  []
    [external_pressure_top_x]
    type = Pressure
    boundary = top              #topset (stress@top-x)
    variable = disp_x
    factor = 1
    function = -8e6
    use_displaced_mesh = false 
  []
    [external_pressure_top_y]
    type = Pressure
    boundary = top               #topset (stress@top-y)
    variable = disp_y
    factor = -8e6
  #  function = i think it will default to 1
    use_displaced_mesh = false
  []
    [external_pressure_tmax]
    type = Pressure
    boundary = dmax            #sideset in max theta (stress @ Side/Y dir.)
    variable = disp_x
    factor = -4E6
   # function = 
    use_displaced_mesh = false
  []
    [external_pressure_tmin]
    type = Pressure
    boundary = dmin            #sideset in min theta (stress @ Side/X dir.)
    variable = disp_y
  #  function = 
    factor = -15e6
    use_displaced_mesh = false
  []
  [cavity_pressure_x]
    type = Pressure
    boundary = injection_area   #injection-rate pf inside borehole (X dir.)
    variable = disp_x
    factor = 1
    postprocessor = constrained_effective_fluid_pressure_at_wellbore
    use_displaced_mesh = false
  []
  [cavity_pressure_y]
    type = Pressure 
    boundary = injection_area   
    variable = disp_y
    factor = 1                   #injection-rate pfactor inside borehole (Y dir.)
    postprocessor = constrained_effective_fluid_pressure_at_wellbore
    use_displaced_mesh = false
  []

  [constant_water_injection]
    type = PorousFlowSink
    boundary = injection_area
    variable = pwater
    fluid_phase = 0
    flux_function = -1.38e-7       #injection-rate value inside borehole (Y dir.)
    use_displaced_mesh = false
  []

[Modules]
  [FluidProperties]
    [true_water]
      type = Water97FluidProperties
    []
    [tabulated_water]
      type = TabulatedFluidProperties
      fp = true_water
      temperature_min = 275
      pressure_max = 1E8
      interpolated_properties = 'density viscosity enthalpy internal_energy'
      fluid_property_file = water97_tabulated_11.csv
    []
   []
[]

[Materials]
  [saturation_calculator]
    type = PorousFlow1PhasePP
    phase0_porepressure = pwater
    capillary_pressure = pc
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'mass_frac_phase0_species0'
  []
  [water_viscosity_density]
    type = PorousFlowSingleComponentFluid
    fp = tabulated_water                                #viscosity and density obtained from EOS
    phase = 0
  []
  [relperm_water]
    type = PorousFlowRelativePermeabilityCorey
    n = 4
    s_res = 0.1
    sum_s_res = 0.2
    phase = 0
  []
  [porosity_mat]
    type = PorousFlowPorosity
    fluid = true
    mechanical = true
    thermal = true
    porosity_zero = 0.1
    reference_temperature = 330
    reference_porepressure = 20E6
    thermal_expansion_coeff = 15E-6 # volumetric
    solid_bulk = 8E9 # unimportant since biot = 1
  []
  [permeability_injection_block]
    type = PorousFlowEmbeddedFracturePermeabilityJB
    block = injection_block
    a = 
    e0 = 
    km = 
    b0 = 
    rad_xy = 
    rad_yz = 
    jf = 
  []
  [permeability_caps]
    type = PorousFlowEmbeddedFracturePermeabilityJB
    block = caps
    a = 
    e0 = 
    km = 
    b0 = 
    rad_xy = 
    rad_yz = 
    jf = 
  []

  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 5E9
    poissons_ratio = 0.0
  []
  [strain]
    type = ComputeSmallStrain
    eigenstrain_names = 'initial_stress'
  []
  [initial_strain]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = '20E6 0 0  0 20E6 0  0 0 20E6'
    eigenstrain_name = initial_stress
  []
  [stress]
    type = ComputeLinearElasticStress
  []

  [effective_fluid_pressure_mat]
    type = PorousFlowEffectiveFluidPressure
  []
  [volumetric_strain]
    type = PorousFlowVolumetricStrain
  []
[]

[Postprocessors]
  [effective_fluid_pressure_at_wellbore]
    type = PointValue
    variable = effective_fluid_pressure
    point = '1 0 0'
    execute_on = timestep_begin
    use_displaced_mesh = false
  []
  [constrained_effective_fluid_pressure_at_wellbore]
    type = FunctionValuePostprocessor
    function = constrain_effective_fluid_pressure
    execute_on = timestep_begin
  []
[]

[Functions]
  [constrain_effective_fluid_pressure]
    type = ParsedFunction
    vars = effective_fluid_pressure_at_wellbore
    vals = effective_fluid_pressure_at_wellbore
    value = 'max(effective_fluid_pressure_at_wellbore, 20E6)'
  []
[]

[Preconditioning]
  active = basic
  [basic]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
  []
  [preferred_but_might_not_be_installed]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  end_time = 1E6 #1E6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1E5  #1E4
    growth_factor = 1.2
    optimal_iterations = 10
  []
  nl_abs_tol = 1E-7
[]

[Outputs]
  exodus = true
  [csv]
  type = CSV
  execute_on = 'initial timestep_end'
  []
[]
