# HM single-phase injection of water into an annular model (representing fluid injection
# from a borehole into salt media) containing orthotropics embedded fractures (from: Numerical 
# BenchMark 2: Zill et. al.(2021): Hydro-mechanical continuum modelling of fluid percolation through rock.)
    

[Mesh]
  [efm]
   type = FileMeshGenerator
   file = annular.inp
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  gravity = '0 0 0'
  biot_coefficient = 1.0
  PorousFlowDictator = dictator
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

[Variables]
  [pwater]
    initial_condition = 1.01325e5       
  []
  [disp_x]
    scaling = 1E-5
  []
  [disp_y]
    scaling = 1E-5
  []
[]

[AuxVariables]
    [darcy_vel_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [darcy_vel_y]
    family = MONOMIAL
    order = CONSTANT
  []
  [swater]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_xx]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_yy]
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
    [permeability]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
    [darcy_vel_x]
    type = PorousFlowDarcyVelocityComponent
    component = x
    variable = darcy_vel_x
    fluid_phase = 0                             # OPTIONAL for single-phase
    execute_on = TIMESTEP_END
  []
  [darcy_vel_y]
    type = PorousFlowDarcyVelocityComponent
    component = y
    variable = darcy_vel_y
    fluid_phase = 0                             # OPTIONAL for single-phase
    execute_on = TIMESTEP_END
  []
  [swater]
    type = PorousFlowPropertyAux
    variable = swater
    property = saturation
    phase = 0
    execute_on = timestep_end
  []
  [stress_xx]
    type = RankTwoScalarAux
    variable = stress_xx
    rank_two_tensor = stress
    scalar_type = MinPrincipal
    point1 = '0 0 0'
    point2 = '0 0 1'
    execute_on = timestep_end
  []
  [stress_yy]
    type = RankTwoScalarAux
    variable = stress_yy
    rank_two_tensor = stress
    scalar_type = MidPrincipal
    point1 = '0 0 0'
    point2 = '0 0 1'
    execute_on = timestep_end
  []
  [porosity]
    type = PorousFlowPropertyAux
    variable = porosity
    property = porosity
    execute_on = timestep_end
  []
    [permeability]
    type = PorousFlowPropertyAux
    variable = permeability
    property = permeability
    execute_on = timestep_end
  []
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
    coupling_type = HydroMechanical
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
    component = 0
    use_displaced_mesh = false
  []
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_x
    component = 0
    use_displaced_mesh = false
  []
  [grad_stress_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = false
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_y
    component = 1
    use_displaced_mesh = false
  []
[]


#[Modules]
#  [FluidProperties]
#    [true_water]
#      type = Water97FluidProperties
#    []
#    [tabulated_water]
#      type = TabulatedFluidProperties
#      fp = true_water
#      temperature_min = 275
#      pressure_max = 1E8
#      interpolated_properties = 'density viscosity enthalpy internal_energy'
#      fluid_property_file = water97_tabulated_11.csv
#    []
#   []
#[]

[Modules]
  [FluidProperties]
    [water]
      type = SimpleFluidProperties
      bulk_modulus = 1e9 #2.27e14
      density0 = 1 #1e3 
      viscosity = 1e-3
      thermal_expansion = 0
    []
  []
[]


[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = 313.15
    use_displaced_mesh = false
  []
  [saturation]
    type = PorousFlow1PhaseP
    porepressure = pwater
    capillary_pressure = pc
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [water_viscosity_density]
    type = PorousFlowSingleComponentFluid
  #  fp = tabulated_water     #viscosity and density obtained from EOS
    fp = water
    phase = 0
  []
    [porosity_mat]
    type = PorousFlowPorosity
    fluid = true
    mechanical = true
    thermal = true
    porosity_zero = 0.01
    reference_temperature = 330
    reference_porepressure = 20E6
    thermal_expansion_coeff = 15E-6 # volumetric
    solid_bulk = 8E9 # unimportant since biot = 1
  []
    [permeability]
     type = PorousFlowOrthotropicEmbeddedFracturePermeability
     alpha =  "1e-2 1e-2 1e-2"    
     eps0 = "1e-5 1e-5 1e-5"   
     km = 1e-18                   # 1e-18
     rad_xy = 0                 # 0.785398
     rad_yz = 0                  # 0.785398
     N = "1 0 0  0 1 0  0 0 1"
  []
  [relperm_water]
    type = PorousFlowRelativePermeabilityCorey
    n = 4.0
    s_res = 0.1
    sum_s_res = 0.2
    phase = 0
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
    initial_stress = '0 0 0   0 0 0   0 0 0'
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


[BCs]
  [fix_u_x]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'left'        
  []
  [fix_u_y]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'bottom'         
  []
    [outerBoundary_compressive_stress_x]
    type = Pressure
    boundary = outer_boundary
    factor = 1e5
    variable = disp_x
    use_displaced_mesh = false
  []
    [outerBoundary_compressive_stress_y]
    type = Pressure
    boundary = outer_boundary
    factor = 1e5
    variable = disp_y
    use_displaced_mesh = false
  []
    [outerBoundary_Pressure]
    type = DirichletBC
    boundary = outer_boundary
    variable = pwater
    value = 1e5
    use_displaced_mesh = false
  []
  [cavity_compressive_stress_x]
    type = Pressure
    boundary = injection_area
    variable = disp_x
    factor = 1e6
    use_displaced_mesh = false
  []
  [cavity_compressive_stress_y]
    type = Pressure
    boundary = injection_area
    variable = disp_y
    factor = 1e6
    use_displaced_mesh = false
  []
    [cavity_Pressure]
    type = DirichletBC
    boundary = injection_area
    variable = pwater
    value = 1e6
    use_displaced_mesh = false
  []
  [constant_water_injection]
    type = PorousFlowSink
    boundary = injection_area
    variable = pwater
    fluid_phase = 0
    flux_function = -1e6
    use_displaced_mesh = false
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
  end_time = 1E3 #1E6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1E3  #1E4
    growth_factor = 1.2
    optimal_iterations = 10
  []
  nl_abs_tol = 1E-12
[]

[Outputs]
  exodus = true
  [csv]
  type = CSV
  execute_on = 'initial timestep_end'
  []
[]









