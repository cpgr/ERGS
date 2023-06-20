# HM single-phase injection of water into an annular model (representing fluid injection
# from a borehole into salt media) containing orthotropic embedded fractures (from: Numerical 
# BenchMark 2: Zill et. al.(2021): Hydro-mechanical continuum modelling of fluid percolation through rock.)

# NOTE: Both strain and Vector MUST ALWAYS be in the same direction for
# permeability to ACTIVATE.    

# NOTE: For single rotation, ONLY one plane rotation is needed. 
# The other plane automatically rotates as well. Rotate about the Z-axis (Right-Hand-Rule) 

[Mesh]
  [efm]
   type = FileMeshGenerator
   file = kamlot.msh
  []
[]


[UserObjects]
 [dictator]
  type = PorousFlowDictator
  porous_flow_vars = 'pwater disp_x disp_y disp_z'
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

  [FluidProperties]
    [water]
      type = SimpleFluidProperties
      bulk_modulus = 2e9
      density0 = 900
      viscosity = 7e-3
      thermal_expansion = 0
    []
  []


[Variables]
  [pwater]
    initial_condition = 1e5       
  []
  [disp_x]
    scaling = 1E-5
  []
  [disp_y]
    scaling = 1E-5
  []
  [disp_z]
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
  [darcy_vel_z]
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
  [randm_rad_XY]
    family = MONOMIAL
    order = CONSTANT
  []
  [rotxy_aux]
    order = FIRST
    family = LAGRANGE
  []
    [rotyz_aux]
    order = FIRST
    family = LAGRANGE
  []
  [fixedxy_aux]
    order = FIRST
    family = LAGRANGE
  []
   [fixedyz_aux]
    order = FIRST
    family = LAGRANGE
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
    [darcy_vel_z]
    type = PorousFlowDarcyVelocityComponent
    component = z
    variable = darcy_vel_z
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
   [randm_rad_XY]
    type = MaterialRealAux
    variable = randm_rad_XY
    property = random_xy_rotation_angle_for_each_element_qp
    execute_on = timestep_end
  []
[]


[ICs]
   [rotxy_aux]
    type = RandomIC
    min = 0    #0
    max = 1.57 #3.1415926535
    legacy_generator = false
    variable = rotxy_aux
   []
    [rotyz_aux]
    type = RandomIC
    min = 0
    max = 1.57 #3.1415926535
    legacy_generator = false
    variable = rotyz_aux
   []
 []


[Kernels]
  [time_derivative]
    type = PorousFlowMassTimeDerivative
    variable = pwater
  []
  [flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    use_displaced_mesh = false
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
  [grad_stress_z]
    type = StressDivergenceTensors
    variable = disp_z
    component = 1
    use_displaced_mesh = false
  []
  [poro_z]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    component = 1
    use_displaced_mesh = false
  []
[]


[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = 293.15
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
    fp = water
    phase = 0
  []
    [porosity_sample]
    type = PorousFlowPorosityConst
    porosity = 0.01
    block = sample
  []
    [porosity_casing]
    type = PorousFlowPorosityConst
    porosity = 0.0
    block = casing
  []
    [porosity_fluid_domain]
    type = PorousFlowPorosityConst
    porosity = 1.0
    block = fluid
  []
    [permeability_sample]
     type = PFOrthoEM
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     alpha =  "1e-2 1e-2 1e-2"    
     eps0 = "1e-5 1e-5 1e-5" 
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-21 
     fix_rad_xy = 0 # 0.785398 #1.5708
     fix_rad_yz = 0 # 0.785398 #1.5708
    block = sample
  []
    [permeability_casing]
     type = PFOrthoEM
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     alpha =  "1e-2 1e-2 1e-2"    
     eps0 = "1e-5 1e-5 1e-5" 
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-23 
     fix_rad_xy = 0 # 0.785398 #1.5708
     fix_rad_yz = 0 # 0.785398 #1.5708
    block = casing
  []
    [permeability_fluid_domain]
     type = PFOrthoEM
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     alpha =  "1e-2 1e-2 1e-2"    
     eps0 = "1e-5 1e-5 1e-5" 
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-7 
     fix_rad_xy = 0 # 0.785398 #1.5708
     fix_rad_yz = 0 # 0.785398 #1.5708
    block = fluid
  []
  [relperm_water]
    type = PorousFlowRelativePermeabilityCorey
    n = 0.0
    s_res = 0.1
    sum_s_res = 0.2
    phase = 0
  []

  [elasticity_tensor_sample]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.5e10
    poissons_ratio = 0.25
    block = sample
  []
  [elasticity_tensor_casing]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.1e11
    poissons_ratio = 0.30
    block = casing
  []
  [elasticity_tensor_fluid_domain]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e3
    poissons_ratio = 0.0
    block = fluid
  []
  [strain]
    type = ComputeSmallStrain
    eigenstrain_names = 'initial_stress'
    block = 'sample casing fluid'
  []
  [initial_strain]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = '0 0 0   0 0 0   0 0 0'
    eigenstrain_name = initial_stress
    block = 'sample casing fluid'
  []
  [stress]
    type = ComputeLinearElasticStress
    block = 'sample casing fluid'
  []

  [effective_fluid_pressure_mat]
    type = PorousFlowEffectiveFluidPressure
    block = 'sample casing fluid'
  []
  [volumetric_strain]
    type = PorousFlowVolumetricStrain
    block = 'sample casing fluid'
  []
[]


[BCs]
  [fix_u_left_x]
   type = DirichletBC
   variable = disp_x
    value = 0
    boundary = 'left'        
  []
  [fix_u_right_y]
   type = DirichletBC
   variable = disp_y
    value = 0
    boundary = 'front'        
  []
    [fix_u_bottom_z]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'bottom'         
  []
    [externalBoundary_compressive_stress_x]    #NOTE: Compressive is +
    type = Pressure
    boundary = right
    function = 12e6
    variable = disp_x
    use_displaced_mesh = false
  []
    [externalBoundary_compressive_stress_y] 
    type = Pressure
    boundary = back
    function = 21e6
    variable = disp_y
    use_displaced_mesh = false
  []
    [externalBoundary_compressive_stress_z] 
    type = Pressure
    boundary = top
    function = 8e6
    variable = disp_z
    use_displaced_mesh = false
  []
  [flux]
    type = PorousFlowSink
    boundary = 'injection_area1'
    variable = pwater
    use_mobility = false
    use_relperm = false
    flux_function = 2.7e-6
  []



 #   [outerBoundary_Pressure]
 #   type = DirichletBC
 #   boundary = outer_boundary
 #  variable = pwater
 #   value = 1e5
 #   use_displaced_mesh = false
 # []
 # [cavity_compressive_stress_x]
 #   type = Pressure
 #   boundary = injection_area
 #   variable = disp_x
 #   function = 1e6
 #   use_displaced_mesh = false
 # []
 # [cavity_compressive_stress_y]
 #   type = Pressure
 #   boundary = injection_area
 #   variable = disp_y
 #   function = 1e6
 #   use_displaced_mesh = false
 # []
 # [cavity_Pressure]
 #   type = DirichletBC
 #   boundary = injection_area
 #   variable = pwater
 #   value = 1e6
 #   use_displaced_mesh = false
 # []
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








