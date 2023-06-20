[Mesh]
  [efm]
   type = FileMeshGenerator
   file = kamlot_arc.msh 
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
   bulk_modulus = 1e9
   density0 = 1 
   viscosity = 1e-3
   thermal_expansion = 0
  []
 []


[Variables]
  [pwater]
    initial_condition = 1e5       
  []
  [disp_x]
    scaling = 1E-7
  []
  [disp_y]
    scaling = 1E-7
  []
  [disp_z]
    scaling = 1E-7
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
#  [time_derivative]
#    type = PorousFlowMassTimeDerivative
#    variable = pwater
#  []
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
    component = 2
    use_displaced_mesh = false
  []
  [poro_z]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    component = 2
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
    [porosity_rocksalt]
    type = PorousFlowPorosityConst
    porosity = 0.01
    block = rocksalt
  []
    [porosity_casing]
    type = PorousFlowPorosityConst
    porosity = 0.0
    block = 'casing'
  []
    [permeability_rocksalt]
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
    block = rocksalt
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
    block = 'casing'
  []
  [relperm_water]
    type = PorousFlowRelativePermeabilityCorey
    n = 0.0
    s_res = 0.1
    sum_s_res = 0.2
    phase = 0
  []

  [elasticity_tensor_rocksalt]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.5e10
    poissons_ratio = 0.25
    block = rocksalt
  []
  [elasticity_tensor_casing]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.1e11
    poissons_ratio = 0.30
    block = 'casing'
  []
  [strain_rocksalt]
    type = ComputeSmallStrain
    eigenstrain_names = 'initial_stress'
    block = 'rocksalt'
  []
  [strain_casing]
    type = ComputeSmallStrain
    eigenstrain_names = 'initial_stress'
    block = 'casing'
  []
  [initial_strain_rocksalt]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = '0 0 0   0 0 0   0 0 0'
    eigenstrain_name = initial_stress
    block = 'rocksalt'
  []
  [initial_strain_casing]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = '0 0 0   0 0 0   0 0 0'
    eigenstrain_name = initial_stress
    block = 'casing'
  []
  [stress_rocksalt]
    type = ComputeLinearElasticStress
    block = 'rocksalt'
  []
  [stress_casing]
    type = ComputeLinearElasticStress
    block = 'casing'
  []

   [effective_fluid_pressure_rocksalt]
    type = PorousFlowEffectiveFluidPressure
    block = 'rocksalt'
  []
   [effective_fluid_pressure_casing]
    type = PorousFlowEffectiveFluidPressure
    block = 'casing'
  []
  [volumetric_strain_rocksalt]
    type = PorousFlowVolumetricStrain
    block = 'rocksalt'
  []
  [volumetric_strain_casing]
    type = PorousFlowVolumetricStrain
    block = 'casing'
  []
[]


[BCs]
  [fix_u_left_x]
   type = DirichletBC
   variable = disp_x
    value = 0
    boundary = 'left'        
  []
  [fix_u_front_y]
   type = DirichletBC
   variable = disp_y
    value = 0
    boundary = 'right'        
  []
    [fix_u_bottom_z]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'bottom'         
  []
    [externalBoundary_compressive_stress_x]    #NOTE: Compressive is +
    type = Pressure
    function = 1e5
    variable = disp_x
    use_displaced_mesh = false
    boundary =  outer_boundary
  []
    [externalBoundary_compressive_stress_y] 
    type = Pressure
    function = 1e5
    variable = disp_y
    use_displaced_mesh = false
    boundary =  outer_boundary
  []
    [externalBoundary_compressive_stress_z] 
    type = Pressure
    function = 1e5
    variable = disp_z
    use_displaced_mesh = false
    boundary = top
  []
    
  [outerBoundary_Pressure] 
    type = DirichletBC
    variable = pwater
    value = 1e5
    use_displaced_mesh = false
    boundary =  'outer_boundary top bottom'
  []
  [cavity_compressive_stress_x]
    type = Pressure
    boundary = 'injection_area injection_area_down'
    variable = disp_x
    function = 1e6
    use_displaced_mesh = false
  []
  [cavity_compressive_stress_y]
    type = Pressure
    boundary = 'injection_area injection_area_down'
    variable = disp_y
    function = 1e6
    use_displaced_mesh = false
  []
  [cavity_Pressure]
    type = DirichletBC
    boundary = 'injection_area injection_area_down'
    variable = pwater
    value = 1e6
    use_displaced_mesh = false
  []
[]


[Functions]
  [my_flux]
  type = ParsedFunction 
  value = 'if(t <= 300,-1e-10,-2.78e-9)' #'if(t <= 300,-1.39e-7,-2.78e-6)'    
  []
[]


[Postprocessors]
   [pwater]
    type = PointValue
    variable = pwater
    point = '0 0 0'
    execute_on = timestep_end
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
  type = Steady
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








