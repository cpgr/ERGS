# HM single-phase injection of water into a simple square model 
# containing planar embedded fractures  (Benchmark 1 -steady-state)
# NOTE: Both strain and Vector MUST ALWAYS be in the same direction for
# permeability to ACTIVATE.

[Mesh]
  [efmcube]
  type = FileMeshGenerator
  file = square2.inp
 # type = GeneratedMeshGenerator
 # dim = 2
 # nx = 10
 # ny = 10
 # xmin = 0
 #  xmax = 1
 #  ymin = 0
 #  ymax = 1
 []
[]


[GlobalParams]
 displacements = 'disp_x disp_y'
  gravity = '0 0 0'
  biot_coefficient = 0.0
  PorousFlowDictator = dictator
[]


[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pwater disp_x disp_y '
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
    initial_condition = 1e5  #1.01325e5 #
  []
  [disp_x]
  initial_condition = 0
    scaling = 1E-5
  []
  [disp_y]
  initial_condition = 0
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
  [stress_xx]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_yy]
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
    [strain_normal]
    family = MONOMIAL
    order = CONSTANT
  []
  [./permeability_xx]
  order = CONSTANT
  family = MONOMIAL
  [../]   
  [./permeability_yy]
  order = CONSTANT
  family = MONOMIAL
  [../]   
  [./permeability_zz]
  order = CONSTANT
  family = MONOMIAL
  [../]    
  [./permeability_xy]
  order = CONSTANT
  family = MONOMIAL
  [../]   
  [./permeability_xz]
  order = CONSTANT
  family = MONOMIAL
  [../]   
  [./permeability_yz]
  order = CONSTANT
  family = MONOMIAL
  [../]
    [density]
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
  [strain_normal]
    type = MaterialRealAux
    variable = strain_normal
    property = fracture_normal_strain_qp
    execute_on = timestep_end
  []
  [./permeability_xx]
  type = EFPComponents
  variable = permeability_xx
  index_i = 0
  index_j = 0
  [../]
  [./permeability_yy]
  type = EFPComponents
  variable = permeability_yy
  index_i = 1
  index_j = 1
  [../]
  [./permeability_zz]
  type = EFPComponents
  variable = permeability_zz
  index_i = 2
  index_j = 2
  [../]
  [./permeability_xy]
  type = EFPComponents
  variable = permeability_xy
  index_i = 0
  index_j = 1
  [../]
  [./permeability_xz]
  type = EFPComponents
  variable = permeability_xz
  index_i = 0
  index_j = 2
  [../]
  [./permeability_yz]
  type = EFPComponents
  variable = permeability_yz
  index_i = 1
  index_j = 2
  [../]
    [density]
    type = PorousFlowPropertyAux
    variable = density
    property = density
    execute_on = timestep_end
  []
[]


[Kernels]
  [flux_water]
    type = PorousFlowAdvectiveFlux
    variable = pwater
    use_displaced_mesh = false
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
    temperature = 293.15 #313.15
    use_displaced_mesh = false
  []
  [saturation_calculator]
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
   [porosity_mat]
    type = PorousFlowPorosity
    fluid = true
    mechanical = true
    thermal = true
    porosity_zero = 0.0
    reference_temperature = 293.15
    reference_porepressure = 20E6
    thermal_expansion_coeff = 15E-6 # volumetric
    solid_bulk = 8E9 # unimportant since biot = 1
  []
    [permeability]
    type = PorousFlowEmbeddedFracturePermeabilityBase
    a =  0.03
    b0 = 0.00
    e0 = 1e-5
    km = 1e-20
    rad_xy = 0 #0.04 #0.785398 #1.5708
    rad_yz = 0 #0.04 #0.785398 #1.5708
    n = "1 0 0"
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
    youngs_modulus = 1E9
    poissons_ratio = 0.3
  []
  [strain]
    type = ComputeSmallStrain
    eigenstrain_names = 'initial_stress'
  []
  [initial_strain]
    type = ComputeEigenstrainFromInitialStress
   initial_stress = '0 0 0  0 0 0  0 0 0'                      
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
  [u_fix_left_x]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary ='inlet'       
  []
  [u_fix_top_y]
    type = DirichletBC
    variable = disp_y
    value = 1e-4
    boundary = 'top'         
  []
  [u_fix_bottom_y]
    type = DirichletBC
   variable = disp_y
   value = 0
    boundary = 'bottom'      
  []
#  [u_fix_right_x]
#    type = DirichletBC
#    variable = disp_x
#    boundary = 'outlet' 
#    value = 0      
#  []
   [pressure_right_x]
    type = DirichletBC
    boundary = 'outlet'              
    variable = pwater
    value = 1e5
    use_displaced_mesh = false 
  []
  [flux_left]
   type = PorousFlowSink
    boundary = inlet
    variable = pwater
    fluid_phase = 0
    flux_function = -1e-10
    use_displaced_mesh = false
  []
[]

[Postprocessors]
   [fracture_normal_strain]
    type = PointValue
    variable = strain_normal
    point = '0 0 0'
    execute_on = timestep_end
  []
   [pwater_pressure]
    type = PointValue
    variable = pwater
    point = '0 0 0'
    execute_on = timestep_end
  [] 
   [permeability]
    type = PointValue
    variable = permeability
    point = '0 0 0'
    execute_on = timestep_end
  []
    [density]
    type = PointValue
    variable = density
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
    dt = 1E4  # 1E3 1E4
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

#[Debug]
#  show_material_props = true
#[]






