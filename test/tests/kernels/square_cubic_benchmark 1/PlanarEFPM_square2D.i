# HM single-phase injection of water into a simple cubic model 
# containing planar embedded fractures 

[Mesh]
  [efmcube]
  #  type = FileMeshGenerator
  #  file = square.inp
  type = GeneratedMeshGenerator
  dim = 2
  nx = 10
  ny = 10
  xmin = 0
 xmax = 1
  ymin = 0
  ymax = 1
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
    initial_condition = 1.01325e5       #atm
  []
  [disp_x]
    scaling = 1E-5
  []
  [disp_y]
    scaling = 1E-5
  []
[]


[AuxVariables]
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
   # fluid_component = 0
  # coupling_type = HydroMechanical
    variable = pwater
    use_displaced_mesh = false
  []
  [vol_strain_rate_water]
    type = PorousFlowMassVolumetricExpansion
   # fluid_component = 0
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
  [temperature]
    type = PorousFlowTemperature
    temperature = 313.15
    use_displaced_mesh = false
  []
    [biot_modulus]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 1
    solid_bulk_compliance = 4.545e-7
    fluid_bulk_modulus = 2.27e14
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
    fp = tabulated_water                                #viscosity and density obtained from EOS
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
    [permeability]
    type = PorousFlowEmbeddedFracturePermeability
    a =  10
    e0 = 10
    km = 1e-19
    rad_xy = 0.785398
    rad_yz = 0.785398
  #  jf = 1
    n = "0 0 1"

   # type = PorousFlowPermeabilityKozenyCarman
   # poroperm_function = kozeny_carman_phi0
   # phi0 = 0.1
   # n = 2
   # m = 2
   # k0 = 1E-12
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


[BCs]
  [u_fix_left_x]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'left'       
  []
  [u_fix_top_y]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'top'         
  []
  [u_fix_bottom_y]
    type = DirichletBC
   variable = disp_y
   value = 0
    boundary = 'bottom'      
  []
 # [u_fix_right_x]
 #   type = DirichletBC
 #   variable = disp_x
 #   value = 0
 #   boundary = 'right'       
 # []
   [pressure_right_x]
    type = DirichletBC
    boundary = 'right'              
    variable = pwater
    value = 1e5
    use_displaced_mesh = false 
  []
  [flux_left]
    type = PorousFlowSink
    boundary = 'left'
    variable = pwater
    fluid_phase = 0
    flux_function = 1e-10    
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









