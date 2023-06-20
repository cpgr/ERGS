[Mesh]                                       
  [efm]
   type = FileMeshGenerator
   file = kamlot3RZ2.msh   #kamlot3_2D.msh  
  []
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


[GlobalParams]
  displacements = 'disp_x disp_y'
  gravity = '0 0 0'
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
    initial_condition = 101325
    block = ' rocksalt fluid injection_area' # NOTE: pwater exist on casing too, but w/o any active kernel because it is block-restricted. See Kernel for info.
  []
  [disp_x]
    scaling = 1E-5
  []
  [disp_y]
    scaling = 1E-5
  []
[]

[AuxVariables]
  [disp_z]
  []
  [effective_fluid_pressure]
    family = MONOMIAL
    order = CONSTANT
    block = ' rocksalt fluid injection_area' 
  []
    [darcy_vel_x]
    family = MONOMIAL
    order = CONSTANT
    block = ' rocksalt fluid injection_area' 
  []
  [darcy_vel_y]
    family = MONOMIAL
    order = CONSTANT
    block = ' rocksalt fluid injection_area' 
  []
  [darcy_vel_z]
    family = MONOMIAL
    order = CONSTANT
    block = ' rocksalt fluid injection_area' 
  []
  [swater]
    family = MONOMIAL
    order = CONSTANT
    block = ' rocksalt fluid injection_area' 
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
  [effective_fluid_pressure]
    type = ParsedAux
    args = 'pwater '
    function = 'pwater'
    variable = effective_fluid_pressure
    block = ' rocksalt fluid injection_area' 
  []
    [darcy_vel_x]
    type = PorousFlowDarcyVelocityComponent
    component = x
    variable = darcy_vel_x
    fluid_phase = 0            
    block = ' rocksalt fluid injection_area'                  
    execute_on = TIMESTEP_END
  []
  [darcy_vel_y]
    type = PorousFlowDarcyVelocityComponent
    component = y
    variable = darcy_vel_y
    fluid_phase = 0          
    block = ' rocksalt fluid injection_area'                    
    execute_on = TIMESTEP_END
  []
    [darcy_vel_z]
    type = PorousFlowDarcyVelocityComponent
    component = z
    variable = darcy_vel_z
    fluid_phase = 0     
    block = ' rocksalt fluid injection_area'                         
    execute_on = TIMESTEP_END
  []
  [swater]
    type = PorousFlowPropertyAux
    variable = swater
    property = saturation
    phase = 0
    block = ' rocksalt fluid injection_area' 
    execute_on = timestep_end
  []
  [stress_xx]
    type = RankTwoScalarAux
    variable = stress_xx
    rank_two_tensor = stress
    scalar_type = RadialStress
    point1 = '0 0 0'
    point2 = '0 0 1'
    execute_on = timestep_end
  []
  [stress_yy]
    type = RankTwoScalarAux
    variable = stress_yy
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
    max = 1.57 
    legacy_generator = false
    variable = rotyz_aux
   []
 []


[Kernels]
  [time_derivative]
    type = PorousFlowMassTimeDerivative
    variable = pwater                     #NOTE: no need to block-restrict, it is automatic because of block-restriction above (i.e., once variables are block restricted, their kernels become inactive)
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
  [poro_x_rocksalt]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_x
    component = 0
    use_displaced_mesh = false
    biot_coefficient = '1'
    block = 'rocksalt'
  []
  [poro_x_fluid]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_x
    component = 0
    use_displaced_mesh = false
    biot_coefficient = '1'
    block = 'fluid injection_area'
  []
#  [poro_x_casing]
#    type = PorousFlowEffectiveStressCoupling
#    variable = disp_x
#    component = 0
#    use_displaced_mesh = false
#    biot_coefficient = '0'
#  []
  [grad_stress_y]
   type = StressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = false
  []
  [poro_y_rocksalt]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_y
    component = 1
    use_displaced_mesh = false
    biot_coefficient = '1'
    block = 'rocksalt'
  []
  [poro_y_fluid]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_y
    component = 0
    use_displaced_mesh = false
    biot_coefficient = '1'
    block = 'fluid injection_area'
  []
#  [poro_y_casing]
#    type = PorousFlowEffectiveStressCoupling
#    variable = disp_y
#    component = 1
#    use_displaced_mesh = false
#    biot_coefficient = '0'
#    block = 'casing'
#  []
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
    block = ' rocksalt fluid injection_area' 
  []
  [massfrac]
    type = PorousFlowMassFraction
    block = ' rocksalt fluid injection_area' 
  []
  [water_viscosity_density]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
    block = ' rocksalt fluid injection_area' 
  [] 
  [relp]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
    block = ' rocksalt fluid injection_area' 
  []
    [porosity_rocksalt]
    type = PorousFlowPorosityConst
    porosity = 0.01
    block = rocksalt
  []
    [porosity_casing]
    type = PorousFlowPorosityConst
    porosity =  0.0                     
    block = 'casing'
  []
    [porosity_fluid]
    type = PorousFlowPorosityConst
    porosity = 1.0
    block = 'fluid injection_area'
  []
    [permeability_rocksalt]
     type = PFOrthoEM
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-21 
     alpha = " 0.0455 0.0455 0.0455"                    
     eps0 =  "-1e-6 0 0  0 -1e-6 0  0 0 -1e-6"          
    block = rocksalt
  []
    [permeability_casing]
     type = PFOrthoEM
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     alpha = " 0.0455 0.0455 0.0455"                    
     eps0 =  "-1e-6 0 0  0 -1e-6 0  0 0 -1e-6"            
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-23 
     fix_rad_xy = 0 # 0.785398 #1.5708
     fix_rad_yz = 0 # 0.785398 #1.5708   
     block = 'casing'
    []
    [permeability_fluid]
     type = PFOrthoEM
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     alpha = " 0.0455 0.0455 0.0455"                    
     eps0 =  "-1e-6 0 0  0 -1e-6 0  0 0 -1e-6"       
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-7
     fix_rad_xy = 0 # 0.785398 #1.5708
     fix_rad_yz = 0 # 0.785398 #1.5708
    block = 'fluid injection_area'
  []
  
  [elasticity_tensor_rocksalt]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus =  2.5e10
    poissons_ratio = 0.25
    block = 'rocksalt'
  []
  [elasticity_tensor_casing]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.1e11
    poissons_ratio = 0.30
   block = 'casing'
  []
  [elasticity_tensor_fluid]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e7
    poissons_ratio = 0.0
    block = 'fluid injection_area'
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
  [effective_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    block = ' rocksalt fluid injection_area' 
  []
  [volumetric_strain]
    type = PorousFlowVolumetricStrain
  []
[]


[BCs]
  [roller_tmax]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'left'                  
  []
  [u_fix_bottom_y]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'bottom'      
  []
    [externalBoundary_compressive_stress_z] 
    type = Pressure
    function =  19e6                                
    variable = disp_y
    use_displaced_mesh = false
    boundary = ' top'
  []
    [externalBoundary_compressive_stress_r] 
    type = Pressure
    function =  4e6                                                   
    variable = disp_x
    use_displaced_mesh = false
    boundary = 'right'
  []

#  [symmetric_faces_no_flow_pressure] 
#    type = DirichletBC
#   variable = pwater
#    value = 0
#    use_displaced_mesh = false
#    boundary = 'left '
#  []
  [external_boundaries_atm_pressure] 
    type = DirichletBC
    variable = pwater
    value = 101325
    use_displaced_mesh = false
    boundary = 'right top bottom' 
  []

  [water_injection]
    type = PorousFlowSink
    boundary = ' injectionR injectionT injectionB'            
    variable = pwater
    flux_function = my_flux2  
    use_displaced_mesh = false
  []
  [injection_compressive_stress_x]
    type = Pressure
    boundary = 'injectionR'
    variable = disp_x
    use_displaced_mesh = false
    postprocessor = constrained_effective_fluid_pressure_at_wellbore
  []
  [injection_compressive_stress_y]
    type = Pressure
    boundary = 'injectionT injectionB'
    variable = disp_y
    use_displaced_mesh = false
    postprocessor = constrained_effective_fluid_pressure_at_wellbore
  []
[]


[Functions]
  [my_flux1]
  type = ParsedFunction 
  value = 'if(t <= 300, -1.39e-7,-2.78e-6)'   #'if(t <= 300,-1.39e-7,-2.78e-6)'        # 0 converge
  []
  [my_flux2]
  type = ParsedFunction 
  value = 'if(t <= 300,-5.5e-7, -2.78e-5)'   #'if(t <= 300, -5.5e-7, -2.78e-5)'
  []
  [pcw_flux1]
    type = PiecewiseLinear
    x = ' 300      2400' 
    y = '-1.39e-7  -2.78e-6'
  []
  [pcw_flux2]
    type = PiecewiseLinear
    x = '300   2400' 
    y = '-5.5e-7 -2.78e-5'     
  []
  [constrain_effective_fluid_pressure]
    type = ParsedFunction
    vars = effective_fluid_pressure_at_wellbore
    vals = effective_fluid_pressure_at_wellbore
    value = 'max(effective_fluid_pressure_at_wellbore, 20e6)'
  []
[]


[Postprocessors]
   [pwater]
    type = PointValue
    variable = pwater
    point = '0 0.05 0'            #'8e-3 0.05 0'
    execute_on = timestep_end
  []
  [effective_fluid_pressure_at_wellbore]
    type = PointValue
    variable = effective_fluid_pressure
    point = '0 0.05 0'            #'8e-3 0.05 0'
    execute_on = timestep_begin
    use_displaced_mesh = false
  []
  [constrained_effective_fluid_pressure_at_wellbore]
    type = FunctionValuePostprocessor
    function = constrain_effective_fluid_pressure
    execute_on = timestep_begin
  []
[]


[Preconditioning]
  active = preferred_but_might_not_be_installed                
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
  end_time = 2400 
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-1               
    optimal_iterations = 10
  []
   nl_abs_tol = 1e-8                                                                             #1e-11/10 converge for flux1/2
[]

[Outputs]
  exodus = true
  [csv]
  type = CSV
  execute_on = 'initial timestep_end'
  []
[]








