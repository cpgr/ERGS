#note: need to retrieve creep_stress and stress as AD properties!!

[Mesh]
  [efm]
   type = FileMeshGenerator
   file = kamlotFlo2.msh #kamlotFlo2.msh  kamlotFlo3.msh #kamlotFlosquare.msh  #  
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
#    initial_condition = 101325        
    block = ' rocksalt fluid' 
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
  [effective_fluid_pressure]
    family = MONOMIAL
    order = CONSTANT
    block = ' rocksalt fluid' 
  []
    [darcy_vel_x]
    family = MONOMIAL
    order = CONSTANT
    block = ' rocksalt fluid' 
  []
  [darcy_vel_y]
    family = MONOMIAL
    order = CONSTANT
    block = ' rocksalt fluid' 
  []
  [darcy_vel_z]
    family = MONOMIAL
    order = CONSTANT
    block = ' rocksalt fluid' 
  []
  [swater]
    family = MONOMIAL
    order = CONSTANT
    block = ' rocksalt fluid' 
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
   [kelvin_creep_strain]
    order = CONSTANT
    family = MONOMIAL
  []
   [smin_ramp]
    order = FIRST
    family = LAGRANGE
  []
   [smax_ramp]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [effective_fluid_pressure]
    type = ParsedAux
    args = 'pwater'
    function = 'pwater'
    variable = effective_fluid_pressure 
    block = 'rocksalt fluid' 
  []
    [darcy_vel_x]
    type = PorousFlowDarcyVelocityComponent
    component = x
    variable = darcy_vel_x
    fluid_phase = 0                             # OPTIONAL for single-phase
    execute_on = TIMESTEP_END 
    block = ' rocksalt fluid' 
  []
  [darcy_vel_y]
    type = PorousFlowDarcyVelocityComponent
    component = y
    variable = darcy_vel_y
    fluid_phase = 0                             # OPTIONAL for single-phase
    execute_on = TIMESTEP_END
    block = ' rocksalt fluid' 
  []
    [darcy_vel_z]
    type = PorousFlowDarcyVelocityComponent
    component = z
    variable = darcy_vel_z
    fluid_phase = 0                             # OPTIONAL for single-phase
    execute_on = TIMESTEP_END
    block = ' rocksalt fluid' 
  []
  [swater]
    type = PorousFlowPropertyAux
    variable = swater
    property = saturation
    phase = 0
    execute_on = timestep_end
    block = ' rocksalt fluid' 
  []
  [stress_xx]
    type = ADRankTwoScalarAux
    variable = stress_xx
    rank_two_tensor = stress
    scalar_type = MinPrincipal
    point1 = '0 0 0'
    point2 = '0 0 1'
    execute_on = timestep_end
  []
  [stress_yy]
    type = ADRankTwoScalarAux
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
 [kelvin_creep_strain] 
    type = ADMaterialRealAux                   #AD
    variable = kelvin_creep_strain
    property = kelvin_creep_strain
    execute_on = 'initial timestep_end'
  []
    [smin_ramp]
    type = FunctionAux
    variable = smin_ramp
    function = smin_fn
  []
   [smax_ramp]
    type = FunctionAux
    variable = smax_ramp
    function = smax_fn
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
  [chamber_pressure]
    type = ConstantIC
    variable = pwater
    value = 1e6
    block = fluid
  []
  [fluid_pressure]
    type = ConstantIC
    variable = pwater
    value = 101325 
    block = rocksalt
  []
 []


[Kernels]
  [time_derivative_water]
    type = PorousFlowMassTimeDerivative
    variable = 'pwater'
  []
  [flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    use_displaced_mesh = false
    variable = 'pwater'
  []
  [grad_stress_x]                                     #
    type = ADStressDivergenceTensors                    #   
    variable = disp_x                                 #
    component = 0
    use_displaced_mesh = true   
  []
  #[poro_x_rocksalt]
  #  type = PorousFlowEffectiveStressCoupling
  #  variable = disp_x
  #  component = 0
  #  use_displaced_mesh = false
  #  biot_coefficient = '1'
  #  block = 'rocksalt'
  #[]
#  [poro_x_casing]
#    type = PorousFlowEffectiveStressCoupling
#    variable = disp_x
#    component = 0
#    use_displaced_mesh = false
#    biot_coefficient =  '0'                      
#    block = 'casing'
#  []
  #[poro_x_fluid]
  #  type = PorousFlowEffectiveStressCoupling
  #  variable = disp_x
  #  component = 0
  #  use_displaced_mesh = false
  #  biot_coefficient = '1'
  #  block = 'fluid'
  #[]
  [grad_stress_y]
   type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = true  
  []
  #[poro_y_rocksalt]
  #  type = PorousFlowEffectiveStressCoupling
  #  variable = disp_y
  #  component = 1
  #  use_displaced_mesh = false
  #  biot_coefficient = '1'
  #  block = 'rocksalt'
  #[]
#  [poro_y_casing]
#    type = PorousFlowEffectiveStressCoupling
#    variable = disp_y
#    component = 1
#    use_displaced_mesh = false
#    biot_coefficient = '0'              
#    block = 'casing '
#  []
  #[poro_y_fluid]
  #  type = PorousFlowEffectiveStressCoupling
  #  variable = disp_y
  #  component = 1
  #  use_displaced_mesh = false
  #  biot_coefficient = '1'
  #  block = 'fluid'
  #[]
  [grad_stress_z]
    type = ADStressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = true      
  []
  #[poro_z_rocksalt]
  #  type = PorousFlowEffectiveStressCoupling
  #  variable = disp_z
  #  component = 2
  #  use_displaced_mesh = false
  #  biot_coefficient = '1'
  #  block = 'rocksalt'
  #[]
#  [poro_z_casing]
#    type = PorousFlowEffectiveStressCoupling
#    variable = disp_z
#    component = 2
#    use_displaced_mesh = false
#    biot_coefficient = '0'                     
#    block = 'casing'
#  []
  #[poro_z_fluid]
  #  type = PorousFlowEffectiveStressCoupling
  #  variable = disp_z
  #  component = 2
  #  use_displaced_mesh = false
  #  biot_coefficient = '1'
  #  block = 'fluid'
  #[]
[]


[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = 293.15
    use_displaced_mesh = false
  []
  [saturation]
    type = PorousFlow1PhaseP
    porepressure = 'pwater'
    capillary_pressure = pc
    block = ' rocksalt fluid' 
  []
  [massfrac]
    type = PorousFlowMassFraction
    block = ' rocksalt fluid' 
  []
  [water_viscosity_density]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
    block = ' rocksalt fluid' 
  []
  [relperm_water]
    type = PorousFlowRelativePermeabilityCorey
    n = 0.0
    s_res = 0.1
    sum_s_res = 0.2
    phase = 0
    block = ' rocksalt fluid' 
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
    block = 'fluid'
  []
    [permeability_rocksalt]
     type = ADPFOrthoEM
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-21 
     alpha = " 0.024 0.024 0.024"                   
     eps0 =  "-3e-5 0 0  0 -3e-5 0  0 0 -3e-5" 
    block = rocksalt
  []
    [permeability_casing]
     type = ADPFOrthoEM
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     alpha = " 0.024 0.024 0.024"                   
     eps0 =  "-3e-5 0 0  0 -3e-5 0  0 0 -3e-5" 
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-23 
     fix_rad_xy = 0 # 0.785398 #1.5708
     fix_rad_yz = 0 # 0.785398 #1.5708   
     block = 'casing'
    []
    [permeability_fluid]
     type = ADPFOrthoEM
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     alpha = " 0.024 0.024 0.024"                   
     eps0 =  "-3e-5 0 0  0 -3e-5 0  0 0 -3e-5" 
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-7
     fix_rad_xy = 0 # 0.785398 #1.5708
     fix_rad_yz = 0 # 0.785398 #1.5708
    block = 'fluid'
  []


  [elasticity_tensor_rocksalt]
    type = ADComputeIsotropicElasticityTensor                 #
    youngs_modulus = 2.5e10                   #2.5e10
    poissons_ratio = 0.25
    block = rocksalt
 []
  [elasticity_tensor_casing]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 2.1e11                     #2.1e11
    poissons_ratio =  0.30 
    block = 'casing'
  []
  [elasticity_tensor_fluid]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e11                  # 1e-18 #2e10   #1e-4 #1e-4                    #1e-9
    poissons_ratio = 0.0                         #0.00
    block = 'fluid'
  []
  
  [radial_return_stress_r]
    type = ADComputeMultipleInelasticStress
    inelastic_models = "L2creep"
    combined_inelastic_strain_weights = '1.0'
    block = 'rocksalt casing fluid'
  [] 
  [L2creep]
    type = ADLubby2CreepBase
    etaM0 = 1.21e8         
    mvM   = -0.327          
    GK0   = 1.88e5         
    mk    = -0.254         
    etaK0 = 4.98e5         
    mvK   = -0.267         
    Stensor = "12 0 0  0 0  0  0 0 0"
    uniaxial_test = true
    smax_ramp = smax_ramp
    smax_function = false
    smax_fix = 12
    smin_function = false
    smin_fix = 0
    smin_ramp = smin_ramp

 #  internal_solve_full_iteration_history = true
 #  use_substep = true
 #   relative_tolerance = 1e-7 #1e-9  1e-7 #
 #   absolute_tolerance = 1e-5 #1e-8  1e-5 #
 #   max_inelastic_increment = 1e-14
    block = 'rocksalt casing fluid'
  []  
 # [non_ad_L2creep]
 #   type = MaterialADConverter
 #   ad_props_in = 'ad_creep_strain'
 #   reg_props_out = 'reg_creep_strain'
 # [] 

 #  [effective_fluid_pressure]
 #   type = PorousFlowEffectiveFluidPressure
 #   block = ' rocksalt fluid'     
 # []
 # [volumetric_strain]
 #   type = PorousFlowVolumetricStrain        
 # []
[]


 [Modules/TensorMechanics/Master]
   [all]
    strain = SMALL
    incremental = true
    add_variables = true
    generate_output = 'stress_yy creep_strain_xx creep_strain_yy creep_strain_zz elastic_strain_yy'  
    use_automatic_differentiation = true
  []
[]


[BCs]
  [symmetric_face_left_fix_u]
   type = DirichletBC
   variable = disp_x
    value = 0
    boundary = 'left'        
  []
  [symmetric_face_front_fix_u]
   type = DirichletBC
   variable = disp_y
    value = 0
    boundary = 'front'        
  []
    [bottom_face_fix_u]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'bottom'         
  []
    [externalBoundary_compressive_stress_x]    #NOTE: Compressive is +
    type = Pressure
    function = 12e6    #load_x #                 
    variable = disp_x
    use_displaced_mesh = false
    boundary = 'right'
  []
    [externalBoundary_compressive_stress_y] 
    type = Pressure
    function = 21e6   #load_y   #                  
    variable = disp_y
    use_displaced_mesh = false
    boundary = 'back'
  []
#    [externalBoundary_compressive_stress_z] 
#    type = Pressure
#    function = 8e6  #load_z  #                   
#    variable = disp_z
#    use_displaced_mesh = false
#    boundary = 'top'
#  []            # if DirichletBC given a value equal 0, it means fixed (pin). If no value, it means roller. ideally, for uniaxial compressive, all sides should be roller,except fixed bottom and top pull.
  [u_top_pull]
    type = Pressure
    variable = disp_y
    boundary = top
    factor = 1                     #-2.5e5
    function = my_flux1 #top_pull
  []

#  [symmetric_faces_no_flow_pressure] 
#    type = DirichletBC
#   variable = pwater
#    value = 0
#    use_displaced_mesh = false
#    boundary = 'left front '
#  []
  [external_boundaries_atm_pressure] 
    type = DirichletBC
    variable = pwater
    value = 101325
    use_displaced_mesh = false
    boundary = 'right back top bottom'
  []


#  [top_water_injection]
#    type = PorousFlowSink
#    boundary = 'fluidT'                    
#    variable = pwater
#    flux_function = top_flux  
#    use_displaced_mesh = false
#  []
  [water_injection]
    type = PorousFlowSink
    boundary = 'fluidFill1 fluidB'                      #fluidFill1 fluidFill2 fluidL fluidR fluidB  
    variable = pwater
    flux_function = my_flux1  
    use_displaced_mesh = false
  []
  [injection_compressive_stress_x]
    type = Pressure
    boundary = 'fluidFill1'  
    variable = disp_x
    use_displaced_mesh = false
    postprocessor = constrained_effective_fluid_pressure_at_wellbore
  []
  [injection_compressive_stress_y]
    type = Pressure
    boundary = 'fluidFill1'
    variable = disp_y
    use_displaced_mesh = false
    postprocessor = constrained_effective_fluid_pressure_at_wellbore
  []
  [injection_compressive_stress_z]
    type = Pressure
    boundary = 'fluidB'
    variable = disp_z
    use_displaced_mesh = false
    postprocessor = constrained_effective_fluid_pressure_at_wellbore
  []
[]


[Functions]
  [my_flux1]
  type = ParsedFunction 
#  value = 'if(t<=300, 0, if(t <= 1140, -5.5e-4, if(t<= 1740, 0, -5.5e-4)))'  #-2.78e-5        
# value = 'if(t <= 300, 0, if(t<=600, -2.78e-5, if(t <= 1140,-5.5e-4, if(t<= 2040, 0,-5.5e-4))))'  
# value = 'if(t <= 300, 0, if(t<=600, -2.78e-5, if(t <= 1440,-5.5e-4, if(t<= 2040, 0,-5.5e-4))))' 
# value = 'if(t <= 300, 0, if(t<=600, 0, if(t <= 1140,-5.5e-4, if(t<= 2040, 0,-5.5e-4))))' 
# value = 'if(t <= 300, 0, if(t<=600, -4.31e-4, if(t <= 1140,-8.63e-3, if(t<= 2040, 0,-8.63e-3))))' 
# value = 'if(t <= 300, 1.43e-19, if(t<=600, -2.78e-5, if(t <= 1140,-5.5e-4, if(t<= 2040, 0,-5.5e-4))))'
  value = 'if(t<=350, 0, if(t <= 1140,-5.004e-4, if(t<= 1740, 0,-5.004e-4)))' 
#  value = 'if(t<=375, 0 , if(t <= 1125,-3e-4, if(t<= 1740, 0,-3e-4)))' 
  []
  [pcw_flux1]
    type = PiecewiseLinear
    x = '300    600       1140    2040    3600'                   # '300   2400'# 
    y = ' 0  -8.688e-4  -0.017375  0   -0.017375'           #'-1.39e-7  -2.78e-6'#   -8.688e-4/-2.78e-5  -0.017375/-5.5e-4
  [] 
  [load_x]
    type = PiecewiseLinear
    x = '0 300' 
    y = '0 12e6' 
  [] 
  [load_y]
    type = PiecewiseLinear
    x = '0 300' 
    y = '0 21e6' 
  []
  [load_z]
    type = PiecewiseLinear
    x = '0 300' 
    y = '0 8e6' 
  []
  [constrain_effective_fluid_pressure]
    type = ParsedFunction
    vars = effective_fluid_pressure_at_wellbore
    vals = effective_fluid_pressure_at_wellbore
    value = 'max(effective_fluid_pressure_at_wellbore, 20E6)'  
  []
  [top_pull]
   type = PiecewiseLinear
   x = ' 0  35 ' #  should be x= '0 1 15' where loading = 0 to 1 and creep-phase=1 to 15.    
   y = ' 0  -8e6'   #-12e6' #  should be y= '0 -2.5e5 -2.5e5', where 0 to -2.5e5 indicates loading phase, -2.5e5 to -2.5e5 indicates creep phase where load is kept constant. 
  []
   [smax_fn]
    type = ParsedFunction 
   value = 'if(t <= 53, 20, 42.0)'
  []
  [smin_fn]
    type = ParsedFunction 
    value = 'if(t <= 53, 20, 20.0)'
  []
[]


[Postprocessors]
   [p_chamber]
    type = PointValue
    variable = pwater
    point = '0 0 0.05'
    execute_on = timestep_end
  []
   [p_sample]
    type = PointValue
    variable = pwater
    point = '0.03 0.03 0.05'
    execute_on = timestep_end
  []
  [effective_fluid_pressure_at_wellbore]
    type = PointValue
    variable = effective_fluid_pressure
    point = ' 0 0 0.05'
    execute_on = timestep_begin
    use_displaced_mesh = false
  []
  [constrained_effective_fluid_pressure_at_wellbore]
    type = FunctionValuePostprocessor
    function = constrain_effective_fluid_pressure
    execute_on = timestep_begin
  []
  [injection_area]
    type = AreaPostprocessor
    boundary = 'fluidFill1 fluidB'
    execute_on = 'INITIAL'
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
  end_time = 2.4e3 #3600 
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-1           
    optimal_iterations = 10
  []
   nl_abs_tol = 1e-11
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






