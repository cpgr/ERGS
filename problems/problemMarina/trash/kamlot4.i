[Mesh]
  [efm]
   type = FileMeshGenerator
   file = kamlotFlo2.msh #kamlotFlo2.msh  kamlotFlo3.msh #kamlotFlosquare.msh  #  
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  gravity = '0 0 0'
  PorousFlowDictator = dictator
  temperature_unit = Celsius
[]

[UserObjects]
 [dictator]
  type = PorousFlowDictator
  porous_flow_vars = 'pwater zi xnacl disp_x disp_y disp_z '
  number_fluid_phases = 2
  number_fluid_components = 3
 []
 [pc]
  type = PorousFlowCapillaryPressureVG
  alpha = 1E-6
  m = 0.6
  []
 [./pcdummy]
    type = PorousFlowCapillaryPressureConst
    pc = 0
 [../]
 [fs]
    type = PorousFlowBrineCO2
    brine_fp = brine
    co2_fp = co2
    capillary_pressure = pc
  []
[]
  

[FluidProperties]
    [co2sw]
      type = CO2FluidProperties
    []
    [co2]
      type = TabulatedBicubicFluidProperties
      fp = co2sw
      save_file = false
    []
    [water]
      type = Water97FluidProperties
    []
    [watertab]
      type = TabulatedBicubicFluidProperties
      fp = water
      save_file = false
    []
    [brine]
      type = BrineFluidProperties
      water_fp = watertab
    []
  []


[Variables]
  [pwater]
#    initial_condition = 101325        
    block = ' rocksalt fluid' 
  []
  [zi]
    initial_condition = 0        
    block = 'rocksalt fluid' 
  []
  [xnacl]
    initial_condition = 0.340        
    block = 'rocksalt fluid' 
  []
  [disp_x]
    scaling = 1e-0 #1E-5  
  []
  [disp_y]
    scaling = 1e-0 #1E-5  
  []
  [disp_z]
    scaling = 1e-0 #1E-5  
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
    block = ' rocksalt fluid'
  []
  [randm_rad_XY]
    family = MONOMIAL
    order = CONSTANT
    block = ' rocksalt fluid'
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
  
 
 [pgas]
    order = CONSTANT
    family = MONOMIAL       
    block = 'rocksalt fluid' 
  []
 [temperature]
    initial_condition = 100
  []
 [unit]
     initial_condition = 1.0
  []

 [r_m]
  order = CONSTANT
  family = MONOMIAL
  initial_condition = 0.00       
    block = 'rocksalt fluid' 
 []
 [aperture]
  order = CONSTANT
  family = MONOMIAL           
    block = 'rocksalt fluid' 
#  initial_from_file_var = aperture
#  initial_from_file_timestep = "LATEST"
 []
 [aperture_old]
  order = CONSTANT
  family = MONOMIAL           
    block = 'rocksalt fluid' 
 []
 [coeff_c]
  order = CONSTANT
  family = MONOMIAL
  initial_condition = 0.0     
    block = 'rocksalt fluid' 
 []
 [dtime_var]
  initial_condition = 0.0000000000001      
    block = 'rocksalt fluid' 
 []

 [xnacl_control]
  initial_condition = 0.367
 []
  [salt_change_reaction]
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
    block = ' rocksalt fluid'
    execute_on = timestep_end
  []
   [randm_rad_XY]
    type = MaterialRealAux
    variable = randm_rad_XY
    property = random_xy_rotation_angle_for_each_element_qp
    block = ' rocksalt fluid'
    execute_on = timestep_end
  []

  [pgas]
    type = PorousFlowPropertyAux
    property = pressure
    phase = 1
    variable = pgas       
    block = 'rocksalt fluid' 
    execute_on = 'TIMESTEP_END'
  []

  [get_coeff_c]
    type = ParsedAux
    variable = coeff_c
    function = if(swater>0.001,aperture*swater*1251*r_m,aperture*1251*r_m)
    args = 'aperture swater r_m' 
    block = 'rocksalt fluid' 
    execute_on = 'timestep_begin'
  []
  [get_r_m]
    type = ParsedAux
    variable = r_m
    function = if(xnacl<=0.3670,0.0,10)
    args = 'xnacl'        
   block = 'rocksalt fluid' 
    execute_on = 'timestep_begin timestep_end'
  []
  [get_time]
    type = FunctionAux
    variable = dtime_var
    function = time_func           
    block = 'rocksalt fluid' 
    execute_on = 'initial linear timestep_end'
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
  [aperture_input]
    type = SBParameterFieldIC
    variable = aperture
    file_name = kField.txt         #ParameterField.txt  #  #gap000_N18S_0MPa_40_50_mesh.txt  #
 #   type = RandomIC
 #   min = 0    #0
 #   max = 1.57 #3.1415926535
 #   legacy_generator = false
 #   variable = aperture
  []
#  [load_yy]
#    type = FunctionIC
#    variable = pwater     #disp_y
#    function = load_y
#     boundary = 'back'
#  []
#  [load_zz]
#    type = FunctionIC
#    variable = pwater     #disp_z
#    function = load_z
#     boundary = 'top'
#  []
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
  [time_derivative_xnacl]
    type = PorousFlowMassTimeDerivative
    variable = xnacl
  []
  [flux_water_xnacl]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    use_displaced_mesh = false
    variable = xnacl 
  []
  [time_derivative_zi]
    type = PorousFlowMassTimeDerivative
    variable = zi
  []
  [flux_water_zi]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    use_displaced_mesh = false
    variable = zi
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
#  [poro_x_casing]
#    type = PorousFlowEffectiveStressCoupling
#    variable = disp_x
#    component = 0
#    use_displaced_mesh = false
#    biot_coefficient =  '0'                      
#    block = 'casing'
#  []
  [poro_x_fluid]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_x
    component = 0
    use_displaced_mesh = false
    biot_coefficient = '1'
    block = 'fluid'
  []
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
#  [poro_y_casing]
#    type = PorousFlowEffectiveStressCoupling
#    variable = disp_y
#    component = 1
#    use_displaced_mesh = false
#    biot_coefficient = '0'              
#    block = 'casing '
#  []
  [poro_y_fluid]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_y
    component = 1
    use_displaced_mesh = false
    biot_coefficient = '1'
    block = 'fluid'
  []
  [grad_stress_z]
    type = StressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = false    
  []
  [poro_z_rocksalt]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    component = 2
    use_displaced_mesh = false
    biot_coefficient = '1'
    block = 'rocksalt'
  []
#  [poro_z_casing]
#    type = PorousFlowEffectiveStressCoupling
#    variable = disp_z
#    component = 2
#    use_displaced_mesh = false
#    biot_coefficient = '0'                     
#    block = 'casing'
#  []
  [poro_z_fluid]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    component = 2
    use_displaced_mesh = false
    biot_coefficient = '1'
    block = 'fluid'
  []
  
  [react2_xnacl]
    type = PorousFlowHeatMassTransfer
    variable = xnacl
    v = xnacl_control
    transfer_coefficient = coeff_c
    save_in = salt_change_reaction
  []
  [react2_pwater]
    type = PorousFlowHeatMassTransfer
    variable = pwater
    v = xnacl_control
    transfer_coefficient = coeff_c
    save_in = salt_change_reaction
  []
  [react2_zi]
    type = PorousFlowHeatMassTransfer
    variable = zi
    v = xnacl_control
    transfer_coefficient = coeff_c
    save_in = salt_change_reaction
  []
[]


[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = temperature
    use_displaced_mesh = false
  []
  [brineCo2Properties]
    type = PorousFlowFluidState
    gas_porepressure = pwater #pgas
    z = zi
    temperature = temperature
    xnacl = xnacl
    capillary_pressure = pc
    fluid_state = fs
    block = 'rocksalt fluid' 
  []
  [relperm_water0]
    type = PorousFlowRelativePermeabilityCorey
    n = 0.0
    s_res = 0.1
    sum_s_res = 0.2
    phase = 0
    block = 'rocksalt fluid' 
  []
  [relperm_water1]
    type = PorousFlowRelativePermeabilityCorey
    n = 0.0
    s_res = 0.1
    sum_s_res = 0.2
    phase = 1
    block = 'rocksalt fluid' 
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
     type = ergsEmbeddedOrthotropicFracturePermeability #PFOrthoEM
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-21 
     alpha = " 0.024 0.024 0.024"                   
     eps0 =  "-3e-5 0 0  0 -3e-5 0  0 0 -3e-5" 
     Aperture = aperture
     block = rocksalt
  []
    [permeability_casing]
     type = PorousFlowEmbeddedOrthotropicFracturePermeability #PFOrthoEM ergsEmbeddedOrthotropicFracturePermeability  #
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     alpha = " 0.024 0.024 0.024"                   
     eps0 =  "-3e-5 0 0  0 -3e-5 0  0 0 -3e-5" 
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-23 
    fix_rad_xy = 0 # 0.785398 #1.5708
     fix_rad_yz = 0 # 0.785398 #1.5708   
  #  Aperture = aperture
     block = 'casing'
    []
    [permeability_fluid]
     type = ergsEmbeddedOrthotropicFracturePermeability #PFOrthoEM
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     alpha = " 0.024 0.024 0.024"                   
     eps0 =  "-3e-5 0 0  0 -3e-5 0  0 0 -3e-5" 
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-7
     fix_rad_xy = 0 # 0.785398 #1.5708
     fix_rad_yz = 0 # 0.785398 #1.5708
     Aperture = aperture
     block = 'fluid'
  []
  
    [embedded_frac_aperture_evolution]
    type = ergsAperture
    satLIQUID = swater
    Xnacl = xnacl
    XEQ = 0.367
    rm = r_m
    Dt = dtime_var
    initial_fracture_aperture_input = aperture
    block = 'fluid rocksalt'
    []

  [elasticity_tensor_rocksalt]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.5e10                   #2.5e10
    poissons_ratio = 0.25
    block = rocksalt
 []
  [elasticity_tensor_casing]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.1e11                     #2.1e11
    poissons_ratio =  0.30 
    block = 'casing'
  []
  [elasticity_tensor_fluid]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e11                  # 1e-18 #2e10   #1e-4 #1e-4                    #1e-9
    poissons_ratio = 0.0                         #0.00
    block = 'fluid'
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
  [stress_rocksalt]
    type = ComputeLinearElasticStress  
  []
   [effective_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    block = ' rocksalt fluid'     
  []
  [volumetric_strain]
    type = PorousFlowVolumetricStrain     
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
    [externalBoundary_compressive_stress_z] 
    type = Pressure
    function = 8e6  #load_z  #                   
    variable = disp_z
    use_displaced_mesh = false
    boundary = 'top'
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
  
  [time_func]
    type = ParsedFunction
    vals = 'get_dtime_var'
    vars = 'dtime_var'
    value = 'dtime_var'
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

  
 [get_dtime_var]
  type = TimestepSize
  execute_on = 'initial linear timestep_end'
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








