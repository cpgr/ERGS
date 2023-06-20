[Mesh]
[gen]
  type = GeneratedMeshGenerator
  dim = 2
  xmin = 5
  xmax = 10000       
  nx = 100                  
  bias_x = 1.01
  ymin = -500
  ymax = 0.0
  ny = 1
  []
#  coord_type = RZ
[]

[FluidProperties]
  [brine]
    type = BrineFluidProperties
  []
  [co2]
    type = CO2FluidProperties
  []
  [water]
    type = Water97FluidProperties
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zi Xnacl temperature disp_x disp_y'       
    number_fluid_phases = 2
    number_fluid_components = 3
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  []
  [fs]
   type = PorousFlowBrineCO2
    brine_fp = brine
    co2_fp = co2 
    capillary_pressure = pc
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  displacements = 'disp_x disp_y'
  gravity = '0.8  0.8  0.8'
  temperature_unit = Celsius
[]

[Variables]
  [pgas]
    initial_condition = 6e6
  []
  [zi]
  []
  [Xnacl]
    initial_condition =  0.3
  []
  [temperature]
    initial_condition = 200 #275.55 #
  []
  [disp_x]
    scaling = 1E-5
  []
  [disp_y]
    scaling = 1E-5
  []
[]

[ICs]
 [zi]
    type = PorousFlowFluidStateIC
    variable = zi
    saturation = 0.45
    gas_porepressure = pgas
    temperature = temperature
    fluid_state = fs
  []
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
  #[aperture]
  #  type = ergsApertureIC
  #  variable = apertureEvol
  #[]
[]

[AuxVariables]
#  [temperature]
#    initial_condition = 100 
#  []
#  [Xnacl]
#    initial_condition = 0.3
#  []
 [sgas]
   order = CONSTANT
   family = MONOMIAL
 []
  [sbrine]
   order = CONSTANT
   family = MONOMIAL
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
 [ControlXnacl]
  initial_condition = 0.27000000
 []
 [ResidualXnacl]
 []
  [TCoeff]
  initial_condition = 5.0
 []
  #[apertureEvol]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[]
[]

[AuxKernels]
  [sgas]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 1
    variable = sgas
    execute_on = 'TIMESTEP_END'
  []
  [sbrine]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = sbrine
    execute_on = 'TIMESTEP_END'
  []
   [randm_rad_XY]
    type = MaterialRealAux
    variable = randm_rad_XY
    property = random_xy_rotation_angle_for_each_element_qp
    execute_on = timestep_end
  []
[]


[Functions]
 [condition]
  type = ParsedFunction
  value = 'xnacl < 1'
  vars = 'xnacl'
  vals = 'XnaclThresholdFlag'      #This is a postprocessor!
 []
[]


[Controls]
 [Xnacl_threshold]
  type = ConditionalFunctionEnableControl
  conditional_function = condition
  disable_objects = 'Kernel::reaction'
  execute_on = 'initial timestep_begin'
 []
[]


[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative  
    variable = pgas
    fluid_component = 0
  []
  [adv0]
    type = PorousFlowAdvectiveFlux
    variable = pgas
    fluid_component = 0
  []
  [mass1]
    type = PorousFlowMassTimeDerivative  
    variable = zi
    fluid_component = 1
  []
  [adv1]
    type = PorousFlowAdvectiveFlux
    variable = zi
    fluid_component = 1
  []
  [mass2]
    type = PorousFlowMassTimeDerivative  
    variable = Xnacl
    fluid_component = 2
  []
  [adv2]
    type = PorousFlowAdvectiveFlux
   variable = Xnacl
    fluid_component = 2
  []
  [energy]
    type = PorousFlowEnergyTimeDerivative       
    variable = temperature
  []
  [heat_adv]
    type = PorousFlowHeatAdvection
    variable = temperature
  []  
  [conduction]
   type = PorousFlowHeatConduction
    variable = temperature
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

 [reaction]
    type = PorousFlowHeatMassTransfer
    variable = Xnacl
    v = ControlXnacl
    transfer_coefficient = TCoeff
    save_in = ResidualXnacl
 []
[]


[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = temperature
  []
  [brineCo2Properties]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = zi
    temperature = temperature
    xnacl = Xnacl
    capillary_pressure = pc
    fluid_state = fs
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.05
  []
  [permeability]
     type = ergsEmbeddedOrthotropicFracturePermeability
     Random_field = true
     rotation_angleXY = rotxy_aux
     rotation_angleYZ = rotyz_aux
     alpha =  "1e-2 1e-2 1e-2"    
     eps0 = "1e-5 1e-5 1e-5" 
     N = "1 0 0  0 1 0  0 0 1"
     km = 1e-18 
     fix_rad_xy = 0 # 0.785398 #1.5708
     fix_rad_yz = 0 # 0.785398 #1.5708
  #   Aperture = apertureEvol
  #  type = PorousFlowPermeabilityConst
  #  permeability = '50e-10 0 0   0 50e-10 0   0 0 50e-10'
  #   type = PorousFlowPermeabilityVP
  #   permeability0 = '50e-10 0 0   0 50e-10 0   0 0 50e-10'
  #   solid_sat = ssolid
  []
  [relperm0]
    type = PorousFlowRelativePermeabilityCorey
    n = 4
    s_res = 0.30
    sum_s_res = 0.35
    phase = 0
  []
  [relperm1]
    type = PorousFlowRelativePermeabilityCorey
    n = 4
    s_res = 0.30
    sum_s_res = 0.35
    phase = 1
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 1000
    density = 2600
  [] 
  [rock_thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2 0 0  0 2 0  0 0 2'
 []
  [embedded_frac_aperture]
    type = ergsAperture
    Xnacl = Xnacl
    XEQ = 0.277
    satLIQUID = 0.18 #sbrine
    rm = 0.277 
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

[DiracKernels]
  [fluid_produce]   
    type = PorousFlowSquarePulsePointSource
    point = '5 -250 0'
    mass_flux = -65           # note: 65kg/s 0.065m/s  1.307e-7
   variable = pgas
 []
[]
#[BCs]
#  [produce_heat]
#    type = PorousFlowSink
#    variable = pgas
#    boundary = left
#    flux_function = 0.065
#    fluid_phase = 0
#    use_enthalpy = true
#  []
#[]

[VectorPostprocessors]
  [vars]
    type = NodalValueSampler
    sort_by = x
    variable = 'pgas zi Xnacl'
    execute_on = 'timestep_end'
    outputs = spatial
  []
  [auxvars]
    type = ElementValueSampler
    sort_by = x
    variable = 'sbrine sgas'
    execute_on = 'timestep_end'
    outputs = spatial
  []
[]

[Postprocessors]
  [Pgas]
    type = PointValue
   point =  '5 0 0'
    variable = pgas
    execute_on = 'initial TIMESTEP_END'
  []
  [sgas]
    type = PointValue
    point =  '5 0 0'
    variable = sgas
    execute_on = 'initial TIMESTEP_END'
  []
  [sbrine]
    type = PointValue
    point =  '5 0 0'
    variable = sbrine
    execute_on = 'initial TIMESTEP_END'
  []
  [Z]
   type = PointValue
    point =  '5 0 0'
    variable = zi
    execute_on = 'initial TIMESTEP_END'
  []
  [Xnacl]
   type = PointValue
    point =  '5 0 0'
    variable = Xnacl
    execute_on = 'initial TIMESTEP_END'
  []
  [XnaclThresholdFlag]
   type = ergsVariableThreshold
    variable = Xnacl
    threshold = 0.27
    execute_on = 'initial TIMESTEP_END'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres bjacobi lu NONZERO'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 1
  end_time = 2e6
  nl_max_its = 25
  l_max_its = 100
  dtmax = 1e5
  nl_abs_tol = 1e-12
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
  []
[]

[Outputs]
  exodus = true
  sync_times = '1e4 1e5 2e6'
 [time]
    type = CSV
  []
  [spatial]
    type = CSV
    sync_only = true
  []
[]

