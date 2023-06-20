#TOUGH2 EOS3/4 example

[Mesh]
#  [efmcube]
#  type = FileMeshGenerator
#  file = ewasg.msh
#  []
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
  coord_type = RZ
[]


[FluidProperties]
  [co2]
    type = CO2FluidProperties
  []
  [brine]
    type = BrineFluidProperties
  []
  [water]
    type = Water97FluidProperties
  []
  [watertab]
    type = TabulatedFluidProperties
    fp = water
    temperature_min = 273.15
    temperature_max = 573.15
    fluid_property_file = water_fluid_properties.csv
    save_file = false
  []
[]


[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zi temperature'
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
  gravity = '0.8 0.8 0.8'
[]


[Variables]
  [pgas]
    initial_condition = 6e6
  []
  [temperature]
    initial_condition = 295.55
  []
  [zi]
 #   initial_condition =  0.21  #0.45
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
[]

[AuxVariables]
  [Xnacl]
    initial_condition = 0.3
  []
 [sgas]
   order = CONSTANT
   family = MONOMIAL
 []
  [sbrine]
   order = CONSTANT
   family = MONOMIAL
  []
[]


[AuxKernels]
  [sgas]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 1
    variable = sgas
  []
  [sbrine]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = sbrine
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
 # [mass2]
 #   type = PorousFlowMassTimeDerivative
 #   variable = Xnacl
 #   fluid_component = 2
 # []
 # [adv2]
 #   type = PorousFlowAdvectiveFlux
 #   variable = Xnacl
 #   fluid_component = 2
 # []
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
[]


#[DiracKernels]
#  [fluid_produce]   
#    type = PorousFlowSquarePulsePointSource
#    point = '0 0 0'
#    mass_flux = -65
#    variable = zi
#  []
#[]


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
    type = PorousFlowPermeabilityConst
    permeability = '50.e-15 0 0   0 50.e-15 0   0 0 50.e-15'
  []
  [relperm0]
    type = PorousFlowRelativePermeabilityCorey
    n = 4
    s_res = 0.1
    sum_s_res = 0.2
    phase = 0
  []
  [relperm1]
    type = PorousFlowRelativePermeabilityCorey
    n = 4
    s_res = 0.1
    sum_s_res = 0.2
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
[]


[BCs]
  [produce_heat]
    type = PorousFlowSink
    variable = temperature
    boundary = right
    flux_function = 65
    fluid_phase = 0
    use_enthalpy = true
  []
[]


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
  []
  [sgas]
    type = PointValue
    point =  '5 0 0'
    variable = sgas
  []
  [sbrine]
    type = PointValue
    point =  '5 0 0'
    variable = sbrine
  []
  [Xnacl]
    type = PointValue
    point =  '5 0 0'
    variable = Xnacl
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
  end_time = 2e6
  nl_max_its = 25
  l_max_its = 100
  dtmax = 1e5
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 100
  []
[]

[Outputs]
  sync_times = '1e4 1e5 2e6'
  [time]
    type = CSV
  []
  [spatial]
    type = CSV
    sync_only = true
  []
[]




