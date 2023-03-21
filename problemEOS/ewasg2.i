#TOUGH2 EWASG example

[Mesh]
[gen]
  type = GeneratedMeshGenerator
  dim = 2
  xmin = 5
  xmax = 10000 #10000
  nx = 100        #100    
  bias_x = 1.01
  ymin = -500
  ymax = 0.0
  ny = 1
  []
  coord_type = RZ
[]


[FluidProperties]
  [brine]
    type = BrineFluidPropertiesBeta
  []
  [co2]
    type = CO2FluidProperties
  []
  [water]
    type = Water97FluidProperties
  []
  [tabulated_water]
    type = TabulatedBicubicFluidProperties
    fp = water
    temperature_max = 900
    temperature_min = 250
    pressure_max = 1E8
    interpolated_properties = 'density viscosity enthalpy internal_energy'
    fluid_property_file = water97_tabulated_11.csv
  []
  [true_co2]
    type = CO2FluidProperties
  []
  [tabulated_co2]
    type = TabulatedBicubicFluidProperties
    fp = true_co2
    temperature_max = 900
    temperature_min = 250
    pressure_max = 1E8
    interpolated_properties = 'density viscosity enthalpy internal_energy'
    fluid_property_file = co2_tabulated_11.csv
  []
[]


[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zi temperature'
    number_fluid_phases = 3
    number_fluid_components = 3
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  []
  [fs]
    type = PorousFlowBrineSaltCO2
    brine_fp = brine
    water_fp = water
    co2_fp = co2         
    capillary_pressure = pc
  []
[]


[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0.8  0.8  0.8'
  temperature_unit = Celsius
[]


[Variables]
  [pgas]
    initial_condition = 6e6
  []
  [temperature]
    initial_condition = 275.55
  []
  [zi]
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
    permeability = '50e-15 0 0   0 50e-15 0   0 0 50e-15'
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
  [relperm2]
    type = PorousFlowRelativePermeabilityConst
    kr = 0
    phase = 2
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


[DiracKernels]
  [fluid_produce]   
    type = PorousFlowSquarePulsePointSource
    point = '5 -250 0'
    mass_flux = -1.3e-4   # -65  #-0.033   # 
    variable = pgas
  []
[]
#[BCs]
#  [produce_heat]
#    type = PorousFlowSink
#    variable = pgas
#    boundary = left
#    flux_function = 65
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
  nl_abs_tol = 1e-20
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 100
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




