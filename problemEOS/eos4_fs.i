#TOUGH2 EOS3/4 example

# NOTE: set xmax = 10000 and dim = 1 or 2 for faster convergence
# The plot of primary variables versus distance from the source are different depending on 
# how the dispersive coefficients are set. It seems the best way is: '0  2.13E-5  0  2.13E-5'
# z is set to 0.01 for the initial gas saturation condition (i.e., sgas = 0.20)

[Mesh]
 [gen]
  type = GeneratedMeshGenerator
  dim = 2
  nx = 500   #500
  xmax = 10000 #1000
  xmin = 1
  ny = 1
  ymax = 4.5
  ymin = 1
  bias_x = 1.01
  []
  coord_type = RZ
[]


[FluidProperties]
  [true_water]
    type = Water97FluidProperties
  []
  [tabulated_water]
    type = TabulatedBicubicFluidProperties
    fp = true_water
    temperature_max = 600
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
    temperature_min = 250
    temperature_max = 600
    pressure_max = 1E8
    interpolated_properties = 'density viscosity enthalpy internal_energy'
    fluid_property_file = co2_tabulated_11.csv
  []
  [air]
    type = IdealGasAirFluidProperties
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas z temperature'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    m = 0.45
    alpha = 1e1 
  []
  [fs]
    type = PorousFlowWaterAir   #PorousFlowWaterNCG # 
    water_fp = true_water
    gas_fp =  air      #true_co2 #                     
    capillary_pressure = pc
  []
[]


[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
  temperature_unit = Celsius
[]


[Variables]
  [pgas]
    initial_condition = 1e5
  []
  [z]
   scaling = 1e4
  []
  [temperature]
    initial_condition = 18 # 291.15     #
  []
[]


[ICs]
  [z]
    type = PorousFlowFluidStateIC
    variable = z
    saturation = 0.2
    gas_porepressure = pgas
    temperature = temperature
    fluid_state = fs
  []
[]


[AuxVariables]
  [sgas]
   order = CONSTANT
   family = MONOMIAL
  []
  [swater]
   order = CONSTANT
   family = MONOMIAL
  []
  [xAir]
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
  [swater]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = swater
  []
  [xAir]
    type = PorousFlowPropertyAux
    variable = xAir
    property = mass_fraction
    phase = 1
    fluid_component = 1
    execute_on = timestep_end
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
   [disp0]
    type = PorousFlowDispersiveFlux
    variable = pgas
    disp_trans = '0 0'
    disp_long = '0 0'
    fluid_component = 0
    use_displaced_mesh = false
  []
  [mass1]
    type = PorousFlowMassTimeDerivative
    variable = z
    fluid_component = 1
  []
  [adv1]
    type = PorousFlowAdvectiveFlux
    variable = z
    fluid_component = 1
  []
   [disp1]
    type = PorousFlowDispersiveFlux
    variable = z
    disp_trans = '0 0'
    disp_long = '0 0'
    fluid_component = 1
    use_displaced_mesh = false
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
  [WaterAirProperties]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = z
    temperature = temperature
    capillary_pressure = pc
    fluid_state = fs
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '20.e-15 0 0   0 20.e-15 0   0 0 20.e-15'
  []
  [relperm0]
    type = PorousFlowRelativePermeabilityVG
    m = 0.45
    phase = 0
  []
  [relperm1]
    type = PorousFlowRelativePermeabilityVG
    m = 0.45
    phase = 1
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.01
  []
  [diffusivity]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '2.13E-5 2.13E-5 0 0'
    tortuosity = '0.25 0.25'
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 800
    density = 2550
  [] 
  [rock_thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2 0 0  0 2 0  0 0 2'
  []
[]


[DiracKernels]
  [heat_source]
    type = PorousFlowSquarePulsePointSource
    point = '1 2.25 0'
    mass_flux = 3e3
    variable = temperature
  []
[]
#[BCs]
#  [flux]
#    type = PorousFlowSink
#    boundary = 'left'
#    variable = temperature
#    use_mobility = false
#    use_relperm = false
#    mass_fraction_component = 1
#    fluid_phase = 0
#    flux_function = -3e3
# []
#[]


[VectorPostprocessors]
  [vars]
    type = NodalValueSampler
    sort_by = id
    variable = 'pgas temperature'
    execute_on = 'timestep_end'
    outputs = spatial
  []
  [auxvars]
    type = ElementValueSampler
    sort_by = id
    variable = 'xAir swater'
    execute_on = 'timestep_end'
    outputs = spatial
  []
[]


[Postprocessors]
  [xAir]
    type = PointValue
    point = '1 1 0'
    variable = xAir
  []
  [Pgas]
    type = PointValue
   point =  '1 1 0'
    variable = pgas
  []
  [T]
    type = PointValue
    point = '1 1 0'
    variable = temperature
  []
  [swater]
    type = PointValue
    point = '1 1 0'
    variable = swater
  []
[]


[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]


[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 1e4
  num_steps = 1e6
  end_time = 3.56E8
  nl_abs_tol = 1e-20
  nl_max_its = 40 #25
  l_max_its = 100
  dtmax = 5e6
  [TimeStepper]
     type = IterationAdaptiveDT
     dt = 10
  []
[]

[Outputs]
  exodus = true
  sync_times = '1e5 1.6e6 2.6e6 2e7 3e8'
  [time]
    type = CSV
  []
  [spatial]
    type = CSV
    sync_only = true
  []
[]




