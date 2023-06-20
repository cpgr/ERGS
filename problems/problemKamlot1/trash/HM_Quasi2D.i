# Quasi-2D Model
# Fluid Property
rho_w = 1173     # kg.m^3
mu_w = 1.252e-3  # Pa.s

phi = 0.26
perm = 2.96e-13     # m^2

# injection rate = 10 m^3 /day
q_water = -10    # -ve means injection
injection_rate2 = ${fparse q_water*rho_w/86400}  # Kg/s
injection_rate3 = ${fparse q_water*rho_w/(86400*6.0)}  # Kg/s/area
injection_rate = ${fparse q_water*rho_w/(2*3.14*0.2*6*86400)}  # Kg/s/area  in the area 2pi*0.2*6
simulation_time = ${fparse 10*86400} # 1000 hours

#Kb =  # drained bulk modulus of porous medium
#Ks =  # Ks =inf then biot_coeff=1 -> incomressible rock bulk modulus of solid
#biot_coeff = ${fparse 1.0 -(Kb/Ks)}

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 100
  bias_x = 1.003
  xmin = 0.2
  xmax = 200
  ymin = 0
  ymax = 1
  ny = 1
  zmin = -1506
  zmax = -1500
  nz = 5
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 -9.81'
  displacements = 'disp_x disp_y disp_z'
  biot_coefficient = 1
[]

[Variables]
  [pwater]
  []
  [disp_x]
    scaling = 1e-6
  []
  [disp_z]
    scaling = 1e-6
  []
[]

[AuxVariables]
  [disp_y]
  []
[]

[ICs]
  [porepressure]
    type = FunctionIC
    function = ppic
    variable = pwater
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowFullySaturatedMassTimeDerivative
    multiply_by_density = false
    coupling_type = HydroMechanical
    variable = pwater
    use_displaced_mesh = false
  []
  [flux]
    type = PorousFlowFullySaturatedDarcyBase
    multiply_by_density = false
    variable = pwater
    use_displaced_mesh = false
  []
  [grad_stress_x]
    type = StressDivergenceTensors
    variable = disp_x
    use_displaced_mesh = false
    component = 0
  []
  [grad_stress_z]
    type = StressDivergenceTensors
    variable = disp_z
    use_displaced_mesh = false
    component = 2
  []
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_x
    use_displaced_mesh = false
    component = 0
  []
  [poro_z]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    use_displaced_mesh = false
    component = 2
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pwater disp_x disp_z'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    m = 0.8
    alpha = 1e-5
  []
[]

[Modules]
  [FluidProperties]
    [water]
      type = SimpleFluidProperties
      bulk_modulus = 2.27e14
      density0 = ${rho_w}
      viscosity = ${mu_w}
      cv = 4149.0
      cp = 4149.0
      porepressure_coefficient = 0.0
      thermal_expansion = 0
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
  [ppss]
    #type = PorousFlow1PhaseP
    type = PorousFlow1PhaseFullySaturated
    porepressure = pwater
    use_displaced_mesh = false
  []
  [massfrac]
    type = PorousFlowMassFraction
    use_displaced_mesh = false
  []
  [water]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
  []
  [porosity_reservoir]
    type = PorousFlowPorosityConst
    porosity = ${phi}
    use_displaced_mesh = false
  []
  [permeability_reservoir]
    type = PorousFlowPermeabilityConst
    permeability = '${perm} 0 0  0 ${perm} 0  0 0 ${perm}'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 0 # unimportant in this fully-saturated situation
    phase = 0
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 10e6
    poissons_ratio = 0.3
    use_displaced_mesh = false
  []
  [strain]
    type = ComputeSmallStrain
    #eigenstrain_names = 'ini_stress'
  []
  # [ini_strain]
  #   type = ComputeEigenstrainFromInitialStress
  #   initial_stress = 'ini_xx 0 0  0 ini_xx 0 0 0 ini_zz'
  #   eigenstrain_name = ini_stress
  # []
  [stress]
    type = ComputeLinearElasticStress
  []
  [eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    use_displaced_mesh = false
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
    use_displaced_mesh = false
  []
[]

[BCs]
  [outer_pressure_fixed]
    type = FunctionDirichletBC
    boundary = right
    variable = pwater
    function = ppic
    use_displaced_mesh = false
  []

  [water_injection]
    type = PorousFlowSink  # 1kg/s/m^2
    boundary = left
    variable = pwater
    fluid_phase = 0
    flux_function = 'min(t/100.0,1)*${injection_rate}'
    use_displaced_mesh = false

  []
  [u_left]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'left'
    use_displaced_mesh = false
  []

  [u_top_bottom]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'back front'
    use_displaced_mesh = false
  []
  [u_right]
    type = Pressure
    variable = disp_x
    boundary = 'right'
    function = ux_total
    use_displaced_mesh = false
  []
[]


[Preconditioning]
  active = 'smp'
  [smp]
    type = SMP
    full = true
    # petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2               1E-11       1E-8        50'
  []
  [mumps]
    type = SMP
    full = true
    # petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    # petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -snes_rtol -snes_atol -snes_max_it'
    # petsc_options_value = 'gmres      lu       mumps                         NONZERO               1E-5       1E-2       50'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = ${simulation_time}
  dtmax = 3600
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.5
  []
[]

[Outputs]
  print_linear_residuals = false
  sync_times = '0 86400  172800  259200  345600  432000  518400  604800  691200
  777600  864000  950400 1036800 1123200 1209600 1296000 1382400 1468800 1555200
  1641600 1728000 1814400 1900800 1987200 2073600 2160000 2246400 2332800 2419200
  2505600 2592000 2678400 2764800 2851200 2937600 3024000 3110400 3196800 3283200 3369600 3456000 3542400'
  perf_graph = false
  exodus = true
  # [csv]
  #   type = CSV
  #   sync_only = true
  # []
[]


[Functions]
  [ppic]
    type = ParsedFunction
    #value = '10e6 - ${rho_w}*9.81*z'
    value = '-11507.1*z'
  []
  [ini_xx]
    type = ParsedFunction
    vars = 'g biot'
    vals = '9.81 1.0'
    value = '0.8*(2650*g*z + biot*(10e6 - ${rho_w}*g*z))'
  []
  [ini_zz]
    type = ParsedFunction
    vars = 'g biot'
    vals = '9.81 1.0'
    value = '2650*g*z + biot*(10e6 - ${rho_w}*g*z)'
  []
  [uz_total]
    type = ParsedFunction
    value = '22229.1*z'
  []
  [ux_total]
    type = ParsedFunction
    value = '15560.4*z'
  []
[]