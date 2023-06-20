[Mesh]
  type = GeneratedMesh
  dim = 2
  nx =40
  ny =50
  xmax = 0.02512
  ymax = 0.0327
[]


[AuxVariables]
 [aperture]
  order = CONSTANT
  family = MONOMIAL
 []
[]


[ICs]
  [./aperture_input]
    type = SBParameterFieldIC
    variable = aperture
    file_name = gap000_N18S_0MPa_40_50_mesh.txt
  [../]
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
  print_linear_residuals = false
  perf_graph = true
  csv = true
  checkpoint = true
  [./out]
    type = Exodus
    execute_on = 'initial timestep_end'
  [../]
[]