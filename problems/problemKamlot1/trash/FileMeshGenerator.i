[Mesh]
  [EFMcube]
    type = FileMeshGenerator
    file = gmsh_simple_cube.inp 
  []
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
  [roller_front]
    type = DirichletBC
    variable = u
    value = 0
    boundary = 'front' 
  []
    [roller_back]
    type = DirichletBC
    variable = u
    value = 0
    boundary = 'back' 
  []
  [pinned_top]
    type = DirichletBC
    variable = u
    value = 0
    boundary = 'top'
  []
    [pinned_bottom]
    type = DirichletBC
    variable = u
    value = 0
    boundary = 'bottom'
  []
  [pinned_L]
    type = DirichletBC
    variable = u
    value = 100
    boundary = 'left'
  []
    [pinned_R]
    type = DirichletBC
    variable = u
    value = 0
    boundary = 'right'
  []
[]

[Executioner]
  type = Steady

  solve_type = 'PJFNK'
[]

[Outputs]
  file_base = EFPMcube
  exodus = true
[]