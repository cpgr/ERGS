[Mesh]
  [annular]
    type = AnnularMeshGenerator
    nr = 10
    rmin = 1.0
    rmax = 10
    growth_r = 1.4
    nt = 4
    dmin = 0
    dmax = 90
  []

  [./extrude]
    type = MeshExtruderGenerator
    input = 'annular'
    num_layers = 3
    extrusion_vector = '0 0 12' 
   # existing_subdomains = '3'
   # layers = '3'
   # new_ids = 2
   # show_info = false
    bottom_sideset = 'bottom'
    top_sideset = 'top'
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
  [roller_tmax]
    type = DirichletBC
    variable = u
    value = 0
    boundary = dmax
  []
  [roller_tmin]
    type = DirichletBC
    variable = u
    value = 0
    boundary = dmin
  []
  [pinned_top_bottom_x]
    type = DirichletBC
    variable = u
    value = 0
    boundary = 'top bottom'
  []
  [pinned_top_bottom_y]
    type = DirichletBC
    variable = u
    value = 0
    boundary = 'top bottom'
  []
[]

[Executioner]
  type = Steady

  solve_type = 'PJFNK'
[]

[Outputs]
  file_base = out_quad_angle
  exodus = true
[]