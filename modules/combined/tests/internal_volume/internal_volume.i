#
# Internal Volume Test
#
# This test is designed to compute the internal volume of a space considering
#   an embedded volume inside.
#
# The mesh is composed of one block (1) with an interior cavity of volume 8.
#   Block 2 sits in the cavity and has a volume of 1.  Thus, the total volume
#   is 7.
# The internal volume is scaled by two and adjusted by negative seven.  Thus,
#   the net result is seven.
#

[Mesh]#Comment
  file = internal_volume.e
  displacements = 'disp_x disp_y disp_z'
[]

[Functions]
  [./step]
    type = PiecewiseLinear
    x = '0. 1. 2. 3.'
    y = '0. 0. 1e-2 0.'
    scale_factor = 0.5
  [../]
[]

[Variables]

  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]

[]

[SolidMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
[]

[BCs]

  [./no_x]
    type = DirichletBC
    variable = disp_x
    boundary = 100
    value = 0.0
  [../]

  [./no_y]
    type = DirichletBC
    variable = disp_y
    boundary = 100
    value = 0.0
  [../]

  [./prescribed_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 100
    function = step
  [../]
[]

[Materials]
  [./stiffStuff]
    type = LinearIsotropicMaterial
    block = 1

    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    youngs_modulus = 1e6
    poissons_ratio = 0.3
    thermal_expansion = 1e-5
    t_ref = 400.
  [../]

  [./stiffStuff2]
    type = LinearIsotropicMaterial
    block = 2

    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    youngs_modulus = 1e6
    poissons_ratio = 0.3
    thermal_expansion = 1e-5
    t_ref = 400.
  [../]
[]

[Executioner]

  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'




  nl_abs_tol = 1e-10

  l_max_its = 20

  start_time = 0.0
  dt = 1.0
  #num_steps = 3
  end_time = 3.0
[]

[Postprocessors]
  [./internalVolume]
    type = InternalVolume
    boundary = 100
    scale_factor = 2
    addition = -7
  [../]

  [./dispZ]
    type = ElementAverageValue
    block = '1 2'
    variable = disp_z
  [../]
[]

[Outputs]
  file_base = out
  exodus = true
  csv = true
  output_on = 'initial timestep_end'
  [./console]
    type = Console
    perf_log = true
    output_on = 'timestep_end failed nonlinear linear'
  [../]
[]
