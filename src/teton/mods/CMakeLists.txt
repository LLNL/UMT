############################
# Embed version info
############################
configure_file (cmake_defines_mod.F90.in
                cmake_defines_mod.F90)

target_sources( teton PRIVATE
  AngleSet_mod.F90
  BoundaryList_mod.F90
  Boundary_mod.F90
  CommSet_mod.F90
  Communicator_mod.F90
  Datastore_mod.F90
  Editor_mod.F90
  Geometry_mod.F90
  GreyAcceleration_mod.F90
  GroupSet_mod.F90
  Material_mod.F90
  MemoryAllocator_mod.F90
  Options_mod.F90
  Profile_mod.F90
  QuadratureData_mod.F90
  QuadratureList_mod.F90
  Quadrature_mod.F90
  RadIntensity_mod.F90
  SetData_mod.F90
  Size_mod.F90
  TimeStepControls_mod.F90
  ZoneSet_mod.F90
  constant_mod.F90
  flags_mod.F90
  io_mod.F90
  iter_control_list_mod.F90
  iter_control_mod.F90
  kind_mod.F90
  radconstant_mod.F90
  system_info_mod.F90
  ${CMAKE_CURRENT_BINARY_DIR}/cmake_defines_mod.F90
  default_iter_controls_mod.F90
)
