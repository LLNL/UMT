! Flags Module:  Flags for communicating with host codes

module flags_mod

use ISO_C_BINDING

private

  INTEGER(C_INT), parameter, public :: &
       convControl_fluxIter = 11, &
       convControl_tempIter = 12, &
       convControl_invalid = 0

  INTEGER(C_INT), parameter, public :: &
       dtControl_radTemp  = 21, &
       dtControl_elecTemp = 22, &
       dtControl_compton  = 23, &
       dtControl_slowConv = 24, &
       dtControl_noConv   = 25, &
       dtControl_minDt    = 26, &
       dtControl_maxDt    = 27, &
       dtControl_none     = 28, &
       dtControl_invalid  = 0

  INTEGER(C_INT), parameter, public :: &
       bcType_none    = 31, &
       bcType_refl    = 32, &
       bcType_shared  = 33, &
       bcType_temp    = 34, &
       bcType_vac     = 35, &
       bcType_fds     = 36, &
       bcType_invalid = 0

  INTEGER(C_INT), parameter, public :: &
       geometry_slab     = 41, &
       geometry_sphere   = 42, &
       geometry_cylinder = 43, &
       geometry_rz       = 44, &
       geometry_xyz      = 45, &
       geometry_invalid  = 0

  INTEGER(C_INT), parameter, public :: &
       bcShape_none    = 51, &
       bcShape_iso     = 52, &
       bcShape_fds     = 53, &
       bcShape_invalid = 0

  INTEGER(C_INT), parameter, public :: &
       comptonType_none                = 61, &
       comptonType_FP                  = 62, &
       comptonType_Boltzmann           = 63, &
       comptonType_Thomson             = 64, &
       comptonType_External_Model      = 65, &
       comptonType_invalid             = 0

end module flags_mod
