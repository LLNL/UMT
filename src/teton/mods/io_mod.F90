! I/O Module:  various input/output options and units.

module io_mod
implicit none

  integer, parameter, public :: nopac = 4, &
                                nin   = 5, &
                                nout  = 6, &
                                neout = 8, &
                                nbout = 9

end module io_mod
