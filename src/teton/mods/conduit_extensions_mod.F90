#include "macros.h"
!=======================================================================
! Extensions to conduit fortran interface
!
! These will be added to future conduit version.
!=======================================================================

module conduit_extensions_mod

  interface
    integer(c_int) function conduit_node_fetch_path_as_int(cnode, path) bind(C)
      use, intrinsic :: iso_c_binding, only : C_INT, C_PTR, C_CHAR
      implicit none
      type(C_PTR), value, intent(IN) :: cnode
      character(kind=C_CHAR), intent(IN) :: path(*)
      integer(kind=C_INT) :: res
    end function conduit_node_fetch_path_as_int
  end interface

contains

!--------------------------------------------------------------------------
  function conduit_node_obj_fetch_path_as_int(obj,path) result(res)
      use iso_c_binding
      use conduit
      use conduit_obj
      implicit none

      class(node) :: obj
      character(*) :: path
      integer(C_INT) :: res
      res = conduit_node_fetch_path_as_int(obj%cnode, trim(path) // C_NULL_CHAR)
  end function conduit_node_obj_fetch_path_as_int

end module conduit_extensions_mod
