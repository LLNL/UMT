#include "macros.h"
!=======================================================================
! Datastore Module
!
! This module provides a centralized location to store or publish data
! that teton wants to share with the C++ layer or other codes.
!
! It relies on the Conduit library, which supports a tree of key/value
! nodes.
!
! The conduit node is encapsulated in a 'datastore' type, so this could be
! refactored to use sidre or another key/value tree in the future.
!
! The conduit node can be accessed directly via
! get_datastore_cnode (Fortran)
! teton_get_datastore_cnode (C/C++)
!
! The get_datastore_cnode functions will automatically create and initialize
! the root conduit node if it is not already created.
!=======================================================================
 
module Datastore_mod 
  use conduit_obj, only : node
  use iso_c_binding, only : c_ptr
  implicit none

  private

  type datastore_type
    type(node) :: root
    logical :: is_initialized = .FALSE.

    contains
      procedure :: save_hdf5
      procedure :: initialize

  end type datastore_type

  type(datastore_type), public :: theDatastore

contains

!=======================================================================
! Create and initialize conduit node, if it doesn't already exist.
!=======================================================================
  subroutine initialize(self)
    use conduit_obj, only : conduit_node_obj_create
    class(datastore_type) :: self

    if (.not. theDatastore%is_initialized) then
        theDatastore%root = conduit_node_obj_create()
        theDatastore%is_initialized = .TRUE.
    endif

  end subroutine

!=======================================================================
! Save datastore contents to hdf5.
! Uses conduit relay library to save.
!=======================================================================
  subroutine save_hdf5(self,path)
    use iso_c_binding, only : c_null_ptr, c_char
    class(datastore_type) :: self
    character(kind=c_char), intent(in) :: path(*)
    
    interface
      subroutine c_conduit_relay_io_save(cnode, path, protocol, coptions) bind(C, name="conduit_relay_io_save")
        use iso_c_binding, only : c_ptr, c_char
        implicit none

        type(c_ptr), value, intent(in) :: cnode
        character(kind=c_char), intent(in) :: path(*)
        character(kind=c_char), intent(in) :: protocol(*)
        type(c_ptr), value, intent(in) :: coptions
      end subroutine c_conduit_relay_io_save
    end interface
    
    call c_conduit_relay_io_save(theDatastore%root%cnode, path, "hdf5", c_null_ptr)

  end subroutine


!=======================================================================
! Get C pointer to conduit node
!=======================================================================
  type(c_ptr) function get_datastore_cptr() result(cnode_cptr) bind(c, name="teton_get_datastore_cptr")
    call theDatastore%initialize() ! Will only initialize node if not already.
    cnode_cptr = theDatastore%root%cnode
  end function get_datastore_cptr

end module Datastore_mod
