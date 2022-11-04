#include "macros.h"
!*****************************************************************************
! Called from the host code to write/read the angle-dependent radiation      *
! intensity (Psi), corner electron temperature, and volume to/from a         *
! conduit tree.                                                              *
!                                                                            *
! This conduit tree can either be saved directly to HDF5 files, or imported  *
! into the SIDRE datastore.                                                  *
!                                                                            *
! This restart is using the SIDRE three-stage restart load model.            *
!                                                                            *
! stage 0: 'prepare for load operation'                                      *
!   Teton should prepare for the host code to load data into conduit         *
! stage 1: 'owned data has been loaded'                                      *
!   Teton can now copy out any data from conduit that it needs.  This is     *
!   typically scalar data or small arrays.                                   *
! stage 2: 'external data has been loaded'                                   *
!   Conduit has the ability to hold pointers to data owned by Teton.         *
!   At stage 2, any pointers to Teton data that were added to conduit have   *
!   now been directly populated, avoiding the need for temporary copies and  *
!   copying data.                                                            *
!*****************************************************************************
! NOTES:                                                                     *
!                                                                            *
! This restart differs from the SILO load/drop restart in that Volume is     *
! being explicitly saved/restored.  It does not need to be recalculated on   *
! the restart.                                                               *
!                                                                            *
! Note from black27:                                                         *
! - I had some difficulty getting bit-perfect restart loads when I tried     *
! recomputing Volume when restarting in Marbl.  This requires a              *
! high order->low order mesh transformation, then updating the teton mesh    *
! node positions, then calling  teton_getvolume().  Essentially a similar    *
! operation as that needed for the SILO restart.  I opted to just restart    *
! the Volume field.                                                          *
!*****************************************************************************


!*****************************************************************************
! Usage:
! A host should needs to perform the following steps for a restart.
!
! Restart load
! ---
! #1 initialize teton
! #2 call teton_conduitcheckpoint_prep_for_load()
! #3 load the non-external data into conduit
! #4 call teton_conduitcheckpoint_data_loaded()
! #5 load the external (data teton owns) into conduit
! #6 call teton_conduitcheckpoint_external_data_loaded()
! #7 call teton_conduitcheckpoint_teardown()
!
! Restart save
! ---
! call teton_conduitcheckpoint_prep_for_save()
! save the data to file
! call teton_conduitcheckpoint_teardown()
!
!*****************************************************************************
module ConduitCheckpoint_mod
   use conduit_obj, only : node
   use, intrinsic :: iso_c_binding, only : c_ptr
   use kind_mod
   implicit none

   ! Conduit tree containing checkpoint data
   type(node) :: root

contains

   type(c_ptr) function get_root_cptr() result(cptr) bind(c,name="teton_conduitcheckpoint_get_cptr")
      use, intrinsic :: iso_c_binding, only : c_associated
      use conduit_obj, only : conduit_node_obj_create
      if (.not. c_associated( root%cnode )) then
         root = conduit_node_obj_create()
      endif
      cptr = root%cnode
   end function

!*****************************************************************************
! Prepare checkpoint node for a load
!*****************************************************************************
   subroutine prepForLoad() bind(c, name="teton_conduitcheckpoint_prep_for_load")
      use, intrinsic :: iso_c_binding, only : c_double, c_size_t, c_associated
      use conduit_obj, only : conduit_node_obj_create

      if (.not. c_associated( root%cnode )) then
         root = conduit_node_obj_create()
      endif

   end subroutine

!*****************************************************************************
! Notify teton that 'conduit owned' data has been read into checkpoint node.
! This is data copied into conduit ( not pointers we've added to conduit that
! point to teton allocated data ).
! 
! Teton will also prepare the node by adding pointers to fields it wants loaded
! in on the next stage.  This implementation expects that Teton's modules
! and fields have been constructed and allocated.  Ie, that teton has been
! initialized already.
!*****************************************************************************
   subroutine dataLoaded() bind(c, name="teton_conduitcheckpoint_data_loaded")
      use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
      use, intrinsic :: iso_fortran_env, only: real32, int32
      use, intrinsic :: iso_c_binding, only : c_size_t
      use Geometry_mod
      use QuadratureList_mod
      use SetData_mod
      use Material_mod
      use Size_mod
      use MemoryAllocator_mod
      use iter_control_mod
      use iter_control_list_mod

      type(node) :: current
      type(SetData), pointer :: set
      integer :: setID
      character(len=64) :: setIDString
      integer(kind=c_size_t) :: numElements
      integer :: maxIts
      real(real32) :: nan
      real(adqt) :: relTol


      type(IterControl), pointer  :: itCon => NULL()

      ! Let's set the fields' we want loaded from the restart to NaN.  This will ensure that if these fields fail to be
      ! populated the compiler will catch this if FPE checking is enabled.
      nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)

      current = root
      TETON_VERIFY(current%fetch_path_as_int32("num_corners") == Size%ncornr, "Teton chkpt load error: Num corners does not match expected number.")
      
      ! Publish a pointer to corner electron temperature field.
      ! Sidre reads in any 'external' entries ( pointers ), but leaves the node as
      !   dtype=empty
      !   value=null
      !
      ! I couldn't find a Conduit Fortran API call that lets you update that entry.
      ! For now just delete the entry and re-add it. -- Aaron
      call current%remove_child("corner_electron_temperature")
      numElements= Size%ncornr
      Mat%tec = nan
      call current%set_path_external_float64_ptr("corner_electron_temperature", Mat%tec, numElements)

      TETON_VERIFY( current%fetch_path_as_int32("num_phase_space_sets") == Quad%nSets, "Teton chkpt load error: Num phase space sets does not match expected number.")
      
      SetLoop: do setID=1,Quad%nSets

         set => getSetData(Quad, SetID)
         write(setIDString,*) setID
         setIDString = adjustl(setIDString)

         current = current%child("phase_space_set_"//setIDString)

         TETON_VERIFY(current%fetch_path_as_int32("num_groups") == set%Groups, "Teton chkpt load error: groups does not match expected number.")
         TETON_VERIFY(current%fetch_path_as_int32("num_corners") == set%nCorner, "Teton chkpt load error: corners does not match expected number.")
         TETON_VERIFY(current%fetch_path_as_int32("num_angles") == set%NumAngles, "Teton chkt load error: angles does not match expected number.")

         numElements = set%Groups * set%nCorner * set%NumAngles
         ! Publish a pointer to PSI for each set.
         ! Same issue as for the corner electron temperatures above, for now just delete the entry and re-add it. -- Aaron
         call current%remove_child("radiation_energy_density")
         set%Psi = nan
         call current%set_path_external_float64_ptr("radiation_energy_density", set%Psi, numElements)

         current = current%parent()

      enddo SetLoop

      current = current%child("geometry")

      call current%remove_child("Volume")
      numElements = Size%ncornr
      Geom%Volume = nan
      call current%set_path_external_float64_ptr("Volume", Geom%Volume, numElements)

      current = current%parent()

      current = current%child("control_parameters")
      current = current%child("iteration")

      itCon => getIterationControl(IterControls,"temperature")
      relTol = current%fetch_path_as_float64("outerTempRelTol")
      maxIts = current%fetch_path_as_int32("outerMaxIt")
      call setControls(itCon, epsilonPoint=relTol, maxNumberOfIterations=maxIts)
      call current%remove_child("outerTempRelTol")
      call current%remove_child("outerMaxIt")

      itCon => getIterationControl(IterControls,"intensity")
      relTol = current%fetch_path_as_float64("outerPhiRelTol")
      call setControls(itCon, epsilonPoint=relTol)
      call current%remove_child("outerPhiRelTol")

      itCon => getIterationControl(IterControls,"grey")
      relTol = current%fetch_path_as_float64("greyRelTol")
      maxIts = current%fetch_path_as_int32("greyMaxIt")
      call setControls(itCon, epsilonPoint=relTol, maxNumberOfIterations=maxIts)
      call current%remove_child("greyRelTol")
      call current%remove_child("greyMaxIt")

      itCon => getIterationControl(IterControls,"incidentFlux")
      relTol = current%fetch_path_as_float64("incidentFluxRelTol")
      maxIts = current%fetch_path_as_int32("incidentFluxMaxIt")
      call setControls(itCon, epsilonPoint=relTol, maxNumberOfIterations=maxIts)
      call current%remove_child("incidentFluxRelTol")
      call current%remove_child("incidentFluxMaxIt")

      itCon => getIterationControl(IterControls,"nonLinear")
      relTol = current%fetch_path_as_float64("innerNLRelTol")
      maxIts = current%fetch_path_as_int32("innerNLMaxIt")
      call setControls(itCon, epsilonPoint=relTol, maxNumberOfIterations=maxIts)
      call current%remove_child("innerNLRelTol")
      call current%remove_child("innerNLMaxIt")

      current = current%parent()
      call current%remove_child("iteration")

      current = current%parent()
      call current%remove_child("control_parameters")

   end subroutine

!*****************************************************************************
! Notify teton that tree of external data has been updated.
!
! Teton doesn't really need to do anything here in the Fortran layer.
!
! In the TetonConduitInterface.cc wrapper around this function, it calls
! teton_rtedit(mCornerTemperature.data());
! in order to get the C++ copy of the corner temps updated.
! 
! Note from black27:
! - I think we can get rid of the extra C++ copy of the corner temps at some
! point.
!*****************************************************************************
   subroutine external_data_loaded() bind(c, name="teton_conduitcheckpoint_external_data_loaded")
   end subroutine

!*****************************************************************************
! Populate checkpoint node with data for save.
!*****************************************************************************
   subroutine prepForSave() bind(c, name="teton_conduitcheckpoint_prep_for_save")
      use, intrinsic :: iso_c_binding, only : c_size_t, c_associated
      use conduit_obj, only : conduit_node_obj_create
      use QuadratureList_mod
      use SetData_mod
      use Material_mod
      use Size_mod
      use MemoryAllocator_mod
      use Geometry_mod
      use iter_control_mod
      use iter_control_list_mod

      type(node) :: current
      type(SetData), pointer :: set
      integer :: setID
      character(len=64) :: setIDString
      integer(kind=c_size_t) :: numElements
      type(IterControl), pointer  :: itCon => NULL()

      if (.not. c_associated( root%cnode )) then
         root = conduit_node_obj_create()
      endif

      current = root

      call current%set_path("num_corners", Size%ncornr )
      call current%set_path("num_phase_space_sets", Quad%nSets )

      ! Add pointer to corner electron temperature field.
      numElements = Size%ncornr
      call current%set_path_external_float64_ptr("corner_electron_temperature", Mat%tec, numElements)

      ! Add pointer to PSI for each set.
      SetLoop: do setID=1,Quad%nSets

         set => getSetData(Quad, SetID)
         write(setIDString,*) setID
         setIDString = adjustl(setIDString)

         current = current%add_child("phase_space_set_"//setIDString)
         call current%set_path("num_groups", set%Groups)
         call current%set_path("num_corners", set%nCorner)
         call current%set_path("num_angles", set%NumAngles)

         numElements = set%Groups * set%nCorner * set%NumAngles
         call current%set_path_external_float64_ptr("radiation_energy_density", set%Psi, numElements)

         current = current%parent()

      enddo SetLoop

      current = current%add_child("geometry")
      numElements = Size%ncornr
      call current%set_path_external_float64_ptr("Volume", Geom%Volume, numElements)

      current = current%parent()

      current = current%add_child("control_parameters")
      current = current%add_child("iteration")
      itCon => getIterationControl(IterControls,"temperature")
      call current%set_path("outerTempRelTol", getEpsilonPoint(itCon))
      call current%set_path("outerMaxIt", getMaxNumberOfIterations(itCon))
      itCon => getIterationControl(IterControls,"intensity")
      call current%set_path("outerPhiRelTol", getEpsilonPoint(itCon))
      itCon => getIterationControl(IterControls,"grey")
      call current%set_path("greyRelTol", getEpsilonPoint(itCon))
      call current%set_path("greyMaxIt", getMaxNumberOfIterations(itCon))
      itCon => getIterationControl(IterControls,"incidentFlux")
      call current%set_path("incidentFluxRelTol", getEpsilonPoint(itCon))
      call current%set_path("incidentFluxMaxIt", getMaxNumberOfIterations(itCon))
      itCon => getIterationControl(IterControls,"nonLinear")
      call current%set_path("innerNLRelTol", getEpsilonPoint(itCon))
      call current%set_path("innerNLMaxIt", getMaxNumberOfIterations(itCon))
      current = current%parent()
      current = current%parent()


   end subroutine

!*****************************************************************************
! Delete the checkpoint node
! 
! Note from black27:
! Conduit doesn't have the 'reset' C function wrapped for Fortran currently.
! That function clears out a conduit tree.  I opted to just put this up in the
! TetonConduitInterface.cc layer.
!*****************************************************************************
   subroutine teardown() bind(c, name="teton_conduitcheckpoint_teardown")
      use conduit_obj, only : conduit_node_obj_create
      ! not implemented yet
   end subroutine

end module ConduitCheckpoint_mod
