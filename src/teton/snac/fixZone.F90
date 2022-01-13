!***********************************************************************
!                       Last Update:  12/2017, PFN                     *
!                                                                      *
!   FIXZONE      - This routine fixes a cycle within a zone.           *
!                                                                      *
!***********************************************************************
   subroutine fixZone(zone, meshCycles, cycleList) 

   use kind_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Arguments

   integer,    intent(in)    :: zone
   integer,    intent(inout) :: meshCycles

   integer,    intent(inout) :: cycleList(Size%ncornr)

!  Local Variables

   integer :: c

!  Add this zone to the cycle list

   do c=1,Geom% numCorner(zone)
     meshCycles            = meshCycles + 1
     cycleList(meshCycles) = Geom% cOffSet(zone) + c
   enddo


   return
   end subroutine fixZone 

