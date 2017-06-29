!***********************************************************************
!                        Version 1:  05/92, PFN                        *
!                                                                      *
!   RSWPMD - Controls transport sweeps for 2D and 3D grids             *
!                                                                      *
!   Input:                                                             *
!                                                                      *
!   Output:  PSIR  - corner angular photon intensity      (E/A/t/ster) *
!            PSIB  - angular photon intensity on boundary              *
!                    elements                             (E/A/t/ster) *
!                                                                      *
!***********************************************************************
   subroutine rswpmd(PSIB, PSIR, PHI, angleLoopTime, intensityIter, tempIter)
 
   use kind_mod
   use constant_mod
   use Size_mod
   use Material_mod
   use Geometry_mod
   use Quadrature_mod
   use ZoneData_mod
   use cudafor

   implicit none

!  Arguments

   real(adqt), intent(inout) :: psib(QuadSet%Groups,Size%nbelem,QuadSet%NumAngles),  &
                                psir(QuadSet%Groups,Size%ncornr,QuadSet%NumAngles),  &
                                Phi(QuadSet%Groups,Size%ncornr), angleLoopTime

   integer, intent(in) :: intensityIter, tempIter ! current flux and temperature iteration from rtmainsn

!  Local

   integer          :: c, c0, nCorner, zone 

   character(len=8) :: ipath

!  Constants

   parameter (ipath='sweep')

!  Add Scattering Source to total
   !$omp parallel do private(Z,nCorner,c0,c)
   do zone=1,Size%nzones
     Z => getZoneData(Geom, zone)
     nCorner   = Z% nCorner
     c0        = Z% c0
     do c=1,nCorner
       Z% STotal(:,c) = Mat%SFixed(:,c0+c) + Size%wtiso*Mat%sigs(:,zone)*Phi(:,c0+c)
     enddo
   enddo

!  Update device-memory copies of STotal, sigt, sigtinv

   call setZones_GPU_STotal(Geom%d_GPU_ZData, Geom%d_ZData)

!  Follow particles through the mesh:

   call snflwxyz(ipath, PSIB, PSIR, PHI, angleLoopTime, intensityIter, tempIter)



   return
   end subroutine rswpmd

