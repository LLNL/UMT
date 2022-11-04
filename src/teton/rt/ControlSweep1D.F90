!***********************************************************************
!                      Last Update:  10/2016, PFN                      *
!                                                                      *
!   ControlSweep1D - Controls transport sweeps for 1D spatial grids.   *
!                                                                      *
!***********************************************************************
   subroutine ControlSweep1D(savePsi) 

   use flags_mod
   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use mpi_param_mod
   use mpif90_mod

   implicit none

!  Arguments

   logical(kind=1), intent(in) :: savePsi

!  Local

   integer    :: Groups
   integer    :: NumAngles
   integer    :: setID
   integer    :: nSets

   real(adqt) :: time1, time2, dtime


   nSets = getNumberOfSets(Quad)

   time1 = MPIWtime()

   if (Size%igeom == geometry_sphere .or. Size%igeom == geometry_cylinder) then

!  Transport sweeps in curvilinear geometry

     if ( savePsi ) then
!$omp parallel do default(none) schedule(static) &
!$omp& shared(nSets, Quad) &
!$omp& private(Groups, NumAngles)

       do setID=1,nSets

!        Post receives for all data
         call InitExchange(setID)

         Groups    = getNumberOfGroups(Quad, setID)
         NumAngles = getNumberOfAngles(Quad, setID)

         call SweepSphereSCBsv(setID, Groups, NumAngles)

       enddo

!$omp end parallel do
     else

!$omp parallel do default(none) schedule(static) &
!$omp& shared(nSets, Quad) &
!$omp& private(Groups, NumAngles)

       do setID=1,nSets

!        Post receives for all data
         call InitExchange(setID)

         Groups    = getNumberOfGroups(Quad, setID)
         NumAngles = getNumberOfAngles(Quad, setID)

         call SweepSphereSCB(setID, Groups, NumAngles)

       enddo

!$omp end parallel do
     endif

   elseif (Size%igeom == geometry_slab) then

!  Transport sweeps in planar geometry

!$omp parallel do default(none) schedule(static) &
!$omp& shared(nSets, Quad, savePsi) &
!$omp& private(Groups, NumAngles)

     do setID=1,nSets

!      Post receives for all data
       call InitExchange(setID)

       Groups    = getNumberOfGroups(Quad, setID)
       NumAngles = getNumberOfAngles(Quad, setID)

       call SweepPlanarSCB(setID, Groups, NumAngles, savePsi)

     enddo

!$omp end parallel do

   endif

   call getPhiTotal

   time2 = MPIWtime()
   dtime = (time2 - time1)/sixty
   Size%SweepTimeCycle = Size%SweepTimeCycle + dtime


   return
   end subroutine ControlSweep1D 

