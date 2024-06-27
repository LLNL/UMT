!***********************************************************************
!                        Last Update:  01/2016, PFN                    *
!                                                                      *
!   DestructMeshData  - Removes all modules that are mesh dependent.   *
!                       This is useful if the mesh is changed during   *
!                       a simulation. Note that the memory within      *
!                       each module is released and then the module    *
!                       itself is deallocated. This allows the module  *
!                       constructors to be reused for the new mesh.    *
!                                                                      *
!***********************************************************************


   subroutine DestructMeshData(nonLTE) BIND(C,NAME="teton_destructmeshdata")

   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use Geometry_mod
   use Material_mod
   use RadIntensity_mod
   use Quadrature_mod
   use QuadratureList_mod
   use Boundary_mod
   use BoundaryList_mod
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   use ComptonControl_mod
#endif
   use Editor_mod
   use TimeStepControls_mod
   use iter_control_list_mod
   use GreyAcceleration_mod
   use AngleSet_mod
   use GroupSet_mod
   use CommSet_mod
   use SetData_mod
   use ZoneSet_mod
   use DataStore_mod


   implicit none

!  Arguments

   logical(C_BOOL), intent(in) :: nonLTE

!  Local

   type(SetData),    pointer :: Set
   type(AngleSet),   pointer :: ASet
   type(GroupSet),   pointer :: GSet
   type(CommSet),    pointer :: CSet

   integer         :: setID
   integer         :: sharedID
   integer         :: nSets
   integer         :: nGTASets
   integer         :: nAngleSets
   integer         :: nGroupSets
   integer         :: nCommSets
   integer         :: nShared

   logical(kind=1) :: GTASet

!  Constants

   nSets      = getNumberOfSets(Quad)
   nAngleSets = getNumberOfAngleSets(Quad)
   nGroupSets = getNumberOfGroupSets(Quad)
   nCommSets  = getNumberOfCommSets(Quad)
   nGTASets   = getNumberOfGTASets(Quad)
   nShared    = getNumberOfShared(RadBoundary)

!  Deallocate Geometry and Radiation modules 

   call destruct(Geom)
   call destruct(Rad)
   deallocate(Rad)
   call Mat%destruct(nonLTE)

!  Deallocate Phase-Spaces Set data

   GTASet = .FALSE.

   do setID=1,nSets
     Set => getSetData(Quad, setID)

     call Set%destruct(GTASet)
   enddo

   GTASet = .TRUE.

   do setID=nSets+1,nSets+nGTASets
     Set => getSetData(Quad, setID)

     call Set%destruct(GTASet)
   enddo

!  Angle Sets

   if (Size% ndim > 1) then

     do setID=1,nAngleSets+nGTASets
       ASet => getAngleSetData(Quad, setID)

         do sharedID=1,nShared
           call destructIncidentTest(ASet, sharedID)
         enddo

     enddo

   endif

   do setID=1,nAngleSets+nGTASets
     ASet => getAngleSetData(Quad, setID)

     call destruct(ASet)
   enddo

!  Group Sets

   do setID=1,nGroupSets
     GSet => getGroupSetData(Quad, setID)

     call destruct(GSet)
   enddo

!  Zone Sets

   call destruct(ZSet)
   deallocate(ZSet)

!  Communication Sets

   do setID=1,nCommSets+nGTASets
     CSet  => getCommSetData(Quad, setID)

     call destruct(CSet)
   enddo

!  Deallocate Quadrature data (One for Sn, One for GTA)

   do setID=1,2
     QuadSet => getQuadrature(Quad, setID)

     call destruct(QuadSet)
   enddo

!  Deallocate Boundary, Compton, Edit and Iteration control modules

   call destruct(RadBoundary)
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   call destruct(Compton)
#endif
   call destruct(RadEdit)
   call destruct(IterControls)
   call destruct(GTA)
   call destruct(Quad)

!  Deallocate parent modules

   deallocate( Geom )
   deallocate( Mat )
   deallocate( Quad )
   deallocate( RadBoundary )
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   deallocate( Compton )
#endif
   deallocate( DtControls )
   deallocate( RadEdit )
   deallocate( IterControls )
   deallocate( GTA )
   deallocate( Size )


   call theDatastore%root%reset()

   return
   end subroutine DestructMeshData 

