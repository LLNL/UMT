!***********************************************************************
!                       Last Update:  10/2016, PFN                     *
!                                                                      *
!   addGreyCorrections - Update scalar and angle-dependent intensities *
!                        with grey corrections.                        *
!                                                                      *
!***********************************************************************
   subroutine addGreyCorrections 

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use RadIntensity_mod
   use GreyAcceleration_mod
   use Communicator_mod
   use QuadratureList_mod
   use BoundaryList_mod
   use Boundary_mod
   use SetData_mod
   use AngleSet_mod
   use CommSet_mod

   implicit none

!  Local

   type(SetData),          pointer  :: Set
   type(AngleSet),         pointer  :: ASet
   type(CommSet),          pointer  :: CSet
   type(Communicator),     pointer  :: CommT
   type(Boundary),         pointer  :: BdyT

   integer    :: zone
   integer    :: nzones
   integer    :: b
   integer    :: b0
   integer    :: c
   integer    :: c0
   integer    :: nCorner
   integer    :: angle
   integer    :: exit_angle
   integer    :: NumAngles
   integer    :: g
   integer    :: g0
   integer    :: Groups
   integer    :: i
   integer    :: setID
   integer    :: nSets
   integer    :: sharedID
   integer    :: nShared
   integer    :: reflID
   integer    :: nReflecting
   integer    :: nBdyElem

   real(adqt) :: wtiso
   real(adqt) :: sumRad
   real(adqt) :: correction
                                                                                       
!  Mesh Constants

   nzones      =  Size% nzones 
   wtiso       =  Size% wtiso
   nShared     =  getNumberOfShared(RadBoundary)
   nReflecting =  getNumberOfReflecting(RadBoundary)
   nSets       =  getNumberOfSets(Quad)

!  Compute the group-dependent corrections

!$omp parallel do default(none) schedule(static) &
!$omp& shared(Size, Geom, Rad, GTA, nzones) &
!$omp& private(c0, nCorner, sumRad, correction) 

   ZoneLoop: do zone=1,nzones

     nCorner = Geom% numCorner(zone)
     c0      = Geom% cOffSet(zone)

     Rad% radEnergy(zone) = zero

     do c=1,nCorner

       sumRad = zero

       do g=1,Size% ngr
         correction            = GTA%GreyCorrection(c0+c)*GTA% Chi(g,c0+c)
         Rad% PhiTotal(g,c0+c) = Rad% PhiTotal(g,c0+c) + correction 
         sumRad                = sumRad + Rad% PhiTotal(g,c0+c)
       enddo

       Rad% radEnergy(zone) = Rad% radEnergy(zone) + Geom% Volume(c0+c)*sumRad

     enddo

   enddo ZoneLoop

!$omp end parallel do

!  Update Set dependent boundary fluxes

!$omp parallel do default(none) schedule(static) &
!$omp& private(NumAngles, g0, Groups, c, nBdyElem, b, b0, exit_angle) &
!$omp& private(Set, ASet, CSet, CommT, BdyT) &
!$omp& shared(nSets, Quad, GTA, RadBoundary, wtiso, nReflecting, nShared)

   SetLoop: do setID=1,nSets

     Set       => getSetData(Quad, setID)
     ASet      => getAngleSetFromSetID(Quad, setID)
     CSet      => getCommSetFromSetID(Quad, setID)

     NumAngles =  Set% NumAngles
     g0        =  Set% g0
     Groups    =  Set% Groups 

!    In spatially decomposed runs, also correct the boundary flux

     SharedLoop: do sharedID=1,nShared

       do angle=1,NumAngles
         CommT => getMessage(CSet, sharedID, angle)
         do i=1,CommT% nSend
           b =  CommT% ListSend(1,i)
           c =  CommT% ListSend(2,i)

           do g=1,Groups
             Set% PsiB(g,b,angle) = Set% PsiB(g,b,angle) + wtiso*  & 
                                    GTA%GreyCorrection(c)*GTA% Chi(g0+g,c)
           enddo
         enddo
       enddo

     enddo SharedLoop

!    Opposing Reflecting Boundaries

     ReflLoop: do reflID=1,nReflecting

       BdyT      => getReflecting(RadBoundary, reflID)
       nBdyElem  =  getNumberOfBdyElements(BdyT)
       b0        =  getFirstBdyElement(BdyT) - 1

       do i=1,ASet% nExit(reflID)
         exit_angle = ASet% ExitAngleList(i,reflID)
         do b=1,nBdyElem
           c = BdyT% BdyToC(b)

           do g=1,Groups
             Set% PsiB(g,b0+b,exit_angle) = Set% PsiB(g,b0+b,exit_angle) + wtiso*  &
                                            GTA%GreyCorrection(c)*GTA% Chi(g0+g,c) 
           enddo
         enddo
       enddo

     enddo ReflLoop

   enddo SetLoop

!$omp end parallel do

 
   return
   end subroutine addGreyCorrections 

