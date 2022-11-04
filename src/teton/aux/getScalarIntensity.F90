!***********************************************************************
!                         Last Update: 10/2016 PFN                     *
!                                                                      *
!    getScalarIntensity -     Called from host to retrieve the         *
!                             scalar radiation intensity (Phi) from    *
!                             the SetData module.                      *
!                                                                      *
!***********************************************************************

   subroutine getScalarIntensity(PhiTotal) &
     BIND(C,NAME="teton_getscalarintensity")


   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use constant_mod
   use Geometry_mod
   use SetData_mod
   use AngleSet_mod
   use QuadratureList_mod

   implicit none

!  Arguments

   real(C_DOUBLE), intent(inout) :: PhiTotal(Size%ngr,Size%ncornr)

!  Local

   type(SetData),    pointer :: Set
   type(AngleSet),   pointer :: ASet

   integer    :: angle
   integer    :: numAngles
   integer    :: c
   integer    :: g
   integer    :: g0
   integer    :: Groups
   integer    :: setID
   integer    :: nSets
   integer    :: zSetID
   integer    :: nZoneSets

   real(adqt) :: quadwt

!  Constants

   nSets     = getNumberOfSets(Quad)
   nZoneSets = getNumberOfZoneSets(Quad)


   PhiTotal(:,:) = zero

!$omp parallel do default(none) schedule(static) &
!$omp& private(Set, ASet, g0, Groups, NumAngles, quadwt) &
!$omp& shared(nZoneSets, nSets, Quad, Geom, PhiTotal)

   ZoneSetLoop: do zSetID=1,nZoneSets

     do setID=1,nSets

       Set       => getSetData(Quad, setID)
       ASet      => getAngleSetFromSetID(Quad, setID)

       NumAngles =  Set% NumAngles
       Groups    =  Set% Groups
       g0        =  Set% g0

       AngleLoop: do Angle=1,NumAngles
         quadwt =  ASet% weight(Angle)

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           do g=1,Groups
             PhiTotal(g0+g,c) = PhiTotal(g0+g,c) + &
                                quadwt*Set% Psi(g,c,Angle)
           enddo
         enddo
       enddo AngleLoop

     enddo

   enddo ZoneSetLoop


   return
   end subroutine getScalarIntensity 

!***********************************************************************
!                         Last Update: 07/2019 PFN                     *
!                                                                      *
!    setScalarIntensity -  This is the same function as above, but     *
!                          sets the internal scalar radiation          *
!                          intensity (Rad% PhiTotal).                  *
!                                                                      *
!***********************************************************************

   subroutine setScalarIntensity


   use kind_mod
   use Size_mod
   use constant_mod
   use Geometry_mod
   use RadIntensity_mod
   use SetData_mod
   use AngleSet_mod
   use QuadratureList_mod

   implicit none


!  Local

   type(SetData),    pointer :: Set
   type(AngleSet),   pointer :: ASet

   integer    :: angle
   integer    :: numAngles
   integer    :: c
   integer    :: g
   integer    :: g0
   integer    :: Groups
   integer    :: setID
   integer    :: nSets
   integer    :: zSetID
   integer    :: nZoneSets

   real(adqt) :: quadwt

!  Constants

   nSets     = getNumberOfSets(Quad)
   nZoneSets = getNumberOfZoneSets(Quad)


   Rad% PhiTotal(:,:) = zero

!$omp parallel do default(none) schedule(static) &
!$omp& private(Set, ASet, g0, Groups, NumAngles, quadwt) &
!$omp& shared(nZoneSets, nSets, Quad, Geom, Rad)

   ZoneSetLoop: do zSetID=1,nZoneSets

     do setID=1,nSets

       Set       => getSetData(Quad, setID)
       ASet      => getAngleSetFromSetID(Quad, setID)

       NumAngles =  Set% NumAngles
       Groups    =  Set% Groups
       g0        =  Set% g0

       AngleLoop: do Angle=1,NumAngles
         quadwt =  ASet% weight(Angle)

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           do g=1,Groups
             Rad% PhiTotal(g0+g,c) = Rad% PhiTotal(g0+g,c) + &
                                     quadwt*Set% Psi(g,c,Angle)
           enddo
         enddo
       enddo AngleLoop

     enddo

   enddo ZoneSetLoop


   return
   end subroutine setScalarIntensity

