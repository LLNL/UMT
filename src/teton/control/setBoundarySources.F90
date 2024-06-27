#include "macros.h"
!***********************************************************************
!                        Last Update:  01/2012, PFN                    *
!                                                                      *
!   RTBDRY - Calculates the radiation intensity incident on the        *
!            problem boundary.  Three types of boundary conditions     *
!            are currently allowed:  vacuum, fds and milne.  A vacuum  *
!            condition yields a zero incoming intensity.  For FDS, the *
!            frequency-dependent intensity is specified by the user    *
!            and is assumed to be isotropic.  For a MILNE condition,   *
!            the temperature is specified and the intensity is a       *
!            Planckian in frequency space and is isotropic.            *
!                                                                      *
!***********************************************************************

   subroutine setBoundarySources(setID) 

   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use Quadrature_mod
   use AngleSet_mod
   use BoundaryList_mod
   use Boundary_mod
   use Profile_mod
   use SetData_mod

   implicit none

!  Arguments

   integer, intent(in)  :: setID

!  Local

   type(SetData),  pointer  :: Set
   type(Boundary), pointer  :: BdyT 
   type(Profile),  pointer  :: ProfID

   integer         :: prof
   integer         :: nSource
   integer         :: nBdyElem
   integer         :: b
   integer         :: b1 
   integer         :: b2
   integer         :: g
   integer         :: g0
   integer         :: angle
   integer         :: NumAngles
   integer         :: Groups

   real(adqt)      :: wtiso
   real(adqt)      :: wtbeam

   logical(kind=1) :: profileOn

!  Some temporary variables used for beam boundary conditions:
   type(AngleSet), pointer :: ASet
   real(adqt)              :: polarAngleValue
   real(adqt)              :: refPolarValue

!  Constants

   Set       => getSetData(Quad, setID)
   ASet      => getAngleSetFromSetID(Quad, setID)

   g0        =  Set% g0
   Groups    =  Set% Groups
   NumAngles =  Set% NumAngles 
   wtiso     =  Size% wtiso
   nSource   =  getNumberOfSource(RadBoundary) 

!  Initialize Boundary Flux

   Set% PsiB(:,:,:) = zero

   ProfileLoop: do prof=1,nSource 

     BdyT   => getSource(RadBoundary, prof)
     ProfID => getProfile(RadBoundary, prof) 
 
     profileOn = ProfID% profileOn 

     CheckStatus: if ( profileOn ) then

       nBdyElem = getNumberOfBdyElements(BdyT)
       b1       = getFirstBdyElement(BdyT)
       b2       = b1 + nBdyElem - 1

!***********************
!  Temperature or FDS  *
!***********************

!      PsiB is per steradian

       if ( ProfID%Isotropic ) then

          do angle=1,NumAngles
            do b=b1,b2
              do g=1,Groups
                Set% PsiB(g,b,angle) = wtiso*ProfID% InterpValues(g0+g)
              enddo
            enddo
          enddo

       else ! Beam along polar axis

          QuadSet => getQuadrature(Quad, Set%quadID)
          TETON_VERIFY(QuadSet%TypeName == 'lobatto', "Beam boundary conditions only work with a Lobatto quadrature set.")

          ! Check whether the beam should be in the +paxis or -paxis direction
          ! Assume all boundary elements are facing the same direction,
          !    so we should just check the first one.
          if ( BdyT%A_bdy(QuadSet%PolarAxis,1) > zero ) then
            ! Outward normal is in +paxis direction,
            !   so we want beam in -paxis direction
            refPolarValue = -one
          else
            ! Outward normal is in -paxis direction,
            !   so we want beam in +paxis direction
            refPolarValue = one
          endif

          do angle=1,NumAngles ! Is it only ever the last two angles?

            polarAngleValue = Set%PolarAngleList(Set%PolarAngle(angle))

            if (abs(polarAngleValue-refPolarValue) > adqtEpsilon ) then
               ! Skip this angle if it's not along the polar axis
               cycle
            else

               wtbeam = 1./ASet%Weight(angle)
               do b=b1,b2
                 do g=1,Groups
                   Set% PsiB(g,b,angle) = wtbeam*ProfID% InterpValues(g0+g)
                 enddo
               enddo
               exit

            endif ! angle along polar axis

          enddo ! loop over angles in set

       endif ! if isotropic

     endif CheckStatus
                                                                                         
   enddo ProfileLoop

 
   return
   end subroutine setBoundarySources 

