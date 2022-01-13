!***********************************************************************
!                        Last Update:  03/2012, PFN                    *
!                                                                      *
!   setRadiationForce   - Calculates the radiation force.              *
!                                                                      *
!***********************************************************************
 
   subroutine setRadiationForce() BIND(C,NAME="teton_setradiationforce")

   USE ISO_C_BINDING
   use kind_mod
   use QuadratureList_mod

   implicit none 

!  Local

   integer         :: setID
   integer         :: nSets
   logical(kind=1) :: Force
 
!***********************************************************************
!  Compute the radiation force on the matter                           *
!***********************************************************************

   nSets = getNumberOfSets(Quad)
   Force = .TRUE.

!$omp parallel do private(setID) schedule(dynamic)

   SetLoop: do setID=1,nSets

     call setRadiationMoments(setID, Force)

   enddo SetLoop

!$omp end parallel do


   return
   end subroutine setRadiationForce 


