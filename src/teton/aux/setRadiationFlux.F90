!***********************************************************************
!                       Last Update:  03/2012, PFN                     *
!                                                                      *
!   setRadiationFlux    - Calculates the zone-average radiation        *
!                         flux vector.                                 * 
!                                                                      *
!***********************************************************************
 
   subroutine setRadiationFlux() BIND(C,NAME="teton_setradiationflux")

   use kind_mod
   use QuadratureList_mod

   implicit none 

!  Local

   integer         :: setID
   integer         :: nSets
   logical(kind=1) :: Force

!  Constants

!***********************************************************************
!  Compute the radiation flux                                          *
!***********************************************************************

   nSets = getNumberOfSets(Quad)
   Force = .FALSE.

!$omp parallel do private(setID) schedule(dynamic)

   SetLoop: do setID=1,nSets

     call setRadiationMoments(setID, Force)

   enddo SetLoop

!$omp end parallel do


   return
   end subroutine setRadiationFlux 


