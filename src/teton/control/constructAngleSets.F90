!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   constructAngleSets - Creates batches of energy groups with the     *
!                        same quadrature set that can be transported   *
!                        together.                                     *
!                                                                      *
!***********************************************************************
 
   subroutine constructAngleSets(QuadDef, gnu) 

   use kind_mod
   use Size_mod
   use QuadratureList_mod

   implicit none

!  Arguments

   integer,    intent(in)       :: QuadDef(6,Size%ngr+1)
   real(adqt), intent(in)       :: gnu(Size%ngr+1)

!  Local

   integer    :: g, gnext, g1, g2, ngr, QuadID
   integer    :: quadType, Order, Type_set, Order_set, Groups
   integer    :: NumAngles, NumMoments, NPolar, NAzimuthal, PolarAxis
   real(adqt) :: GrpBnds(Size%ngr+1)

!  Find the number of angle sets and the number of groups in them 

   QuadID     = 1 
   Type_set   = QuadDef(1,1) 
   Order_set  = QuadDef(2,1) 
   Groups     = 0 
   NumMoments = 1 
   ngr        = Size%ngr
   g1         = 1

   GroupLoop: do g=1,ngr+1

     Groups = Groups + 1

     gnext = g + 1
     if (g < ngr) then
       quadType = QuadDef(1,gnext)
       Order    = QuadDef(2,gnext)
     else
!  This forces the last set to be completed
       quadType  = 0 
       Order = 0 
     endif

     if ( quadType == Type_set   .and.   &
             Order == Order_set ) then 

     else

!  Set the current quadrature definition

       NPolar     = QuadDef(3,g)
       NAzimuthal = QuadDef(4,g)
       PolarAxis  = QuadDef(5,g)
       NumAngles  = QuadDef(6,g)
       g2         = g + 1

       if (g < ngr+1) then
         GrpBnds(1:Groups+1) = Gnu(g1:g2)
       else
         GrpBnds(1) = Gnu(1) 
         GrpBnds(2) = Gnu(ngr+1)
       endif

       call setQuadrature(Quad,          &
                          QuadID,        &
                          Groups,        &
                          NumAngles,     &
                          NumMoments,    &
                          Order_set,     &
                          NPolar,        &
                          NAzimuthal,    &
                          PolarAxis,     &
                          Type_set,      &
                          GrpBnds)

!  Start a new batch

       if (g < ngr+1) then
         QuadID    = QuadID + 1
         Groups    = 0 
         Type_set  = QuadDef(1,gnext) 
         Order_set = QuadDef(2,gnext) 
         g1        = g + 1
       endif

     endif

   enddo GroupLoop


   return
   end subroutine constructAngleSets 

