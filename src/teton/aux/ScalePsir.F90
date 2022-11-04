!***********************************************************************
!                       Last Update:  03/2012, MAL                     *
!                                                                      *
!   ScalePsir - Initially scaled the intensity, given a relative       *
!               change in radiation energy density (isotropic).        *
!               Now just adds the isotropic change to the intensity.   *
!                                                                      *
!   Input:      Phi      - scalar flux             (E/A/t)             *
!               deltaPhi - change in scalar flux   (E/A/t)             *
!                                                                      *
!   Output:     psi  - radiation intensity         (E/A/t/sterad)      *
!                                                                      *
!***********************************************************************
 
   subroutine ScalePsir(deltaPhi) BIND(C,NAME="teton_scalepsir")

   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use QuadratureList_mod
   use SetData_mod

   implicit none 

!  Arguments

   real(C_DOUBLE), intent(in)    :: deltaPhi(Size%ngr, Size%ncornr)

!  Local

   type(SetData), pointer  :: Set

   integer    :: angle
   integer    :: c
   integer    :: g 
   integer    :: g0
   integer    :: nCorner
   integer    :: Groups
   integer    :: NumAngles
   integer    :: setID
   integer    :: nSets

!  Constants

   nSets   = getNumberOfSets(Quad)
   nCorner = Size% ncornr

   SetLoop: do setID=1,nSets

     Set       => getSetData(Quad, setID)

     Groups    =  Set% Groups
     NumAngles =  Set% NumAngles
     g0        =  Set% g0    

     do angle = 1,NumAngles
       do c = 1,nCorner 
         do g=1,Groups
            Set% Psi(g,c,angle) = Set% Psi(g,c,angle) +   &
                                  Size%wtiso*deltaPhi(g0+g,c)
         enddo
       enddo
     enddo

   enddo SetLoop


   return
   end subroutine ScalePsir
