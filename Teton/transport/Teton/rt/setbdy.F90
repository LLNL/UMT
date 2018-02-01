!***********************************************************************
!                        Version 1:  12/94, PFN                        *
!                                                                      *
!   SETBDY - Sets boundary fluxes PSIB to specified incident boundary  *
!            fluxes.                                                   *
!                                                                      *
!   Output:  PSIB  - boundary intensities                      (E/A/t) *
!                                                                      *
!***********************************************************************

   subroutine setbdyD(anglebatch, Angles, psir,  PSIB, streamid) 

     ! sets the boundary, a batch of angles at a time.

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use ProfileList_mod
   use BoundaryList_mod
   use Boundary_mod
   
   use cudafor

   implicit none

!  Arguments

   integer, intent(in) :: anglebatch ! number of angles to process (usually number of angles in a bin)

   integer, device, intent(in)    :: Angles(anglebatch)

   real(adqt), device, intent(in)    :: psir(Size%ngr,Size%ncornr,anglebatch)

   real(adqt), device, intent(inout) :: psib(Size%ngr,Size%nbelem,anglebatch)

   integer(kind=cuda_stream_kind), intent(in) :: streamid



!  Local

   integer :: i,ia,ib,ic,nangSN, mm, ig
   integer :: nReflecting, nVacuum, nSource, nShared 
   integer :: nBdyElem, b0, b1, b2, profID

   integer :: Groups

   integer, device, pointer :: BdyToC(:)
!   integer :: tid, nth, iabeg, iaend
!   integer, external :: omp_get_thread_num, omp_get_num_threads

!  Constants

   Groups = Size%ngr

   !nangSN      = Size%nangSN
   nSource     = getNumberOfSource(RadBoundary)
   nReflecting = getNumberOfReflecting(RadBoundary)
   nVacuum     = getNumberOfVacuum(RadBoundary)
   nShared     = getNumberOfShared(RadBoundary)



!  Reflecting
   
   do i=1,nReflecting
      Bdy      => getReflecting(RadBoundary, i)
      nBdyElem =  getNumberOfBdyElements(Bdy)
      b0       =  getFirstBdyElement(Bdy) - 1

      ! might need this:
      BdyToC => Bdy%d_BdyToC

      !$cuf kernel do(3) <<< (*,*), (*,*), stream=streamid >>>            
      do mm=1,anglebatch
         do ib=1,nBdyElem
            do ig=1,Groups
               ic = BdyToC(ib)
               psib(ig,b0+ib,mm) = psir(ig,ic,mm)
            enddo
         enddo
      enddo

   enddo

!  Vacuum

   do i=1,nVacuum
     Bdy      => getVacuum(RadBoundary, i)
     nBdyElem =  getNumberOfBdyElements(Bdy)
     b1       =  getFirstBdyElement(Bdy)
     b2       =  b1 + nBdyElem - 1
     !$cuf kernel do(3) <<< (*,*), (*,*), stream=streamid >>>
     do mm=1,anglebatch
        do ib=b1,b2
           do ig=1,Groups
              psib(ig,ib,mm) = zero
           enddo
        enddo
     enddo
   enddo

!  Sources
                                                                                                    
   do i=1,nSource
     Bdy => getSource(RadBoundary, i)
     nBdyElem = getNumberOfBdyElements(Bdy)
     b1       = getFirstBdyElement(Bdy)
     b2       = b1 + nBdyElem - 1
     profID   = getProfileID(Bdy)
     do mm=1,anglebatch
        ia=angles(mm)
        do ib=b1,b2
           psib(:,ib,mm) = SourceProfiles% Psi_Inc(:,ia,profID)
        enddo
     enddo
   enddo

!  Shared
                                                                                                    
   do i=1,nShared
     Bdy      => getShared(RadBoundary, i)
     nBdyElem =  getNumberOfBdyElements(Bdy)
     b0       =  getFirstBdyElement(Bdy) - 1

      BdyToC => Bdy%d_BdyToC

      !$cuf kernel do(3) <<< (*,*), (*,*), stream=streamid >>>            
      do mm=1,anglebatch
         do ib=1,nBdyElem
            do ig=1,Groups
               ic = BdyToC(ib)
               psib(ig,b0+ib,mm) = psir(ig,ic,mm)
            enddo
         enddo
      enddo

   enddo


   return
 end subroutine setbdyD




   subroutine setbdy(psir, PSIB) 

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use ProfileList_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none

!  Arguments

   real(adqt), intent(in)    :: psir(Size%ngr,Size%ncornr,Size%nangSN) 

   real(adqt), intent(inout) :: psib(Size%ngr,Size%nbelem,Size%nangSN)

!  Local

   integer :: i,ia,ib,ic,nangSN
   integer :: nReflecting, nVacuum, nSource, nShared 
   integer :: nBdyElem, b0, b1, b2, profID

   integer :: tid, nth, iabeg, iaend
   integer, external :: omp_get_thread_num, omp_get_num_threads

!  Constants

   nangSN      = Size%nangSN
   nSource     = getNumberOfSource(RadBoundary)
   nReflecting = getNumberOfReflecting(RadBoundary)
   nVacuum     = getNumberOfVacuum(RadBoundary)
   nShared     = getNumberOfShared(RadBoundary)

!$omp parallel private(i,Bdy,nBdyElem,b0,ib,ic,ia)

!  Reflecting

   do i=1,nReflecting
     Bdy      => getReflecting(RadBoundary, i)
     nBdyElem =  getNumberOfBdyElements(Bdy)
     b0       =  getFirstBdyElement(Bdy) - 1
!$omp do
     do ib=1,nBdyElem
       ic = Bdy% BdyToC(ib)
       do ia=1,nangSN
         psib(:,b0+ib,ia) = psir(:,ic,ia)
       enddo
     enddo
   enddo

!  Vacuum

   do i=1,nVacuum
     Bdy      => getVacuum(RadBoundary, i)
     nBdyElem =  getNumberOfBdyElements(Bdy)
     b1       =  getFirstBdyElement(Bdy)
     b2       =  b1 + nBdyElem - 1
!$omp do
     do ib=b1,b2
       do ia=1,nangSN
          psib(:,ib,ia) = zero
       enddo
     enddo
   enddo

!  Sources
                                                                                                    
   do i=1,nSource
     Bdy => getSource(RadBoundary, i)
     nBdyElem = getNumberOfBdyElements(Bdy)
     b1       = getFirstBdyElement(Bdy)
     b2       = b1 + nBdyElem - 1
     profID   = getProfileID(Bdy)
     do ia=1,nangSN
!$omp do
       do ib=b1,b2
         psib(:,ib,ia) = SourceProfiles% Psi_Inc(:,ia,profID)
       enddo
     enddo
   enddo

!  Shared
                                                                                                    
   do i=1,nShared
     Bdy      => getShared(RadBoundary, i)
     nBdyElem =  getNumberOfBdyElements(Bdy)
     b0       =  getFirstBdyElement(Bdy) - 1
!$omp do
     do ib=1,nBdyElem
       ic = Bdy% BdyToC(ib)
       do ia=1, nangSN  
         psib(:,b0+ib,ia) = psir(:,ic,ia)
       enddo
     enddo
   enddo
!$omp end parallel

   return
   end subroutine setbdy

