!***********************************************************************
!                        Version 1:  09/96, PFN                        *
!                                                                      *
!   SNREFLECT - This routine, called by SNFLWXYZ and SNFLWRZA,         *
!               computes the boundary flux (PSIB) for angle m, on      *
!               reflecting boundaries for which angle m is incoming.   *
!                                                                      *
!   Input:                                                             *
!                                                                      *
!   Output:                                                            *
!                                                                      *
!***********************************************************************
  module snreflect

  contains
    
    subroutine snreflectD(anglebatch, Angles, d_psib, pinned_psib, streamid) 

      use kind_mod
      use constant_mod
      use Size_mod
      use Quadrature_mod
      use BoundaryList_mod
      use Boundary_mod
      use cudafor
      
      implicit none
      
      !  Arguments
      
      integer,  intent(in)    :: anglebatch!, groups, nbelem
      
      integer, device, intent(in)    :: Angles(anglebatch)
      
      real(adqt), device, intent(out) :: d_psib(QuadSet%Groups,Size%nbelem, anglebatch)
      real(adqt), device, intent(in) :: pinned_psib(QuadSet%Groups,Size%nbelem, QuadSet%NumAngles) ! zerocopy

!      integer, intent(in) :: mm1

      integer(kind=cuda_stream_kind), intent(in) :: streamid

      
      !  Local Variables

      integer, device, pointer :: d_ReflAngle(:)


      integer    :: i, ib, Mref, Minc, mm, ig, m0
      
      integer    :: nBdyElem, b1, b2, Groups
      

      integer      :: nReflecting, set ! num reflecting boundaries, set=QuadSet%QuadID
      
      nReflecting = getNumberOfReflecting(RadBoundary)
      
      set  =  QuadSet% QuadID
      
      Groups = QuadSet% Groups
      

      !allocate( d_ReflAngle(QuadSet%NumAngles) )


      !  Loop over reflecting-boundary sets:
      
      ReflectingLoop: do i=1,nReflecting

         ! Bdy       => getReflecting(RadBoundary, i)
         ! nBdyElem  =  getNumberOfBdyElements(Bdy)
         ! b1        =  getFirstBdyElement(Bdy)
         ! b2        =  b1 + nBdyElem - 1

         ! Mref      =  getReflectedAngle(Bdy, Minc)

         Bdy => RadBoundary%iBoundary(i)
         nBdyElem = Bdy%NumBdyElem
         b1 = Bdy%BdyElem1
         b2        =  b1 + nBdyElem - 1

         d_ReflAngle => RadBoundary%iBoundary(i)%iRef(set)% d_ReflAngle

         ! loop over angles in anglebatch: (THIS CAN ACTUALLY BE DONE IN PARALLEL. IN FACT SHOULD DO EVEN ON HOST)
         !$cuf kernel do(3) <<< (*,*), (*,*), stream=streamid >>>            
         do mm=1, anglebatch
            do ib=b1,b2
               do ig=1,Groups         
                  Minc = Angles(mm)
                  
                  Mref = d_ReflAngle(Minc)

                  if (Mref > 0) then

                     d_psib(ig,ib,mm) = pinned_psib(ig,ib,Mref)
                  endif
               enddo
            enddo
         enddo 

      enddo ReflectingLoop


      return
    end subroutine snreflectD


    ! attributes(global) subroutine snreflectD(anglebatch, Angles, PSIB, nReflecting, set) 

    !   use kind_mod
    !   use constant_mod
    !   use Size_mod
    !   use Quadrature_mod
    !   use BoundaryList_mod
    !   use Boundary_mod
      
    !   implicit none
      
    !   !  Arguments
      
    !   integer, value,  intent(in)    :: anglebatch, groups, ncornr, nbelem
      
    !   integer, device, intent(in)    :: Angles(anglebatch)
      
    !   real(adqt),device, intent(inout) :: psib(QuadSet%Groups,Size%nbelem, anglebatch) 

    !   integer, value, intent(in)     :: nReflecting, set ! num reflecting boundaries, set=QuadSet%QuadID
      
    !   !  Local Variables

    !   integer    :: i, ib, Mref, Minc
      
    !   integer    :: nBdyElem, b1, b2
      

      
      
    !   ! loop over angles in anglebatch:
    !   do mm=blockIdx%x, anglebatch, gridDim%x
    !      Minc = Angles(mm)

    !      !  Loop over reflecting-boundary sets:
         
    !      ReflectingLoop: do i=1,nReflecting

    !         ! Bdy       => getReflecting(RadBoundary, i)
    !         ! nBdyElem  =  getNumberOfBdyElements(Bdy)
    !         ! b1        =  getFirstBdyElement(Bdy)
    !         ! b2        =  b1 + nBdyElem - 1

    !         ! Mref      =  getReflectedAngle(Bdy, Minc)

    !         d_Bdy => RadBoundary%d_iBoundary(i)
    !         nBdyElem = d_Bdy%NumBdyElem
    !         b1 = d_Bdy%BdyElem1
    !         b2        =  b1 + nBdyElem - 1

    !         Mref = d_Bdy%iRef(set)% d_ReflAngle(Minc)

    !         if (Mref > 0) then

    !            do ib=b1,b2
    !               do ig=1,Groups
    !                  psib(ig,ib,Minc) = psib(ig,ib,Mref)
    !               enddo
    !            enddo

    !         endif

    !      enddo ReflectingLoop

    !   enddo


    !   return
    ! end subroutine snreflectD




   subroutine snreflect(Minc, PSIB) 

   use kind_mod
   use constant_mod
   use Size_mod
   use Quadrature_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none

!  Arguments

   integer,    intent(in)    :: Minc

   real(adqt), intent(inout) :: psib(QuadSet%Groups,Size%nbelem,  &
                                     QuadSet%NumAngles) 

!  Local Variables

   integer    :: i, ib, Mref

   integer    :: nReflecting, nBdyElem, b1, b2

!  Constants

   nReflecting = getNumberOfReflecting(RadBoundary)

!  Loop over reflecting-boundary sets:
 
   ReflectingLoop: do i=1,nReflecting

     Bdy       => getReflecting(RadBoundary, i)
     nBdyElem  =  getNumberOfBdyElements(Bdy)
     b1        =  getFirstBdyElement(Bdy)
     b2        =  b1 + nBdyElem - 1

     Mref      =  getReflectedAngle(Bdy, Minc)

     if (Mref > 0) then

       do ib=b1,b2
         psib(:,ib,Minc) = psib(:,ib,Mref)
       enddo

     endif
 
   enddo ReflectingLoop 

 
   return
   end subroutine snreflect


 end module snreflect
