!=======================================================================
! construct dynamic memory 
!=======================================================================
   subroutine constructDynMemory(setID, maxPerPlane)

   use kind_mod
   use Size_mod
   use constant_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod
   use Options_mod
   use MemoryAllocator_mod

   implicit none

!  Arguments 

   integer, intent(in)      :: setID
   integer, intent(in)      :: maxPerPlane

!  Local

   type(SetData),  pointer  :: Set
   type(AngleSet), pointer  :: ASet
   type(SweepSet), pointer  :: Swp

   integer                  :: totalCycles
   integer                  :: dom
   integer                  :: nHyperDomains
   integer                  :: sweepVersion

!  Allocate Memory 

   if (Size% ndim == 1) then
     return
   endif

   Set  => getSetData(Quad, setID)
   ASet => getAngleSetFromSetID(Quad, setID)

   totalCycles   = max(ASet% totalCycles, 1) 
   nHyperDomains = getNumberOfHyperDomains(Quad,1)
   sweepVersion  = Options% getSweepVersion()

   call Allocator%allocate(Size%usePinnedMemory, Set%label, "cyclePsi", Set% cyclePsi, Set% Groups, totalCycles)

!  These are onlu used in the GPU sweep
   if (Size% useGPU) then

     if ( sweepVersion == 1 ) then

       do dom=1,nHyperDomains
         Swp => Set% SweepPtr(dom)
         call Allocator%allocate(Size%usePinnedMemory, Set%label, "Q", Swp% Q, Set% Groups, Size% maxCorner, maxPerPlane)
         call Allocator%allocate(Size%usePinnedMemory, Set%label, "S", Swp% S, Set% Groups, Size% maxCorner, maxPerPlane)
       enddo

     do dom=1,nHyperDomains
       Swp => Set% SweepPtr(dom)
       Swp% Q(:,:,:) = zero
       Swp% S(:,:,:) = zero
     enddo

     endif

     call Allocator%allocate(Size%usePinnedMemory, Set%label, "PsiInt", Set% PsiInt, Set% Groups, ASet% maxInterface, Set% numAngles)

     Set% PsiInt(:,:,:) = zero
   endif


   return

   end subroutine constructDynMemory 


   subroutine initCyclePsi(setID)

   use kind_mod
   use Size_mod
   use constant_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod

   implicit none

!  Arguments 

   integer, intent(in)      :: setID

!  Local

   type(SetData),  pointer  :: Set
   type(AngleSet), pointer  :: ASet

   integer                  :: m
   integer                  :: mCycle
   integer                  :: c
   integer                  :: offSet
   integer                  :: angle

!  Allocate Memory 

   Set  => getSetData(Quad, setID)
   ASet => getAngleSetFromSetID(Quad, setID)

!  Initialize

   do angle=1,ASet% numAngles

     offSet    = ASet% cycleOffSet(angle)

     do m=1,ASet% numCycles(angle)

       mCycle  = offSet + m
       c       = ASet% cycleList(mCycle)

       Set% cyclePsi(:,mCycle) = Set% Psi(:,c,angle)
     enddo

   enddo


   return

   end subroutine initCyclePsi

!=======================================================================
! initialize using cycle Psi 
!=======================================================================
   subroutine initFromCycleList(setID, angle, Groups, Psi1)

   use kind_mod
   use Size_mod
   use constant_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod

   implicit none

!  Arguments 

   integer,    intent(in)    :: setID
   integer,    intent(in)    :: angle
   integer,    intent(in)    :: Groups

   real(adqt), intent(inout) :: Psi1(Groups,Size%ncornr+Size%nbelem)

!  Local

   type(SetData),  pointer  :: Set
   type(AngleSet), pointer  :: ASet

   integer                  :: m
   integer                  :: mCycle
   integer                  :: c
   integer                  :: offSet

!  Allocate Memory 

   Set    => getSetData(Quad, setID)
   ASet   => getAngleSetFromSetID(Quad, setID)

   offSet =  ASet% cycleOffSet(angle)

!  Initialize

   do m=1,ASet% numCycles(angle)

     mCycle = offSet + m
     c      = ASet% cycleList(mCycle)

     Psi1(:,c) = Set% cyclePsi(:,mCycle)

   enddo


   return

   end subroutine initFromCycleList 

!=======================================================================
! update cyclePsi 
!=======================================================================
   subroutine updateCycleList(setID, angle, Groups, Psi1)

   use kind_mod
   use Size_mod
   use constant_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod

   implicit none

!  Arguments 

   integer,    intent(in)   :: setID
   integer,    intent(in)   :: angle
   integer,    intent(in)   :: Groups

   real(adqt), intent(in)   :: Psi1(Groups,Size%ncornr+Size%nbelem)

!  Local

   type(SetData),  pointer  :: Set
   type(AngleSet), pointer  :: ASet

   integer                  :: m
   integer                  :: mCycle
   integer                  :: c
   integer                  :: offSet

!  Allocate Memory 

   Set    => getSetData(Quad, setID)
   ASet   => getAngleSetFromSetID(Quad, setID)

   offSet =  ASet% cycleOffSet(angle)

!  Initialize

   do m=1,ASet% numCycles(angle)
     mCycle = offSet + m
     c      = ASet% cycleList(mCycle)

     Set% cyclePsi(:,mCycle) = Psi1(:,c)
   enddo


   return

   end subroutine updateCycleList

