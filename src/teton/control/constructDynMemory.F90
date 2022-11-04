!=======================================================================
! construct dynamic memory 
!=======================================================================
   subroutine constructDynMemory(setID, maxZonesPerPlane)

   use kind_mod
   use Size_mod
   use constant_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod
   use MemoryAllocator_mod

   implicit none

!  Arguments 

   integer, intent(in)      :: setID
   integer, intent(in)      :: maxZonesPerPlane

!  Local

   type(SetData),  pointer  :: Set
   type(AngleSet), pointer  :: ASet

   integer                  :: totalCycles

!  Allocate Memory 

   if (Size% ndim == 1) then
     return
   endif

   Set  => getSetData(Quad, setID)
   ASet => getAngleSetFromSetID(Quad, setID)

   totalCycles = max(ASet% totalCycles, 1) 

   call Allocator%allocate(Size%usePinnedMemory, Set%label, "cyclePsi", Set% cyclePsi, Set% Groups, totalCycles)

!  These are onlu used in the GPU sweep
   if (Size% useGPU) then
     call Allocator%allocate(Size%usePinnedMemory, Set%label, "Q",        Set% Q,Set% Groups, Size% maxCorner, maxZonesPerPlane)
     call Allocator%allocate(Size%usePinnedMemory, Set%label, "S",        Set% S,Set% Groups, Size% maxCorner, maxZonesPerPlane)

     Set% Q(:,:,:) = zero
     Set% S(:,:,:) = zero
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

