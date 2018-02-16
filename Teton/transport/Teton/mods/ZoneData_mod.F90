! ZoneData Module:  Contains data structures for zone information 

module ZoneData_mod 

  use kind_mod
  use constant_mod
  use cudafor

  private

! public interfaces

  public constructZone, constructGPUZone, constructZones_SoA, setZones_GPU_mesh, setZones_GPU_STotal
                                                                                 
  type, public :: ZoneData 

     integer              :: nCorner            ! number of corners
     integer              :: nCFaces            ! number of corner faces 
     integer              :: c0                 ! global corner ID of first corner
     integer              :: nFaces             ! number of zone faces
     real(adqt)           :: VolumeZone         ! zone volume
     real(adqt)           :: EnergyDensityOld   ! old energy density


     real(adqt), contiguous, pointer  :: Area(:)            ! corner area
     real(adqt), contiguous, pointer  :: Radius(:,:,:)      ! average radius

     real(adqt), pinned, allocatable  :: Volume(:)          ! corner volume
     real(adqt), pinned, allocatable  :: VolumeOld(:)       ! old corner volume
     !real(adqt), pinned, allocatable  :: volumeRatio(:)
     real(adqt), pinned, allocatable  :: Sigt(:)            ! total opacity
     real(adqt), pinned, allocatable  :: SigtInv(:)         ! reciprocal of total opacity
     real(adqt), pinned, allocatable  :: STotal(:,:)        ! fixed + scattering source
     !real(adqt), pinned, allocatable  :: STime(:,:,:)       ! time-dependent source
     real(adqt), pinned, allocatable  :: A_fp(:,:,:)        ! outward normals on corner faces 
     real(adqt), pinned, allocatable  :: A_ez(:,:,:)        !
     integer,    pinned, allocatable  :: Connect(:,:,:)     ! nearest neighbor connectivity 

  end type ZoneData 

  type(ZoneData), pointer, public :: Z
!$OMP threadprivate(Z)

  type, public :: GPU_ZoneData

     integer              :: nCorner            ! number of corners
     integer              :: nCFaces            ! number of corner faces 
     integer              :: c0                 ! global corner ID of first corner
     integer              :: nFaces             ! number of zone faces


     real(adqt), device, allocatable  :: Volume(:)          ! corner volume
     real(adqt), device, allocatable  :: volumeRatio(:)
     real(adqt), device, allocatable  :: Sigt(:)            ! total opacity
     real(adqt), device, allocatable  :: SigtInv(:)         ! reciprocal of total opacity
     real(adqt), device, allocatable  :: STotal(:,:)        ! fixed + scattering source
     real(adqt), device, allocatable  :: A_fp(:,:,:)        ! outward normals on corner faces 
     real(adqt), device, allocatable  :: A_ez(:,:,:)        !
     integer,    device, allocatable  :: Connect_reorder(:,:,:)     ! nearest neighbor connectivity 

  end type GPU_ZoneData



  type, public :: ZoneData_SoA

     ! integer, device, allocatable :: nCorner(:)
     ! integer, device, allocatable :: nCFaces(:)
     ! integer, device, allocatable :: c0(:)

     ! real(adqt), device, allocatable  :: Volume(:,:)        ! corner volume
     ! real(adqt), device, allocatable  :: volumeRatio(:,:)        ! old corner volume just needed to scale by volume changes
     ! real(adqt), device, allocatable  :: Sigt(:,:)          ! total opacity
     ! real(adqt), device, allocatable  :: SigtInv(:,:)       ! reciprocal of total opacity
     ! real(adqt), device, allocatable  :: STotal(:,:,:)      ! fixed + scattering source
!     real(adqt), device, allocatable  :: STime(:,:,:)     ! time dependent source new layout:(ig,c0+c,zone)
     real(adqt), pinned, allocatable  :: STime(:,:,:)     ! time dependent source new layout:(ig,c0+c,zone)

     !real(adqt), device, allocatable  :: A_fp(:,:,:,:)      ! outward normals on corner faces 
     !real(adqt), device, allocatable  :: A_ez(:,:,:,:)      !
     !integer,    device, allocatable  :: Connect(:,:,:,:)   ! nearest neighbor connectivity 

     ! create CPU memory area to hold these stored dot products:
     real(adqt), pinned, allocatable :: omega_A_fp(:,:,:,:) ! size: mF*mC*nZ*nA
     real(adqt), pinned, allocatable :: omega_A_ez(:,:,:,:) ! size: mF*mC*nZ*nA
     !integer, device, allocatable :: Connect_reorder(:,:,:,:) ! mF,mC,3,nZ
     

  end type ZoneData_SoA

  interface constructZone
    module procedure ZoneData_ctor
  end interface

  interface constructGPUZone
    module procedure GPU_ZoneData_ctor
  end interface

  interface constructZones_SoA
    module procedure ZoneData_SoA_ctor
  end interface

  interface setZones_GPU_mesh
    module procedure ZoneData_GPU_init_mesh
  end interface


  interface setZones_GPU_STotal
    module procedure ZoneData_GPU_update_STotal
  end interface


!  interface destruct
!    module procedure ZoneData_dtor
!  end interface

contains

!=======================================================================
! construct interface
!=======================================================================
                                                                                   
  subroutine ZoneData_ctor(self,        &
                           nCorner,     &
                           nCFaces,     &
                           corner0,     &
                           nFaces,      &
                           Connect)

    use Size_mod

    implicit none

!   Passed variables

    type(ZoneData), intent(inout)    :: self

    integer, intent(in)              :: nCorner
    integer, intent(in)              :: nCFaces       
    integer, intent(in)              :: corner0
    integer, intent(in)              :: nFaces
    integer, intent(in)              :: Connect(3,Size%maxcf,Size%maxCorner) 

!   Local

    integer          :: cID, i 

!   Set Properties

    self% nCorner = nCorner 
    self% nCFaces = nCFaces
    self% c0      = corner0
    self% nFaces  = nFaces

    allocate( self % VolumeOld(self% nCorner) )

    allocate( self % Volume(self% nCorner) )
    !allocate( self % volumeRatio(self% nCorner) )
    allocate( self % Sigt(Size% ngr) )
    allocate( self % SigtInv(Size% ngr) )
    allocate( self % A_fp(Size% ndim,self% nCFaces,self% nCorner) )
    allocate( self % A_ez(Size% ndim,self% nCFaces,self% nCorner) )
    allocate( self % Connect(3,self% nCFaces,self% nCorner) )
    allocate( self % STotal(Size% ngr,self% nCorner) )
    !allocate( self % STime(Size% ngr,self% nCorner,Size% nangSN) )

    if (Size%ndim == 2) then
      allocate( self % Area(self% nCorner) )
      allocate( self % Radius(2,2,self% nCorner) )
    endif

    do cID=1,self% nCorner
      do i=1,self% nCFaces
        self % Connect(1,i,cID) = Connect(1,i,cID)
        self % Connect(2,i,cID) = Connect(2,i,cID)
        self % Connect(3,i,cID) = Connect(3,i,cID) - corner0
      enddo
    enddo

    return

  end subroutine ZoneData_ctor



  subroutine GPU_ZoneData_ctor(self,        &
                           nCorner,     &
                           nCFaces,     &
                           corner0,     &
                           nFaces,      &
                           Connect)

    use Size_mod

    implicit none

!   Passed variables

    type(GPU_ZoneData), intent(inout)    :: self

    integer, intent(in)              :: nCorner
    integer, intent(in)              :: nCFaces       
    integer, intent(in)              :: corner0
    integer, intent(in)              :: nFaces
    integer, intent(in)              :: Connect(3,Size%maxcf,Size%maxCorner) 

!   Local

    integer          :: cID, i 

!   Set Properties

    self% nCorner = nCorner 
    self% nCFaces = nCFaces
    self% c0      = corner0
    self% nFaces  = nFaces

    allocate( self % Volume(self% nCorner) )
    allocate( self % volumeRatio(self% nCorner) )
    allocate( self % Sigt(Size% ngr) )
    allocate( self % SigtInv(Size% ngr) ) ! this not needed.
    allocate( self % A_fp(Size% ndim,self% nCFaces,self% nCorner) )
    allocate( self % A_ez(Size% ndim,self% nCFaces,self% nCorner) )
    !allocate( self % Connect(3,self% nCFaces,self% nCorner) )
    allocate( self % Connect_reorder(self%nCFaces, self%nCorner,3) )
    allocate( self % STotal(Size% ngr,self% nCorner) )


    ! think this is where the problem of too many tiny transfers comes in

    ! do cID=1,self% nCorner
    !   do i=1,self% nCFaces
    !     ! Different order than ZData, better to load coalesced into shared memory 
    !     self % Connect_reorder(i,cID,1) = Connect(1,i,cID)
    !     self % Connect_reorder(i,cID,2) = Connect(2,i,cID)
    !     self % Connect_reorder(i,cID,3) = Connect(3,i,cID) - corner0
    !   enddo
    ! enddo

    return

  end subroutine GPU_ZoneData_ctor




  subroutine ZoneData_SoA_ctor(self)

    use Size_mod

    implicit none

!   Passed variables

    type(ZoneData_SoA), intent(inout)    :: self

!   Set Properties

    !allocate( self % nCorner(Size% nzones))
    !allocate( self % nCFaces(Size% nzones))
    !allocate( self % c0(Size% nzones))

    !allocate( self % Volume(Size% maxCorner ,Size% nzones) )
    !allocate( self % volumeRatio(Size% maxCorner ,Size% nzones) )
    !allocate( self % Sigt(Size% ngr,Size% nzones) )
    !allocate( self % SigtInv(Size% ngr,Size% nzones) )
    !allocate( self % A_fp(Size% ndim,Size% maxcf,Size% maxCorner,Size% nzones) )
    !allocate( self % A_ez(Size% ndim,Size% maxcf,Size% maxCorner,Size% nzones) )
    

    allocate( self% omega_A_fp(Size% maxcf, Size% maxCorner,Size% nzones, Size% nangSN) )
    allocate( self% omega_A_ez(Size% maxcf, Size% maxCorner,Size% nzones, Size% nangSN) )


    !allocate( self % Connect(3,Size% maxcf,Size% maxCorner,Size% nzones) )
    !allocate( self % Connect_reorder(Size% maxcf, Size% maxCorner,3,Size% nzones) )
    
    !allocate( self % STotal(Size% ngr, Size% maxCorner, Size% nzones) )
    allocate( self % STime(Size% ngr, Size%ncornr, Size% nangSN ) )

    return

  end subroutine ZoneData_SoA_ctor



!   attributes(global) subroutine ZoneData_SoA_update_STotal_kernel(self,ZData,nzones,maxCorner,&
!        ngr,ndim,maxcf,nangsn)
    
!     ! this routine sets STotal, sigt, sigtinv on the GPU. 
!     ! It uses zero copy to populate the data on the GPU.

!     use Size_mod

!     implicit none

! !   Passed variables

!     type(ZoneData_SoA), device, intent(inout)    :: self
!     type(ZoneData),     device, intent(in)       :: ZData(nzones)

!     integer, value, intent(in) :: nzones,maxCorner,ngr,ndim,maxcf,nangsn

! !   Local

!     integer :: i,zone,c,ic,id,mic,mnd,nCorner, angle

!     i = threadIdx%x
!     zone = (blockIdx%y-1)*blockDim%y + threadIdx%y

!     if (zone <= nzones) then

!       !nCorner = ZData(zone)% nCorner
!       nCorner = self % nCorner(zone)

!       if (i <= ngr) then
!         self % Sigt(i,zone) = ZData(zone)% Sigt(i)
!         self % SigtInv(i,zone) = ZData(zone)% SigtInv(i)
!         do c=1,nCorner
!           self % STotal(i,c,zone) = ZData(zone)% STotal(i,c)
!         enddo
!       endif
  
!     endif
!   end subroutine ZoneData_SoA_update_STotal_kernel



!   attributes(global) &
!   subroutine ZoneData_SoA_init_mesh_kernel(self,ZData,nzones,maxCorner,ngr,ndim,maxcf,nangsn)

!     ! This routine sets up the mesh description on the GPU.

!     use Size_mod

!     implicit none

! !   Passed variables

!     type(ZoneData_SoA), device, intent(inout)    :: self
!     type(ZoneData),     device, intent(in)       :: ZData(nzones)

!     integer, value, intent(in) :: nzones,maxCorner,ngr,ndim,maxcf,nangsn

! !   Local

!     integer :: i,zone,c,ic,id,mic,mnd,nCorner, angle

!     i = threadIdx%x
!     zone = (blockIdx%y-1)*blockDim%y + threadIdx%y

!     if (zone <= nzones) then

!       nCorner = ZData(zone)% nCorner
!       self % nCorner(zone) = nCorner

!       self % nCFaces(zone) = ZData(zone)% nCFaces
!       self % c0(zone) = ZData(zone)% c0

!       if (i <= nCorner) then
!         self % Volume(i,zone) = ZData(zone)% Volume(i)
!         self % volumeRatio(i,zone) = ZData(zone)% VolumeOld(i)/ZData(zone)% Volume(i) 
!         ! try below later:
!         !self % volumeRatio(i+c0) = ZData(zone)% volumeRatio(i) ! only needed for first iteration in the set.
!       endif

!       ! if (i <= ngr) then
!       !   self % Sigt(i,zone) = ZData(zone)% Sigt(i)
!       !   self % SigtInv(i,zone) = ZData(zone)% SigtInv(i)
!       !   do c=1,nCorner
!       !     self % STotal(i,c,zone) = ZData(zone)% STotal(i,c)
!       !   enddo
!       ! endif
  
!       mnd = max(ndim,3)
!       mic = mnd * maxcf

!       if (i <= mic) then
!         ic = (i-1)/mnd + 1 ! split thread block x-dimension loop into two dims
!         id = i - ((ic-1)*mnd) ! remainder

!         do c=1,nCorner
!           if (id <= ndim) then
!             self % A_fp(id,ic,c,zone) = ZData(zone)% A_fp(id,ic,c)
!             self % A_ez(id,ic,c,zone) = ZData(zone)% A_ez(id,ic,c)
!           endif
!           if (id <= 3) then
!             self % Connect(id,ic,c,zone) = ZData(zone)% Connect(id,ic,c)
!           endif
!         enddo
!       endif

!     endif
!   end subroutine ZoneData_SoA_init_mesh_kernel


!   subroutine ZoneData_SoA_init_mesh(self,ZData)

!     use Size_mod

!     implicit none

! !   Passed variables

!     type(ZoneData_SoA), device, intent(inout)    :: self
!     type(ZoneData),     device, intent(in)       :: ZData(:) ! (nzones)

! !   Local

!     type(dim3) :: threads,blocks

! !   Function

!     !threads = dim3(max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf),16,1)
!     threads = dim3(max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf),4,1)
!     ! Number of threads as implemented will not give correct results when > 1024
!     if (max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf) .GT. 1024) then
!        print *,"ERROR: requested threads.x = ",  max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf)
!     endif
!     blocks  = dim3(1,(Size%nzones+threads%y-1)/threads%y,1)

!     call ZoneData_SoA_init_mesh_kernel<<<blocks,threads>>>(self,ZData,Size%nzones, &
!                                                       Size%maxCorner,Size%ngr,&
!                                                       Size%ndim,Size%maxcf,Size%nangsn)


!   end subroutine ZoneData_SoA_init_mesh


!   subroutine ZoneData_SoA_update_STotal(self,ZData)

!     use Size_mod

!     implicit none

! !   Passed variables

!     type(ZoneData_SoA), device, intent(inout)    :: self
!     type(ZoneData),     device, intent(in)       :: ZData(:) ! (nzones)

! !   Local

!     type(dim3) :: threads,blocks

! !   Function

!     !threads = dim3(max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf),16,1)
!     threads = dim3(max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf),4,1)
!     ! Number of threads as implemented will not give correct results when > 1024
!     if (max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf) .GT. 1024) then
!        print *,"ERROR: requested threads.x = ",  max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf)
!     endif
!     blocks  = dim3(1,(Size%nzones+threads%y-1)/threads%y,1)

!     call ZoneData_SoA_update_STotal_kernel<<<blocks,threads>>>(self,ZData,Size%nzones, &
!                                                       Size%maxCorner,Size%ngr,&
!                                                       Size%ndim,Size%maxcf,Size%nangsn)


!   end subroutine ZoneData_SoA_update_STotal



  attributes(global) subroutine ZoneData_GPU_update_STotal_kernel(self,ZData,nzones,maxCorner,&
       ngr,ndim,maxcf,nangsn)
    
    ! this routine sets STotal, sigt, sigtinv on the GPU. 
    ! It uses zero copy to populate the data on the GPU.

    use Size_mod

    implicit none

!   Passed variables

    type(GPU_ZoneData), device, intent(inout)    :: self(nzones)
    type(ZoneData),     device, intent(in)       :: ZData(nzones)

    integer, value, intent(in) :: nzones,maxCorner,ngr,ndim,maxcf,nangsn

!   Local

    integer :: i,zone,c,ic,id,mic,mnd,nCorner, angle

    i = threadIdx%x
    zone = (blockIdx%y-1)*blockDim%y + threadIdx%y

    if (zone <= nzones) then

      nCorner = self(zone) % nCorner

      if (i <= ngr) then
        self(zone) % Sigt(i) = ZData(zone)% Sigt(i)
        self(zone) % SigtInv(i) = ZData(zone)% SigtInv(i)
        do c=1,nCorner
          self(zone) % STotal(i,c) = ZData(zone)% STotal(i,c)
        enddo
      endif
  
    endif
  end subroutine ZoneData_GPU_update_STotal_kernel



  attributes(global) &
  subroutine ZoneData_GPU_init_mesh_kernel(self,ZData,nzones,maxCorner,ngr,ndim,maxcf,nangsn)

    ! This routine sets up the mesh description on the GPU.

    use Size_mod

    implicit none

!   Passed variables

    type(GPU_ZoneData), device, intent(inout)    :: self(nzones)
    type(ZoneData),     device, intent(in)       :: ZData(nzones)

    integer, value, intent(in) :: nzones,maxCorner,ngr,ndim,maxcf,nangsn

!   Local

    integer :: i,zone,c,icf,idim,mic,mnd,nCorner, angle

    i = threadIdx%x
    zone = (blockIdx%y-1)*blockDim%y + threadIdx%y

    if (zone <= nzones) then

      nCorner = ZData(zone)% nCorner
      
      if (i <= nCorner) then
        self(zone) % Volume(i) = ZData(zone)% Volume(i)
        self(zone) % volumeRatio(i) = ZData(zone)% VolumeOld(i)/ZData(zone)% Volume(i) 
        ! try below later:
        !self % volumeRatio(i+c0) = ZData(zone)% volumeRatio(i) ! only needed for first iteration in the set.
      endif

      mnd = max(ndim,3)
      mic = mnd * maxcf

      if (i <= mic) then
        icf = (i-1)/mnd + 1 ! split thread block x-dimension loop into two dims
        idim = i - ((icf-1)*mnd) ! remainder

        do c=1,nCorner
          if (idim <= ndim) then
            self(zone) % A_fp(idim,icf,c) = ZData(zone)% A_fp(idim,icf,c)
            self(zone) % A_ez(idim,icf,c) = ZData(zone)% A_ez(idim,icf,c)
          endif
          if (idim <= 3) then
             self(zone) % Connect_reorder(icf,c,idim) = ZData(zone)% Connect(idim,icf,c)
          endif
        enddo
      endif

    endif
  end subroutine ZoneData_GPU_init_mesh_kernel




  subroutine ZoneData_GPU_init_mesh(self,ZData)

    use Size_mod

    implicit none

!   Passed variables

    type(GPU_ZoneData), device, intent(inout)    :: self(:)
    type(ZoneData),     device, intent(in)       :: ZData(:) ! (nzones)

!   Local

    type(dim3) :: threads,blocks

!   Function

    !threads = dim3(max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf),16,1)
    threads = dim3(max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf),4,1)
    ! Number of threads as implemented will not give correct results when > 1024
    if (max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf) .GT. 1024) then
       print *,"ERROR: requested threads.x = ",  max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf)
    endif
    blocks  = dim3(1,(Size%nzones+threads%y-1)/threads%y,1)

    call ZoneData_GPU_init_mesh_kernel<<<blocks,threads>>>(self,ZData,Size%nzones, &
                                                      Size%maxCorner,Size%ngr,&
                                                      Size%ndim,Size%maxcf,Size%nangsn)


  end subroutine ZoneData_GPU_init_mesh





  subroutine ZoneData_GPU_update_STotal(self,ZData)

    use Size_mod

    implicit none

!   Passed variables

    type(GPU_ZoneData), device, intent(inout)    :: self(:)
    type(ZoneData),     device, intent(in)       :: ZData(:) ! (nzones)

!   Local

    type(dim3) :: threads,blocks

!   Function

    !threads = dim3(max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf),16,1)
    threads = dim3(max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf),4,1)
    ! Number of threads as implemented will not give correct results when > 1024
    if (max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf) .GT. 1024) then
       print *,"ERROR: requested threads.x = ",  max(Size%maxCorner,Size%ngr,max(Size%ndim,3)*Size%maxcf)
    endif
    blocks  = dim3(1,(Size%nzones+threads%y-1)/threads%y,1)

    call ZoneData_GPU_update_STotal_kernel<<<blocks,threads>>>(self,ZData,Size%nzones, &
                                                      Size%maxCorner,Size%ngr,&
                                                      Size%ndim,Size%maxcf,Size%nangsn)


  end subroutine ZoneData_GPU_update_STotal


!=======================================================================
! destruct interface
!=======================================================================
                                                                                    
!  subroutine ZoneData_dtor(self, Ndim)


!    implicit none

!   Passed variables
                                                                                     

!    return

!  end subroutine ZoneData_dtor


end module ZoneData_mod

