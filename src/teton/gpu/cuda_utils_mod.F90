! Wraps several C utility functions in CudaStepBoltz.cu
! These are currently being used in UpdateMaterialCoupling, so
! are being provided to the FORTRAN code.
module cuda_utils_mod

#if defined (TETON_ENABLE_CUDA)
use iso_c_binding

contains

   subroutine fcheapSync()
      implicit none

      interface 
         subroutine cheapSync () bind(c)
            use iso_c_binding
         end subroutine cheapSync
      end interface

      call cheapSync()
   end subroutine fcheapSync

   subroutine fcheapSetStream( i )
      implicit none
      integer :: i

      interface 
         subroutine cheapSetStream ( i ) bind(c)
            use iso_c_binding
            integer (c_int) :: i
         end subroutine cheapSetStream
      end interface

      call cheapSetStream(i)
   end subroutine fcheapSetStream

   subroutine fcheapSyncStream( i )
      implicit none
      integer :: i

      interface
         subroutine cheapSyncStream ( i ) bind(c)
         use iso_c_binding
       
         integer (c_int) :: i

         end subroutine cheapSyncStream
      end interface

      call cheapSyncStream(i)

   end subroutine fcheapSyncStream

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   subroutine getBCSolverMemEstimate(num_bytes )
      use ComptonControl_mod
      use Size_mod
      implicit none

      integer (c_int) :: ngr
      integer (c_int) :: tabG2Gs
      integer (c_int) :: tabTaus
      integer (c_int) :: zoneBatchSize
      integer (c_int) :: maxCornersPerZone

      integer (c_size_t) :: num_bytes

      interface 
         subroutine getrequiredgpumemory ( ngr, tabG2Gs, tabTaus, zoneBatchSize, maxCornersPerZone, num_bytes) bind (c)
            use iso_c_binding
       
            integer (c_int) :: ngr
            integer (c_int) :: tabG2Gs
            integer (c_int) :: tabTaus
            integer (c_int) :: zoneBatchSize
            integer (c_int) :: maxCornersPerZone
            integer (c_size_t) :: num_bytes

         end subroutine getrequiredgpumemory
      end interface

      num_bytes = 0

      if ( getUseBoltzmann(Compton) .AND. Size%useCUDASolver ) then
         ngr = Size%ngr
         tabG2Gs = Size%nBCITabG2Gs
         tabTaus = Size%nBCITabTaus
         zoneBatchSize = Size%zoneBatchSize
         maxCornersPerZone = Size%maxCorner

         call getrequiredgpumemory( ngr, tabG2Gs, tabTaus, zoneBatchSize, maxCornersPerZone, num_bytes )
      endif

      end subroutine getBCSolverMemEstimate

#endif

   subroutine fallocateGpuMemory( ngr, tabG2Gs, tabTaus, zoneBatchSize, maxCornersPerZone )
      implicit none

      integer (c_int) :: ngr
      integer (c_int) :: tabG2Gs
      integer (c_int) :: tabTaus
      integer (c_int) :: zoneBatchSize
      integer (c_int) :: maxCornersPerZone

      interface 
         subroutine allocateGpuMemory ( ngr, tabG2Gs, tabTaus, zoneBatchSize, maxCornersPerZone) bind (c)
            use iso_c_binding
       
            integer (c_int) :: ngr
            integer (c_int) :: tabG2Gs
            integer (c_int) :: tabTaus
            integer (c_int) :: zoneBatchSize
            integer (c_int) :: maxCornersPerZone

         end subroutine allocateGpuMemory
      end interface

      call allocateGpuMemory( ngr, tabG2Gs, tabTaus, zoneBatchSize, maxCornersPerZone )

      end subroutine fallocateGpuMemory

   subroutine ffreeGpuMemory()
      implicit none

      interface 
         subroutine freeGpuMemory () bind (c)
            use iso_c_binding
         end subroutine freeGpuMemory
      end interface

      call freeGpuMemory()
   end subroutine ffreeGpuMemory
#endif

end module cuda_utils_mod
