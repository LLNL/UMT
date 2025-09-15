#include "macros.h"
!=======================================================================
! Functions for getting/setting runtime options
!
! The newer C++ conduit interface has an initialize() method that will
! get called early on.
!
! For the older interface there is no central
! initialize function, so we can't assume the datastore or options
! have been initialized before any of the below functions are called.
! To protect against this each call will check if the options defaults
! have been initialized yet and if not, will do so.  This can be
! improved upon later.
!
! A. Black 4/15/2025
!=======================================================================
 
module Options_mod 
  use Datastore_mod, only : theDatastore
  use, intrinsic :: iso_fortran_env, only : int32 
  use, intrinsic :: iso_c_binding, only : c_bool, c_int

#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
  use OMPUtilities_mod, only : get_gpu_processor_count
#endif

  implicit none

  private

  type :: options_type
    logical :: is_initialized = .FALSE.

  contains
    procedure :: getNumOmpMaxThreads
    procedure :: getNumDeviceProcessors
    procedure :: getVerbose
    procedure :: getSweepVersion
    procedure :: getSweepMaxHyperDomains
    procedure :: getGTAMaxHyperDomains
    procedure :: getMinGroupSetSize
    procedure :: initialize
    procedure :: isRankVerbose
    procedure :: isRankZeroVerbose
  end type options_type

  type(options_type), public :: Options

contains

!=======================================================================
! Populate options with default values for any entries that do not yet
! exist.
!=======================================================================
  subroutine initialize(self)
#if defined(TETON_ENABLE_OPENMP)
    use omp_lib
#endif
    use cmake_defines_mod, only : install_prefix, version, git_sha1, system_type, cxx_compiler, fortran_compiler, &
                                  omp_device_team_thread_limit, min_groupset_size, &
                                  max_num_sweep_hyperdomains, max_num_gta_hyperdomains

    class(options_type) :: self
    logical*4 :: temp
    integer :: val

    if ( self%is_initialized ) then
      if ( self%isRankVerbose() > 0 ) then
         print *, "Teton: attempt to initialize default run time options more than once."
      endif
    endif
 
    call theDatastore%initialize()

    temp = theDatastore%root%has_path("options/sweep/kernel/version")
    if (temp) then
      val = theDatastore%root%fetch_path_as_int32("options/sweep/kernel/version")
      ! 0 = use default sweep version, which is 1=zone sweep
      if (val == 0) then
        call theDatastore%root%set_path("options/sweep/kernel/version", 1)
      endif
    else
      call theDatastore%root%set_path("options/sweep/kernel/version", 1)
    endif

    ! Default to current max number of threads value set in OpenMP runtime
    temp = theDatastore%root%has_path("options/concurrency/omp_cpu_max_threads")
    if (.NOT. temp) then
#if defined(TETON_ENABLE_OPENMP)
      val = omp_get_max_threads()
#else
      val = 1
#endif
      call theDatastore%root%set_path("options/concurrency/omp_cpu_max_threads", val)
    endif

    ! Default to the number of processors advertised by the GPU.
    temp = theDatastore%root%has_path("options/concurrency/omp_device_num_processors")
    if (.NOT. temp) then
#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
      val = get_gpu_processor_count()
#else
    ! Default to 1 if no GPU, as we rarely want to use more than one thread team on a CPU.
      val = 1
#endif
      call theDatastore%root%set_path("options/concurrency/omp_device_num_processors", val)
    endif

    ! Default for the max number of threads to use in a device thread team.  We have no use cases, as yet, to not use
    ! the max available (1024) for GPUs.
    temp = theDatastore%root%has_path("options/concurrency/omp_device_team_thread_limit")
    if (.NOT. temp) then
      call theDatastore%root%set_path("options/concurrency/omp_device_team_thread_limit", omp_device_team_thread_limit)
    endif

    ! Build information.
    ! These are populated in the cmake_defines_mod.F90
    call theDatastore%root%set_path("build_meta_data/cmake_install_prefix",install_prefix)
    call theDatastore%root%set_path("build_meta_data/project_version", version)
    call theDatastore%root%set_path("build_meta_data/git_sha1", git_sha1)
    call theDatastore%root%set_path("build_meta_data/cmake_system", system_type)
    call theDatastore%root%set_path("build_meta_data/cmake_cxx_compiler", cxx_compiler)
    call theDatastore%root%set_path("build_meta_data/cmake_fortran_compiler", fortran_compiler)

    ! Default to not verbose, if this was not set already
    temp = theDatastore%root%has_path("options/verbose")
    if (.NOT. temp) then
      call theDatastore%root%set_path("options/verbose", 0)
    endif
      
    ! Default to the sweep max number of hyperdomains from CMake, if this was not set already.
    temp = theDatastore%root%has_path("options/sweep/max_hyperdomains")
    if (.NOT. temp) then
      call theDatastore%root%set_path("options/sweep/max_hyperdomains", max_num_sweep_hyperdomains)
    endif

    ! Default to the gta max number of hyperdomains from CMake, if this was not set already.
    temp = theDatastore%root%has_path("options/gta/max_hyperdomains")
    if (.NOT. temp) then
      call theDatastore%root%set_path("options/gta/max_hyperdomains", max_num_gta_hyperdomains)
    endif

    ! Default the minimum group set size to value that CMake sets in
    ! cmake_defines_mod.F90.in.  CMake currently determines this depending on
    ! the value of CMAKE_<CUDA/HIP>_ARCHITECTURES passed in at build time.
    temp = theDatastore%root%has_path("options/sweep/min_group_set_size")
    if (.NOT. temp) then
      call theDatastore%root%set_path("options/sweep/min_group_set_size", min_groupset_size)
    endif

      
    self%is_initialized = .TRUE.
  end subroutine initialize

!=======================================================================
! Accessor functions for getting options.
!=======================================================================
  integer(kind=int32) function getNumOmpMaxThreads(self) result(numOmpMaxThreads)
    class(options_type) :: self

    if (.NOT. self%is_initialized) then
       call self%initialize()
    endif

    numOmpMaxThreads = theDatastore%root%fetch_path_as_int32("options/concurrency/omp_cpu_max_threads")
    TETON_VERIFY(numOmpMaxThreads > 0, "getNumOmpMaxThreads options/concurrency/omp_cpu_max_threads came back as zero, minimum should be one.")

  end function

  integer(kind=int32) function getNumDeviceProcessors(self) result(numDeviceProcessors)
    class(options_type) :: self

    if (.NOT. self%is_initialized) then
       call self%initialize()
    endif

    numDeviceProcessors = theDatastore%root%fetch_path_as_int32("options/concurrency/omp_device_num_processors")
    TETON_VERIFY(numDeviceProcessors > 0, "getNumDeviceProcessors options/concurrency/omp_device_num_processors came back as zero, minimum should be one.")

  end function

!***********************************************************************
!    getVerbose - Get verbosity level.                                
!
!    verbose=x0 - no verbose output
!    verbose=0x - rank 0 at verbose level x
!    verbose=1x - all ranks at verbose level x
!***********************************************************************
   integer(kind=int32) function getVerbose(self) result(level)
     class(options_type) :: self

     if (.NOT. self%is_initialized) then
        call self%initialize()
     endif

     level = theDatastore%root%fetch_path_as_int32("options/verbose")
   end function getVerbose

!***********************************************************************


!***********************************************************************
!    isRankVerbose - returns the local verbosity level of this rank
!***********************************************************************
   integer function isRankVerbose(self) result(verbose)
      use mpi

      class(options_type) :: self
      integer :: verbose_level
      integer :: rank, ierr

      call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
      verbose_level = self%getVerbose()

      ! verbose_level is a two digit number.
      ! The first digit determines whether only rank 0 is verbose
      ! The second digit determines the verbosity level.
      if (rank == 0 .OR. verbose_level > 10) then
         verbose = mod(verbose_level,10)
      else
         verbose = 0
      endif

      return

   end function isRankVerbose

!***********************************************************************
!    isRankZeroVerbose - returns the verbosity level of rank 0
!      This is useful for cases where all ranks need to do something
!      because of rank 0 being verbose.
!***********************************************************************
   integer function isRankZeroVerbose(self) result(verbose)
      use mpi
      class(options_type) :: self

      integer :: verbose_level

      verbose_level = self%getVerbose()
      verbose = mod(verbose_level,10)

      return

   end function isRankZeroVerbose

!***********************************************************************
!    setVerbose - Set verbosity level.  Accessible from C.             
!
!    verbose=x0 - no verbose output
!    verbose=0x - rank 0 at verbose level x
!    verbose=1x - all ranks at verbose level x
!
!    Note:
!    The older interface prior to teton conduit interface does not
!    have a initialize() starting function, so we can't assume the
!    options have been initialized prior to this call.
!***********************************************************************
   subroutine set_verbose_c(level) BIND(C,NAME="teton_setverbose")
      integer(kind=C_INT), intent(in) :: level

      if (.NOT. Options%is_initialized) then
         call Options%initialize()
      endif

      call theDatastore%root%set_path("options/verbose", level)

      if ( Options%isRankVerbose() > 0 ) then
         print *, "Teton: verbose output enabled."
      endif

      return
   end subroutine set_verbose_c

!***********************************************************************
!    getSweepVersion - Return the sweep implementation version to use
!    from the input.
!
!    The default sweep version is the zone sweep.
!***********************************************************************
  integer(kind=c_int) function getSweepVersion(self) result(sweepVersion)
    class(options_type) :: self

    if (.NOT. self%is_initialized) then
       call self%initialize()
    endif

    sweepVersion = theDatastore%root%fetch_path_as_int32("options/sweep/kernel/version")

    return
  end function getSweepVersion

!***********************************************************************
!    getSweepMaxHyperDomains - Return the max number of sweep hyper-domains to
!    use.
!***********************************************************************
  integer(kind=c_int) function getSweepMaxHyperDomains(self) result(hyperdomains)
    class(options_type) :: self

    if (.NOT. self%is_initialized) then
       call self%initialize()
    endif

    hyperdomains = theDatastore%root%fetch_path_as_int32("options/sweep/max_hyperdomains")
    TETON_VERIFY(hyperdomains > 0, "Getting sweep max # hyper-domains: the number of hyper-domains must be > 0.")

    return
  end function getSweepMaxHyperDomains

!***********************************************************************
!    getGTAMaxHyperDomains - Return the max number of GTA hyper-domains to
!    use.
!***********************************************************************
  integer(kind=c_int) function getGTAMaxHyperDomains(self) result(hyperdomains)
    class(options_type) :: self

    if (.NOT. self%is_initialized) then
       call self%initialize()
    endif

    hyperdomains = theDatastore%root%fetch_path_as_int32("options/gta/max_hyperdomains")
    TETON_VERIFY(hyperdomains > 0, "Getting GTA max hyper-domains: the number of hyper-domains must be > 0.")

    return
  end function getGTAMaxHyperDomains

!***********************************************************************
!    getMinGroupSetSize - Returns the minimum size allowed for group set.
!
!    This value is used when the code decomposes the groups into phase-space sets.
!***********************************************************************
  integer(kind=c_int) function getMinGroupSetSize(self) result(mingroupsetsize)
    class(options_type) :: self

    if (.NOT. self%is_initialized) then
       call self%initialize()
    endif

    mingroupsetsize = theDatastore%root%fetch_path_as_int32("options/sweep/min_group_set_size")
    TETON_VERIFY(mingroupsetsize > 0, "Getting the min group set size: the min size must be > 0.")

    return
  end function getMinGroupSetSize

!***********************************************************************
!    initialize, C accessible version
!    
!    Over time the code that initializes the defaults should be moved
!    to the C++ code, as the TetonConduitInterface is already
!    initializing multiple defaults.
!***********************************************************************
   subroutine initialize_c() BIND(C,NAME="teton_initialize_defaults")
      call Options%initialize()
      return
   end subroutine initialize_c

end module Options_mod
