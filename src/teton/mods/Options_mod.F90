#include "macros.h"
!=======================================================================
! Functions for getting/setting runtime options
!=======================================================================
 
module Options_mod 
  use Datastore_mod, only : theDatastore
  use, intrinsic :: iso_fortran_env, only: int32 
  use, intrinsic :: iso_c_binding, only : c_bool, c_int
  implicit none

  private

  type :: options_type
  contains
    procedure :: getNumOmpMaxThreads
    procedure :: initialize
    procedure :: check
    procedure :: getVerbose
    procedure :: isRankVerbose
    procedure :: isRankZeroVerbose
    procedure :: setVerbose
    procedure :: getSweepVersion
    procedure :: getSweepNumHyperDomains
    procedure :: getGTANumHyperDomains
  end type options_type

  public :: setVerboseOldCVersion

  type(options_type), public :: Options

contains

!=======================================================================
! Populate options with default values.
! TODO - Move any appropriate defaults to
! Radar C++ API( for general TRT options )
! instead of here in the Fortran.
!=======================================================================
  subroutine initialize(self)
#if defined(TETON_ENABLE_OPENMP)
    use omp_lib
#endif
    use cmake_defines_mod

    class(options_type) :: self

    logical*4 :: temp
    integer :: omp_cpu_max_threads

    ! Default to current max number of threads value set in OpenMP runtime
    temp = theDatastore%root%has_path("options/concurrency/omp_cpu_max_threads")
    if (.NOT. temp) then
#if defined(TETON_ENABLE_OPENMP)
      omp_cpu_max_threads = omp_get_max_threads()
#else
      omp_cpu_max_threads = 1
#endif
      call theDatastore%root%set_path("options/concurrency/omp_cpu_max_threads", omp_cpu_max_threads)
    endif

    ! Build information.
    ! These are populated in the cmake_defines_mod.F90
    call theDatastore%root%set_path("build_meta_data/cmake_install_prefix",install_prefix)
    call theDatastore%root%set_path("build_meta_data/project_version", version)
    call theDatastore%root%set_path("build_meta_data/git_sha1", git_sha1)
    call theDatastore%root%set_path("build_meta_data/cmake_system", system_type)
    call theDatastore%root%set_path("build_meta_data/cmake_cxx_compiler", cxx_compiler)
    call theDatastore%root%set_path("build_meta_data/cmake_fortran_compiler", fortran_compiler)

    ! Default to not verbose
    temp = theDatastore%root%has_path("options/verbose")
    if (.NOT. temp) then
      call theDatastore%root%set_path("options/verbose", 0)
    endif
      
    call self%check()
  end subroutine

!=======================================================================
! Validate the option values.
! TODO - Implement this in
! Radar C++ API( for general TRT options )
!
! Leave in Fortran until then, so we can catch any input errors early.
! Also need some assert macros in the C++...
!=======================================================================
  subroutine check(self)
    class(options_type) :: self
    logical*4 :: temp

    ! Retrieve value in local 4 byte logical to workaround buffer overflow in
    ! conduit API.  Will be fixed in future conduit version.
    temp = theDatastore%root%has_path("options/concurrency")
    TETON_VERIFY(temp, "Options tree is missing 'concurrency'.")
    temp = theDatastore%root%has_path("options/concurrency/omp_cpu_max_threads")
    TETON_VERIFY(temp, "Options is missing concurrency/omp_cpu_max_threads.")

    temp = theDatastore%root%has_path("options/verbose")
    TETON_VERIFY(temp, "Options tree is missing verbose.")

  end subroutine
!=======================================================================
! Accessor functions for getting options.
!=======================================================================
  integer(kind=int32) function getNumOmpMaxThreads(self) result(numOmpMaxThreads)
    class(options_type) :: self
    numOmpMaxThreads = theDatastore%root%fetch_path_as_int32("options/concurrency/omp_cpu_max_threads")
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
!    setVerbose - Set verbosity level.
!    
!    verbose=x0 - no verbose output
!    verbose=0x - rank 0 at verbose level x
!    verbose=1x - all ranks at verbose level x
!***********************************************************************
   subroutine setVerbose(self, level)
      class(options_type) :: self
      integer(kind=C_INT), intent(in) :: level
      integer :: rank, ierr

      ! TODO - Remove this when we have a proper teton initialize() function
      call theDatastore%initialize()

      call theDatastore%root%set_path("options/verbose", level)

      if ( self%isRankVerbose() > 0 ) then
         print *, "Teton: verbose output enabled."
      endif

      return
   end subroutine setVerbose

!***********************************************************************
!    setVerbose - Set verbosity level.  Old version that doesn't use   *
!    type bound procedure.                                             *
!***********************************************************************
   subroutine setVerboseOldCVersion(level) BIND(C,NAME="teton_setverbose")
      integer(kind=C_INT), intent(in) :: level

      ! TODO - Remove this when we have a proper teton initialize() function
      call theDatastore%initialize()

      call Options%setVerbose(level)

      return
   end subroutine setVerboseOldCVersion

!***********************************************************************
!    setSweepVersion - Set the sweep implementation version to use
!    from the input.
!***********************************************************************
  subroutine setSweepVersion(sweepversion) BIND(C,NAME="teton_setsweepversion")
    integer(kind=C_INT), intent(in) :: sweepversion

    ! 0 = zone sweep
    ! 1 = corner sweep
    ! Do not set this to anything else, getSweepVersion will handle the zone sweep case.
    if ( sweepversion == 1 ) then
       call theDatastore%initialize()

       call theDatastore%root%set_path("options/sweep/kernel/version", sweepversion)
    endif

    return
  end subroutine setSweepVersion

!***********************************************************************
!    getSweepVersion - Return the sweep implementation version to use
!    from the input.
!
!    The default sweep version is the zone sweep.
!***********************************************************************
  integer(kind=c_int) function getSweepVersion(self) result(sweepVersion)
    class(options_type) :: self
    logical*4 :: temp

    temp = theDatastore%root%has_path("options/sweep/kernel/version")
    if (.NOT. temp) then
      sweepVersion = 0
    else
      sweepVersion = theDataStore%root%fetch_path_as_int32("options/sweep/kernel/version")
    endif

    return
  end function

!***********************************************************************
!    setSweepNumHyperDomains - Set the sweep number of hyper-domains 
!    from the input.
!***********************************************************************
  subroutine setSweepNumHyperDomains(hyperdomains) BIND(C,NAME="teton_setsweepnumhyperdomains")
    integer(kind=C_INT), intent(in) :: hyperdomains

    ! 0 - Let Teton automatically set this value. (default)
    ! >= 1 - Set to this number of hyper-domains.

    TETON_VERIFY(hyperdomains >= 0, "Setting sweep hyper-domains: the number of hyper-domains must be >= 0.")

    call theDatastore%initialize()

    call theDatastore%root%set_path("options/sweep/sn/numhyperdomains", hyperdomains)

    return
  end subroutine setSweepNumHyperDomains

!***********************************************************************
!    getSweepHyperDomains - Return the number of sweep hyper-domains to
!    use.
!***********************************************************************
  integer(kind=c_int) function getSweepNumHyperDomains(self) result(hyperdomains)
    class(options_type) :: self
    logical*4 :: temp

    temp = theDatastore%root%has_path("options/sweep/sn/numhyperdomains")
    if (.NOT. temp) then
!     Default to 0, which indicates user did not provide a setting. Teton
!     will set this to a non-zero value in ConstructPhaseSpaceSets.F90 
      hyperdomains = 0
    else
      hyperdomains = theDataStore%root%fetch_path_as_int32("options/sweep/sn/numhyperdomains")
      TETON_VERIFY(hyperdomains >= 0, "Getting sweep hyper-domains: the number of hyper-domains must be >= 0.")
    endif

    return
  end function

!***********************************************************************
!    setGTANumHyperDomains - Set the new GTA number of hyper-domains 
!    from the input.
!***********************************************************************
  subroutine setGTANumHyperDomains(hyperdomains) BIND(C,NAME="teton_setgtanumhyperdomains")
    integer(kind=C_INT), intent(in) :: hyperdomains

    ! 0 - Let Teton automatically set this value. (default)
    ! >= 1 - Set to this number of hyper-domains.

    TETON_VERIFY(hyperdomains >= 0, "Setting new GTA hyper-domains: the number of hyper-domains must be >= 0.")

    call theDatastore%initialize()

    call theDatastore%root%set_path("options/sweep/gta/numhyperdomains", hyperdomains)

    return
  end subroutine setGTANumHyperDomains

!***********************************************************************
!    getGTANumHyperDomains - Return the number of new GTA hyper-domains to
!    use.
!***********************************************************************
  integer(kind=c_int) function getGTANumHyperDomains(self) result(hyperdomains)
    class(options_type) :: self
    logical*4 :: temp

    temp = theDatastore%root%has_path("options/sweep/gta/numhyperdomains")
    if (.NOT. temp) then
!     Default to 0, which indicates user did not provide a setting. Teton
!     will set this to a non-zero value in ConstructPhaseSpaceSets.F90 
      hyperdomains = 0
    else
      hyperdomains = theDataStore%root%fetch_path_as_int32("options/sweep/gta/numhyperdomains")
      TETON_VERIFY(hyperdomains >= 0, "Getting new GTA hyper-domains: the number of hyper-domains must be >= 0.")
    endif

    return
  end function

end module Options_mod
