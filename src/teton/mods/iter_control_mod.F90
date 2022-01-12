module iter_control_mod

#include "macros.h"
use kind_mod
use constant_mod

!=======================================================================
!                       Version 1.1: 05/2018, maginot1
!-----------------------------------------------------------------------
! Iteration Control
!   This class contains an iteration control entry
!
! epsilonPoint    pointwise relative error convergence criterion
! localError
! globalError
! maxIter         maximum number of iterations allowed to be taken
! nIter           number of iterations actually taken before an iteration was stopped
! nTotIter        total number of iterations over 
! zoneOfMax       zone number of maximum relative pointwise error
! processOfMax    process that took the highest number of iterations or that had the highest error
!
!=======================================================================

private

! public interfaces
  public construct 
  public setControls 
  public resetNumberOfIterations
  public setNumberOfIterations
  public setMaxNumberOfIterations
  public setGlobalMaxIterationsTaken   !
  public setZoneOfMax 
  public setProcessOfMax
  public setLocalError 
  public setGlobalError 
  public destruct
  public getEpsilonPoint 
  public getLocalError 
  public getGlobalError
  public getMaxNumberOfIterations
  public getGlobalMaxIterationsTaken !
  public getNumberOfIterations
  public getTotalNumberOfIterations 
  public getZoneOfMax
  public getProcessOfMax

  type, public :: IterControl
    private
    real(adqt)    :: epsilonPoint   ! a target convergence tolerance
    real(adqt)    :: localError
    real(adqt)    :: globalError
    integer       :: maxIter        ! maximum number of iterations allowed
    integer       :: globalMaxIterTaken ! highest number of iterations taken over all zones/processes (nonlinear Iteration control)
    integer       :: nIter
    integer       :: nTotIter
    integer       :: zoneOfMax
    integer       :: processOfMax
  end type IterControl

  interface construct
    module procedure iter_control_ctor
  end interface

  interface setControls
    module procedure iter_control_set
  end interface

  interface resetNumberOfIterations
    module procedure iter_control_reset_nIter
  end interface

  interface setNumberOfIterations
    module procedure iter_control_set_nIter
  end interface

  interface setMaxNumberOfIterations
    module procedure iter_control_set_maxIter
  end interface

  interface setGlobalMaxIterationsTaken
    module procedure iter_control_set_globalMaxIterTaken
  end interface
  
  interface setZoneOfMax
    module procedure iter_control_set_zoneOfMax
  end interface

  interface setProcessOfMax
    module procedure iter_control_set_processOfMax
  end interface

  interface setLocalError
    module procedure iter_control_set_localError
  end interface

  interface setGlobalError
    module procedure iter_control_set_globalError
  end interface

  interface destruct
    module procedure iter_control_dtor
  end interface

  interface getEpsilonPoint
    module procedure iter_control_get_epsilonPoint
  end interface

  interface getMaxNumberOfIterations
    module procedure iter_control_get_maxIter
  end interface 

  interface getGlobalMaxIterationsTaken
    module procedure iter_control_get_globalMaxIterTaken
  end interface 

  interface getNumberOfIterations
    module procedure iter_control_get_nIter
  end interface

  interface getTotalNumberOfIterations
    module procedure iter_control_get_nTotIter
  end interface

  interface getZoneOfMax
    module procedure iter_control_get_zoneOfMax
  end interface

  interface getProcessOfMax
    module procedure iter_control_get_ProcessOfMax
  end interface

  interface getLocalError
    module procedure iter_control_get_localError
  end interface

  interface getGlobalError
    module procedure iter_control_get_globalError
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================

  subroutine iter_control_ctor(self,epsilonPoint,localError, &
       globalError,maxNumberOfIterations, globalMaxIterTaken, &
       numberOfIterations, &
       totalNumberOfIterations,zoneOfMaximum,processOfMaximum)

!    Construct the iteration control object
!      epsilonPoint              rel. point error convergence criterion
!      localError                maximum pointwise relative error for process
!      globalError               global maximum pointwise relative error
!      maxNumberOfIterations     maximum number of iterations
!      numberOfIterations        number of iterations
!      totalNumberOfIterations   total number of iterations
!      zoneOfMaximum             zone of maximum pointwise error
!      processOfMaximum          process with the maximum pointwise error

!    variable declarations
     implicit none

!    passed variables
     type(IterControl),    intent(inout) :: self
     real(adqt), optional, intent(in)    :: epsilonPoint
     real(adqt), optional, intent(in)    :: localError
     real(adqt), optional, intent(in)    :: globalError
     integer,    optional, intent(in)    :: maxNumberOfIterations
     integer,    optional, intent(in)    :: globalMaxIterTaken
     integer,    optional, intent(in)    :: numberOfIterations
     integer,    optional, intent(in)    :: totalNumberOfIterations
     integer,    optional, intent(in)    :: zoneOfMaximum
     integer,    optional, intent(in)    :: processOfMaximum

!    construct the iteration control object

     if (present(epsilonPoint)) then
        self % epsilonPoint = epsilonPoint
     else
        self % epsilonPoint = 1.0e-4_adqt
     endif

     if (present(localError)) then
        self % localError = localError
     else
        self % localError = 0.0_adqt
     endif

     if (present(globalError)) then
        self % globalError = globalError
     else
        self % globalError = 0.0_adqt
     endif

     if (present(maxNumberOfIterations)) then
        self % maxIter = maxNumberOfIterations
     else
        self % maxIter = 10
     endif

     if (present(numberOfIterations)) then
        self % nIter = numberOfIterations
     else
        self % nIter = 0
     endif

     if (present(totalNumberOfIterations)) then
        self % nTotIter = totalNumberOfIterations
     else
        self % nTotIter = 0
     endif

     if (present(zoneOfMaximum)) then
        self % zoneOfMax = zoneOfMaximum
     else
        self % zoneOfMax = 0
     endif

     if (present(processOfMaximum)) then
        self % processOfMax = processOfMaximum
     else
        self % processOfMax = 0
     endif

     if (present(globalMaxIterTaken)) then
        self % globalMaxIterTaken = globalMaxIterTaken
     else
        self % globalMaxIterTaken = 0
     endif

!    assertions
     tetonAssert(self%epsilonPoint>zero,"Invalid iter control ctor")
     tetonAssert(self%localError>=zero,"Invalid iter control ctor")
     tetonAssert(self%globalError>=zero,"Invalid iter control ctor")
     tetonAssert(self%maxIter>zero,"Invalid iter control ctor")
     tetonAssert(self%globalMaxIterTaken>=zero,"Invalid iter control ctor")
     tetonAssert(self%nIter>=zero,"Invalid iter control ctor")
     tetonAssert(self%nTotIter>=zero,"Invalid iter control ctor")
     tetonAssert(self%zoneOfMax>=zero,"Invalid iter control ctor")
     tetonAssert(self%processOfMax>=zero,"Invalid iter control ctor")

     return
  end subroutine iter_control_ctor

!=======================================================================
! setControls interface
!=======================================================================

  subroutine iter_control_set(self,epsilonPoint,maxNumberOfIterations)

!    Set the iteration controls in the iteration control object
!      epsilonPoint            rel. point error convergence criterion
!      maxNumberOfIterations   maximum number of iterations

!    variable declarations
     implicit none

!    passed variables
     type(IterControl),    intent(inout) :: self
     real(adqt), optional, intent(in)    :: epsilonPoint
     integer,    optional, intent(in)    :: maxNumberOfIterations

!    set the iteration controls of the iteration control object

     if (present(epsilonPoint)) then
!       Don't allow values that are too small.
        if ( epsilonPoint > 1.0e-14_adqt ) then
           self % epsilonPoint = epsilonPoint
        else
           self % epsilonPoint = 1.0e-14_adqt
        endif
     endif

     if (present(maxNumberOfIterations)) then
        self % maxIter = maxNumberOfIterations
     endif

     return
  end subroutine iter_control_set

!=======================================================================
! resetNumberOfIterations interface
!=======================================================================

  subroutine iter_control_reset_nIter(self)

!    Reset the number of iteration in the iteration control object

!    variable declarations
     implicit none

!    passed variables
     type(IterControl), intent(inout) :: self

!    reset the number of iterations
     self % nIter = 0
     self % nTotIter = 0

!    assertions
     tetonAssert(self%nIter==0,"Invalid iter control reset")
     tetonAssert(self%nTotIter==0,"Invalid iter control reset")

     return
  end subroutine iter_control_reset_nIter


!=======================================================================
! setNumberOfIterations interface
!=======================================================================

  subroutine iter_control_set_nIter(self, nIter)

!    Set the number of iterations in the iteration control object

!    variable declarations
     implicit none

!    passed variables
     type(IterControl), intent(inout) :: self
     integer,           intent(in)    :: nIter

!    assertions
     tetonAssert(nIter>=0,"Invalid number of iterations")

!    reset the number of iterations
     self % nIter = nIter
     self % nTotIter = self % nTotIter + nIter

!    assertions
     tetonAssert(self%nIter>=0,"Invalid number of iterations")
     tetonAssert(self%nTotIter>=0,"Invalid number of iterations")

     return
  end subroutine iter_control_set_nIter

!=======================================================================
! setMaxNumberOfIterations interface
!=======================================================================

  subroutine iter_control_set_maxIter(self, maxIter)

!    Set the maximum number of iterations in the iteration control object

!    variable declarations
     implicit none

!    passed variables
     type(IterControl), intent(inout) :: self
     integer,           intent(in)    :: maxIter

!    assertions
     tetonAssert(maxIter>=0, "Invalid maximum number of iterations")

!    reset the number of iterations
     self % maxIter = maxIter

!    assertions
     tetonAssert(self%maxIter>=0, "Invalid number of iterations")

     return
   end subroutine iter_control_set_maxIter

!=======================================================================
! setGlobalMaxIterationsTaken interface
!=======================================================================

  subroutine iter_control_set_globalMaxIterTaken(self, maxTaken)

!    Set the maximum number of iterations in the iteration control object

!    variable declarations
     implicit none

!    passed variables
     type(IterControl), intent(inout) :: self
     integer,           intent(in)    :: maxTaken

!    assertions
     tetonAssert(maxTaken>=0, "Invalid maximum number of iterations taken")

!    reset the number of iterations
     self % globalMaxIterTaken = maxTaken

!    assertions
     tetonAssert(self%globalMaxIterTaken>=0, "Invalid maximum number of iterations taken")

     return
   end subroutine iter_control_set_globalMaxIterTaken

!=======================================================================
! setZoneOfMax interface
!=======================================================================
                                                                                      
  subroutine iter_control_set_zoneOfMax(self, zoneOfMax)
                                                                                      
!    Set the zone number of the maximum relative pointwise error
!     zoneOfMax   zone of maximum relative pointwise error
                                                                                      
!    variable declarations
     implicit none
                                                                                      
!    passed variables
     type(IterControl), intent(inout) :: self
     integer,           intent(in)    :: zoneOfMax 
                                                                                      
!    assertions
     tetonAssert(zoneOfMax>=0,"Invalid zone number")
                                                                                      
!    reset the zone with maximum relative error 
     self % zoneOfMax = zoneOfMax 
                                                                                      
     return
  end subroutine iter_control_set_zoneOfMax

!=======================================================================
! setProcessOfMax interface
!=======================================================================
                                                                                      
  subroutine iter_control_set_processOfMax(self, processOfMax)
                                                                                      
!    Set the process number that has the maximum relative pointwise error
!    processOfMax   process of maximum relative pointwise error
                                                                                      
!    variable declarations
     implicit none
 
!    passed variables
     type(IterControl), intent(inout) :: self
     integer,           intent(in)    :: processOfMax
 
!    assertions
     tetonAssert(processOfMax>=0,"Invalid process number")
 
!    reset the process with maximum relative error
     self % processOfMax = processOfMax
 
!    assertions
     tetonAssert(self%processOfMax>=0,"Invalid process number")
 
     return
  end subroutine iter_control_set_processOfMax

!=======================================================================
! setLocalError interface
!=======================================================================
                                                                                       
  subroutine iter_control_set_localError(self, localError)
                                                                                       
!    Set the maximum relative pointwise error
!    localError   maximum relative pointwise error for process
   
!    variable declarations
     implicit none
   
!    passed variables
     type(IterControl), intent(inout) :: self
     real(adqt),        intent(in)    :: localError
   
!    assertions
     tetonAssert(localError>=0,"Invalid local error")
   
!    reset local maximum relative error
     self % localError = localError 
   
!    assertions
     tetonAssert(self%localError>=0,"Invalid local error")
   
     return
  end subroutine iter_control_set_localError

!=======================================================================
! setGlobalError interface
!=======================================================================
          
  subroutine iter_control_set_globalError(self, globalError)
                                                    
!    Set the maximum relative pointwise error
!    globalError   global maximum relative pointwise error
   
!    variable declarations
     implicit none
   
!    passed variables
     type(IterControl), intent(inout) :: self
     real(adqt),        intent(in)    :: globalError
   
!    assertions
     tetonAssert(globalError>=0,"Invalid global error")
   
!    reset global maximum relative error
     self % globalError = globalError
   
!    assertions
     tetonAssert(self%globalError>=0,"Invalid global error")
   
     return
  end subroutine iter_control_set_globalError

!=======================================================================
! destruct interface
!=======================================================================

  subroutine iter_control_dtor(self)

!    Destruct the iteration control object

!    variable declarations
     implicit none

!    passed variables
     type(IterControl), intent(inout) :: self

!    destruct the iteration control
     self % epsilonPoint = zero
     self % localError   = zero
     self % globalError  = zero
     self % maxIter      = 0
     self % nIter        = 0
     self % nTotIter     = 0
     self % zoneOfMax    = 0
     self % processOfMax = 0

     return
  end subroutine iter_control_dtor

!=======================================================================
! external data access routines
!=======================================================================

!-----------------------------------------------------------------------
  function iter_control_get_epsilonPoint(self) result(epsilonPoint)

!    Return the relative pointwise error convergence criterion
!      epsilonPoint   relative pointwise error convergence criterion

!    variable declarations
     implicit none

!    passed variables
     type(IterControl), intent(in) :: self
     real(adqt)                    :: epsilonPoint

    epsilonPoint = self % epsilonPoint

!   assertions
    tetonAssert(epsilonPoint==self%epsilonPoint,"Invalid data access")

    return
  end function iter_control_get_epsilonPoint

!-----------------------------------------------------------------------
  function iter_control_get_maxIter(self) result(maxIter)

!   Return the maximum number of iterations
!     maxIter   maximum number of iterations

!   variable declarations
    implicit none

!   passed variables
    type(IterControl), intent(in) :: self
    integer                       :: maxIter

    maxIter = self%maxIter

!   assertions
    tetonAssert(maxIter==self%maxIter,"Invalid data access")

    return
  end function iter_control_get_maxIter

!-----------------------------------------------------------------------
  function iter_control_get_globalMaxIterTaken(self) result(maxIter)

!   Return the maximum number of iterations taken globally

!   variable declarations
    implicit none

!   passed variables
    type(IterControl), intent(in) :: self
    integer                       :: maxIter

    maxIter = self%globalMaxIterTaken

!   assertions
    tetonAssert(maxIter==self%globalMaxIterTaken,"Invalid data access")

    return
  end function iter_control_get_globalMaxIterTaken

!-----------------------------------------------------------------------
  function iter_control_get_nIter(self) result(nIter)

!   Return the required number of iterations
!     nIter   number of iterations

!   variable declarations
    implicit none

!   passed variables
    type(IterControl), intent(in) :: self
    integer                       :: nIter

    nIter = self%nIter

!   assertions
    tetonAssert(nIter==self%nIter,"Invalid data access")

    return
  end function iter_control_get_nIter

!-----------------------------------------------------------------------
  function iter_control_get_nTotIter(self) result(nTotIter)

!   Return the total number of iterations required
!     nTotIter   total number of iterations

!   variable declarations
    implicit none

!   passed variables
    type(IterControl), intent(in) :: self
    integer                       :: nTotIter

    nTotIter = self%nTotIter

!   assertions
    tetonAssert(nTotIter==self%nTotIter,"Invalid data access")

    return
  end function iter_control_get_nTotIter

!-----------------------------------------------------------------------
  function iter_control_get_zoneOfMax(self) result(zoneOfMax)

!   Return the zone number of the maximum relative pointwise error
!     zoneOfMax   zone of maximum relative pointwise error

!   variable declarations
    implicit none

!   passed variables
    type(IterControl), intent(in) :: self
    integer                       :: zoneOfMax

    zoneOfMax = self%zoneOfMax

!   assertions
    tetonAssert(zoneOfMax==self%zoneOfMax,"Invalid data access")

    return
  end function iter_control_get_zoneOfMax

!-----------------------------------------------------------------------
  function iter_control_get_processOfMax(self) result(processOfMax)
                                                                                       
!   Return the process number of the maximum relative pointwise error
!     processOfMax   zone of maximum relative pointwise error
                                                                                       
!   variable declarations
    implicit none
                                                                                       
!   passed variables
    type(IterControl), intent(in) :: self
    integer                       :: processOfMax
                                                                                       
    processOfMax = self%processOfMax
                                                                                       
!   assertions
    tetonAssert(processOfMax==self%processOfMax,"Invalid data access")
                                                                                       
    return
  end function iter_control_get_processOfMax

!-----------------------------------------------------------------------
  function iter_control_get_localError(self) result(localError)
                                                                                       
!   Return the maximum relative pointwise error for process
!     localError   maximum relative pointwise error
                                                                                       
!   variable declarations
    implicit none
                                                                                       
!   passed variables
    type(IterControl), intent(in) :: self
    real(adqt)                    :: localError 
                                                                                       
    localError = self%localError
                                                                                       
!   assertions
    tetonAssert(localError==self%localError,"Invalid data access")
                                                                                       
    return
  end function iter_control_get_localError

!-----------------------------------------------------------------------
  function iter_control_get_globalError(self) result(globalError)
                                                                                       
!   Return the global maximum relative pointwise error
!     globalError   maximum relative pointwise error
                                                                                       
!   variable declarations
    implicit none
                                                                                       
!   passed variables
    type(IterControl), intent(in) :: self
    real(adqt)                    :: globalError
                                          
    globalError = self%globalError
                           
!   assertions
    tetonAssert(globalError==self%globalError,"Invalid data access")
                                 
    return
  end function iter_control_get_globalError

end module iter_control_mod
