module EmissionUtils

  use iso_c_binding
  implicit none

  interface
    subroutine integrateBlackBodyGroups(T, k, Bnorm, numGroups, groupBounds, B ) &
                                    bind(C, name="NBB_integrateBlackBodyGroups")
       import
       real(c_double), VALUE, intent(in) :: T
       real(c_double), VALUE, intent(in) :: k
       real(c_double), VALUE, intent(in) :: Bnorm
       integer(c_int), VALUE, intent(in) :: numGroups
       real(c_double), intent(in) :: groupBounds(numGroups)
       real(c_double), intent(out) :: B(numGroups)
    end subroutine integrateBlackBodyGroups

    subroutine integrateBlackBodyAndDerivGroups(T, k, Bnorm, numGroups, groupBounds, B, dBdT ) &
                                          bind(C, name="NBB_integrateBlackBodyAndDerivGroups")
       import
       real(c_double), VALUE, intent(in) :: T
       real(c_double), VALUE, intent(in) :: k
       real(c_double), VALUE, intent(in) :: Bnorm
       integer(c_int), VALUE, intent(in) :: numGroups
       real(c_double), intent(in) :: groupBounds(numGroups)
       real(c_double), intent(out) :: B(numGroups)
       real(c_double), intent(out) :: dBdT(numGroups)
    end subroutine integrateBlackBodyAndDerivGroups

  end interface

end module EmissionUtils
