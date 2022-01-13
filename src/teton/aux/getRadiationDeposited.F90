!***********************************************************************
!                         Version 0: 04/06 PFN                         *
!                                                                      *
!    getRadiationDeposited -  Called from host to get the              *
!                             radiation deposited in the matter        *
!                             during the radiation step and the        *
!                             radiation tempewrature.                  *
!                                                                      *
!***********************************************************************

   subroutine getRadiationDeposited(zone, radEnergyDeposited,  &
                                    radiationTemperature) &
                                    BIND(C,NAME="teton_getradiationdeposited")

   USE ISO_C_BINDING
   use kind_mod
   use Material_mod

   implicit none

!  Arguments

   integer(C_INT),   intent(in)     :: zone
   real(C_DOUBLE),   intent(out)    :: radEnergyDeposited 
   real(C_DOUBLE),   intent(out)    :: radiationTemperature


   radEnergyDeposited   = Mat% denez(zone)
   radiationTemperature = Mat% trz(zone)


   return
   end subroutine getRadiationDeposited 
