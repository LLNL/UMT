!*******************************************************************************
!  PrintEnergies - A subroutine that can be called for debugging purposes.     *
!    This subroutine prints out sums of varables of interest.                  *
!    Useful for comparing output from two versions of Teton.                   *
!*******************************************************************************

subroutine PrintEnergies(routine)

   use constant_mod
   use radconstant_mod
   use Geometry_mod
   use Material_mod
   use RadIntensity_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod
   use GroupSet_mod
   use io_mod
   use Size_mod
   use mpi_param_mod
   use mpif90_mod
   use Options_mod

!  variable declarations
   implicit none

!  passed variables
   character(*), intent(in) :: routine

   integer, parameter :: minimumVerbosityLevel = 3

   integer    :: c, c0, nCorner, zone, nzones
   integer    :: setID, nSets, groupSetID, nGroupSets
   integer    :: Angle, NumAngles 

   type(SetData) ,pointer  :: Set
   type(AngleSet),pointer  :: ASet
   type(GroupSet),pointer  :: GSet

   real(adqt) :: Rho, Cve, Vol
   real(adqt) :: geometryFactor
   real(adqt) :: quadwt

   real(adqt) :: radEnergy_psi
   real(adqt) :: radEnergy_phiTotal
   ! The difference between deltaEMat_denez and deltaEMat_tec gives a sense of how much temperature flooring occurred
   real(adqt) :: deltaEMat_denez ! Energy deposited computed from denez.  This is what gets returned to the host code
   real(adqt) :: deltaEMat_tec ! Change in material energy computed from the change in the corner temperature
   real(adqt) :: sourceEnergy

   ! All ranks need to do work, even if rank 0 is the only one printing
   if (Options%isRankZeroVerbose() < minimumVerbosityLevel) then
     return
   endif

   nzones     = Size% nzones
   nSets      = getNumberOfSets(Quad)
   nGroupSets = getNumberOfGroupSets(Quad)

   radEnergy_psi = zero
   radEnergy_phiTotal = zero
   deltaEMat_denez = zero
   deltaEMat_tec = zero
   sourceEnergy = zero

   ZoneLoop: do zone=1,nzones

     nCorner = Geom% numCorner(zone)
     c0      = Geom% cOffSet(zone)

     Rho = Mat%Rho(zone)
     Cve = Mat%Cve(zone)

     CornerLoop: do c=c0+1,c0+nCorner

       Vol = Geom%Volume(c)

       SetLoop: do setID=1,nSets

         Set       => Quad% SetDataPtr(setID) 
         ASet      => Quad% AngSetPtr(Set% angleSetID)
         NumAngles =  Set% NumAngles

         AngleLoop: do Angle=1,NumAngles

           quadwt = ASet% Weight(Angle)

           ! The SUM is over groups
           radEnergy_psi = radEnergy_psi + quadwt*Vol*SUM(Set%Psi(:,c,Angle))

         enddo AngleLoop

       enddo SetLoop

       GroupSetLoop: do groupSetID=1,nGroupSets

         GSet      => getGroupSetFromSetID(Quad, groupSetID)
         ! The SUM is over groups
         sourceEnergy = sourceEnergy + Vol*SUM(GSet%STotal(:,c))

       enddo GroupSetLoop

       ! The SUM is over groups
       radEnergy_phiTotal = radEnergy_phiTotal + Vol*SUM(Rad% PhiTotal(:,c))
       deltaEMat_tec = deltaEMat_tec + Rho*Cve*(Mat%Tec(c)-Mat%Tecn(c))*Vol

     enddo CornerLoop

     deltaEMat_denez = deltaEMat_denez + Rho*Geom%VolumeZone(zone)*Mat%denez(zone)

   enddo ZoneLoop

   geometryFactor       = getGeometryFactor(Size)
   radEnergy_psi        = geometryFactor*radEnergy_psi/speed_light
   radEnergy_phiTotal   = geometryFactor*radEnergy_phiTotal/speed_light
   deltaEMat_denez      = geometryFactor*deltaEMat_denez
   deltaEMat_tec        = geometryFactor*deltaEMat_tec
   sourceEnergy         = geometryFactor*sourceEnergy

   ! Print out local values for each rank if verbosity level is appropriate:
   if (Options%isRankVerbose() >= minimumVerbosityLevel) then
     ! E24.14E3 --> 24 total places, 14 digits after decimal, 3 spots for exponent?
     print *, "============================================================"
     print '(a,a,a,i6)', " in file/location: ", routine
     print '(a,i6,a,1pE24.14E3)', " Rank ", Size%myRankInGroup, " Radiation Energy computed from psi: ", radEnergy_psi
     print '(a,i6,a,1pE24.14E3)', " Rank ", Size%myRankInGroup, " Radiation Energy computed from Rad%phiTotal: ", radEnergy_phiTotal
     print '(a,i6,a,1pE24.14E3)', " Rank ", Size%myRankInGroup, " Material Energy Change computed from denez (what gets returned to host code): ", deltaEMat_denez
     print '(a,i6,a,1pE24.14E3)', " Rank ", Size%myRankInGroup, " Material Energy Change computed from change in Tec (Teton internal only): ", deltaEMat_tec
     print '(a,i6,a,1pE24.14E3)', " Rank ", Size%myRankInGroup, " Source Energy computed from GSet%STotal: ", sourceEnergy
     print *, "============================================================"
   endif

   call MPIAllReduce(radEnergy_psi, "sum", MY_COMM_GROUP)
   call MPIAllReduce(radEnergy_phiTotal, "sum", MY_COMM_GROUP)
   call MPIAllReduce(deltaEMat_denez, "sum", MY_COMM_GROUP)
   call MPIAllReduce(deltaEMat_tec, "sum", MY_COMM_GROUP)
   call MPIAllReduce(sourceEnergy, "sum", MY_COMM_GROUP)

   ! Print out local values for each rank if verbosity level is appropriate:
   if (Size%myRankInGroup == 0) then
     ! 1pE24.14E3 --> 24 total places, 14 digits after decimal, 3 spots for exponent?
     print *, "============================================================"
     print '(a,a,a,i6)', " in file/location: ", routine
     print '(a,1pE24.14E3)', " GLOBAL Radiation Energy computed from psi: ", radEnergy_psi
     print '(a,1pE24.14E3)', " GLOBAL Radiation Energy computed from Rad%phiTotal: ", radEnergy_phiTotal
     print '(a,1pE24.14E3)', " GLOBAL Material Energy Change computed from denez (what gets returned to host code): ", deltaEMat_denez
     print '(a,1pE24.14E3)', " GLOBAL Material Energy Change computed from change in Tec (Teton internal only): ", deltaEMat_tec
     print '(a,1pE24.14E3)', " GLOBAL Source Energy computed from GSet%STotal: ", sourceEnergy
     print *, "============================================================"
   endif

   return

end subroutine PrintEnergies
