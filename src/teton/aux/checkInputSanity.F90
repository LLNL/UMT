!***********************************************************************
!                        Created:  04/2018, PGM                        *
!                                                                      *
!     checkInputSanity: Have teton check all hydro code inputs for     *
!                       sanity. For example,positive densities and     *
!                       volumes.                                       *
!                                                                      *
!***********************************************************************
 
subroutine checkInputSanity(killOnBad,   &
                            complaintLevel, &
                            numCatsToCheck,   &
                            arrayOfCatsToCheck, &                  
                            numBadInputCategories) &
                            BIND(C,NAME="teton_checkinputsanity")

   use ISO_C_BINDING
   use kind_mod
   use constant_mod
   use radconstant_mod  ! provides sigmaWarn value
   use io_mod           ! provides nout aka file number of stdout
   use Size_mod
   use Geometry_mod
   use Material_mod

   use ieee_arithmetic
   use Datastore_mod, only : theDatastore
   use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                             stdout=>output_unit, &
                                             stderr=>error_unit

   implicit none
   ! Input / Output arguments
   logical(C_BOOL), intent(in)    :: killOnBad
   integer(C_INT),  intent(in)    :: complaintLevel
   integer(C_INT),  intent(in)    :: numCatsToCheck
   integer(C_INT),  intent(in)    :: arrayOfCatsToCheck(numCatsToCheck)
   integer(C_INT),  intent(inout) :: numBadInputCategories

   ! different numbers may be more appropriate, but these are physically meaninful
   real(adqt) :: CVEMIN    = adqtTiny
   real(adqt) :: VOLUMEMIN = adqtTiny
   real(adqt) :: SIGMIN    = adqtTiny
   real(adqt) :: SIGMAX    = sigmaWarn
   real(adqt) :: RHOMIN    = adqtTiny
   real(adqt) :: NEZMIN    = adqtTiny
   real(adqt) :: TEMPMIN
   
   ! Indices
   integer         :: g, zone, c, c0, cats

   ! number of corners in a particular zone
   integer        :: nCorner
   
   ! sizes of things we want to check
   integer         :: nGroups
   integer         :: nZones
   integer         :: nCornersTotal

   ! number of entries of a given category
   integer         :: numBadEntries
   ! total number of values in that category
   integer         :: numEntries
   
   ! Write out category names in this buffer
   character(len=20)        :: categoryName

   logical(kind=1) :: isBadZone = .false.
   
   TEMPMIN = Size% tfloor*0.999999_adqt

   ! Check datastore for optional values:
   CVEMIN    = theDatastore%fetchIfExists("options/iteration/sanitizer/min_cve", CVEMIN)
   VOLUMEMIN = theDatastore%fetchIfExists("options/iteration/sanitizer/min_volume", VOLUMEMIN)
   SIGMIN    = theDatastore%fetchIfExists("options/iteration/sanitizer/min_opacity", SIGMIN)
   SIGMAX    = theDatastore%fetchIfExists("options/iteration/sanitizer/max_opacity", SIGMAX)
   RHOMIN    = theDatastore%fetchIfExists("options/iteration/sanitizer/min_density", RHOMIN)
   NEZMIN    = theDatastore%fetchIfExists("options/iteration/sanitizer/min_electron_density", NEZMIN)
   TEMPMIN   = theDatastore%fetchIfExists("options/iteration/sanitizer/min_temperature", TEMPMIN)
   
   !  Constants
   nGroups        = Size% ngr
   nZones         = Size% nzones
   nCornersTotal  = Size% ncornr

   ! Category format string
100 FORMAT('Rank :', I6,' Bad ',A20,' Data. ',I10,' Entries out of' ,I12,' were bad')

   ! Zone and group-dependent format string
200 FORMAT('Rank: ',I6,' Zone: ',I8,' Group: ',I4,1X,A20,' is bad: ',ES20.8)

   ! Zone material property
300 FORMAT('Rank: ',I6,' Zone: ',I8,1X,A20,' is bad: ',ES20.8)

   ! Corner material property
400 FORMAT('Rank: ',I6,' Zone: ',I8,' Corner: ',I2,1X,A20,' is bad: ',ES20.8)
   
   ! assume all data is good to start
   numBadInputCategories = 0

   do cats=1,numCatsToCheck
      select case ( arrayOfCatsToCheck(cats) )
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
      case (1)
         ! **************************************
         ! sigs
         ! **************************************
         numBadEntries = 0
         numEntries = nZones * nGroups
         categoryName = "SigmaS"
         do zone=1, nZones
            isBadZone = .false.

            do g=1, nGroups
               if( .not. ieee_is_finite(Mat%sigs(g,zone)) .or. Mat%sigs(g,zone) < zero .or. Mat%sigs(g,zone) > SIGMAX) then
                  numBadEntries = numBadEntries + 1
                  if( complaintLevel > 2 ) then
                     call printDetailedBadZoneInfo(zone,categoryName)
                     exit
                  else if( complaintLevel > 1 ) then
                     isBadZone = .true.
                     WRITE(nout , 200) Size%myRankInGroup,zone,g,categoryName,Mat%sigs(g,zone)
                  endif
               endif
            enddo

            if (isBadZone) then
               call printBadZoneInfo(zone,categoryName)
            endif
         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) Size%myRankInGroup,categoryName,numBadEntries,numEntries
            endif
         endif
      case (2)
      ! **************************************
      ! siga
      ! **************************************
         numBadEntries = 0
         numEntries = nZones * nGroups
         categoryName = "SigmaA"
         do zone=1, nZones
            isBadZone = .false.

            do g=1, nGroups
               if( .not. ieee_is_finite(Mat%sigs(g,zone)) .or. Mat%siga(g,zone) < SIGMIN .or. Mat%siga(g,zone) > SIGMAX) then
                  numBadEntries = numBadEntries + 1
                  if( complaintLevel > 2 ) then
                     call printDetailedBadZoneInfo(zone,categoryName)
                     exit
                  else if( complaintLevel > 1 ) then
                     isBadZone = .true.
                     WRITE(nout , 200) Size%myRankInGroup,zone,g,categoryName,Mat%siga(g,zone)
                  endif
               endif
            enddo

            if (isBadZone) then
               call printBadZoneInfo(zone,categoryName)
            endif
         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) Size%myRankInGroup,categoryName, numBadEntries, numEntries
            endif
         endif
      case (3)
      ! **************************************
      ! cve
      ! **************************************
         numBadEntries = 0
         numEntries = nZones
         categoryName = "Eff. Specific Heat"
         do zone=1, nZones
            if( .not. ieee_is_finite(Mat%cve(zone)) .or. Mat%cve(zone) < CVEMIN) then
               numBadEntries = numBadEntries + 1
               if( complaintLevel > 2 ) then
                  call printDetailedBadZoneInfo(zone,categoryName)
                  exit
               else if( complaintLevel > 1 ) then
                  WRITE(nout , 300) Size%myRankInGroup,zone, categoryName, Mat% cve(zone)
                  call printBadZoneInfo(zone,categoryName)
               endif
            endif
         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) Size%myRankInGroup,categoryName, numBadEntries, numEntries
            endif
         endif
      case (4)
      ! **************************************
      ! rho
      ! **************************************
         numBadEntries = 0
         numEntries = nZones
         categoryName = "Density"
         do zone=1, nZones
            if( .not. ieee_is_finite(Mat%rho(zone)) .or. Mat%rho(zone) < RHOMIN) then
               numBadEntries = numBadEntries + 1
               if( complaintLevel > 2 ) then
                  call printDetailedBadZoneInfo(zone,categoryName)
                  exit
               else if( complaintLevel > 1 ) then
                  WRITE(nout , 300) Size%myRankInGroup,zone, categoryName, Mat% rho(zone)
                  call printBadZoneInfo(zone,categoryName)
               endif
            endif
         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) Size%myRankInGroup,categoryName, numBadEntries, numEntries
            endif
         endif
      case (5)
      ! **************************************
      ! tez
      ! **************************************
         numBadEntries = 0
         numEntries = nZones
         categoryName = "Electron Temp"
         do zone=1, nZones
            if( .not. ieee_is_finite(Mat%tez(zone)) .or. Mat%tez(zone) < TEMPMIN) then
               numBadEntries = numBadEntries + 1
               if( complaintLevel > 1 ) then
                  WRITE(nout , 300) Size%myRankInGroup,zone,categoryName,Mat% tez(zone)
                  call printBadZoneInfo(zone,categoryName)
               endif
            endif
         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) Size%myRankInGroup,categoryName, numBadEntries, numEntries
            endif
         endif
      case (6)
      ! **************************************
      ! nez
      ! **************************************
         numBadEntries = 0
         numEntries = nZones
         categoryName = "Electron Density"
         do zone=1, nZones
            if( .not. ieee_is_finite(Mat%nez(zone)) .or. Mat%nez(zone) < NEZMIN) then
               numBadEntries = numBadEntries + 1
               if( complaintLevel > 2 ) then
                  call printDetailedBadZoneInfo(zone,categoryName)
                  exit
               else if( complaintLevel > 1 ) then
                  WRITE(nout , 300) Size%myRankInGroup,zone,categoryName,Mat% nez(zone)
                  call printBadZoneInfo(zone,categoryName)
               endif
            endif
         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) Size%myRankInGroup,categoryName, numBadEntries, numEntries
            endif
         endif
#endif
      case (7)
      ! **************************************
      ! Corner volumes
      ! **************************************
         numBadEntries = 0
         numEntries = nCornersTotal
         categoryName = "Corner Volume"
         do zone=1, nZones
            isBadZone = .false.
            nCorner = Geom% numCorner(zone) 
            c0      = Geom% cOffSet(zone) 
            do c=1,nCorner
               if( .not. ieee_is_finite(Geom%Volume(c0+c)) .or. Geom%Volume(c0+c) < VOLUMEMIN) then
                  isBadZone = .true.
                  numBadEntries = numBadEntries + 1
                  if( complaintLevel > 2 ) then
                     call printDetailedBadZoneInfo(zone,categoryName)
                     exit
                  else if( complaintLevel > 1 ) then
                     WRITE(nout , 400) Size%myRankInGroup,zone,c,categoryName,Geom% Volume(c0+c)
                  endif
               endif
            enddo
            if (isBadZone) then
               call printBadZoneInfo(zone,categoryName)
            endif
         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) Size%myRankInGroup,categoryName, numBadEntries, numEntries
            endif
         endif
      end select
   enddo
   if(killOnBad .AND. (numBadInputCategories > 0)) then
      flush(stdout)
      call f90fatal("Bad Teton input data found")
   endif
   
   return
 end subroutine checkInputSanity

subroutine printBadZoneInfo(zone, categoryName)

  use io_mod           ! provides nout aka file number of stdout
  use Size_mod
  use Geometry_mod
  use Material_mod

  integer, intent(in)            :: zone
  character(len=20), intent(in)  :: categoryName

500 FORMAT('Rank:', I5, ', Bad ', A, ' data was found in zone ', I8, ', centered at:', 1ES15.7, ' with rho = ', ES10.3, ', Te =', ES10.3, ', and Tr =', ES10.3)
600 FORMAT('Rank:', I5, ', Bad ', A, ' data was found in zone ', I8, ', centered at:', 2ES15.7, ' with rho = ', ES10.3, ', Te =', ES10.3, ', and Tr =', ES10.3)
700 FORMAT('Rank:', I5, ', Bad ', A, ' data was found in zone ', I8, ', centered at:', 3ES15.7, ' with rho = ', ES10.3, ', Te =', ES10.3, ', and Tr =', ES10.3)
  if (Size%ndim == 1) then
      WRITE(nout, 500) Size%myRankInGroup, TRIM(categoryName), zone, getZoneCenter(Geom,zone), &
                       Mat%rho(zone), Mat%Tez(zone), Mat%Trz(zone)
  else if (Size%ndim == 2) then
      WRITE(nout, 600) Size%myRankInGroup, TRIM(categoryName), zone, getZoneCenter(Geom,zone), &
                       Mat%rho(zone), Mat%Tez(zone), Mat%Trz(zone)
  else
      WRITE(nout, 700) Size%myRankInGroup, TRIM(categoryName), zone, getZoneCenter(Geom,zone), &
                       Mat%rho(zone), Mat%Tez(zone), Mat%Trz(zone)
  endif

  return
end subroutine printBadZoneInfo

subroutine printDetailedBadZoneInfo(zone, categoryName)

  use io_mod           ! provides nout aka file number of stdout
  use Size_mod
  use Geometry_mod
  use Material_mod

  integer, intent(in)            :: zone
  character(len=20), intent(in)  :: categoryName

  integer :: c0, nCorner, c

  call printBadZoneInfo(zone, categoryName)

  nCorner = Geom% numCorner(zone)
  c0      = Geom% cOffSet(zone)

  write(nout, *) 'cve = ', Mat%cve(zone), ', nez = ', Mat%nez(zone)
  write(nout, *) 'SigmaA = ', Mat%sigA(:,zone)
  write(nout, *) 'SigmaS = ', Mat%sigS(:,zone)
  write(nout, *) 'Corner volumes = ', Geom%Volume(c0+1:c0+nCorner)
  write(nout, *) 'Corner temperatures = ', Mat%Tec(c0+1:c0+nCorner)
  write(nout, *) 'Zone vertices: '
  do c=1,nCorner
    write(nout, *) Geom%px(:,c0+c)
  enddo


  return
end subroutine printDetailedBadZoneInfo

