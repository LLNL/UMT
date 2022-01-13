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
   use ZoneData_mod

   implicit none
   ! Input / Output arguments
   logical(C_BOOL), intent(in)    :: killOnBad
   integer(C_INT),  intent(in)    :: complaintLevel
   integer(C_INT),  intent(in)    :: numCatsToCheck
   integer(C_INT),  intent(in)    :: arrayOfCatsToCheck(numCatsToCheck)
   integer(C_INT),  intent(inout) :: numBadInputCategories

   !  Local pointers
   type(ZoneData), pointer        :: ZonePtr => NULL()

   ! different numbers may be more appropriate, but these are physically meaninful
   ! Note that opacities are now clipped at 1e50 in teton_setopacity
   real(adqt), parameter          :: CVEMIN    = adqtEpsilon
   real(adqt), parameter          :: VOLUMEMIN = adqtEpsilon
! TODO: Make SIGMIN/sigmaWarn tunable by the host codes? See TETON-131.
   real(adqt), parameter          :: SIGMIN    = adqtEpsilon
   real(adqt), parameter          :: RHOMIN    = adqtEpsilon
   real(adqt), parameter          :: NEZMIN    = adqtEpsilon
   real(adqt)                     :: TEMPMIN   
   
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

   ! maybe this is faster than comparing to an integer every time?
   logical(C_BOOL) :: printEvery
   printEvery = (2 == complaintLevel)
   
   TEMPMIN = Size% tfloor
   
   !  Constants
   nGroups        = Size% ngr
   nZones     = Size% nzones
   nCornersTotal   = Size% ncornr

   
   ! Category format string
100 FORMAT('Bad ',A20,' Data. ',I10, ' Entries out of' , I12, ' were bad')

   ! sigma format string
200 FORMAT('Zone: ',I8, ' Group: ',I4,' Sigma',A20,' is bad: ',ES20.8)

   ! Zone material property
300 FORMAT('Zone: ' , I8, 1X , A20, ' is bad: ',ES20.8)

   ! Corner material property
400 FORMAT('Zone: ', I8, 'Corner: ',I2,1X,A20, ' is bad: ', ES20.8)
   
   ! assume all data is good to start
   numBadInputCategories = 0


   do cats=1,numCatsToCheck
      select case ( arrayOfCatsToCheck(cats) )
      case(1)
         ! **************************************
         ! sigs
         ! **************************************
         numBadEntries = 0
         numEntries = nZones * nGroups
         categoryName = "SigmaS"
         do zone=1, nZones

            do g=1, nGroups
               if( Mat% sigs(g,zone) < SIGMIN .or. Mat% sigs(g,zone) > sigmaWarn) then
                  numBadEntries = numBadEntries + 1
                  if( printEvery ) then
                     WRITE(nout , 200) zone,g,categoryName,Mat%sigs(g,zone)
                  endif
               endif
            enddo

         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) categoryName, numBadEntries, numEntries
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

            do g=1, nGroups
               if( Mat% siga(g,zone) < SIGMIN .or. Mat% siga(g,zone) > sigmaWarn) then
                  numBadEntries = numBadEntries + 1
                  if( printEvery ) then
                     WRITE(nout , 200) zone,g,categoryName,Mat%siga(g,zone)
                  endif
               endif
            enddo

         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) categoryName, numBadEntries, numEntries
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
            if( Mat% cve(zone) < CVEMIN) then
               numBadEntries = numBadEntries + 1
               if( printEvery ) then
                  WRITE(nout , 300) zone, categoryName, Mat% cve(zone)
               endif
            endif
         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) categoryName, numBadEntries, numEntries
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
            if( Mat% rho(zone) < RHOMIN) then
               numBadEntries = numBadEntries + 1
               if( printEvery ) then
                  WRITE(nout , 300) zone, categoryName, Mat% rho(zone)
               endif
            endif
         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) categoryName, numBadEntries, numEntries
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
            if( Mat% tez(zone) < TEMPMIN) then
               numBadEntries = numBadEntries + 1
               if( printEvery ) then
                  WRITE(nout , 300) zone,categoryName,Mat% tez(zone)
               endif
            endif
         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) categoryName, numBadEntries, numEntries
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
            if( Mat% nez(zone) < NEZMIN) then
               numBadEntries = numBadEntries + 1
               if( printEvery ) then
                  WRITE(nout , 300) zone,categoryName,Mat% nez(zone)
               endif
            endif
         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) categoryName, numBadEntries, numEntries
            endif
         endif
      case (7)
      ! **************************************
      ! Corner volumes
      ! **************************************
         numBadEntries = 0
         numEntries = nCornersTotal
         categoryName = "Corner Volume"
         do zone=1, nZones
            ZonePtr => getZoneData(Geom, zone)
            nCorner = ZonePtr% nCorner
            c0      = ZonePtr% c0
            do c=1,nCorner
               if( Geom% Volume(c0+c)  < VOLUMEMIN) then
                  numBadEntries = numBadEntries + 1
                  if( printEvery ) then
                     WRITE(nout , 400) zone,c,categoryName,Geom% Volume(c0+c) 
                  endif
               endif
            enddo
         enddo
         if( numBadEntries > 0) then
            numBadInputCategories = numBadInputCategories + 1
            if (complaintLevel == 1) then
               WRITE( nout , 100) categoryName, numBadEntries, numEntries
            endif
         endif
      end select
   enddo
   if(killOnBad .AND. (numBadInputCategories > 0)) then
      call f90fatal("Bad Teton input data found")
   endif
   
   return
 end subroutine checkInputSanity



