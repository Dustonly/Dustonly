!================================================================
! $Id: src_dust.f90,v 0 2018/06/29 Faust $
!================================================================
!                                                               !
!   MUSCAT : MUltiScale Chemistry Transport Code                !
!                                                               !
!   Institute for Tropospheric Research                         !
!   Permoserstr. 15, D-04303 Leipzig, Germany                   !
!                                                               !
!---------------------------------------------------------------!
!   Contact: Ralf Wolke (wolke@tropos.de)                       !
!================================================================
!
MODULE src_dust
!---------------------------------------------------------------------
! Description:
! New module for dust emissions in muscat
! The aim is to simplify the dust Code
! create a simple structur: 1 Modul containing a organizes routine, logical subroutines called in organize
! remove outdated dependencies and unused input
! rename unlogical named vars
! translate into more modern Fortran
!---------------------------------------------------------------------
!
! Current Code Owner:Institut für Troposphärenforschung e.V. Leipzig, Matthias Faust
! email:  faust@tropos.de
! History:
! Version      Date       Name
! ----------   ---------- ----
! V0           2018-06-29 Start Coding
! V0.1         2018-07-12 add read_ascii
! V0.2         2018_07_16 add copy2block
! V0.3         2018_07_17 add read_nc for Vegetation
!
!
! Code Description:
! Language: Fortran 90.
! Software Standards: COSMO Standards for Source Code Development
!---------------------------------------------------------------------


  ! Modules MUSCAT
  USE  mo_dust
  ! USE  dust_org

#ifndef OFFLINE
  USE  partition
  USE  sub_geo
  USE  mo_ctm , ONLY: nalloc         & ! Size of Allocated Arrays
                     ,SurfRef        & ! Resolution of Landuse Data
                     ,SurfLevel

  USE sub_block, ONLY: tracer_name   & ! name of tracers
                      ,nt              ! number of tracer

  ! Modules COSMO
  USE data_modelconfig, ONLY: &
    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction
    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    startlon_tot, & ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
    startlat_tot    ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)

  USE data_parallel, ONLY:  &
    my_cart_id,             &
    icomm_cart,             & ! communicator for the virtual cartesian topology
    imp_reals,              & ! datatypes for MPI
    num_compute               ! number of compute PEs

  ! Modules External
  USE  mo_mpi
#endif



  ! Global Variables


  ! Local Variables
  PRIVATE

  ! subroutines for calling outside the module
  PUBLIC :: organize_dust,quick_nc,quick_ascii




  ! subroutines in this module
  CONTAINS
  ! organize_dust
  ! read_ascii
  ! copy2block

  !+ organize dust
  !---------------------------------------------------------------------
  SUBROUTINE organize_dust(yaction,subdomain,flux)
  !---------------------------------------------------------------------
  ! Description:
  !   This subroutine organize the dust emission sheme in muscat
  !
  !---------------------------------------------------------------------

    ! Modules
#ifndef OFFLINE
    USE data_runcontrol,    ONLY : hstop
#else
    USE offline_org
#endif

    IMPLICIT NONE

    !- Begin of header -----------------------------------------------------
    ! Parameter:
    CHARACTER(LEN=*), INTENT(IN)            :: &
      yaction ! action to be performed

    TYPE (rectangle),  INTENT(IN)  :: &
      subdomain

    REAL(8), OPTIONAL, INTENT(INOUT)        :: &
        flux(ntz,subdomain%nty,subdomain%ntx,nt)

    ! Variables vor Initialization
    INTEGER        :: &
      ib1,            &
      js,             &
      igx0, igx1,     &
      igy0, igy1,     &
      ix0,  ix1,      &
      iy0,  iy1,      &
      i,filenum!,      & ! loop
      ! ierr              ! error code for mpi

    INTEGER :: ifind  ! muscat funktion /INIT/ifind.f90
    EXTERNAL   ifind  ! search indices for species


    CHARACTER(20)  :: &
      string
    CHARACTER(120) :: &
      ifile(10)
    CHARACTER(120) :: &
      filename

    INTEGER ::    &
      ifile_num,  &
      ifile_dim(10)


    INTEGER  ::  &
    dimveg,      & ! time dimension of vegetation  = 12 if monthly, = nSimuDays if daily
    dim            ! dimensions of input data



    REAL(8), POINTER   ::  &
      copy2d(:,:),         & ! Pointer to copy data into block struckture
      copy3d(:,:,:)


    INTEGER        :: ierr   ! Error code
    CHARACTER(120) :: yerr   ! Error message



    !- End of header -----------------------------------------------------


    !---------------------------------------------------------------------
    ! Subroutine Body
    !---------------------------------------------------------------------

    IF (lddebug) PRINT*, 'Enter organize_dust, yaction=',yaction

    ierr = 0

    ! ------------------------------------
    ! +-+-+- Section 1 Init + Input -+-+-+
    ! ------------------------------------
    IF (yaction == "init") THEN

      ifile_num = 0

      ! +-+-+- Sec 1.1 Check -+-+-+
      ! Check consistency of namelist settings


      ! dust_scheme need right values
      IF (dust_scheme < 0 .OR. dust_scheme > 2) THEN
        ierr = 100001
        yerr = 'wrong value for dust_scheme'
        PRINT*,'ERROR    src_dust "init" '
        PRINT*,'         #',ierr
        PRINT*,'         ',yerr
        STOP
      END IF

      ! threshold_scheme need right values
      IF (threshold_scheme < 0 .OR. threshold_scheme > 2) THEN
        ierr = 100001
        yerr = 'wrong value for threshold_scheme'
        PRINT*,'ERROR    src_dust "init" '
        PRINT*,'         #',ierr
        PRINT*,'         ',yerr
        STOP
      END IF

      ! soiltypeFile always necessary
      IF (TRIM(soiltypeFile) == 'without') THEN
        ierr = 100002
        yerr = 'SoilTypeFile is missing'
        PRINT*,'ERROR    src_dust "init" '
        PRINT*,'         #',ierr
        PRINT*,'         ',yerr
        STOP
      ELSE
        ifile_num = ifile_num + 1
        ifile(ifile_num) = 'soil'
        IF (soilmaptype == 1) ifile_dim(ifile_num) = 1
        IF (soilmaptype == 2) ifile_dim(ifile_num) = 3
      END IF


      IF (soilmaptype == 1) THEN
        nmode = 4
      ELSEIF (soilmaptype == 2) THEN
        nmode = 3
      ELSE
        ierr = 1000023
        yerr = 'wrong value for soilmaptype'
        PRINT*,'ERROR    src_dust "init" '
        PRINT*,'         #',ierr
        PRINT*,'         ',yerr
        STOP
      END IF


      ! psrcType need right values
      IF (psrcType < 0 .OR. psrcType > 2) THEN
        ierr = 100003
        yerr = 'wrong value for psrcType'
        PRINT*,'ERROR    src_dust "init" '
        PRINT*,'         #',ierr
        PRINT*,'         ',yerr
        STOP
      END IF

      ! psrcFile always necessary
      IF (psrcType > 0) THEN
        IF (TRIM(psrcFile) == 'without') THEN
          ierr = 100004
          yerr = 'psrcFile is missing'
          PRINT*,'ERROR    src_dust "init" '
          PRINT*,'         #',ierr
          PRINT*,'         ',yerr
          STOP
        ELSE
          ifile_num = ifile_num + 1
          ifile(ifile_num) = 'source'
          ifile_dim(ifile_num) = 1
        END IF
      END IF

      ! z0File maybe necessary
      IF (lwithz0) THEN
        IF (TRIM(z0File) == 'without') THEN
          ierr = 100005
          yerr = 'z0File is missing'
          PRINT*,'ERROR    src_dust "init" '
          PRINT*,'         #',ierr
          PRINT*,'         ',yerr
          STOP
        ELSE
          ifile_num = ifile_num + 1
          ifile(ifile_num) = 'z0'
          ifile_dim(ifile_num) = 1
        END IF
      ELSE
        z0File = 'without'
      ENDIF

      ! veg_scheme ==2 need biomes as input
      IF (veg_scheme == 2 .AND. .NOT. lwithbiom) THEN
        PRINT*,'WARNING  src_dust "init" '
        PRINT*,'         veg_scheme == 2 needs biomes as input but,'
        PRINT*,'         lwithbiom = .FALSE.'
        PRINT*,' '
        PRINT*,'         CHANGING SETTINGS'
        PRINT*,'         lwithbiom = .TRUE.'
        PRINT*,' '
        lwithbiom = .TRUE.
      END IF

      ! biomFile maybe necessary
      IF (lwithbiom) THEN
        IF (TRIM(biomeFile) == 'without') THEN
          ierr = 100006
          yerr = 'biomFile is missing'
          PRINT*,'ERROR    src_dust "init" '
          PRINT*,'         #',ierr
          PRINT*,'         ',yerr
          STOP
        ELSE
          ifile_num = ifile_num + 1
          ifile(ifile_num) = 'biom'
          ifile_dim(ifile_num) = 1
        END IF
      ELSE
        biomeFile = 'without'
      ENDIF

      ! veg_scheme need right values
      IF (veg_scheme < 0 .OR. veg_scheme > 2) THEN
        ierr = 100007
        yerr = 'wrong value for veg_scheme'
        PRINT*,'ERROR    src_dust "init" '
        PRINT*,'         #',ierr
        PRINT*,'         ',yerr
        STOP
      END IF

      ! any vegFile maybe necessary
      IF (veg_scheme > 0) THEN
        IF (TRIM(vegmonFile) == 'without' .AND. TRIM(vegdayFile) == 'without') THEN
          ierr = 100008
          yerr = 'vegmonFile or vegdayFile is missing'
          PRINT*,'ERROR    src_dust "init" '
          PRINT*,'         #',ierr
          PRINT*,'         ',yerr
          STOP
        END IF
        ! if vegmonFile and vegdayFile are specified then vegdayFile is prefered
        IF (TRIM(vegmonFile) /= 'without' .AND. TRIM(vegdayFile) /= 'without') THEN
          vegmonFile = 'without'
        ENDIF
        ! set flag for daily vegetation
        lvegdaily = .FALSE.
        IF (TRIM(vegdayFile) /= 'without') lvegdaily = .TRUE.
      ELSE
        vegmonFile = 'without'
        vegdayFile = 'without'
      ENDIF

      ! Init of vegitation time dimension
      !   climatological or daily data
      IF (veg_scheme > 0 ) THEN
        ! If there is no daily data monthly data will be used
        dimveg = 12
        ! when daily data is available calc max number of days in the simulation
        IF (lvegdaily) dimveg = CEILING(hstop/24) + 1 ! hardcoding at this point my be helpfull for some tests 366 370!1096!1858!366!
      ELSE
        ! if no vegetation scheme is used then dimveg = 1
        dimveg = 1
      END IF

      IF (lvegdaily) THEN
        ifile_num = ifile_num + 1
        ifile(ifile_num) = 'vegday'
        ifile_dim(ifile_num) = dimveg
      ELSEIF( .NOT. lvegdaily .AND. TRIM(vegmonFile) /= 'without') THEN
        ifile_num = ifile_num + 1
        ifile(ifile_num) = 'vegmon'
        ifile_dim(ifile_num) = dimveg
      END IF


      ! set flag for minimum vegetation file
      lvegmin = .FALSE.
      IF (TRIM(vegminFile) /= 'without') lvegmin = .TRUE.
      IF (lvegmin) THEN
        ifile_num = ifile_num + 1
        ifile(ifile_num) = 'vegmin'
        ifile_dim(ifile_num) = 1
      END IF

      ! soil moisture need input stream
      ! at the moment soil moisture only for the offline version
#ifdef OFFLINE
      IF (moist_scheme > 0) THEN
        IF (TRIM(moistFile) == 'without') THEN
          ierr = 100008
          yerr = 'moistFile is missing'
          PRINT*,'ERROR    src_dust "init" '
          PRINT*,'         #',ierr
          PRINT*,'         ',yerr
          STOP
        ELSE
          ifile_num = ifile_num + 1
          ifile(ifile_num) = 'moist'
          ifile_dim(ifile_num) = lasttstep+1-firsttstep
        END IF
      ELSE
        moistFile = 'without'
      END IF
#else
      moist_scheme=0
#endif



      ! +-+-+- Sec 1.2 Init -+-+-+

      ! - Init for feedback?
      ! Set Indices of Dust Tracer
      ! necessary for AOD Calc?
      DO js=1,DustBins
#ifndef OFFLINE
        string = TRIM(ADJUSTL(DustName(js)))
        DustInd(js)  = ifind(nt, string, tracer_name)   ! ifind -> muscat funktion /INIT/ifind.f90
#else
        DustInd(js) = js
#endif
        IF (DustInd(js) <= 0)  THEN
          WRITE(*,8010)  string
          STOP  'Dust_Init: Error in Input Data !!'
        END IF
      END DO


8010  format(1x,50('+')//    &
      'Dust_Init: Dust name ',a20,' not included as tracer!' / 1x,50('+'))

      ! Define corners of subdomain?
      IF (SurfLevel >= 0)  THEN
        SurfRef = 2**SurfLevel
        igx0 = domain%igx0 * SurfRef
        igy0 = domain%igy0 * SurfRef
        igx1 = domain%igx1 * SurfRef
        igy1 = domain%igy1 * SurfRef
      ELSE IF (SurfLevel < 0)  THEN
        SurfRef = 2**(-SurfLevel)
        igx0 = domain%igx0 / SurfRef
        igy0 = domain%igy0 / SurfRef
        igx1 = domain%igx1 / SurfRef
        igy1 = domain%igy1 / SurfRef
      END IF




      ! - Init of 'dust' type 'dust_subdomain', see data_dust

      ! allocate datatype 'dust' as array with nb number of blocks
      ALLOCATE(dust(nb))
      ! ALLOCATE(dust_ini(nb))
      ! ALLOCATE(dust_flux(nb))


      ! Bolck strukture of muscat, loop over blocks
      DO ibLoc=1,nbLoc
#ifndef OFFLINE
        ib1 = LocGlob(ibLoc)
#else
        ib1 = 1
#endif
        ! allocate arrays in dust for each block
        ALLOCATE(dust(ib1)%soilprop(decomp(ib1)%iy0+1:decomp(ib1)%iy1,      &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,1:SoilFrac,1:SoilNumb))
        ALLOCATE(dust(ib1)%lai(decomp(ib1)%iy0+1:decomp(ib1)%iy1,           &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,1:SoilFrac,1:dimveg))
        ALLOCATE(dust(ib1)%vegmin(decomp(ib1)%iy0+1:decomp(ib1)%iy1,        &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,1))
        ALLOCATE(dust(ib1)%alpha(decomp(ib1)%iy0+1:decomp(ib1)%iy1,         &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,1:SoilFrac))
        ALLOCATE(dust(ib1)%c_eff(decomp(ib1)%iy0+1:decomp(ib1)%iy1,         &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,1:SoilFrac))
        ALLOCATE(dust(ib1)%umin2(decomp(ib1)%iy0+1:decomp(ib1)%iy1,         &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,1:SoilFrac))
        ALLOCATE(dust(ib1)%lai_eff(decomp(ib1)%iy0+1:decomp(ib1)%iy1,       &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,1:SoilFrac,1:dimveg))
        ALLOCATE(dust(ib1)%w_str(decomp(ib1)%iy0+1:decomp(ib1)%iy1,         &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1))
        ALLOCATE(dust(ib1)%d_emis(decomp(ib1)%iy0+1:decomp(ib1)%iy1,        &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,1:nt))

        ! Allocate dust_ini
        ALLOCATE(dust(ib1)%biome(decomp(ib1)%iy0+1:decomp(ib1)%iy1,      &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1))
        ALLOCATE(dust(ib1)%cult(decomp(ib1)%iy0+1:decomp(ib1)%iy1,      &
                decomp(ib1)%ix0+1:decomp(ib1)%ix1))
        ALLOCATE(dust(ib1)%veg(decomp(ib1)%iy0+1:decomp(ib1)%iy1,       &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,dimveg))
        ALLOCATE(dust(ib1)%vmoist(decomp(ib1)%iy0+1:decomp(ib1)%iy1,       &
                decomp(ib1)%ix0+1:decomp(ib1)%ix1,1:lasttstep+1-firsttstep))
        ALLOCATE(dust(ib1)%vegmin2(decomp(ib1)%iy0+1:decomp(ib1)%iy1,    &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1))

        ! Allocate dust_flux
        ALLOCATE(dust(ib1)%soiltype(decomp(ib1)%iy0+1:decomp(ib1)%iy1,  &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1))
        ALLOCATE(dust(ib1)%z0(decomp(ib1)%iy0+1:decomp(ib1)%iy1,        &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1))
        ALLOCATE(dust(ib1)%source(decomp(ib1)%iy0+1:decomp(ib1)%iy1,   &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1))
        ALLOCATE(dust(ib1)%alpha2(decomp(ib1)%iy0+1:decomp(ib1)%iy1,    &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1))
        ALLOCATE(dust(ib1)%feff(decomp(ib1)%iy0+1:decomp(ib1)%iy1,     &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,dimveg))
        ALLOCATE(dust(ib1)%veff(decomp(ib1)%iy0+1:decomp(ib1)%iy1,     &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,dimveg))
        ALLOCATE(dust(ib1)%mfac(decomp(ib1)%iy0+1:decomp(ib1)%iy1,     &
                decomp(ib1)%ix0+1:decomp(ib1)%ix1))
        ALLOCATE(dust(ib1)%d_emis(decomp(ib1)%iy0+1:decomp(ib1)%iy1,   &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,1:nt))


        ALLOCATE(dust(ib1)%soilmap(decomp(ib1)%iy0+1:decomp(ib1)%iy1,   &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,nmode))

        ALLOCATE(ustar(decomp(ib1)%iy0+1:decomp(ib1)%iy1,   &
                decomp(ib1)%ix0+1:decomp(ib1)%ix1))


        dust(ib1)%biome(:,:)=0.
        dust(ib1)%cult(:,:)=0.
        dust(ib1)%veg(:,:,:)=0.
        dust(ib1)%vmoist(:,:,:)=0.
        dust(ib1)%vegmin2(:,:)=0.
        dust(ib1)%soiltype(:,:)=0.
        dust(ib1)%z0(:,:)=0.001 !cm
        dust(ib1)%source(:,:)=0.
        dust(ib1)%alpha2(:,:)=0.
        dust(ib1)%feff(:,:,:)=1.
        dust(ib1)%veff(:,:,:)=1.
        dust(ib1)%mfac(:,:)=1.
        dust(ib1)%d_emis(:,:,:)=0.

        dust(ib1)%soilmap = 0.

        ustar = 0

      END DO


      ! +-+-+- Sec 1.3 Input -+-+-+
      DO i = 1, ifile_num
        print*,ifile(i)
        AllOCATE (read_input(igy0+1:igy1,igx0+1:igx1,ifile_dim(i)))

        ! only prozess #0 open files
        IF (my_cart_id == 0) THEN
          CALL read_infile(ifile(i),read_input,ierr,yerr)


          IF (ierr /= 0) THEN
          print*, ierr
            ierr = 100009
            PRINT*,'ERROR    src_dust "init" '
            PRINT*,'         #',ierr
            PRINT*,'         ERROR reading',TRIM(ifile(i))
            PRINT*,'         ',yerr
            STOP
          END IF

        END IF !(my_cart_id == 0)

#ifndef OFFLINE
        ! distribut input to all px if necessary
        IF (num_compute > 1) THEN
          CALL MPI_BCAST(read_input,size(read_input),imp_reals,0,icomm_cart,ierr)
        END IF
#endif
        ! copy to muscat block structur
        DO ibLoc=1,nbLoc
#ifndef OFFLINE
          ib1 = LocGlob(ibLoc)
#else
          ib1 = 1
#endif
          IF (TRIM(ifile(i)) == 'soil' .AND. soilmaptype == 1) copy2d => dust(ib1)%soiltype
          IF (TRIM(ifile(i)) == 'soil' .AND. soilmaptype == 2) copy3d => dust(ib1)%soilmap
          IF (TRIM(ifile(i)) == 'source' ) copy2d => dust(ib1)%source
          IF (TRIM(ifile(i)) == 'z0')      copy2d => dust(ib1)%z0
          ! IF (filenum == 4) copy2d => dust(ib1)%biome
          IF (TRIM(ifile(i)) == 'vegday' ) copy3d => dust(ib1)%veg
          ! IF (filenum == 6) copy3d => dust(ib1)%veg
          ! IF (filenum == 7) copy2d => dust(ib1)%vegmin2
          IF (TRIM(ifile(i)) == 'moist') copy3d => dust(ib1)%vmoist
          !
          IF (ifile_dim(i) == 1) CALL copy2block(decomp(ib1),2,read_input,too2d=copy2d)
          IF (ifile_dim(i)  > 1) CALL copy2block(decomp(ib1),3,read_input,too3d=copy3d)



        END DO

        DEALLOCATE(read_input)
      END DO ! i = 1, ifile_num



      ! +-+-+- Sec 1.4 physical init -+-+-+
      ! the physical init is sperated in different subroutines
      !   this will make it easy to include other schemes later

      ! block loop
      DO ibLoc=1,nbLoc
#ifndef OFFLINE
        ib1 = LocGlob(ibLoc)
#else
        ib1 = 1
#endif

        ! +-+-+- Sec 1.4.1 thresold -+-+-+
        CALL init_thresold()

        ! +-+-+- Sec 1.4.2 Soilmap -+-+-+
        IF (soilmaptype ==  1) THEN
          CALL init_soilmap(decomp(ib1))
        END IF


        ! +-+-+- Sec 1.4.3 Moisture -+-+-+
        IF (moist_scheme == 1) THEN
          CALL fecan('init',decomp(ib1),ntstep)
        END IF

        ! +-+-+- Sec 1.4.4 Preferential Sources -+-+-+
        IF (psrcType > 0) THEN
          CAll init_psrc(decomp(ib1))
        END IF

        ! +-+-+- Sec 1.4.5 Alpha -+-+-+
        CALL init_alpha(decomp(ib1),2)


        ! +-+-+- Sec 1.4.6 Dust Flux -+-+-+
        IF (dust_scheme == 1) THEN
          ! init of the dust emission sheme by Tegen et al. 2002
          CALL init_tegen(decomp(ib1),ndays=dimveg)!ierr,yerr)
        ELSEIF (dust_scheme == 2) THEN
          ! init of the dust emission sheme by Tegen et al. 2002
          CALL tegen02('init',decomp(ib1))!ierr,yerr)
        END IF


        ! +-+-+- Sec 1.4.6 Vegetation -+-+-+
        IF (veg_scheme == 1) THEN
          CALL okin_vegetation(decomp(ib1),dimveg)
        ELSE IF (veg_scheme == 2) THEN
          CALL linear_vegetation(decomp(ib1),dimveg)
        END IF


        ! +-+-+- Sec 1.4.7 Surface Roughness -+-+-+
        IF (lwithz0) THEN
          ! Drag partition
          CALL roughness(decomp(ib1),dimveg)
        END IF



      END DO


    ! +-+-+- Sec 1.5 clean up -+-+-+

    ! DEALLOCATE(dust_it)

    call quick_nc('z0',var2d=dust(1)%z0)

    ! ------------------------------------
    ! +-+-+- Section 2 Dust flux calculation -+-+-+
    ! ------------------------------------
    ELSEIF (yaction == "calc") THEN

      CALL get_ustar(decomp(ib1))

      IF (moist_scheme == 1) THEN
        CALL fecan('calc',decomp(ib1),ntstep)
      END IF



      IF (dust_scheme == 1) THEN
        CALL emission_tegen(subdomain,flux)
      END IF

      IF (dust_scheme == 2) THEN
        ! init of the dust emission sheme by Tegen et al. 2002
        CALL tegen02(yaction,subdomain)!ierr,yerr)
      END IF
      ! STOP 'TESTING'

    END IF ! (yaction == "***")

    IF (lddebug) PRINT*, 'Leave organize_dust, yaction=',yaction,', ierr=',ierr,''//NEW_LINE('')

  END SUBROUTINE organize_dust

  !+ init_thresold
  !---------------------------------------------------------------------
  SUBROUTINE init_thresold()
  !---------------------------------------------------------------------
  ! Description:

  ! psrcType == 1
  ! Preferential Sources = Potential lakes
  !   If a grid box is identified as dust source then the soiltype switches to
  !   best emission properties (100% silt).
  !
  !
  !--------------------------------------------------------------------

    USE mo_dust
    USE tegen02_param
#ifdef OFFLINE
    USE offline_org
#endif

    IMPLICIT NONE

    INTEGER :: &
      n

    REAL(8) :: &
      dp,      & ! Current diameter
      dmy_B,   & ! Reynolds number, dummy for thresold friction velocity
      dmy_K   ! Reynolds number, dummy for thresold friction velocity


    IF (lddebug) PRINT*, 'Enter init_thresold'

    ! +-+-+- Sec 1 define particle diameter -+-+-+

    ! define particle diameter
    dp_meter(1) = Dmin
    DO n = 2, nclass
      dp_meter(n) = dp_meter(n-1) * EXP(Dstep)
    END DO


    ! +-+-+- Sec 2  calculat thresold friction velocity -+-+-+

    IF (threshold_scheme == 0) THEN
      Uth = 0.2 ! [ms-1] minimum value
      PRINT*, '*** WARNING! ***'
      PRINT*, 'threshold_scheme == 0'
      PRINT*, 'threshold friction velocity is set to minimum (0.2 ms-1)'
      PRINT*,''
    END IF

    ! particle size loop, calculation of Uth for every particle size
    DO n = 1, nclass
      dp = dp_meter(n)

      IF (threshold_scheme == 1) THEN

        dmy_K = SQRT(rop * g * dp / roa) * SQRT(10000. + 0.006 /(rop * g * dp ** 2.5))   ! Marticorena 95 eq(4) but in [m]
        dmy_B = a_rnolds * (dp ** x_rnolds) + b_rnolds                                   ! Marticorena 95 eq(5)
        IF (dmy_B < 10) THEN
          Uth(n) = 0.0013 * dmy_K / SQRT(1.928 * (dmy_B ** 0.092) - 1.)                  ! Marticorena 95 eq(6) but in [m]
        ELSE
          Uth(n) = 0.001207 * dmy_K * ( 1. -0.0858 * EXP(-0.0617 * (dmy_B - 10.)) )      ! Marticorena 95 eq(7) but in [m]
        END IF

      ELSEIF (threshold_scheme == 2) THEN
        Uth(n) = SQRT(0.0123 * (rop/roa * g *dp + 3.e-4/(roa*dp)) )
      END IF
    END DO ! n = 1, nclass

    IF (lddebug) PRINT*, 'Leave init_thresold',''//NEW_LINE('')

  END SUBROUTINE init_thresold

  !+ init_soilmap
  !---------------------------------------------------------------------
  SUBROUTINE init_soilmap(subdomain)
  !---------------------------------------------------------------------
  ! Description:

  !--------------------------------------------------------------------

    USE mo_dust
    USE dust_tegen_data
#ifdef OFFLINE
    USE offline_org
#endif

    IMPLICIT NONE

    TYPE(rectangle), INTENT(IN) :: subdomain

    INTEGER :: &
      i,j

    REAL(8), POINTER ::  &
      soiltype(:,:),     &
      soilmap(:,:,:)

    soiltype => dust(subdomain%ib)%soiltype(:,:)
    soilmap => dust(subdomain%ib)%soilmap(:,:,:)

    IF (lddebug) PRINT*, 'Enter init_soilmap'

    ! start lon-lat-loop
    DO i=1,subdomain%ntx
      DO j=1,subdomain%nty

        IF (soiltype(j,i) > 0.0) THEN
          soilmap(j,i,1) = solspe(soiltype(j,i),3)
          soilmap(j,i,2) = solspe(soiltype(j,i),6)
          soilmap(j,i,3) = solspe(soiltype(j,i),9)
          soilmap(j,i,4) = solspe(soiltype(j,i),12)
        ELSE
          soilmap(j,i,:) = 0.0
        END IF

      END DO
    END DO
    ! end lon-lat-loop

    IF (lddebug) PRINT*, 'Leave init_soilmap',''//NEW_LINE('')

  END SUBROUTINE init_soilmap

  !+ init_alpha
  !---------------------------------------------------------------------
  SUBROUTINE init_alpha(subdomain,alpha_type)
  !---------------------------------------------------------------------
  ! Description:

  !--------------------------------------------------------------------

    USE mo_dust
    USE dust_tegen_data
#ifdef OFFLINE
    USE offline_org
#endif

    IMPLICIT NONE

    TYPE(rectangle), INTENT(IN) :: subdomain

    INTEGER, INTENT(IN) :: alpha_type

    INTEGER :: &
      i,j

    REAL(8), POINTER ::  &
      alpha(:,:), &
      soiltype(:,:), &
      soilmap(:,:,:)

    alpha    => dust(subdomain%ib)%alpha2(:,:)
    soiltype => dust(subdomain%ib)%soiltype(:,:)
    soilmap  => dust(subdomain%ib)%soilmap(:,:,:)

    IF (lddebug) PRINT*, 'Enter init_alpha'

    ! start lon-lat-loop
    DO i=1,subdomain%ntx
      DO j=1,subdomain%nty

        !alpha_type == 1 : lookup table
        IF (alpha_type == 1) THEN

          ! if the soil type is in the lookup table
          IF (soiltype(j,i) > 0 .AND. soiltype(j,i) < nats) THEN
            alpha(j,i) = solspe(soiltype(j,i),13)
          ELSE
            alpha(j,i) = 0.
          END IF

        ! alpha_type == 2 : calc from fraction
        ELSEIF (alpha_type == 2) THEN
          IF (soilmaptype == 1) THEN
            alpha(j,i) = soilmap(j,i,1) * 1.E-7 &
                  + soilmap(j,i,2) * 1.E-6 &
                  + soilmap(j,i,3) * 1.E-5 &
                  + soilmap(j,i,4) * 1.E-6
            IF (soilmap(j,i,4) > 0.45) THEN
              alpha(j,i) = soilmap(j,i,1) * 1.E-7 &
                    + soilmap(j,i,2) * 1.E-6 &
                    + soilmap(j,i,3) * 1.E-5 &
                    + soilmap(j,i,4) * 1.E-7
            END IF
          ELSEIF (soilmaptype == 2) THEN
            alpha(j,i) = soilmap(j,i,1) * 1.E-6 &
                  + soilmap(j,i,2) * 1.E-5 &
                  + soilmap(j,i,3) * 1.E-6
            IF (soilmap(j,i,3) > 0.45) THEN
              alpha(j,i) = soilmap(j,i,1) * 1.E-6 &
                    + soilmap(j,i,2) * 1.E-5 &
                    + soilmap(j,i,3) * 1.E-7
            END IF
          ENDIF

        ! alpha_type == 3 : clay fraction eq from Marticorena
        !ELSEIF (alpha_type == 3) THEN

        END IF

      END DO
    END DO
    ! end lon-lat-loop

    ! all alpha values are in [cm-1] but for for dust_scheme == 2 we need alpha in [m-1]
    IF (dust_scheme == 2)  THEN
      alpha = alpha*100. !  [cm-1] = 100 [m-1]
    END IF

  !  call quick_nc('alpha',var2d=alpha)
  !  call quick_nc('soil',var3d=soilmap)

  ! call quick_ascii('alpha',alpha)
  ! call quick_ascii('sand',soilmap(:,:,1),pmin=0.,pmax=1.)
  ! call quick_ascii('silt',soilmap(:,:,2),pmin=0.,pmax=1.)
  ! call quick_ascii('clay',soilmap(:,:,3),pmin=0.,pmax=1.)

  IF (lddebug) PRINT*, 'Leave init_alpha',''//NEW_LINE('')

  END SUBROUTINE init_alpha


  !+ init_psrc
  !---------------------------------------------------------------------
  SUBROUTINE init_psrc(subdomain)
  !---------------------------------------------------------------------
  ! Description:

  ! psrcType == 1
  ! Preferential Sources = Potential lakes
  !   If a grid box is identified as dust source then the soiltype switches to
  !   best emission properties (100% silt).
  !
  !
  !--------------------------------------------------------------------

    USE mo_dust
    USE dust_tegen_data
#ifdef OFFLINE
    USE offline_org
#endif

    IMPLICIT NONE

    TYPE(rectangle), INTENT(IN) :: subdomain

    INTEGER :: &
      i,j

    REAL(8), POINTER ::  &
      psrc(:,:),     &
      soilmap(:,:,:)

    psrc => dust(subdomain%ib)%source(:,:)
    soilmap => dust(subdomain%ib)%soilmap(:,:,:)

    IF (lddebug) PRINT*, 'Enter init_psrc'

    ! start lon-lat-loop
    DO i=1,subdomain%ntx
      DO j=1,subdomain%nty

        IF (psrcType == 1) THEN ! preferential source scheme by Tegen02
          IF (psrc(j,i) > 0.5 .AND. sum(soilmap(j,i,:)) > 0.0) THEN
            soilmap(j,i,:) = 0.0
            soilmap(j,i,nmode-1) = 1.
          END IF
        ELSEIF (psrcType == 2) THEN ! MSG source scheme by schepanski08
          IF (psrc(j,i) >= 2 .AND. sum(soilmap(j,i,:)) > 0.0) THEN
            soilmap(j,i,:) = 0.0
            soilmap(j,i,nmode-1) = 1.
          END IF
        END IF

      END DO
    END DO
    ! end lon-lat-loop
    ! call quick_nc('source',var2d=psrc)
    ! call quick_ascii('source',psrc,pmin=0.,pmax=1.)

    IF (lddebug) PRINT*, 'Leave init_psrc',''//NEW_LINE('')

  END SUBROUTINE init_psrc



  !+ init_tegen
  !---------------------------------------------------------------------
  SUBROUTINE tegen02(yaction,subdomain)
  !---------------------------------------------------------------------
  ! Description:
  !   This subroutine performes the initialization for
  !   the dust emisson scheme by Tegen et al. 2002
  !   https://doi.org/10.1029/2001JD000963
  !
  !   The Physics used by Tegen based on the paper of
  !   Marticorena and Bergametti 1995
  !   https://doi.org/10.1029/95JD00690
  !
  ! The code based on Tegen et al. 2002
  !--------------------------------------------------------------------

    ! Modules
    USE mo_dust
    USE tegen02_param
    USE dust_tegen_data
#ifdef OFFLINE
    USE offline_org
#endif


    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)            :: &
      yaction ! action to be performed

    TYPE(rectangle), INTENT(IN) :: subdomain


    INTEGER :: &
      i,j,n,m, & ! loops
      tnow

    REAL(8) :: &
      dp,      & ! Current diameter
      dM,      &   ! mass size distribution
      stot,    &    ! total basal surface
      uwind,   &
      vwind,   &
      tot_wind,   &
      time_start, &
      time_now, &
      uthp, &
      dmy_R, &
      s_rel,  &
      dflux

    REAL(8) :: &
      fluxtot (ntrace), &
      fluxbin (ntrace)

    real    ::  T1,T2

    REAL(8), POINTER :: source(:,:)
    REAL(8), POINTER :: soilmap(:,:,:)
    REAL(8), POINTER :: feff(:,:,:)
    REAL(8), POINTER :: veff(:,:,:)
    REAL(8), POINTER :: mfac(:,:)
    REAL(8), POINTER :: z0(:,:)
    REAL(8), POINTER :: alpha(:,:)
    REAL(8), POINTER :: DustEmis(:,:,:)

#ifndef OFFLINE
    REAL(8), POINTER :: EmiRate(:,:,:,:)
#endif

    source   => dust(subdomain%ib)%source(:,:)
    soilmap  => dust(subdomain%ib)%soilmap(:,:,:)
    feff     => dust(subdomain%ib)%feff(:,:,:)
    veff     => dust(subdomain%ib)%veff(:,:,:)
    mfac     => dust(subdomain%ib)%mfac(:,:)
    z0       => dust(subdomain%ib)%z0(:,:)
    alpha    => dust(subdomain%ib)%alpha2(:,:)
    DustEmis => dust(subdomain%ib)%d_emis(:,:,:)


#ifndef OFFLINE
    EmiRate  => block(subdomain%ib)%EmiRate(:,:,:,:)
#endif

    IF (lddebug) PRINT*, 'Enter tegen02, yaction=',yaction

    IF (yaction == 'init') THEN
      call cpu_time(T1)


      ALLOCATE(median_dp(nmode))

      median_dp(nmode)   = median_dp_clay
      median_dp(nmode-1) = median_dp_silt
      median_dp(nmode-2) = median_dp_sand
      IF (nmode == 4) median_dp(nmode-3) = median_dp_csand

      ALLOCATE(srel_map(subdomain%nty,subdomain%ntx,nclass))

      ! +-+-+- Sec xx srel calculation -+-+-+
      ! start lon-lat-loop
      DO i = 1,subdomain%ntx
        DO j = 1,subdomain%nty

          stot = 0.
          DO n = 1, nclass
            dp = dp_meter(n)

            ! calc soil mass size distribution Marticorena 95 eq(29)
            dM = 0.
            DO m = 1, nmode
              dM = dM + soilmap(j,i,m)/(SQRT(2. * pi) * LOG(sigma)) &
                      * EXP((LOG(dp) - LOG(median_dp(m)))**2. / (-2. * LOG(sigma)**2.))
            END DO

            ! size distribution of basal surface Marticorena 95 eq(30)
            srel_map(j,i,n) = dM / (2./3. * rop * dp) * Dstep

            ! total basal surface Marticorena 95 eq(31)
            stot = stot + srel_map(j,i,n)

          END DO ! n = 1, nclass


          DO n = 1, nclass
            dp = dp_meter(n)

            IF (stot > 0.) THEN
              srel_map(j,i,n) = srel_map(j,i,n)/stot
            ELSE
              srel_map(j,i,n) = 0.
            END IF

          END DO ! n = 1, nclass


        END DO ! j=1,subdomain%nty
      END DO ! i=1,subdomain%ntx


      call cpu_time(T2)

      print*, 'init time:',T2-T1

    ELSEIF (yaction == 'calc') THEN
      ! +-+-+- Sec 1 Set the actually date -+-+-+

      ! the drag partition (feff) has a dependency on time
      IF (veg_scheme > 0) THEN
        IF (lvegdaily) THEN
          ! find the exact day and set "tnow" with the number of the actually day
          READ(ydate_ini(9:10),*) time_start
          time_start=time_start + hstart
          time_now=time_start+ntstep*dt/3600.
          tnow=time_now/24 + 1 ! convert to integer
        ELSE ! if the drag partition has monthly values "tnow" is the number of the month
          READ(StartDate,'(4x,i2)') tnow
        END IF
      ELSE
        tnow = 1
      END IF

      DO i=1,subdomain%ntx
        DO j=1,subdomain%nty
          call cpu_time(T1)

          ! +-+-+- Sec 2 update of the meteorological variables -+-+-+

          IF(feff(j,i,tnow) > 0. .AND. ustar(j,i) > 0. ) THEN

            m = 1 ! index of dust bin
            dflux = 0.
            fluxbin = 0.

            DO n = 1, nclass
              dp = dp_meter(n)

              ! calculate the recent threshold friction velocity
              ! from the particle thr. fric. velo.,
              ! the drag partition from soil roughness and vegetation
              ! and the moisture factor
              uthp = uth(n)/feff(j,i,tnow) * mfac(j,i)

              s_rel = srel_map(j,i,n)

              IF (m > DustBins) m = DustBins
              IF (m <= DustBins) THEN
                IF (dp > dustbin_top(m)) m = m+1
              END IF


              IF (ustar(j,i) > uthp) THEN
                dmy_R = uthp/ustar(j,i)
              ELSE
                dmy_R = 1
              END IF


              ! Horizontal dust flux
              dflux = roa/g * ustar(j,i)**3 * (1+dmy_R) * (1-dmy_R**2) * s_rel
              ! END IF

              ! Vertical dust flux
              dflux = dflux * alpha(j,i)

              ! bin-wise integration
              IF (m <= DustBins) THEN
                fluxbin(m) = fluxbin(m)+dflux
              END IF


            END DO ! n = 1, nclass


          ENDIF

          ! DO n=1,ntrace
          !   fluxtot(n) = fluxtot(n) + fluxbin(n)
          ! END DO

          DO n=1,ntrace
            ! fluxtot: g/cm2/sec --> kg/m2/sec
            !

            ! Mask Effective area determined by preferential source fraction:
            ! only for psrcType = 2
            ! IF (psrcType == 1) THEN
            !   fluxbin(n) = fluxbin(n) * source(j,i)
            ! END IF

            ! Mask Effective area determined by vegetation fraction:
            ! only for veg_scheme = 2
            IF (veg_scheme == 2) THEN
              fluxbin(n) = fluxbin(n) * veff(j,i,tnow)
            END IF

            ! write output in [g m-2 s-1]
            DustEmis(j,i,n) = fluxbin(n) * 1.E3
          END DO

          call cpu_time(T2)
          !print*, 'calc time step:',T2-T1
        END DO ! j
      END DO ! i
      !call quick_ascii('Dust',SUM(DustEmis,dim=3))
    END IF ! yaction

    IF (lddebug) PRINT*, 'Leave tegen02, yaction=',yaction,''//NEW_LINE('')

  END SUBROUTINE tegen02

  !+ init_tegen
  !---------------------------------------------------------------------
  SUBROUTINE init_tegen(subdomain,ndays)
  !---------------------------------------------------------------------
  ! Description:
  !   This subroutine performes the initialization for
  !   the dust emisson scheme by Tegen et al. 2002
  !   https://doi.org/10.1029/2001JD000963
  !
  !   The Physics used by Tegen based on the paper of
  !   Marticorena and Bergametti 1995
  !   https://doi.org/10.1029/95JD00690
  !
  ! The code based on Tegen et al. 2002
  !--------------------------------------------------------------------

    ! Modules
    USE mo_dust
    USE dust_tegen_param
    USE dust_tegen_data
#ifdef OFFLINE
    USE offline_org
#endif


    IMPLICIT NONE

    TYPE(rectangle), INTENT(IN) :: subdomain
    INTEGER, OPTIONAL, INTENT(IN) :: ndays

    INTEGER    :: &
      nn,ns,kk,nm,      & ! loops
      nd,nsi,np,        &
      i,j,t,            &
      isoiltype

    REAL (8)   :: &
      dp,               & ! particle size inside the loop
      Stotal,           & ! total surface area of particles
      StotalV,          & ! total surface volume of particles
      su,               & ! surface area of particles
      suV,              & ! surface volume of particles
      su_loc,           & ! local su step
      su_locV!,          & ! local suV step

    ! 1D Arrays
    REAL (8)   :: &
      ! Uth(Nclass),              & ! threshold friction velocity
      ! Uth_bod(Nclass),          & ! threshold friction velocity
      gransize(nclass),         & !granulometric class size
      dlastj (nclass),          & ! ???
      dkmin(nclass),            & !min original size class dkmin
      dkmax(nclass),            & !max original size class dkmax
      su_class(nclass),         & !surface occupied by each granulometric class
      su_classV(nclass),        & !volume occupied by each granulometric class
      utest(nats)

    ! ! 2D Arrays
    ! REAL(8) ::                  &
    !   srel(nats,nclass),        & !
    !   srelV(nats,nclass),       & !
    !   su_srelV(nats,nclass)!,    & !

    REAL(8) ::          &
      dmy_B,            & ! Dummy for B Factor (Marticorena 95)
      dmy_K,            & ! Dummy for K Factor (Marticorena 95)
      dmy_xk,           & ! Dummy for Marticorena 95 eq. (29)
      dmy_xl,           & ! Dummy for Marticorena 95 eq. (29)
      dmy_xm,           & ! Dummy for Marticorena 95 eq. (29)
      dmy_xn!,           & ! Dummy for Marticorena 95 eq. (29)

    REAL(8), POINTER ::  &
      soiltype(:,:),     &
      source(:,:),       &
      z0(:,:),           &
      alpha(:,:),        &
      feff(:,:,:)

    REAL(8), ALLOCATABLE :: printvar(:,:,:,:)

    soiltype => dust(subdomain%ib)%soiltype(:,:)
    source => dust(subdomain%ib)%source(:,:)
    z0 => dust(subdomain%ib)%z0(:,:)
    alpha => dust(subdomain%ib)%alpha2(:,:)
    feff => dust(subdomain%ib)%feff(:,:,:)
    !  => dust_ini(ib1)%biome
    !  => dust_ini(ib1)%veg
    !  => dust_ini(ib1)%veg
    !  => dust_ini(ib1)%vegmin

    IF (lddebug) PRINT*, 'Enter init_tegen'

    ! +-+-+- Sec 1 init of threshold friction velocity Uth -+-+-+
    ! Marticorena and Bergametti 1995 Sec 2.2, eq. (3)-(7)
    ! https://doi.org/10.1029/95JD00690

    gransize(:) = 0.
    Uth(:)  = 0.
    su_srelV(:,:) =  0.
    srelV(:,:) = 0.
    srel(:,:)= 0.

    nn = 0
    dp = Dmin


    Uth = Uth * 100 ! ms-1 to cms-1

    ! +-+-+- Sec 2 init of surface saltation flux -+-+-+
    ! Marticorena and Bergametti 1995 Sec 2.3.1, see eq. ()

    dkmin=0.
    dkmax=0.

    ! soiltype loop, calculation surface flux for every possible soiltype
    DO ns = 1,nats !ns (soil type)
      Stotal = 0.
      StotalV = 0.
      dp = Dmin
      kk = 0
      su_class(:) = 0.
      su_classV(:) = 0.
      utest(:) = 0.


      ! particle size loop
      DO WHILE (dp <= Dmax+1E-5)
        kk = kk + 1
        su = 0.
        suV =0.

        ! loop of soiltype modes
        ! (coarse sand, med/fine sand, silt, clay)
        DO nm = 1, nmode
          nd  = ((nm - 1) *3 ) + 1 ! index of solspe that point to mean diameter
          nsi = nd + 1             ! index of solspe that point to sigma
          np  = nd + 2             ! index of solspe that point to the fraction of sand/silt/clay

          IF (solspe(ns,nd) == 0.) THEN
            su_loc = 0.
          !  su_locV=0.
          ELSE
            !Marticorena 95, eq. (29)
            dmy_xk = solspe(ns,np)/(sqrt(2.* pi2)*log(solspe(ns,nsi)))
            dmy_xl = ((log(dp)-log(solspe(ns,nd)))**2)/                 &
                 (2.*(log(solspe(ns,nsi)))**2)
            dmy_xm = dmy_xk * exp(-dmy_xl)
            !Marticorena 95, eq. (30)
            dmy_xn =  rop*(2./3.)*(dp/2.) !surface
            su_loc = (dmy_xm*Dstep/dmy_xn)
            su_locV = (dmy_xm*Dstep)
          END IF
          su  = su  + su_loc
          suV = suV + su_locV
        END DO !nmode

        su_class(kk) = su
        su_classV(kk) = suV
        Stotal = Stotal + su
        StotalV = StotalV + suV
        dlastj(kk)=dp*5000.
        dkmin(kk)=dp
        dp = dp * exp(Dstep)
        dkmax(kk)=dp
      END DO ! (dp <= Dmax+1E-5)

      DO nn = 1,nclass
        IF (Stotal.eq.0.)THEN
          srel(ns,nn) = 0.
          srelV(ns,nn) = 0.
        ELSE
          !Marticorena 95, eq. (33)
          srel(ns,nn) = su_class(nn)/Stotal
          srelV(ns,nn) = su_classV(nn)/StotalV
          utest(ns)=utest(ns)+srelV(ns,nn)
          su_srelV(ns,nn)=utest(ns)
          ! if (ns == 27 .and. nn == 99) print*, 'utest',ns,nn,utest(ns),su_srelV(ns,nn)!,srelV(ns,nn),su_classV(nn),StotalV,su_class(nn),Stotal
        END IF
      END DO !nn=1,nclass
    END DO !ns (soil type)

    ! +-+-+- Sec 3 Prepare the flux calculation -+-+-+
    ! start lon-lat-loop
    DO i=1,subdomain%ntx
      DO j=1,subdomain%nty

        ! +- Sec 3.1 Selection of potential dust sources areas
        IF (psrcType == 1) THEN  ! IT02
          ! Preferential Sources = Potential lakes
          !   If a grid box is identified as dust source then the soiltype switches to
          !   best emission properties (soiltype=10), to avoid an underestimation.
          IF(source(j,i) > 0.5) THEN
            soiltype(j,i) = 10
            IF (z0(j,i) <= 0.) z0(j,i) = 0.001 ! [cm]
          ENDIF

        ELSEIF (psrcType == 2) THEN ! MSG
          ! Preferential Sources = Potential lakes
          IF(soiltype(j,i) < 13) soiltype(j,i) = 9
          IF(source(j,i) >= 2) THEN
            IF (soiltype(j,i) == 9) soiltype(j,i) = 10
            IF (z0(j,i) <= 0.) z0(j,i) = 0.001 ! [cm]
          ENDIF

        ENDIF



        isoiltype = int(soiltype(j,i)) ! index of soiltype

        ! change everthing that not fits the soiltype table to 9 = ICE
        IF(isoiltype < 1 .or. isoiltype > nats) isoiltype=9

        ! init alpha from lookup table
        !alpha(j,i) = solspe(isoiltype,nmode*3+1)

        ! avoid z0 = 0. at any place
        IF (z0(j,i) == 0.0) z0(j,i)=1.E-9






    !
    !     !---------------------------------------------------------------------------------------
    !     !       Factors according to dsf increase seen in data **
    !     !---------------------------------------------------------------------------------------
    !    umin2(j,i,1)=umin
    !    IF (sp(j,i,1,3).EQ.2.) THEN      !rivm !
    !
    !     ! 0.96(rx1.5) 0.89(sx2) 0.72(rx3) 0.89(rx2)   rx1.5 for umin=18: 0.93
    !     !---------------------------------------------------------------------------------------
    !      IF ((sp(j,i,1,5).EQ.21..OR.sp(j,i,1,5).EQ.27)       &
    !         .OR.(sp(j,i,1,5).EQ.13.OR.sp(j,i,1,5).EQ.14))    &
    !          umin2(j,i,1) = umin*0.93
    !         ! 0.97(rx1.5) 0.95(sx2) 0.73(rx3) 0.93(rx2)rx1.5 for umin=18: 0.99
    !         !---------------------------------------------------------------------------------------
    !      IF ((sp(j,i,1,5).EQ.19.OR.sp(j,i,1,5).EQ.20))       &
    !           umin2(j,i,1) = umin*0.99
    !    END IF !cult=2
    !
    !     !  crop 0.5 (rivm/sage x10), 0.58 (rivmx5) 0.68 rivm x 3.5, 0.77 rivm x 2.5, 0.86 rivm x 2
    !     !---------------------------------------------------------------------------------------
    !    IF (sp(j,i,1,3).eq.1.) THEN      !rivm !
    !      IF ((sp(j,i,1,5).EQ.21.OR.sp(j,i,1,5).EQ.27)     &
    !      .OR.(sp(j,i,1,5).EQ.13.OR.sp(j,i,1,5).EQ.14))    &
    !       umin2(j,i,1)=umin*0.73                         !&
    !       !             umin2(j,i,1)=umin*0.5
    !    END IF !cult=1

      END DO
    END DO
    ! end lon-lat-loop


    IF (lddebug) PRINT*, 'Leave init_tegen',''//NEW_LINE('')

  END SUBROUTINE init_tegen


  !+ emission_tegen
  !---------------------------------------------------------------------
  SUBROUTINE emission_tegen(subdomain,flux)
  !---------------------------------------------------------------------
  ! Description:
  !   This subroutine performes the dust flux calculation for
  !   the dust emisson scheme by Tegen et al. 2002
  !   https://doi.org/10.1029/2001JD000963
  !
  !   The Physics used by Tegen based on the paper of
  !   Marticorena and Bergametti 1995
  !   https://doi.org/10.1029/95JD00690
  !
  ! The code based on Tegen et al. 2002
  !--------------------------------------------------------------------

    ! Modules
    USE mo_dust
    USE dust_tegen_param
    USE dust_tegen_data
#ifndef OFFLINE
    USE partition
    USE sub_block
    USE sub_geo
    USE sub_met
    USE mo_ctm
    USE mo_gas, ONLY: ConvPart
    ! USE dust_org
    USE data_io,  ONLY : ydate_ini
    USE data_runcontrol,    ONLY : hstart, ntstep
    USE data_modelconfig, ONLY: dt
    USE data_fields, ONLY:  ustar_fv
#else
    USE offline_org
#endif

    IMPLICIT NONE

    TYPE (rectangle), INTENT(IN)    :: &
      subdomain

    REAL(8),          INTENT(INOUT) :: &
      Flux(ntz,subdomain%nty,subdomain%ntx,nt)

    ! Internal variables
    INTEGER :: it
    REAL(8) :: del,tnm,sum

    INTEGER ::   &
      i,j,nn,    &
      k,kk,kkk,  &
      kfirst,    &
      kkmin,     &
      tnow,      &
      i_s1,      &
      i_s11!,     &



    REAL (8)   :: &
      dp,               & ! particle size inside the loop
      dpd,              & ! particle size inside the loop
      uthp,             & ! threshold friction velocity for each particle diam.
      fdp1,fdp2,        & !
      flux_umean,       & ! dust flux with in terms of the mean Ustar thrsld
      flux_diam,        & ! vertical dust flux,
      Ustar_var,        & !
      dbstart,          & !minimum dust bin for sandblasting
      en_kin,           & !kinetic energy of dust
      dlast,            &
      cultfac


    REAL(8) :: uwind, vwind,van
    !REAL(8) :: mfac             !factor due to soil moisture [Fecan, F. et al., 1999]
    REAL(8) :: FDust(subdomain%nty,subdomain%ntx,ntrace)
    REAL(8) :: time_start,time_now

    REAL(8) :: &
      ustar_cm
    REAL(8) ::                                &
      dpk (ntrace),                           & !dpk
      dbmin (ntrace),                         & !bin size limit
      dbmax (ntrace),                         & !bin size limit
      fluxtot (ntrace),                       & !total dust flux at 06am and 06 pm fluxtot06
      fluxtyp (nclass),                       & !
      fluxbin (ntrace)!,                       & !


    ! Pointer
#ifndef OFFLINE
    REAL(8), POINTER :: dxK(:,:)
    REAL(8), POINTER :: dyK(:,:)
    REAL(8), POINTER :: dz(:,:,:)
    REAL(8), POINTER :: usur(:,:)
    REAL(8), POINTER :: vsur(:,:)
    REAL(8), POINTER :: qrsur(:,:)
    REAL(8), POINTER :: rhosur(:,:)
#endif

    REAL(8), POINTER :: soiltype(:,:)
    REAL(8), POINTER :: source(:,:)
    REAL(8), POINTER :: alpha(:,:)
    REAL(8), POINTER :: feff(:,:,:)
    REAL(8), POINTER :: veff(:,:,:)
    REAL(8), POINTER :: mfac(:,:)
    REAL(8), POINTER :: z0(:,:)
    ! REAL(8), POINTER :: umin2(:,:)
    REAL(8), PARAMETER :: umin2=umin
    REAL(8), POINTER :: w_str(:,:)
    REAL(8), POINTER :: DustEmis(:,:,:), EmiRate(:,:,:,:)




    IF (nDust == 0) RETURN
    IF (DustMod <= 0) RETURN


    soiltype => dust(subdomain%ib)%soiltype(:,:)
    source   => dust(subdomain%ib)%source(:,:)
    alpha    => dust(subdomain%ib)%alpha2(:,:)
    feff     => dust(subdomain%ib)%feff(:,:,:)
    veff     => dust(subdomain%ib)%veff(:,:,:)
    mfac     => dust(subdomain%ib)%mfac(:,:)
    z0       => dust(subdomain%ib)%z0(:,:)
    ! lai_eff => dust(subdomain%ib)%lai_eff(:,:,:,:)
    ! umin2    => dust(subdomain%ib)%umin2(:,:,1)
    !
    ! w_str   => dust(subdomain%ib)%w_str(:,:,1)
    !
    DustEmis => dust(subdomain%ib)%d_emis(:,:,:)
#ifndef OFFLINE
    EmiRate  => block(subdomain%ib)%EmiRate(:,:,:,:)
#endif

IF (lddebug) PRINT*, 'Enter emission_tegen'

    ! +-+-+- Sec 1 Set the actually date -+-+-+

    ! the drag partition (feff) has a dependency on time
    IF (veg_scheme > 0) THEN
      IF (lvegdaily) THEN
        ! find the exact day and set "tnow" with the number of the actually day
        READ(ydate_ini(9:10),*) time_start
        time_start=time_start + hstart
        time_now=time_start+ntstep*dt/3600.
        tnow=time_now/24 + 1 ! convert to integer
      ELSE ! if the drag partition has monthly values "tnow" is the number of the month
        READ(StartDate,'(4x,i2)') tnow
      END IF
    ELSE
      tnow = 1
    END IF



    DO i=1,subdomain%ntx
      DO j=1,subdomain%nty
        ! +-+-+- Sec 2 update of the meteorological variables -+-+-+

        !---  zero setting
        dbmin(:)=0.E0
        dbmax(:)=0.E0
        fluxtot(:) = 0.E0
        fluxtyp(:)=0.0
        fluxbin(:)=0.0


        i_s1 = int(soiltype(j,i))  ! soiltype value of grid cell

        IF (i_s1 == 0) i_s1 = 9    ! set "0" (mostly 0->Water) to 9->Ice

        i_s11=i_s1                 ! ??? i_s11 strong emitting soiltype for optimisation??

        IF (i_s1 == 10 .OR. i_s1 == 12) i_s11 = 11  ! i_s1=10  ->  Potential Lakes,
                                                    ! i_s1=12  ->  Potential Lakes Australia,
                                                    ! i_s11=11 ->  Potential Lakes (Clay)

        ! IF((i_s1.EQ.10.OR.i_s1.EQ.12).AND.van.GT.10)i_s11=11

        IF (i_s1 == 13) i_s11 = 14  ! i_s1=13  ->  Potential Lakes BODELE,
                                    ! i_s11=14 ->  Potential Lakes BODELE (Clay)








        ! Friction velocity of the wind (ustar)
        ustar_cm = ustar(j,i)*100.



        ! +-+-+- Sec 3 Flux calculation -+-+-+
        IF (feff(j,i,tnow) > 0.) THEN
          IF (ustar_cm > 0 .AND. ustar_cm > umin2/feff(j,i,tnow) ) THEN
            kk = 0
            dp = Dmin
            DO WHILE (dp <= Dmax+1E-5)
              kk = kk+1

              ! original Tegen Code
              ! ! Is this reduction necessary (MF)?
              Uthp=uth(kk)*umin2/umin*u1fac !reduce threshold for cultivated soils
              ! drag coeff
              Uthp=Uthp/feff(j,i,tnow)
              ! ! moist
              ! Uthp=Uthp*mfac(j,i,tnow)
              ! Marticorena:
              fdp1 = 1.+(Uthp/ustar_cm)
              fdp2 = 1.-(Uthp/ustar_cm)**2.
              ! fdp1 = (1.-(Uthp * Ustar))
              ! fdp2 = (1.+(Uthp * Ustar))**2.

              ! ! Shao:
              ! fdp1 = (1.-(Uthp/(feff(j,i,tnow) * Ustar))**2)
              ! fdp2 = 1.

              IF (fdp1 <= 0 .OR. fdp2 <= 0) THEN
                flux_umean = 0.
              ELSE
                flux_umean = srel(i_s1,kk) * fdp1 * fdp2 * cd * ustar_cm**3 *alpha(j,i)
                flux_diam = flux_umean

                Ustar_var = Ustar_min

                ! ! taking subgrid-scale variations of Ustar into account
                ! DO WHILE(Ustar_var <= Ustar_max+1E-5)
                !   flux_diam = flux_diam + flux_umean * 2*(Ustar_var/(Ustar**2)) * &
                !               exp(-(Ustar_var/Ustar)**2) * (Ustar_var*exp(Ustar_step)-Ustar_var)
                !   Ustar_var=Ustar_var*exp(Ustar_step)
                ! END DO !Ustar_var

                ! calculate kinetic energy of saltating particle fraction en_kin
                !   and minimum dust particle size for emitted dust (Alfaro et al JGR 1999)

                !#############################################################################
                !V_start
                !   This is only done for diatomite fields (Bodele) (Tegen et al ACP 2006).
   	            !   In other regions, all sizes are mobilised.
   	            !----------------------------------------------------------------------
                dbstart=dp
                IF(i_s1 == 13) THEN
                  en_kin = 0.5*1./24.*rop_bod*pi2*(dp**3)*(95.*van)**2  !g*cm2/s2

                  IF (en_kin > 0.) dbstart = -1./0.53*log(en_kin*20/14.9)*0.0001

                  IF(dbstart < dmin) dbstart = dmin
                ELSE
                  dbstart = dmin        !all sizes mobilised
                END IF

                !V_end
                !#############################################################################

                IF (dbstart >= dp) THEN
                  fluxtyp(kk)=fluxtyp(kk)+flux_diam
                ELSE
                  ! loop over dislocated dust particle sizes
                  dpd=dmin
                  kkk=0
                  kfirst=0
                  DO WHILE(dpd <= dp+1e-5)
                    kkk=kkk+1

                    IF (dpd >= dbstart) THEN
                      IF(kfirst == 0) kkmin=kkk
                      kfirst=1

                      ! scaling with relative contribution of dust size  fraction
                      IF (kk > kkmin) THEN
                        fluxtyp(kkk) = fluxtyp(kkk) +flux_diam               &
                                      *srelV(i_s11,kkk)/((su_srelV(i_s11,kk) &
                                      -su_srelV(i_s11,kkmin)))
                      END IF ! (kk > kkmin)
                    END IF ! (dpd >= dbstart)
                    dpd=dpd*exp(dstep)
                  END DO ! dpd
                  ! end of saltation loop
                END IF ! (dbstart >= dp)
              END IF ! (fdp1 <= 0 .OR. fdp2 <= 0)

              dp = dp * exp(Dstep)
            END DO ! WHILE (dp <= Dmax+1E-5)

            ! +-+-+- Sec 4 Section assign fluxes to bins -+-+-+

            dp=dmin
            dlast=dmin
            nn=1
            kk=0

            DO WHILE (dp <= dmax+1e-5)
              kk=kk+1
              IF (nn <= ntrace) THEN
                fluxbin(nn) = fluxbin(nn)+fluxtyp(kk)
                IF(mod(kk,nbin).EQ.0) THEN
                  dbmax(nn)=dp*10000.*0.5  !radius in um
                  dbmin(nn)=dlast*10000.*0.5           !
                  dpk(nn)=sqrt(dbmax(nn)*dbmin(nn))
                  nn=nn+1
                  dlast=dp
                END IF
              END IF
              dp = dp * exp(Dstep)
            END DO !dp

            DO nn=1,ntrace
              fluxtot(nn) = fluxtot(nn) + fluxbin(nn)
            END DO

          END IF   ! (Ustar > 0 .AND. Ustar > umin2/feff(j,i,tnow) )
        END IF   ! (feff(j,i,tnow) > 0.)



        ! clfc1=0
        cultfac=1.

        DO nn=1,ntrace
          ! fluxtot: g/cm2/sec --> kg/m2/sec
          ! MASK: Effective area determined by cultfac/snow
          fdust(j,i,nn) = fluxtot(nn) * 10000 / 1000          &
                          *cultfac!*(1.-snow365(j,i))

          ! Mask Effective area determined by preferential source fraction:
          ! only for psrcType = 2
          IF (psrcType == 2) THEN
            fdust(j,i,nn) = fdust(j,i,nn) * source(j,i)
          END IF

          ! Mask Effective area determined by vegetation fraction:
          ! only for veg_scheme = 2
          IF (veg_scheme == 2) THEN
            fdust(j,i,nn) = fdust(j,i,nn) * veff(j,i,tnow)
          END IF


          ! ! MASK: Soil moisture threshold, using w0 !MF don't use this
          ! IF(qrsur(j,i).GE.w0) THEN
          !   fdust(j,i,nn)=0.
          ! END IF

        END DO ! nn=1,ntrace

        !---  transformation to chemistry units  (RW)
        FDust(j,i,:) = FDust(j,i,:) * 1.E3         ! kg/m2/s ==> g/m2/s
        ! FDust(j,i,:) = FDust(j,i,:) * ConvPart     ! chemistry units (nradm=1): g/m2/s ==> g/m2/s * mol2part


        !---  add fluxes to right hand side
        DO nn=1,DustBins

          DustEmis(j,i,DustInd(nn)) = FDust(j,i,nn)
#ifndef OFFLINE
          Flux(1,j,i,DustInd(nn))   = Flux(1,j,i,DustInd(nn)) + FDust(j,i,nn)/dz(1,j,i)
          !---  summarize biogenic and total emission rates
          EmiRate(EmiIndBio,j,i,DustInd(nn)) = EmiRate(EmiIndBio,j,i,DustInd(nn)) + FDust(j,i,nn)
          EmiRate(EmiIndSum,j,i,DustInd(nn)) = EmiRate(EmiIndSum,j,i,DustInd(nn)) + FDust(j,i,nn)
#endif
        END DO
      END DO
    END DO

    IF (lddebug) PRINT*, 'Leave emission_tegen',''//NEW_LINE('')

  END SUBROUTINE emission_tegen

  !+ get_ustar
  !---------------------------------------------------------------------
  SUBROUTINE get_ustar(subdomain)
  !---------------------------------------------------------------------
  ! Description:
  ! This subroutine provides the current friction velocity (ustar).
  !
  ! Methode 1:
  !   Calculate ustar with the wind profile method using the 10 m wind
  !   or the wind of the lowes model layer and the roughness length.
  !
  ! Methode 2:
  !   Use pre-calculated data of ustar.
  !--------------------------------------------------------------------

#ifdef OFFLINE
    USE offline_org
#endif

    IMPLICIT NONE

    TYPE(rectangle), INTENT(IN) :: subdomain

    INTEGER :: &
      i,j

    REAL(8) :: &
      uwind, &
      vwind, &
      tot_wind

    REAL(8), POINTER :: z0(:,:)
#ifndef OFFLINE
    REAL(8), POINTER :: dxK(:,:)
    REAL(8), POINTER :: dyK(:,:)
    REAL(8), POINTER :: dz(:,:,:)
    REAL(8), POINTER :: usur(:,:)
    REAL(8), POINTER :: vsur(:,:)
    REAL(8), POINTER :: qrsur(:,:)
    REAL(8), POINTER :: rhosur(:,:)
#endif

    REAL(8), PARAMETER :: &
      VK = 0.4

    z0       => dust(subdomain%ib)%z0(:,:)
#ifndef OFFLINE
    dxK     => geo  (subdomain%ib)%dxK(:,:)
    dyK     => geo  (subdomain%ib)%dyK(:,:)
    dz      => geo  (subdomain%ib)%dz(:,:,:)
    rhosur  => meteo(subdomain%ib)%rho(1,:,:,ScalCur)
    qrsur   => meteo(subdomain%ib)%QRSur(:,:,ScalCur)
    usur    => meteo(subdomain%ib)%u(1,:,:,WindCur)
    vsur    => meteo(subdomain%ib)%v(1,:,:,WindCur)
#endif

    IF (lddebug) PRINT*, 'Enter get_ustar'

    IF (fricvelo_scheme == 1) THEN
      DO i = 1,subdomain%ntx
        DO j = 1,subdomain%nty

#ifndef OFFLINE
          !---  flux initialisations
          uwind = usur(j,i+1)/dyK(j,i+1)+usur(j,i)/dyK(j,i)
          uwind = 0.5E0 * uwind / dz(1,j,i)
          vwind = vsur(j+1,i)/dxK(j+1,i)+vsur(j,i)/dxK(j,i)
          vwind = 0.5E0 * vwind / dz(1,j,i)
          tot_wind = SQRT(uwind**2+vwind**2)/rhosur(j,i)
#else
          tot_wind = SQRT(u(j,i,ntstep)**2+v(j,i,ntstep)**2)
#endif


          ustar(j,i) = (VK * tot_wind )/(log( dz(1,j,i)/(z0(j,i)/100.))) ! [m/s]
        END DO
      END DO
    ELSEIF (fricvelo_scheme == 2) THEN
      ustar(:,:) = ust(:,:,ntstep)

    END IF

    IF (lddebug) PRINT*, 'Leave get_ustar',''//NEW_LINE('')

  END SUBROUTINE get_ustar

  !+ okin_vegetation
  !---------------------------------------------------------------------
  SUBROUTINE okin_vegetation(subdomain,dimveg)
  !---------------------------------------------------------------------
  ! Description:
  ! Okin aprroach (G. Okin, 2008, DOI:10.1029/2007JF000758 )
  !
  ! Reduction of dust emission caused by vegetation.
  ! The fractional plant cover (FCOVER) (dust%veg(:,:,:))
  ! decreases the drag partition (dust%feff). A reduced
  ! drag partition leads to a decrease in the horizontal dust flux.
  !--------------------------------------------------------------------
#ifdef OFFLINE
    USE offline_org
#endif

    IMPLICIT NONE

    TYPE(rectangle), INTENT(IN) :: subdomain

    INTEGER, INTENT(IN) :: dimveg

    INTEGER :: &
      vegnow,  &  ! time loop
      i,j         ! loops

    REAL(8), PARAMETER :: pi = 3.14159265359

    REAL(8) :: &
      gap, gap1, gap2, &  ! gap size between plants
      gapheight,       &  ! funktion of gap size / plant height
      SSR                 ! shear stress ratio

    REAL(8), POINTER ::  &
      veg(:,:,:),        &
      vegmin(:,:),       &
      feff(:,:,:)!,       &


    veg     => dust(subdomain%ib)%veg(:,:,:)
    vegmin  => dust(subdomain%ib)%vegmin2(:,:)
    feff    => dust(subdomain%ib)%feff(:,:,:)

    IF (lddebug) PRINT*, 'Enter okin_vegetation'

    ! start lon-lat-loop
    DO i=1,subdomain%ntx
      DO j=1,subdomain%nty
        DO vegnow=1,dimveg

         ! -- gap size/height --
         ! Formula of okin depends on the fraction of
         !   gap between the plants(L) and plant height (h) --> L/h.

         ! -- scale up Fcover a bit becaus of the offset FCOVER_MIN = 0.1 !TEST maybe dump
         ! veg(j,i,vegnow) = 10./9. * veg(j,i,vegnow) - 1./9.
         IF (lvegmin) veg(j,i,vegnow) = 1./(1-vegmin(j,i)) * veg(j,i,vegnow) - 1./(1-vegmin(j,i))*(vegmin(j,i))
         IF(veg(j,i,vegnow) < 0.) veg(j,i,vegnow)=0.
         IF(veg(j,i,vegnow) > 1.) veg(j,i,vegnow)=1.

         ! -- paramet for gap size (Matthias Faust) --
         ! Uniform plants with r = 0.05 m are distributed as grid over the whole Area.
         ! The total plant area is equal to the fractional plant cover (FCOVER).
         ! The gap between the plants is calc as the mean of column gap and diagonal gap.

         ! cover should be bigger than one plant
         IF (veg(j,i,vegnow) >= pi*0.05**2 ) THEN
           gap1=1./(sqrt((veg(j,i,vegnow))/(pi*0.05**2))-1.) -2.*0.05
           gap2=sqrt( 2.*(1./(sqrt((veg(j,i,vegnow))/(pi*0.05**2))-1.))**2) -2.*0.05
           gap = (gap1+gap2)/2.

           ! -- paramet for plant height (Matthias Faust) --
           ! Okin is only calculated on agricultural areas.
           ! Most common crop is wheat. At vegetation maximum, wheat is approx 1m height
           !   and the fractional plant cover is nearly 1.
           ! Because of this, FCOVER is used as plant height.

           gapheight = gap/(veg(j,i,vegnow) + 1.E-8) ! + 1.E-8 to avoid division by zero

           ! -- paramet of the SSR (Okin)
           ! There is a constant c1 with the best fit value of c1=4.8 by okin.
           ! In the aprroach there is a minimum SSR. For Okins experiments the best fit of
           !   the best fit of SSR0 is 0.32. However I (MF) assume for dense vegetated
           !   crop land SSR0 = 0.
           ! With this asumptions the Formula for the SSR is:

           SSR = 1. - 1./((1./4.8)*gapheight + 1.)

           ! the drag partition is defined as sqrt(SSR)
           feff(j,i,vegnow) = SQRT(SSR)


         ELSE ! if cover is smaller than one plant
           feff(j,i,vegnow) = 1
         ENDIF

         ! This should not happen, but just in case
         IF(feff(j,i,vegnow) < 0.) feff(j,i,vegnow)=0.
         IF(feff(j,i,vegnow) > 1.) feff(j,i,vegnow)=1.

        END DO
      END DO
    END DO
    ! end lon-lat-loop

    IF (lddebug) PRINT*, 'Leave okin_vegetation',''//NEW_LINE('')

  END SUBROUTINE okin_vegetation

  !+ linear_vegetation
  !---------------------------------------------------------------------
  SUBROUTINE linear_vegetation(subdomain,dimveg)
  !---------------------------------------------------------------------
  ! Description:
  ! Based on the fractional vegetation cover (dust%veg(:,:,:))
  ! an effective fraction for dust emission is calculated (dust%veff).
  ! veff=1 -> full emission, veff=0 -> no emission
  ! veff reached its minimum when the threshold vegetation cover is reached.
  ! The threshold is defined as veg_lim in dust_tegen_param (data_dust.f90)
  ! Theoretically, a dependency on cultivation class and biome is included.
  ! At the moment the cultivation class is constantly set cult = 3 (???)
  ! For the biomes, the model distinguishes between shrubland and others.
  !
  ! The code based on Tegen et al. 2002
  ! https://doi.org/10.1029/2001JD000963
  !--------------------------------------------------------------------
    USE dust_tegen_param
    USE dust_tegen_data

#ifdef OFFLINE
    USE offline_org
#endif

    IMPLICIT NONE

    TYPE(rectangle), INTENT(IN) :: subdomain

    INTEGER, INTENT(IN) :: dimveg

    INTEGER :: &
      vegnow,  &  ! time loop
      i,j,     &  ! loops
      idust,   &  ! flag for dust active biome
      vflag,   &  ! flag for vegetation above threshold
      vgtp        ! index of biome

    REAL(8) ::  &
      veg_max (subdomain%nty,subdomain%ntx),  & !maximum veg
      veg_min (subdomain%nty,subdomain%ntx)     !minimum veg

    REAL(8), POINTER ::  &
      veg(:,:,:),        &
      cult(:,:),         &
      biome(:,:),        &
      veff(:,:,:)!,       &

    veg   => dust(subdomain%ib)%veg(:,:,:)
    cult  => dust(subdomain%ib)%cult(:,:)
    biome => dust(subdomain%ib)%biome(:,:)
    veff  => dust(subdomain%ib)%veff(:,:,:)

    IF (lddebug) PRINT*, 'Enter linear_vegetation'

    ! start lon-lat-loop
    DO i=1,subdomain%ntx
      DO j=1,subdomain%nty

        ! Test on the vegetation type
        cult(j,i)=3.
        IF(cult(j,i) < 0. .OR. cult(j,i) > 27.) cult(j,i)=3.

        idust = 0
        vgtp = int(biome(j,i))
        IF (vgtp > 0) idust = active(vgtp)
        IF(cult(j,i) <= 2 .AND. biome(j,i) >= 1 .AND. biome(j,i) <= 27) idust=1


        IF (idust == 1) THEN

          ! search for max LAI:
          veg_max(j,i)=0.
          veg_min(j,i)=1.
          DO vegnow=1,dimveg
            IF(veg(j,i,vegnow) > veg_max(j,i)) veg_max(j,i) = veg(j,i,vegnow)
            IF(veg(j,i,vegnow) < veg_min(j,i)) veg_min(j,i) = veg(j,i,vegnow)
            veff(j,i,vegnow)=0.
          END DO

          ! check if any monthly LAI is above threshold (limit=veg_lim2)
          vflag=0
          DO vegnow=1,dimveg
            IF (veg(j,i,vegnow) > veg_lim2) vflag=1
          END DO

          IF(cult(j,i) <= 2 ) vflag=0

          !---------------------------------------------------------------------------------------
          !--  Calculate effective surface for VEG <veg_lim (as proxy for
          !--  veg. cover), shrubby vegetation is determined by max
          !--  annual LAI, grassy by monthly LAI. Relative shrub area is
          !--  set for active biomes in b2000_veg.data
          !---------------------------------------------------------------------------------------
          DO vegnow=1,dimveg
            ! biome 13 Tropical xerophytic shrubland  14 Temperate xerophytic shrubland
            IF((biome(j,i) >= 13.) .AND. (biome(j,i) <= 14.) .AND. (vflag == 0)) THEN
             veff(j,i,vegnow)=1.-(veg(j,i,vegnow))*1./veg_lim
            ELSE
             veff(j,i,vegnow)= 1.-( veg(j,i,vegnow) * (1.-shrub(int(biome(j,i))))  &
                               +veg_max(j,i)*shrub(int(biome(j,i))) )*1./veg_lim

            END IF


            ! veff(j,i,vegnow)=1.-(veg(j,i,vegnow))/veg_lim  ! grass only (IT)
            ! veff(j,i,vegnow)=1.-(veg_max(j,i))/veg_lim     !shrub only (IT)


            ! veff(j,i,vegnow)=(veg(j,i,vegnow)-1)**4        ! linear dependensy of veg (MF)

            IF(veff(j,i,vegnow) <= 0.) veff(j,i,vegnow) = 0.
            IF(veff(j,i,vegnow) >= 1.) veff(j,i,vegnow) = 1.

          END DO
        END IF
      END DO
    END DO
    ! end lon-lat-loop

    IF (lddebug) PRINT*, 'Leave linear_vegetation',''//NEW_LINE('')

  END SUBROUTINE linear_vegetation


  !+ roughness
  !---------------------------------------------------------------------
  SUBROUTINE roughness(subdomain,dimsimu)
  !---------------------------------------------------------------------
  ! Description:
  !
  ! Reduction of dust emission caused by roughness elements.
  ! A rough surface (dust%z0) decreases the drag partition
  ! (dust%feff). A reduced drag partition leads to a decrease
  ! in the horizontal dust flux.
  !
  ! The Physics based on the paper of Marticorena and Bergametti 1995
  ! https://doi.org/10.1029/95JD00690
  !
  ! The code based on Tegen et al. 2002
  ! https://doi.org/10.1029/2001JD000963
  !--------------------------------------------------------------------

    USE dust_tegen_param
    USE dust_tegen_data

#ifdef OFFLINE
    USE offline_org
#endif

    IMPLICIT NONE

    TYPE(rectangle), INTENT(IN) :: subdomain

    INTEGER, INTENT(IN) :: dimsimu

    INTEGER :: &
      dnow,    &  ! time loop
      i,j!,     &  ! loops


    REAL(8) ::  &
      local_feff,        &  ! feff inside the loop
      z0s,               &  ! small scale roughness length
      AAA,               &
      BB,                &
      CCC,               &  ! dummys
      DDD,               &
      EE,                &
      FF


    REAL(8), POINTER ::  &
      z0(:,:),           &
      feff(:,:,:)!,       &

    z0   => dust(subdomain%ib)%z0(:,:)
    feff => dust(subdomain%ib)%feff(:,:,:)

    IF (lddebug) PRINT*, 'Enter roughness'

    ! z0s  = dp / 30.
    z0s  = 0.001 !! en cm, these Marticorena p.85
    ! d1   = 0.  !100
    local_feff = 0.

    ! start lon-lat-loop
    DO i=1,subdomain%ntx
      DO j=1,subdomain%nty
      ! z0 and efficient fraction feff
      ! partition of energy between the surface and the elements of rugosity, these pp 111-112

      IF (z0(j,i) <= 0.) THEN
       z0(j,i) = 0.
       local_feff = 1.
      ELSE
       AAA = log(z0(j,i)/z0s)
       BB = log(aeff*(xeff/z0s)**0.8)
       CCC = 1.- AAA/BB

       !         * partition between Z01 and Z02 *
       ! IF (d1 == 0) THEN
       !   FF = 1.
       ! ELSE
       !   DDD = log(z0(j,i)/z0(j,i))   !Z02/Z01
       !   EE = log(aeff * (d1/z0(j,i))**0.8)
       !   FF = 1.- DDD/EE
       ! END IF

        FF = 1.

        ! fraction efficace totale
        local_feff = FF*CCC


        IF (local_feff < 0.) local_feff = 0.
        IF (local_feff > 1.) local_feff = 1.
      END IF
      feff(j,i,:) = local_feff
      END DO
    END DO
    ! end lon-lat-loop
    call quick_nc('feff',var2d=feff(:,:,1))

    IF (lddebug) PRINT*, 'Leave roughness'//NEW_LINE('')

  END SUBROUTINE roughness


  !+ roughness
  !---------------------------------------------------------------------
  SUBROUTINE fecan(yaction,subdomain,time_now)
  !---------------------------------------------------------------------
  ! Description:
  !
  !--------------------------------------------------------------------

    USE dust_tegen_param
    USE dust_tegen_data
    USE mo_dust

#ifdef OFFLINE
    USE offline_org
#endif

    IMPLICIT NONE

    TYPE(rectangle), INTENT(IN) :: subdomain

    INTEGER, INTENT(IN) :: &
      time_now

    CHARACTER(LEN=*), INTENT(IN)            :: &
      yaction ! action to be performed

    INTEGER :: &
      i,j      ! loops

    REAL(8) :: &
      clay ,   &  ! soil clay contant [%]
      moist       ! gravimeric soil moisture [%]


    REAL(8), POINTER ::  &
      vmoist(:,:,:),           &
      mfac(:,:),       &
      w_str(:,:),       &
      soilmap(:,:,:)!,     &

    vmoist   => dust(subdomain%ib)%vmoist(:,:,:)
    mfac => dust(subdomain%ib)%mfac(:,:)
    w_str => dust(subdomain%ib)%w_str(:,:)
    soilmap => dust(subdomain%ib)%soilmap(:,:,:)

    IF (lddebug) PRINT*, 'Enter fecan, yaction=',yaction

    IF (yaction == 'init') THEN
      ! call quick_nc('vmoist',var2dtime=vmoist)

    ! start lon-lat-loop
      DO i=1,subdomain%ntx
        DO j=1,subdomain%nty
          IF (soilmaptype == 1) THEN
            clay = soilmap(j,i,4) * 100.
          ELSE IF (soilmaptype == 2) THEN
            clay = soilmap(j,i,3) * 100.
          END IF

          ! Calculation of the threshold soil moisture (w')
          w_str(j,i) = 0.0014*clay**2 + 0.17*clay

        END DO
      END DO

      ! end lon-lat-loop

    ELSEIF (yaction == 'calc') THEN ! yaction == 'init'

    ! start lon-lat-loop
      DO i=1,subdomain%ntx
        DO j=1,subdomain%nty

          ! calculate gravimeric soil moisture from the volumeric soil moisture
          ! moist_g = moist_v * rho_water / rho_soil * 100%   -> [kg_water/m3_soil / kg_soil/m3_soil] = [%]
          moist = vmoist(j,i,time_now) * 1000/2650 * 100


          IF (moist <= w_str(j,i)) THEN
            mfac(j,i) = 1
          ELSE
            mfac(j,i) = (1 + 1.21 * ( moist - w_str(j,i))**0.68 )**0.5
          END IF

        END DO
      END DO

    END IF ! yaction == 'calc'

    IF (lddebug) PRINT*, 'Leave fecan, yaction=',yaction,''//NEW_LINE('')

  END SUBROUTINE fecan

  !+ read_ascii
  !---------------------------------------------------------------------
  SUBROUTINE read_ascii(filename,outvar)
  !---------------------------------------------------------------------
  ! Description:
  !   This subroutine read ascii input for the dust emissions
  !
  !   Input file needs a minimum header of 2 lines
  !   the lines BEGIN_DATA and VARNAME need to be
  !   directly above the values
  !
  !   Minimium header:
  !     BEGIN_DATA
  !     VARNAME
  !     0.0
  !     0.0
  !     ...
  !
  !
  !   Additional header lines can specify the pol location and the model Domain
  !   (like the 'old' muscat input files)
  !   If this additional information is given in the input file
  !   then the valuses have to be consistent with the model domain
  !
  !   Extended header:
  !     #--------------------------------------------------------------
  !     # pollon pollat
  !     #--------------------------------------------------------------
  !     POL
  !     -170.000000 40.000000
  !     #
  !     #--------------------------------------------------------------
  !     # i startlon_tot dlon
  !     # j startlat_tot dlat
  !     #--------------------------------------------------------------
  !     KOORD
  !     42 -2.1 0.1
  !     42 -2.1 0.1
  !     #--------------------------------------------------------------
  !     #                 Daten
  !     #--------------------------------------------------------------
  !     BEGIN_DATA
  !     SOIL
  !     0.0
  !     0.0
  !
  !   The header can include even more information,
  !   only POL and KOORD are requested as keywords.
  !   The rest will be ignored
  !
  !---------------------------------------------------------------------
    ! USE data_modelconfig, ONLY: &
    !   ie_tot,       & ! number of grid points in zonal direction
    !   je_tot,       & ! number of grid points in meridional direction
    !   pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    !   pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
    !   dlon,         & ! grid point distance in zonal direction (in degrees)
    !   dlat,         & ! grid point distance in meridional direction (in degrees)
    !   startlon_tot, & ! transformed longitude of the lower left grid point
    !                   ! of the total domain (in degrees, E>0)
    !   startlat_tot    ! transformed latitude of the lower left grid point
    !                   ! of the total domain (in degrees, N>0)

    IMPLICIT NONE

    ! patameters
    CHARACTER(*), INTENT(IN) :: &
      filename

    REAL(8), INTENT(INOUT)   :: &
      outvar(:,:,:)

    ! local vars
    INTEGER        :: &
      unit,           &   ! unit number for file
      ios,            &   ! error code for open and closing files
      eof,            &   ! error code when eof is reached
      i,j,t,n,        &   ! loop
      nlines,         &
      nhlines

    REAL(8)           :: &
      read_ie_tot,       & ! read form input: number of grid points in zonal direction
      read_je_tot,       & ! read form input: number of grid points in meridional direction
      read_pollon,       & ! read form input: longitude of the rotated north pole (in degrees, E>0)
      read_pollat,       & ! read form input: latitude of the rotated north pole (in degrees, N>0)
      read_dlon,         & ! read form input: grid point distance in zonal direction (in degrees)
      read_dlat,         & ! read form input: grid point distance in meridional direction (in degrees)
      read_startlon_tot, & ! read form input: transformed longitude of the lower left grid point
                           ! of the total domain (in degrees, E>0)
      read_startlat_tot    ! read form input: transformed latitude of the lower left grid point
                           ! of the total domain (in degrees, N>0)

    CHARACTER(120)  :: &
      hline                ! single line of header
    CHARACTER(120), ALLOCATABLE :: &
      header(:)            ! header of input file

    IF (lddebug) PRINT*, 'Enter read_ascii'

    ! Print short status
    PRINT*,'read_ascii ', filename

    ! - 1) Open file
    OPEN(newunit=unit, file=filename, iostat=ios,status="old", action="read")
    IF ( ios /= 0 ) stop "Error opening file "

    ! - 2) Read header
    nlines = 0
    nhlines = 0
    !Read line by line and count number of lines in the file
    DO
      ! read single line
      READ(unit,*,iostat=eof) hline
      ! leave loop at eof
      IF ( eof /= 0 ) EXIT
      ! count number of lines in the file
      nlines = nlines+1
      ! search for the keyword 'BEGIN_DATA' which marked the end of the header
      IF (TRIM(hline) == 'BEGIN_DATA') nhlines=nlines+1
    END DO

    ! ERROR if 'BEGIN_DATA' was not found
    IF (nhlines == 0) THEN
      PRINT*,'ERROR:   Reading file: ',filename
      PRINT*,'         Unexpected EOF'
      PRINT*,'         Tag BEGIN_DATA missing in header'
      STOP 'ERROR reading ascii input for dust emission'
    END IF

    ! Allocate header with n+1 lines
    AllOCATE (header(nhlines))

    ! jump back to the top of the file
    REWIND(unit)

    ! read  header
    READ(unit,'(A)') (header(i),i=1,nhlines)

    ! Check header information
    !   if POL and KOORD infomation is stored in the header
    !   check if this data is correct (abort if not)
    !   if POL or KOORD is missing, print a warning (no abort)

    ! first search for POL
    IF (ANY(header == 'POL')) THEN
      DO i=1,nhlines
        IF (header(i) == 'POL') THEN
          ! read in the data
          READ (header(i+1),*) read_pollon, read_pollat
          !checking
          IF (read_pollon /= pollon .OR. read_pollat /= pollat) THEN
            PRINT*,'ERROR:   Wrong information in file: ',filename
            PRINT*,'         Wrong values at POL tag'
            STOP 'ERROR reading ascii input for dust emission'
          END IF
          !leave loop
          EXIT
        END IF
      END DO
    ELSE
      PRINT*,'WARNING:   Missing information in file: ',filename
      PRINT*,'           No POL tag was found'
    END IF !ANY(header == 'POL'))

    ! search for KOORD
    IF (ANY(header == 'KOORD')) THEN
      DO i=1,nhlines
        IF (header(i) == 'KOORD') THEN
          ! read in the data
          READ (header(i+1),*) read_ie_tot, read_startlon_tot, read_dlon
          READ (header(i+2),*) read_je_tot, read_startlat_tot, read_dlat
          !checking
          IF (read_ie_tot /= ie_tot .OR. read_je_tot /= je_tot                         .OR. &
              read_startlon_tot /= startlon_tot .OR. read_startlat_tot /= startlat_tot .OR. &
              read_dlon /= dlon .OR. read_dlat /= dlat                                  ) THEN
            PRINT*,'ERROR:   Wrong information in file: ',filename
            PRINT*,'         Wrong values at KOORD tag'
            STOP 'ERROR reading ascii input for dust emission'
          END IF
          !leave loop
          EXIT
        END IF
      END DO
    ELSE
      PRINT*,'WARNING:   Missing information in file: ',filename
      PRINT*,'           No KOORD tag was found'
    END IF !ANY(header == 'KOORD'))


    ! - 3) Read Data

    ! check size
    ! not enough values
    IF (nlines-nhlines < ie_tot*je_tot) THEN
      PRINT*,'ERROR:   Reading file: ',filename
      PRINT*,'         not enough values'
      STOP 'ERROR reading ascii input for dust emission'
    ! Too much values, + 1 because last line may be 'END_DATA'
    ELSEIF (nlines-nhlines > ie_tot*je_tot +1) THEN
      PRINT*,'ERROR:   Reading file: ',filename
      PRINT*,'         too much values'
      STOP 'ERROR reading ascii input for dust emission'
    END IF


    ! read data
    READ(unit,*) (((outvar(j,i,t),t=1,size(outvar,3)),i=1,ie_tot),j=1,je_tot)


    DEALLOCATE(header)

    CLOSE(unit)
    IF ( ios /= 0 ) STOP "Error closing file"

    IF (lddebug) PRINT*, 'Leave read_ascii',''//NEW_LINE('')

  END SUBROUTINE read_ascii



  SUBROUTINE read_infile(infile,outvar,ierror,yerrmsg)

    ! Modules
    USE netcdf
    USE mo_dust
#ifndef OFFLINE
    USE data_io,  ONLY : ydate_ini    ! start of the forecast
#endif
    ! USE data_modelconfig,   ONLY :   &
    !   ie_tot,                  & ! number of grid points in zonal direction
    !   je_tot                     ! number of grid points in meridional direction


    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)    :: &
      infile        ! name of field

    REAL(8),        INTENT(INOUT) :: &
      outvar(:,:,:)        ! Output var

    INTEGER,        INTENT(OUT)   :: &
      ierror

    CHARACTER(*),   INTENT(OUT)   :: &
      yerrmsg       ! error message


    ! local Vars

    CHARACTER(15)     :: &
      varname(3)       ! name of variable that shoud be read

    CHARACTER(120)     :: &
      filename       ! name of variable that shoud be read

    INTEGER            :: &
      ntimes         ! number of possible days in the model run
    LOGICAL           :: &
      timecheck     ! check if time fits model run


    INTEGER  :: &
      i,j,t,iv, &  ! loop
      ary_size, &  ! array size
      varnum,   &  !
      istart,   &  ! index of the start date in the var file
      idate,    &  ! ydate read into integer
      istat,    &  ! local error code
      ncID,     &  ! id of nc files
      varID,    &  ! id of the  var
      timeID,   &  ! id of time var
      dimID,    &  ! id of a dimension
      dimlen       ! length of the above dimension

    REAL  :: &
      var_scale, &    ! index of the start date in the var file
      var_offset    ! index of the start date in the var file

    INTEGER, ALLOCATABLE :: &
      times(:)

    REAL, ALLOCATABLE :: &
      var_read(:,:,:)


    CHARACTER (len = 10) :: &
      lon_names(4), & ! common names for the longitudes
      lat_names(4)    ! common names for the latitudes


      lon_names(1)='x'
      lon_names(2)='lon'
      lon_names(3)='rlon'
      lon_names(4)='longitude'

      lat_names(1)='y'
      lat_names(2)='lat'
      lat_names(3)='rlat'
      lat_names(4)='latitude'

    ! start subroutine
    ! ---------------------------------------------------------


    IF (lddebug) PRINT*, 'Enter read_infile'

    ! Definitions
    varnum = 1
    ntimes = 0
    istart = 0

    IF (infile == 'soil') THEN
      filename = TRIM(soiltypeFile)
      IF (soilmaptype == 1) THEN
        varname = 'soiltype'
      ELSEIF (soilmaptype == 2) THEN
        varnum = 3
        varname(1)='sand'
        varname(2)='silt'
        varname(3)='clay'
      END IF
    ELSEIF (infile == 'source') THEN
      filename = TRIM(psrcFile)
	    varname = 'source'
    ELSEIF (infile == 'z0') THEN
      filename = TRIM(z0File)
      varname = 'z0'
    ELSEIF (infile == 'moist') THEN
      filename = TRIM(moistFile)
      varname='swvl1'
      ntimes = SIZE(outvar(1,1,:))
      timecheck = .FALSE.
    ELSEIF (infile == 'vegday') THEN
      filename = TRIM(vegdayFile)
      varname='FCOVER'
      ntimes = SIZE(outvar(1,1,:))
      timecheck = .True.
    END IF

    IF (filename(LEN(TRIM(filename))-2:) /= '.nc') THEN
      CALL read_ascii(TRIM(filename),outvar)
      RETURN
    ENDIF

    ! Print short status
    PRINT*,'read_nc ', TRIM(filename)

    ! open the file
    istat = nf90_open(filename, nf90_nowrite, ncid)
    IF (istat /= nf90_noerr) THEN
      ierror  = 10001
      yerrmsg = TRIM(nf90_strerror(istat))
      RETURN
    ENDIF

    ! get id of lat dimension
    DO i = 1, size(lat_names)
      istat = nf90_inq_dimid(ncID, TRIM(lat_names(i)), dimID)
      IF (istat == nf90_noerr) EXIT
    END DO
    IF (istat /= nf90_noerr) THEN
      ierror  = 10002
      PRINT*,''
      PRINT*,'ERROR'
      PRINT*,'  No latitude dimension was found in ',infile
      PRINT*,'  possible names for the latitude dimension:'
      PRINT*,'  ',lat_names
      PRINT*, ''
      yerrmsg = TRIM(nf90_strerror(istat))
      RETURN
    ENDIF

    ! read the length of lat dimension
    istat = nf90_inquire_dimension(ncID, dimID, len = dimlen)
    IF (istat /= nf90_noerr) THEN
      ierror  = 10003
      yerrmsg = TRIM(nf90_strerror(istat))
      RETURN
    ENDIF

    ! Check the Size
    IF (dimlen /= je_tot) THEN
      ierror = 10004
      yerrmsg = 'Error reading '//TRIM(infile)//' file: wrong lat dimension'
      RETURN
    END IF

    ! get id of lon dimension
    DO i = 1, size(lon_names)
      istat = nf90_inq_dimid(ncID, TRIM(lon_names(i)), dimID)
      IF (istat == nf90_noerr) EXIT
    END DO
    IF (istat /= nf90_noerr) THEN
      ierror  = 10005
      PRINT*,''
      PRINT*,'ERROR'
      PRINT*,'  No longitude dimension was found in ',infile
      PRINT*,'  possible names for the longitude dimension:'
      PRINT*,'  ',lon_names
      PRINT*, ''
      yerrmsg = TRIM(nf90_strerror(istat))
      RETURN
    ENDIF

     ! read the length of lon dimension
     istat = nf90_inquire_dimension(ncID, dimID, len = dimlen)
     IF (istat /= nf90_noerr) THEN
       ierror  = 10006
       yerrmsg = TRIM(nf90_strerror(istat))
       RETURN
     ENDIF

     ! Check the Size
     IF (dimlen /= ie_tot) THEN
       ierror = 10007
       yerrmsg = 'Error reading '//TRIM(infile)//' file: wrong lon dimension'
       RETURN
     END IF

    IF (ntimes > 0) THEN
      ! get id of time dimension
      istat = nf90_inq_dimid(ncID, 'time', dimID)
      IF (istat /= nf90_noerr) THEN
        ierror  = 10008
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF

      ! read the length of the time dimension
      istat = nf90_inquire_dimension(ncID, dimID, len = dimlen)
      IF (istat /= nf90_noerr) THEN
        ierror  = 10009
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF

      ! allocate the var ( times ) that hold the available dates in the nc file
      istat = 0
      ALLOCATE (times(dimlen), STAT=istat)
      IF (istat /= 0) THEN
        ierror = 10010
        yerrmsg = 'allocation of times failed'
        RETURN
      ENDIF

      ! get the id of the time var
      istat = nf90_inq_varid(ncID, 'time', timeID)
      IF (istat /= nf90_noerr) THEN
        ierror  = 10005
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF

      ! get the dates in the nc file
      istat = nf90_get_var(ncid, timeID, times)
      IF (istat /= nf90_noerr) THEN
        ierror  = 10011
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF

      ! transform date string to an integer
      READ(ydate_ini(1:8),*) idate

      ! istart -> fist read time step
      istart = 1

      ! match the dates of the model run and the var
      IF (timecheck) THEN
        IF (ANY(times(:) == idate)) THEN
          ! find index of istat
          DO WHILE(times(istart) /= idate)
            istart = istart + 1
          END DO
          IF (istart + ntimes - 1 > size(times)) THEN
            ierror  = 10012
            yerrmsg = 'Error reading '//TRIM(infile)//' file: not enough dates in the file'
            RETURN
          END IF
        ELSE
          ierror  = 10013
          yerrmsg = 'Error reading '//TRIM(infile)//' file: model start date is not in the file'
          RETURN
        END IF
      END IF
    END IF
    ! Allocate var to read
    istat = 0

    ary_size = ntimes
    IF ( ntimes == 0 ) ary_size = 1
    IF ( varnum  > 0 ) ary_size = varnum
    IF ( ntimes  > 0 ) ary_size = ntimes

    ALLOCATE (var_read(ie_tot,je_tot,ary_size), STAT=istat)
    IF (istat /= 0) THEN
      ierror = 10014
      yerrmsg = 'allocation of var_read failed'
      RETURN
    ENDIF

    ! loop for all vars
    DO iv = 1, varnum
      ! get the id of the var
      istat = nf90_inq_varid(ncID, varname(iv) , varID)
      IF (istat /= nf90_noerr) THEN
        ierror  = 10015
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF

      ! get value of the scaling factor
      istat = nf90_get_att(ncID, varID, 'scale_factor',var_scale)
      IF (istat /= nf90_noerr) THEN
        var_scale = 1
      ENDIF
      istat = nf90_get_att(ncID, varID, 'scale',var_scale)
      IF (istat /= nf90_noerr) THEN
        var_scale = 1
      ENDIF

      ! get value of the offset
      istat = nf90_get_att(ncID, varID, 'add_offset',var_offset)
      IF (istat /= nf90_noerr) THEN
        ! ierror  = 10016
        var_offset = 0
        ! yerrmsg = TRIM(nf90_strerror(istat))
        !RETURN
      ENDIF

      ! get the var
      !IF (istart < 1) istart = 1
      istat = nf90_get_var(ncid, varID, var_read(:,:,iv),start= (/1,1,istart/),count=(/ie_tot,je_tot,ntimes/))
      IF (istat /= nf90_noerr) THEN
        ierror  = 10017
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF

      ! wirte the var into outvar
      DO i=1,ie_tot
        DO j=1,je_tot
          IF (varnum == 1) THEN
            outvar(j,i,:)=var_read(i,j,:)*var_scale+var_offset
          ELSEIF (varnum > 1) THEN
            outvar(j,i,iv)=var_read(i,j,iv)*var_scale+var_offset
          END IF
        END DO
      END DO

    END DO

    DEALLOCATE (var_read)
    IF (ntimes > 0) DEALLOCATE (times)

    IF (lddebug) PRINT*, 'Leave read_infile, ierr=',ierror,''//NEW_LINE('')

  END SUBROUTINE read_infile


  !+ copy2block
  !---------------------------------------------------------------------
  SUBROUTINE copy2block(subdomain,dim,from, too2d,too3d)
  !---------------------------------------------------------------------
  ! Description:
  !   This subroutine copy the input data into
  !     the block struckture of muscat
  !
  !   This subroutine is basicly the subroutine 'soilprop_copy'
  !     only a few changes where made.
  !   The subroutine copys the data from the var 'from(je_tot,ie_tot,:)'
  !     too the var 'too(je_sub,ie_sub,:)'
  !     for 2d (if dim=2) and 3d (if dim=3)
  !----------------------------------------------------------------

#ifdef OFFLINE
    USE offline_org
#endif

    IMPLICIT NONE

    ! patameters
    TYPE (rectangle) subdomain      ! actual subdomain (INTENT(IN))

    INTEGER, INTENT(IN)   :: &
      dim                           ! dimension of copy var (dim = 2 .or. 3 )

    REAL(8), INTENT(IN)   :: &
      from(:,:,:)                   ! var that will be copied

    REAL(8), OPTIONAL, INTENT(INOUT)   :: &
      too2d(:,:)                    ! var that will store the copy in block structur (2d)

    REAL(8), OPTIONAL, INTENT(INOUT)   :: &
      too3d(:,:,:)                  ! var that will store the copy in block structur (3d)


    INTEGER :: i,j,k,jx,jy,in,jn,RefLoc  ! loops and
    INTEGER :: igx0,igy0,SurfLoc         ! helping integers

    IF (lddebug) PRINT*, 'Enter copy2block'

 !----------------------------------------------------------------
 !---  Fine Grid
 !----------------------------------------------------------------
    RefLoc  =  2.e0**abs(subdomain%refine-SurfLevel)
    SurfLoc =  2.e0**SurfLevel

    IF (subdomain%refine>=SurfLevel) THEN
      in=0
      DO i=subdomain%ix0+1,subdomain%ix1
        in=in+1
        jx=(i+ReFloc-1)/RefLoc
        jn=0
        DO j=subdomain%iy0+1,subdomain%iy1
          jn=jn+1
          jy=(j+ReFloc-1)/RefLoc
          IF (dim == 2) too2d(jn,in) = from(jy,jx,1)
          IF (dim == 3) too3d(jn,in,:) = from(jy,jx,:)
        END DO
      END DO
 !----------------------------------------------------------------
 !---  Coarser Grid
 !----------------------------------------------------------------
    ELSE
      igx0 = SurfLoc * subdomain%igx0
      in=0
      DO i=subdomain%ix0+1,subdomain%ix1
        in=in+1
        igy0 = SurfLoc * subdomain%igy0
        jn=0
        DO j=subdomain%iy0+1,subdomain%iy1
          jn=jn+1
          DO jx=igx0+1,igx0+RefLoc
            DO jy=igy0+1,igy0+RefLoc
              IF (dim == 2) too2d(jn,in) = from(jy,jx,1)
              IF (dim == 3) too3d(jn,in,:) = from(jy,jx,:)
            END DO
          END DO
          igy0 = igy0 + RefLoc
         END DO

         igx0 = igx0 + RefLoc
      END DO
    END IF

    IF (lddebug) PRINT*, 'Leave copy2block',''//NEW_LINE('')

  END SUBROUTINE copy2block


  SUBROUTINE quick_nc(name,var2d,var2dtime,var3d,var3dtime)
  !---------------------------------------------------------------------
  ! Description:
  !   quick output of data fields in netcdf format.
  !   Mainly for debugging.
  !
  ! Usage:
  !   call quick_nc('name',vartype=printvar)
  !   name definens the filename and the variable name in the nc file
  !   possible variable types
  !   var2d     : a two dimensional field                         -> printvar(:,:)
  !   var2dtime : a two dimensional field with a time dimension   -> printvar(:,:,:)
  !   var3d     : a three dimensional field                       -> printvar(:,:,:)
  !   var3dtime : a three dimensional field with a time dimension -> printvar(:,:,:,:)
  !
  !----------------------------------------------------------------

    USE mo_dust
    USE netcdf
    USE offline_org

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: name

    REAL(8), OPTIONAL,  INTENT(IN) :: &
      var2d    (:,:)  ,  &
      var2dtime(:,:,:),  &
      var3d    (:,:,:),  &
      var3dtime(:,:,:,:)

    REAL(8), ALLOCATABLE :: &
      print2d    (:,:)  ,  &
      print2dtime(:,:,:),  &
      print3d    (:,:,:),  &
      print3dtime(:,:,:,:)

    CHARACTER(30) :: &
      fname

    INTEGER :: &
      i,j, &
      istat,&
      ncID,&
      xID,&
      yID,&
      zID,&
      tID,&
      varID

    INTEGER, ALLOCATABLE :: &
      varshape(:)

    IF (lddebug) PRINT*, 'Enter quick_nc'

    xID=0
    yID=0
    zID=0
    tID=0

    ! Define filename
    fname=name//'.nc'

    ! create file
    IF (present(var2d) .OR. present(var3d) .OR. present(var2dtime) .OR. present(var3dtime) ) THEN
     print*,"create quick_nc file: ",fname
     istat=nf90_create(TRIM(fname), NF90_SHARE, ncID)
     if(istat /= nf90_NoErr) print*,'quick_nc fail #1 ',nf90_strerror(istat)

    ELSE
     PRINT*, 'quick_nc data field is missing, return'
     RETURN
    end if

    ! get dimension of the print variable
    IF (present(var2d)) THEN
      ALLOCATE(varshape(2))
      varshape = shape(var2d)
      AllOCATE(print2d(varshape(2),varshape(1)))
      DO i = 1, varshape(2)
        DO j = 1, varshape(1)
          print2d(i,j) = var2d(j,i)
        END DO
      END DO
    ELSEIF(present(var2dtime)) THEN
      ALLOCATE(varshape(3))
      varshape = shape(var2dtime)
      AllOCATE(print2dtime(varshape(2),varshape(1),varshape(3)))
      DO i = 1, varshape(2)
        DO j = 1, varshape(1)
          print2dtime(i,j,:) = var2dtime(j,i,:)
        END DO
      END DO
    ELSEIF (present(var3d)) THEN
      ALLOCATE(varshape(3))
      varshape = shape(var3d)
      AllOCATE(print3d(varshape(2),varshape(1),varshape(3)))
      DO i = 1, varshape(2)
        DO j = 1, varshape(1)
          print3d(i,j,:) = var3d(j,i,:)
        END DO
      END DO
    ELSEIF (present(var3dtime)) THEN
      ALLOCATE(varshape(4))
      varshape = shape(var3dtime)
      AllOCATE(print3dtime(varshape(2),varshape(1),varshape(3),varshape(4)))
      DO i = 1, varshape(2)
        DO j = 1, varshape(1)
          print3dtime(i,j,:,:) = var3dtime(j,i,:,:)
        END DO
      END DO
    END IF


    ! Define dimensions

    istat=nf90_def_dim(ncID, 'x', varshape(2), xID)
    if(istat /= nf90_NoErr) print*,'quick_nc fail #2 ',nf90_strerror(istat)
    istat=nf90_def_dim(ncID, 'y', varshape(1), yID)
    if(istat /= nf90_NoErr) print*,'quick_nc fail #3 ',nf90_strerror(istat)

    IF (present(var3d) .OR. present(var3dtime)) THEN
      istat=nf90_def_dim(ncID, 'z', varshape(3), zID)
      if(istat /= nf90_NoErr) print*,'quick_nc fail #4 ',nf90_strerror(istat)
    END IF

    IF (present(var2dtime)) THEN
      istat=nf90_def_dim(ncID, 'time', varshape(3), tID)
      if(istat /= nf90_NoErr) print*,'quick_nc fail #5 ',nf90_strerror(istat)
    END IF

    IF (present(var3dtime)) THEN
      istat=nf90_def_dim(ncID, 'time', varshape(4), tID)
      if(istat /= nf90_NoErr) print*,'quick_nc fail #6 ',nf90_strerror(istat)
    END IF

    ! Define variables
    IF (present(var2d)) THEN
      istat=nf90_def_var(ncID, name, NF90_FLOAT, (/xID,yID/), varID)
      if(istat /= nf90_NoErr) print*,'quick_nc fail #7 ',nf90_strerror(istat)
    ELSEIF (present(var2dtime)) THEN
      istat=nf90_def_var(ncID, name, NF90_FLOAT, (/xID,yID,tID/), varID)
      if(istat /= nf90_NoErr) print*,'quick_nc fail #8 ',nf90_strerror(istat)
    ELSEIF (present(var3d)) THEN
      istat=nf90_def_var(ncID, name, NF90_FLOAT, (/xID,yID,zID/), varID)
      if(istat /= nf90_NoErr) print*,'quick_nc fail #9 ',nf90_strerror(istat)
    ELSEIF (present(var3dtime)) THEN
      istat=nf90_def_var(ncID, name, NF90_FLOAT, (/xID,yID,zID,tID/), varID)
      if(istat /= nf90_NoErr) print*,'quick_nc fail #10 ',nf90_strerror(istat)
    END IF

    istat=nf90_enddef(ncID)

    ! wirte variable to file
    IF (present(var2d)) THEN
      istat=nf90_put_var(ncID, varID,print2d(:,:))
      if(istat /= nf90_NoErr) print*,'quick_nc fail #11 ',nf90_strerror(istat)
    ELSEIF (present(var2dtime)) THEN
      istat=nf90_put_var(ncID, varID,print2dtime(:,:,:))
      if(istat /= nf90_NoErr) print*,'quick_nc fail #12 ',nf90_strerror(istat)
    ELSEIF (present(var3d)) THEN
      istat=nf90_put_var(ncID, varID,print3d(:,:,:))
      if(istat /= nf90_NoErr) print*,'quick_nc fail #13 ',nf90_strerror(istat)
    ELSEIF (present(var3dtime)) THEN
      istat=nf90_put_var(ncID, varID,print3dtime(:,:,:,:))
      if(istat /= nf90_NoErr) print*,'quick_nc fail #14 ',nf90_strerror(istat)
    END IF

    istat=nf90_close(ncID)
    if(istat /= nf90_NoErr) print*,'quick_nc fail #15 ',nf90_strerror(istat)

    DEALLOCATE(varshape)
    IF (present(var2d)) THEN
      DEAllOCATE(print2d)
    ELSEIF(present(var2dtime)) THEN
      DEAllOCATE(print2dtime)
    ELSEIF (present(var3d)) THEN
      DEAllOCATE(print3d)
    ELSEIF (present(var3dtime)) THEN
      DEAllOCATE(print3dtime)
    END IF

    IF (lddebug) PRINT*, 'Leave quick_nc, ierr=',istat,''//NEW_LINE('')

  END SUBROUTINE quick_nc

  SUBROUTINE quick_ascii(name,var,pmin,pmax)
  !---------------------------------------------------------------------
  ! Description:
  !   quick output of data fields with ascii symols in the terminal
  !   Mainly for debugging.
  !
  ! Usage:
  !   call quick_ascii('name',var(:,:))
  !   name definens variable name
  !   var is a 2d variable(j,i)
  !
  !   Optional:
  !     pmin min value of the plot
  !     pmax max value of the plot
  !
  !----------------------------------------------------------------

    USE mo_dust
    USE netcdf
    USE offline_org

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: name

    REAL(8),  INTENT(IN) :: &
      var    (:,:)

    REAL, OPTIONAL, INTENT(IN) :: &
      pmin, &
      pmax

    INTEGER :: &
      i,j,b,    & ! loops
      shift_x,  &
      shift_y,  &
      i_screen, & ! size of the terminal
      j_screen, & ! size of the terminal
      i_plot,   & ! size of the plot
      j_plot,   & ! size of the plot
      i_fact,   & ! grid boxes per plot cell
      j_fact,   & ! grid boxes per plot cell
      funit,    &
      pos

    REAL(8), ALLOCATABLE :: &
      printvar(:,:)

    REAL(8) ::     &
      plotmax,     &
      plotmin,     &
      plotbins(12)

    CHARACTER(1) :: &
      symbols(12),  &
      add_sym

    CHARACTER(9) :: &
      str_bin

    CHARACTER(1000), ALLOCATABLE :: &
      plotstrings(:)

    IF (lddebug) PRINT*, 'Enter quick_ascii'

    ! get screen size
    CALL EXECUTE_COMMAND_LINE('tput cols > screen.tmp')
    CALL EXECUTE_COMMAND_LINE('tput lines >> screen.tmp')
    OPEN(newunit = funit, file = 'screen.tmp')
    READ(funit,*) i_screen
    READ(funit,*) j_screen
    CLOSE(funit)
    CALL EXECUTE_COMMAND_LINE('rm screen.tmp')

    i_screen = i_screen/2 - 2 - 10
    j_screen = j_screen - 3



    i_fact = max(1,ceiling(float(ie_tot)/float(i_screen)))
    j_fact = max(1,ceiling(float(je_tot)/float(j_screen)))


    IF (i_fact > j_fact) j_fact = i_fact
    IF (i_fact < j_fact) i_fact = j_fact


    i_plot = ie_tot/i_fact
    j_plot = je_tot/j_fact



    AllOCATE(printvar(j_plot,i_plot))
    printvar = 0.

    DO i=1, i_plot
      DO j=1, j_plot

        DO shift_x=1,i_fact
          DO shift_y=1, j_fact
            printvar(j,i) = printvar(j,i) + var((j-1)*j_fact+shift_y,(i-1)*i_fact+shift_x)
          END DO
        END DO

        printvar(j,i) = printvar(j,i)/(i_fact*j_fact)

      END DO
    END DO

    IF (present(pmax)) THEN
      plotmax=pmax
    ELSE
      plotmax=maxval(printvar)
    END IF

    IF (present(pmin)) THEN
      plotmin=pmin
    ELSE
      plotmin=minval(printvar)
    END IF


    symbols( 1)='.'
    symbols( 2)=','
    symbols( 3)='-'
    symbols( 4)='~'
    symbols( 5)=':'
    symbols( 6)=';'
    symbols( 7)='!'
    symbols( 8)='='
    symbols( 9)='*'
    symbols(10)='#'
    symbols(11)='$'
    symbols(12)='@'


    DO i=1, 12
      plotbins(i) = plotmin + i*(plotmax-plotmin)/13
    END DO

    ALLOCATE(plotstrings(j_plot+3))
    plotstrings=''
    plotstrings(3:j_plot+2) = '|'





    DO i=1,i_plot*2+1
      plotstrings(2) = TRIM(plotstrings(2))//'-'
      plotstrings(j_plot+3) = TRIM(plotstrings(j_plot+3))//'-'
    END DO

    DO j=1,j_plot
      DO i=1,i_plot
        DO b=1,12
          IF (printvar(j,i) < plotbins(1)) THEN
            add_sym = 'x'
          ELSEIF (printvar(j,i) > plotbins(b)) THEN
            add_sym = symbols(b)
          END IF

        END DO

        plotstrings(j+2)=TRIM(plotstrings(j+2))//' '//add_sym
      END DO
      plotstrings(j+2)=TRIM(plotstrings(j+2))//' |'

      pos = scan(plotstrings(j+2),'x')
      DO WHILE (pos > 0)
        plotstrings(j+2)(pos:pos)=' '
        pos = scan(plotstrings(j+2),'x')
      END DO

    END DO


    ! legend
    DO i=1,12
      WRITE(str_bin, '(ES8.2)') plotbins(i)
      pos = 15
      plotstrings(pos-i)=TRIM(plotstrings(pos-i))//' '//symbols(i) // ' > '//str_bin
    END DO

    ! min, mean, med, max
    WRITE(str_bin, '(ES8.2)') minval(printvar)
    plotstrings(1) = 'Minimum: '//str_bin
    WRITE(str_bin, '(ES8.2)') maxval(printvar)
    plotstrings(1) = TRIM(plotstrings(1))//'    Maximum: '//str_bin


    plotstrings(j_plot+3)(2:2+len(name))=name
    plotstrings(j_plot+3)=' '//TRIM(plotstrings(j_plot+3))
    plotstrings(2)=' '//TRIM(plotstrings(2))
    plotstrings(1)=' '//TRIM(plotstrings(1))


    ! IF (present(doflush)) THEN
    !   IF (doflush) CALL FLUSH()
    ! END IF
    DO j=j_plot+3,1,-1
      print*,TRIM(plotstrings(j))
    END DO


    IF (lddebug) PRINT*, 'Leave quick_ascii',''//NEW_LINE('')




  END SUBROUTINE


END MODULE src_dust
