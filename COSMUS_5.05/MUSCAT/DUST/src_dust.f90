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
  PUBLIC :: organize_dust




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
      SoilFile
    CHARACTER(120) :: &
      filename

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

    ierr = 0

    ! ------------------------------------
    ! +-+-+- Section 1 Init + Input -+-+-+
    ! ------------------------------------
    IF (yaction == "init") THEN

      ! +-+-+- Sec 1.1 Check -+-+-+
      ! Check consistency of namelist settings


      ! dust_scheme need right values
      IF (dust_scheme < 0 .OR. dust_scheme > 1) THEN
        ierr = 100001
        yerr = 'wrong value for dust_scheme'
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
      IF (TRIM(psrcFile) == 'without') THEN
        ierr = 100004
        yerr = 'psrcFile is missing'
        PRINT*,'ERROR    src_dust "init" '
        PRINT*,'         #',ierr
        PRINT*,'         ',yerr
        STOP
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
        IF (TRIM(vegmonFile) /= 'without' .AND. TRIM(vegdayFile) /= 'without') vegmonFile = 'without'
        ! set flag for daily vegetation
        lvegdaily = .FALSE.
        IF (TRIM(vegdayFile) /= 'without') lvegdaily = .TRUE.
      ELSE
        vegmonFile = 'without'
        vegdayFile = 'without'
      ENDIF

      ! set flag for minimum vegetation file
      lvegmin = .FALSE.
      IF (TRIM(vegminFile) /= 'without') lvegmin = .TRUE.

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
        END IF
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
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,1:SoilFrac))
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
                decomp(ib1)%ix0+1:decomp(ib1)%ix1,dimveg))
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
                decomp(ib1)%ix0+1:decomp(ib1)%ix1,dimveg))
        ALLOCATE(dust(ib1)%d_emis(decomp(ib1)%iy0+1:decomp(ib1)%iy1,   &
                 decomp(ib1)%ix0+1:decomp(ib1)%ix1,1:nt))

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
        dust(ib1)%mfac(:,:,:)=1.
        dust(ib1)%d_emis(:,:,:)=0.

      END DO



      ! +-+-+- Sec 1.3 Input -+-+-+

      ! - Loop over all possible input files
      DO filenum = 1,8

        ! dim = 2 mostly
        dim = 2

        ! specify filename,
        ! soiltype (2d)
        IF (filenum == 1) filename = soiltypeFile

        ! source (2d)
        IF (filenum == 2) filename = psrcFile

        ! z0 (2d)
        IF (filenum == 3) filename = z0File

        ! biome (2d)
        IF (filenum == 4) filename = biomeFile

        ! vegmon (3d)
        IF (filenum == 5) filename = vegmonFile
        IF (filenum == 5) dim = 3

        ! vegday (3d)
        IF (filenum == 6) filename = vegdayFile
        IF (filenum == 6) dim = 3

        ! vegmin (2d)
        IF (filenum == 7) filename = vegminFile

        ! vegday (3d)
        IF (filenum == 8) filename = moistFile
        IF (filenum == 8) dim = 3

        ! if (filename == without) nothing happen
        IF (TRIM(filename) /= 'without') THEN

          ! allocate var to store the input
          IF (dim == 2) AllOCATE (read_input(igy0+1:igy1,igx0+1:igx1,1))
          IF (dim == 3) AllOCATE (read_input(igy0+1:igy1,igx0+1:igx1,dimveg))

          ! only prozess #0 open files
          IF (my_cart_id == 0) THEN
            ! decide ascii or netcdf
            IF (filename(LEN(TRIM(filename))-2:) == '.nc') THEN
              IF (filenum == 1) STOP 'nc input not supported yet for soiltype yet' ! CALL read_nc(TRIM(filename),'varname',read_input,1     ,.FALSE.,ierr,yerr)
              IF (filenum == 2) STOP 'nc input not supported yet for psrc yet'     ! CALL read_nc(TRIM(filename),'varname',read_input,1     ,.FALSE.,ierr,yerr)
              IF (filenum == 3) STOP 'nc input not supported yet for z0 yet'       ! CALL read_nc(TRIM(filename),'varname',read_input,1     ,.FALSE.,ierr,yerr)
              IF (filenum == 4) STOP 'nc input not supported yet for biome'    ! CALL read_nc(TRIM(filename),'varname',read_input,1     ,.FALSE.,ierr,yerr)
              IF (filenum == 5) CALL read_nc(TRIM(filename),'FCOVER' ,read_input,dimveg,.FALSE.,ierr,yerr)
              IF (filenum == 6) CALL read_nc(TRIM(filename),'FCOVER' ,read_input,dimveg,.TRUE. ,ierr,yerr)
              IF (filenum == 7) CALL read_nc(TRIM(filename),'FCOVER' ,read_input,1     ,.FALSE.,ierr,yerr)
              IF (filenum == 8) CALL read_nc(TRIM(filename),'swvl1' ,read_input,dimveg     ,.FALSE.,ierr,yerr)

              IF (ierr /= 0) THEN
                print*, ierr
                ierr = 100009
                PRINT*,'ERROR    src_dust "init" '
                PRINT*,'         #',ierr
                PRINT*,'         ERROR reading',TRIM(filename)
                PRINT*,'         ',yerr
                STOP
              END IF
            ELSE
              CALL read_ascii(TRIM(filename),read_input)
            END IF
          END IF ! (my_cart_id == 0)

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
            IF (filenum == 1) copy2d => dust(ib1)%soiltype
            IF (filenum == 2) copy2d => dust(ib1)%source
            IF (filenum == 3) copy2d => dust(ib1)%z0
            IF (filenum == 4) copy2d => dust(ib1)%biome
            IF (filenum == 5) copy3d => dust(ib1)%veg
            IF (filenum == 6) copy3d => dust(ib1)%veg
            IF (filenum == 7) copy2d => dust(ib1)%vegmin2
            IF (filenum == 8) copy3d => dust(ib1)%vmoist

            IF (dim == 2) CALL copy2block(decomp(ib1),dim,read_input,too2d=copy2d)
            IF (dim == 3) CALL copy2block(decomp(ib1),dim,read_input,too3d=copy3d)


          END DO

          DEALLOCATE(read_input)
        END IF ! (TRIM(filename) /= 'without')

      END DO ! (filenum = 1,7)


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


        ! +-+-+- Sec 1.4.1 dust flux -+-+-+

        IF (dust_scheme == 1) THEN
          ! init of the dust emission sheme by Tegen et al. 2002
          ! CALL init_tegen(decomp(ib1))!ierr,yerr)
          CALL init_tegen(decomp(ib1),ndays=dimveg)!ierr,yerr)
        END IF


        ! +-+-+- Sec 1.4.2 vegetation init -+-+-+


        IF (veg_scheme == 1) THEN
          CALL okin_vegetation(decomp(ib1),dimveg)
        ELSE IF (veg_scheme == 2) THEN
          CALL linear_vegetation(decomp(ib1),dimveg)
        END IF


        ! +-+-+- Sec 1.4.2 surface roughness -+-+-+
        IF (lwithz0) THEN
          ! Drag partition
          CALL roughness(decomp(ib1),dimveg)
        END IF


        ! +-+-+- Sec 1.4.3 moisture -+-+-+
        IF (moist_scheme == 1) THEN
          print*, 'call fecan'
          CALL fecan(decomp(ib1),dimveg)
        END IF

      END DO


    ! +-+-+- Sec 1.5 clean up -+-+-+

    ! DEALLOCATE(dust_it)

    ! ------------------------------------
    ! +-+-+- Section 2 Dust flux calculation -+-+-+
    ! ------------------------------------
    ELSEIF (yaction == "calc") THEN

      CALL emission_tegen(subdomain,flux)
      ! STOP 'TESTING'

    END IF ! (yaction == "***")


    ! STOP 'TESTING'
  END SUBROUTINE organize_dust

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
    ! particle size loop, calculation of Uth for every particle size
    DO WHILE(dp.lE.Dmax + 1E-05)
      nn = nn + 1
      gransize(nn) = dp
      dmy_B = a_rnolds * (dp ** x_rnolds) + b_rnolds
      dmy_K = SQRT(rop * g * dp / roa) * SQRT(1. + 0.006 /(rop * g * dp ** 2.5))
      IF (dmy_B < 10) THEN
        Uth(nn) = 0.129 * dmy_K / SQRT(1.928 * (dmy_B ** 0.092) - 1.)
      ELSE
        Uth(nn) = 0.129 * dmy_K * ( 1. -0.0858 * EXP(-0.0617 * (dmy_B - 10.)) )
      END IF
      dp = dp * EXP(Dstep)
    END DO


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
        IF (psrcType == 0) THEN  ! IT02
          ! Preferential Sources = Potential lakes
          !   If a grid box is identified as dust source then the soiltype switches to
          !   best emission properties (soiltype=10), to avoid an underestimation.
          IF(source(j,i) > 0.5) THEN
            soiltype(j,i) = 10
            IF (z0(j,i) <= 0.) z0(j,i) = 0.001 ! [cm]
          ENDIF

        ELSEIF (psrcType == 1) THEN ! MSG
          ! Preferential Sources = Potential lakes
          IF(soiltype(j,i) < 13) soiltype(j,i) = 9
          IF(source(j,i) >= 2) THEN
            IF (soiltype(j,i) == 9) soiltype(j,i) = 10
            IF (z0(j,i) <= 0.) z0(j,i) = 0.001 ! [cm]
          ENDIF

        ELSEIF (psrcType == 2) THEN ! AC Dust
          ! For the agricultural dust, the source map provide the fraction of cropland
          ! in the grid box. Every grid box with a cropland fraction > 0 is allowed to
          ! emit dust. The total dust emission of the grid box is scaled late with the
          ! cropland fraction.
          IF(source(j,i) > 0.) THEN
            IF (z0(j,i) <= 0.) z0(j,i) = 0.001 ! [cm]
          ELSE
            soiltype(j,i) = 0.
          ENDIF
        ENDIF



        isoiltype = int(soiltype(j,i)) ! index of soiltype

        ! change everthing that not fits the soiltype table to 9 = ICE
        IF(isoiltype < 1 .or. isoiltype > nats) isoiltype=9

        ! init alpha from lookup table
        alpha(j,i) = solspe(isoiltype,nmode*3+1)

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



    ! !--------TEST NC OUTPUT     only activate when needed
    !        AllOCATE(printvar(subdomain%ntx,subdomain%nty,1,ndays))
    !
    !        do i=1,subdomain%ntx
    !          do j=1,subdomain%nty
    !
    !              ! do t=1,ndays
    !              !   ! ! if (lai(j,i,1,t) /= 0. .and. lai(j,i,1,t) < 1) then
    !              !   ! if (lai(j,i,1,t) /= maxval(lai(:,:,1,:))) then
    !              !   !   if (sp(j,i,1,2) == 0.) then
    !              !   !     printvar(i,j,1,t)=99.
    !              !   !   else
    !              !   !     printvar(i,j,1,t)=lai(j,i,1,t)!lai_eff(j,i,1,t)
    !              !   !   end if
    !              !   ! ELSEIF(lai(j,i,1,t) == 0.) then
    !              !   !   printvar(i,j,1,t)=-99.
    !              !   ! end if
    !              ! printvar(i,j,1,t)=feff(j,i,t)
    !              ! ! printvar(i,j,1,t)=lai_eff(j,i,1,t)
    !              ! ! printvar(i,j,1,t)=lai(j,i,1,t)!vegmin(j,i,1)
    !              ! end do
    !
    !
    !
    !              ! printvar(i,j,1,1)=vegmin(j,i,1)
    !
    !              ! soiltype
    !              printvar(i,j,1,1)=soiltype(j,i)
    !
    !              if (i==20 .and. j==20) printvar(i,j,1,1)=50
    !
    !
    !              ! !biom
    !              ! printvar(i,j,1,1)=sp(j,i,1,5)
    !
    !              !cosmo lai
    !              ! printvar(i,j,1,1)=newlai(i,j)
    !
    !               !cosmo z0
    !               ! printvar(i,j,1,1)=sp (j,i,1,4)
    !
    !               !cult
    !               ! printvar(i,j,1,1)=sp (j,i,1,3)
    !
    !          end do
    !        end do
    !
    !        ! print*, 'pvar 20/20',printvar(20,20,1,1)
    !        ! print*, 'pvar 19/21',printvar(19,21,1,1)
    !       if (subdomain%ib == 1) then
    !
    !          ! call quick_nc(0,'efflaipw4.nc','EFFLAI',printvar(:,:,:,:),ie_tot,je_tot,1,ndays,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !          ! call quick_nc(0,'z0.nc','Z0',printvar(:,:,:,:),ie_tot,je_tot,1,1,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !          ! call quick_nc(0,'cult.nc','CULT',printvar(:,:,:,:),ie_tot,je_tot,1,1,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !         !  call quick_nc(0,'biom.nc','biom',printvar(:,:,:,:),ie_tot,je_tot,1,1,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !         ! call quick_nc(0,'efflai.nc','efflai',printvar(:,:,:,:),ie_tot,je_tot,1,ndays,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !         ! call quick_nc(0,'vegmin.nc','vegmin',printvar(:,:,:,:),ie_tot,je_tot,1,1,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !                 ! call quick_nc(0,'lai.nc','lai',printvar(:,:,:,:),ie_tot,je_tot,1,ndays,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !
    !          call quick_nc(0,'soiltype.nc','soiltype',printvar(:,:,:,:),ie_tot,je_tot,1,1,subdomain%igx0, &
    !          subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !       ELSE
    !         CALL SLEEP(1)
    !       end if
    !
    !       do i=0, num_compute
    !         if (i==subdomain%ib) then
    !            ! call quick_nc(1,'efflaipw4.nc','EFFLAI',printvar(:,:,:,:),ie_tot,je_tot,1,ndays,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb
    !             ! call quick_nc(1,'z0.nc','Z0',printvar(:,:,:,:),ie_tot,je_tot,1,1,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !            ! call quick_nc(1,'cult.nc','CULT',printvar(:,:,:,:),ie_tot,je_tot,1,1,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !           !  call quick_nc(1,'biom.nc','biom',printvar(:,:,:,:),ie_tot,je_tot,1,1,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !           ! call quick_nc(1,'efflai.nc','efflai',printvar(:,:,:,:),ie_tot,je_tot,1,ndays,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !           ! call quick_nc(1,'vegmin.nc','vegmin',printvar(:,:,:,:),ie_tot,je_tot,1,1,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !           ! call quick_nc(1,'lai.nc','lai',printvar(:,:,:,:),ie_tot,je_tot,1,ndays,subdomain%igx0,subdomain%igx1,subdomain%igy0,subdomain%igy1,subdomain%ib,nb)
    !           call quick_nc(1,'soiltype.nc','soiltype',printvar(:,:,:,:),ie_tot,je_tot,1,1,subdomain%igx0+1,subdomain%igx1+1, &
    !           subdomain%igy0+1,subdomain%igy1+1,subdomain%ib,nb)
    !         ELSE
    !           CALL SLEEP(1)
    !         end if
    !       end do

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
    USE mo_gas, ONLY: mol2part, nradm
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
      ustar
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
    REAL(8), POINTER :: mfac(:,:,:)
    REAL(8), POINTER :: z0(:,:)
    ! REAL(8), POINTER :: umin2(:,:)
    REAL(8), PARAMETER :: umin2=umin
    REAL(8), POINTER :: w_str(:,:)
    REAL(8), POINTER :: DustEmis(:,:,:), EmiRate(:,:,:,:)




    IF (nDust == 0) RETURN
    IF (DustMod <= 0) RETURN


#ifndef OFFLINE
    dxK     => geo  (subdomain%ib)%dxK(:,:)
    dyK     => geo  (subdomain%ib)%dyK(:,:)
    dz      => geo  (subdomain%ib)%dz(:,:,:)
    rhosur  => meteo(subdomain%ib)%rho(1,:,:,ScalCur)
    qrsur   => meteo(subdomain%ib)%QRSur(:,:,ScalCur)
    usur    => meteo(subdomain%ib)%u(1,:,:,WindCur)
    vsur    => meteo(subdomain%ib)%v(1,:,:,WindCur)
#endif

    soiltype => dust(subdomain%ib)%soiltype(:,:)
    source   => dust(subdomain%ib)%source(:,:)
    alpha    => dust(subdomain%ib)%alpha2(:,:)
    feff     => dust(subdomain%ib)%feff(:,:,:)
    veff     => dust(subdomain%ib)%veff(:,:,:)
    mfac     => dust(subdomain%ib)%mfac(:,:,:)
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
        ! print*, 'dxK',dyK(j,i),'dyK',dyK(j,i),'dz',dz(1,j,i)

        ! +-+-+- Sec 2 update of the meteorological variables -+-+-+

#ifndef OFFLINE
        !---  flux initialisations
        uwind = usur(j,i+1)/dyK(j,i+1)+usur(j,i)/dyK(j,i)
        uwind = 0.5E0 * uwind / dz(1,j,i)
        vwind = vsur(j+1,i)/dxK(j+1,i)+vsur(j,i)/dxK(j,i)
        vwind = 0.5E0 * vwind / dz(1,j,i)

        van = SQRT(uwind**2+vwind**2)/rhosur(j,i)
#else
        van = SQRT(u(j,i,ntstep)**2+v(j,i,ntstep)**2)
#endif


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

        IF(feff(j,i,tnow) <= 0.) THEN
         ustar = 0.
        ELSE
         ustar = (VK * van *100.)/(log(0.5E0 * 100 * dz(1,j,i)/z0(j,i))) !!cm/s
         ! print*, i,j,ustar,dz(1,j,i)/z0(j,i),ustar_fv(i,j)
        END IF  !! IF(feff(j,i,tnow).LE.0.)


        ! +-+-+- Sec 3 Flux calculation -+-+-+
        ! print*,tnow,shape(feff)
        IF (feff(j,i,tnow) > 0.) THEN
          IF (Ustar > 0 .AND. Ustar > umin2/feff(j,i,tnow) ) THEN
            kk = 0
            dp = Dmin
            DO WHILE (dp <= Dmax+1E-5)
              kk = kk+1

              ! original Tegen Code
              ! ! Is this reduction necessary (MF)?
              Uthp=uth(kk)*umin2/umin*u1fac !reduce threshold for cultivated soils
              ! drag coeff
              ! Uthp=Uthp/feff(j,i,tnow)
              ! ! moist
              ! Uthp=Uthp*mfac(j,i,tnow)
              ! Marticorena:
              fdp1 = (1.-(Uthp/(feff(j,i,tnow) * Ustar)))
              fdp2 = (1.+(Uthp/(feff(j,i,tnow) * Ustar)))**2.
              ! fdp1 = (1.-(Uthp * Ustar))
              ! fdp2 = (1.+(Uthp * Ustar))**2.

              ! ! Shao:
              ! fdp1 = (1.-(Uthp/(feff(j,i,tnow) * Ustar))**2)
              ! fdp2 = 1.

              IF (fdp1 <= 0 .OR. fdp2 <= 0) THEN
                flux_umean = 0.
              ELSE
                flux_umean = srel(i_s1,kk) * fdp1 * fdp2 * cd * Ustar**3 *alpha(j,i)
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

                        ! if (fluxtyp(kkk) /= fluxtyp(kkk) .or. fluxtyp(kkk) > 10.) THEN
                        !   print*, 'su_srelV',i_s11,kk,kkmin,srelV(i_s11,kkk),su_srelV(i_s11,kk),-su_srelV(i_s11,kkmin)
                        !   STOP
                        ! end if
                        ! if (fluxtyp(kkk) == 0. .and. i_s1 >0 .and. i_s1 /= 8 .and. i_s1 /= 9)print*,'ftyp', i_s1,i_s11,srelV(i_s11,kkk),su_srelV(i_s11,kk),su_srelV(i_s11,kkmin)

                      END IF ! (kk > kkmin)
                    END IF ! (dpd >= dbstart)
                    dpd=dpd*exp(dstep)
                  END DO !dpd
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
              ! if (fluxtot(nn) == 0. .and. i_s1 >0 .and. i_s1 /= 8 .and. i_s1 /= 9)print*,'ftot', i_s1,fluxtot(nn),fluxbin(nn)
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
        IF (nradm == 1)  THEN                      ! chemistry units: g/m2/s ==> g/m2/s * mol2part
          FDust(j,i,:) = FDust(j,i,:) * mol2part
        END IF

        ! IF (FDust(j,i,1) /= FDust(j,i,1) ) print*,'Fdust', i,j,FDust(j,i,1)

        !---  add fluxes to right hand side
        DO nn=1,DustBins
          Flux(1,j,i,DustInd(nn))   = Flux(1,j,i,DustInd(nn)) + FDust(j,i,nn)/dz(1,j,i)
          DustEmis(j,i,DustInd(nn)) = FDust(j,i,nn)

#ifndef OFFLINE
          !---  summarize biogenic and total emission rates
          EmiRate(EmiIndBio,j,i,DustInd(nn)) = EmiRate(EmiIndBio,j,i,DustInd(nn)) + FDust(j,i,nn)
          EmiRate(EmiIndSum,j,i,DustInd(nn)) = EmiRate(EmiIndSum,j,i,DustInd(nn)) + FDust(j,i,nn)
#endif
        END DO
      END DO
    END DO
  END SUBROUTINE emission_tegen

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
        feff = FF*CCC


        IF (local_feff < 0.) local_feff = 0.
        IF (local_feff > 1.) local_feff = 1.
      END IF
      feff(j,i,:) = local_feff
      END DO
    END DO
    ! end lon-lat-loop

  END SUBROUTINE roughness


  !+ roughness
  !---------------------------------------------------------------------
  SUBROUTINE fecan(subdomain,dimsimu)
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
    USE mo_dust

#ifdef OFFLINE
    USE offline_org
#endif

    IMPLICIT NONE

    TYPE(rectangle), INTENT(IN) :: subdomain

    INTEGER, INTENT(IN) :: dimsimu

    INTEGER :: &
      dnow,    &  ! time loop
      i,j,t,     &  ! loops
      isoiltype


    REAL(8) ::  &
      w_str,        &  ! feff inside the loop
      z0s,               &  ! small scale roughness length
      AAA,               &
      BB,                &
      CCC,               &  ! dummys
      DDD,               &
      EE,                &
      FF


    REAL(8), POINTER ::  &
      vmoist(:,:,:),           &
      mfac(:,:,:),       &
        soiltype(:,:)!,     &

    vmoist   => dust(subdomain%ib)%vmoist(:,:,:)
    mfac => dust(subdomain%ib)%mfac(:,:,:)
    soiltype => dust(subdomain%ib)%soiltype(:,:)



    ! start lon-lat-loop
    DO i=1,subdomain%ntx
      DO j=1,subdomain%nty
        isoiltype = int(soiltype(j,i))
        w_str = 0.0014*(solspe(isoiltype,nmode*3)*100)**2 + 0.17*(solspe(isoiltype,nmode*3)*100)
        DO dnow=1,dimsimu
            !---------------------------------------------------------------------------------------
            !       Calculation of the threshold soil moisture (w')  [Fecan, F. et al., 1999]
            !---------------------------------------------------------------------------------------

        !       !          W0   = 0.99! solspe(isoiltype,nmode*3+2)
        !    feff = 0.
        !
        ! ! calculation of mfac: ratio between threshold friction velocities wet/dry
        ! ! [Fecan, F. et al., 1999]

          IF ((vmoist(j,i,dnow)*100 - w_str) > 0.0) THEN
            mfac(j,i,dnow) = (1 + 1.21 * ( vmoist(j,i,dnow)*100 - w_str)**0.68 )**0.5
          ELSE
            mfac(j,i,dnow) = 1000.
          END IF
            ! mfac(j,i,dnow)=mfac(j,i,dnow)/2
            ! mfac=1.
          ! print*, 'mf',mfac(j,i,dnow),vmoist(j,i,dnow)*100,w_str
        END DO
      END DO
    END DO

    ! end lon-lat-loop

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

  END SUBROUTINE read_ascii





  SUBROUTINE read_nc(infile,varname,outvar,ndays,timecheck,ierror,yerrmsg)

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

    ! Parameter
    CHARACTER(*), INTENT(IN)    :: &
      infile        ! filename
    CHARACTER(*),   INTENT(IN)    :: &
      varname       ! name of variable that shoud be read

    INTEGER,        INTENT(IN)    :: &
      ndays         ! number of possible days in the model run
    LOGICAL,        INTENT(IN)    :: &
      timecheck     ! check if time fits model run

    REAL(8),        INTENT(INOUT) :: &
      outvar(:,:,:)        ! Output var

    INTEGER,        INTENT(OUT)   :: &
      ierror

    CHARACTER(*),   INTENT(OUT)   :: &
      yerrmsg       ! error message


    ! local Vars
    INTEGER  :: &
      i,j,t,    &  ! loop
      istart,   &  ! index of the start date in the var file
      idate,    &  ! ydate read into integer
      istat,    &  ! local error code
      ncID,     &  ! id of nc files
      timeID,   &  ! id of time var
      varID,    &  ! id of the  var
      dimID,    &  ! id of a dimension
      dimlen       ! length of the above dimension

    REAL  :: &
      var_scale, &    ! index of the start date in the var file
      var_offset    ! index of the start date in the var file

    INTEGER, ALLOCATABLE :: &
      times(:)

    REAL, ALLOCATABLE :: &
      var_read(:,:,:)



    ! start subroutine
    ! ---------------------------------------------------------

    ! Print short status
    PRINT*,'read_nc ', infile

    ! open the file
    istat = nf90_open(infile, nf90_nowrite, ncid)
    IF (istat /= nf90_noerr) THEN
      ierror  = 10001
      yerrmsg = TRIM(nf90_strerror(istat))
      RETURN
    ENDIF

    ! get id of lat dimension
    istat = nf90_inq_dimid(ncID, 'y', dimID)
    IF (istat /= nf90_noerr) THEN
      ierror  = 10002
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
      yerrmsg = 'Error reading '//TRIM(varname)//' file: wrong lat dimension'
      RETURN
    END IF

     ! get id of lon dimension
     istat = nf90_inq_dimid(ncID, 'x', dimID)
     IF (istat /= nf90_noerr) THEN
       ierror  = 10005
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
       yerrmsg = 'Error reading '//TRIM(varname)//' file: wrong lon dimension'
       RETURN
     END IF

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
        IF (istart + ndays - 1 > size(times)) THEN
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

    ! Allocate var to read
    istat = 0
    ALLOCATE (var_read(ie_tot,je_tot,ndays), STAT=istat)
    IF (istat /= 0) THEN
      ierror = 10014
      yerrmsg = 'allocation of var_read failed'
      RETURN
    ENDIF


    ! get the id of the var
    istat = nf90_inq_varid(ncID, varname , varID)
    IF (istat /= nf90_noerr) THEN
      ierror  = 10015
      yerrmsg = TRIM(nf90_strerror(istat))
      RETURN
    ENDIF

    ! get value of the scaling factor
    istat = nf90_get_att(ncID, varID, 'scale_factor',var_scale)
    IF (istat /= nf90_noerr) THEN
      ! ierror  = 10016
      var_scale = 1
      ! yerrmsg = TRIM(nf90_strerror(istat))
      !RETURN
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
    istat = nf90_get_var(ncid, varID, var_read,start= (/1,1,istart/),count=(/ie_tot,je_tot,ndays/))
    IF (istat /= nf90_noerr) THEN
      ierror  = 10017
      yerrmsg = TRIM(nf90_strerror(istat))
      RETURN
    ENDIF


    ! AllOCATE(outvar(je_tot,ie_tot,ndays))

    ! wirte the var into outvar
    DO i=1,ie_tot
      DO j=1,je_tot
        DO t=1, ndays
          ! digital values to physical values with the scaling factor
          outvar(j,i,t)=var_read(i,j,t)*var_scale+var_offset
          ! physical max of lai is 7 everthing abov is a missing value, in MUSCAT set to 0.
          ! IF (SG_lai(j,i,1,t) > 7.) SG_lai(j,i,1,t) = 0.
        END DO
      END DO
    END DO

    ! print*,'vegfile inside',outvar(30,30,1)

    DEALLOCATE (times)
    DEALLOCATE (var_read)
    ! DEALLOCATE (outvar)
    RETURN
  END SUBROUTINE read_nc



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

  !----------------------------------------------------------------
  END SUBROUTINE copy2block


  SUBROUTINE quick_nc(prep,name,vname,var,xe,ye,ze,te,xmin,xmax,ymin,ymax,isub,nblocs)

    USE mo_dust
    USE netcdf

#ifdef OFFLINE
    USE offline_org
#endif
    ! USE data_parallel, ONLY: my_cart_id

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: name
    CHARACTER(*), INTENT(IN) :: vname

    REAL(8),INTENT(IN) :: var(:,:,:,:)

    INTEGER, INTENT(IN) :: prep, xe, ye, ze, te, xmin,xmax,ymin,ymax,isub,nblocs

    LOGICAL :: firstbloc

    INTEGER :: i,istat,ncID,xID,yID,zID,tID,varID

    xID=0
    yID=0
    zID=0
    tID=0


    if(prep == 0) then
      print*,my_cart_id,isub,"create nc"
      istat=nf90_create(name, NF90_SHARE, ncID)
      if (xe >1) istat=nf90_def_dim(ncID, 'x', xe, xID)
      if (ye >1) istat=nf90_def_dim(ncID, 'y', ye, yID)
      if (ze >1) istat=nf90_def_dim(ncID, 'z', ze, zID)
      if (te >1) istat=nf90_def_dim(ncID, 't', te, tID)

      IF (xID /=0 .and. yID/=0 .and. zID == 0 .and. tID==0) &
        istat=nf90_def_var(ncID, vname, NF90_FLOAT, (/xID,yID/), varID)
      IF (xID /=0 .and. yID/=0 .and. zID /= 0 .and. tID==0) &
        istat=nf90_def_var(ncID, vname, NF90_FLOAT, (/xID,yID,zID/), varID)
      IF (xID /=0 .and. yID/=0 .and. zID /= 0 .and. tID/=0) &
        istat=nf90_def_var(ncID, vname, NF90_FLOAT, (/xID,yID,zID,tID/), varID)
      IF (xID /=0 .and. yID/=0 .and. zID == 0 .and. tID/=0) &
        istat=nf90_def_var(ncID, vname, NF90_FLOAT, (/xID,yID,tID/), varID)

      istat=nf90_enddef(ncID)

      istat=nf90_close(ncID)

    else
      do i=1,nblocs
        if (i==isub) then
          print*,my_cart_id,isub,'wirte data',xmin,ymin,xmax,ymax,xe,ye

           istat=nf90_open(name, NF90_WRITE, ncID)
           if(istat /= nf90_NoErr) print*,'fail 1',nf90_strerror(istat)
           istat=nf90_inq_varid(ncID, vname, varID)
           if(istat /= nf90_NoErr) print*,'fail 2',nf90_strerror(istat)
           IF (xe /=1 .and. ye/=1 .and. ze == 1 .and. te==1) &
             istat=nf90_put_var(ncID, varID,var(:,:,1,1),    &
                           start = (/ xmin, ymin/), &
                           count = (/ xmax-xmin, ymax-ymin/))
             if(istat /= nf90_NoErr) print*,'fail 3',nf90_strerror(istat)
           IF (xe /=1 .and. ye/=1 .and. ze == 1 .and. te/=1) &
             istat=nf90_put_var(ncID, varID,var(:,:,1,:),    &
                           start = (/ xmin, ymin,1/), &
                           count = (/ xmax-xmin, ymax-ymin,te/))
             if(istat /= nf90_NoErr) print*,'fail 4',nf90_strerror(istat)
          istat=nf90_close(ncID)
          if(istat /= nf90_NoErr) print*,'fail 5',nf90_strerror(istat)
        end if
      end do
    end if

    ! print*,var

  END SUBROUTINE quick_nc


END MODULE src_dust
