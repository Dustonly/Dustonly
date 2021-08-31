!================================================================
! $Id: mo_dust.f90,v 1.5 2007/01/09 17:09:39 wolke Exp $
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
!======================================================================================
MODULE mo_dust
#ifndef OFFLINE
  USE src_aerosol,   ONLY: DustMod ! Flag for Setting Soil Data
#endif
  !======================================================================================
  ! Description:
  ! ------------
  !  This module contains the namelist values for the dust emisson scheme.
  !======================================================================================
  !
  !---------------------------------------------------------------------
  !-- control parameters
  LOGICAL ::    &
    old_dust,      & ! if true use the original dust routines
                     ! if false use new module src_dust
    lwithz0,       & ! =false without z0, =true with z0
    lwithbiom        ! =false without biomes, =true with biomes

  INTEGER ::      &
    nDust,        & ! Flag for Dust Calculations
    dust_scheme,  & ! 1=Tegen02
    veg_scheme,   & ! =0 no vegitation; =1 Okin scheme; =2 linear Tegen
    moist_scheme, &
    fricvelo_scheme, & ! 1 calc from wind and z0;  2 pre-calculated
    psrcType,     & ! Flag for type of potential dust source
                    ! 0 : off, 1 : psrc, 2 : msgsrc
    soilmaptype,      & ! 1 : solspe table, 2 : soilgrids
    threshold_scheme   ! 0 : Marticorena, 1 : Shao  !

  REAL :: &
    moistinc, &    ! time increment of soil moisture
    SSR_min


  !-- Files with soil data
  CHARACTER(120) ::     &
    soiltypeFile,       & ! Filename of Soil Type Data
    psrcFile,           & ! Filename of preferential Dust Sources
    cultFile,           & ! Filename of Cultivation Class
    vegmonFile,         & ! Filename of monthly vegitation cover
    vegdayFile,         & ! Filename of daily vegetation cover
    vegminFile,         & ! Filename of min vegetation cover MF
    z0File,             & ! Filename of Roughness Length
    biomeFile,          & ! Filename of Vegetation Cover/Type Data
    moistFile


  ! description of dust particles for external use
  INTEGER, PARAMETER :: DustBins = 5  !8  ! number of dust particle fractions
  INTEGER :: DustInd(DustBins)            ! indices of dust particles

  REAL :: dustbin_top(DustBins)
  DATA dustbin_top(1) /1.E-6/,  &
       dustbin_top(2) /3.E-6/,  &
       dustbin_top(3) /9.E-6/,  &
       dustbin_top(4) /26.E-6/, &
       dustbin_top(5) /80.E-6/

  CHARACTER(20) :: DustName(DustBins)
  DATA  DustName(1) /'DP_01'/,     &
        DustName(2) /'DP_03'/,     &
        DustName(3) /'DP_09'/,     &
        DustName(4) /'DP_26'/,     &
        DustName(5) /'DP_80'/!,     &
        ! DustName(6) /'DP_240'/,    &
        ! DustName(7) /'DP_720'/,    &
        ! DustName(8) /'DP_2200'/

  ! TYPE dust_fx
  ! REAL(8), POINTER ::   &
  !   soiltype(:,:),      & ! soil properties
  !   z0  (:,:),          & ! roughness length
  !   source  (:,:),      & ! preferential dust source
  !   alpha   (:,:),      & ! ratio horiz/vertical flux
  !   feff    (:,:,:),    & ! drag partition
  !   veff    (:,:,:),    & ! effective vegetation
  !   d_emis  (:,:,:)       ! dust emission
  ! END TYPE dust_fx
  ! TYPE (dust_fx), ALLOCATABLE, TARGET :: dust_flux(:)

  TYPE dust_subdomain
  REAL(8), POINTER ::     &
    soilprop (:,:,:,:),   & !soil properties
    soilmap(:,:,:),       & ! sand, silt, clay map
    lai (:,:,:,:),        & !leafe area index
    vegmin (:,:,:),       & !minimum of vegetation
    alpha (:,:,:),        & !ratio horiz/vertical flux
    c_eff (:,:,:),        & !fraction efficace
    lai_eff (:,:,:,:),    & !effective surface for dust deflation from LAI condition
    w_str (:,:),        & !threshold soil moisture w' (Fecan, F. et al., 1999)
    umin2(:,:,:),         &
    d_emis(:,:,:),         & !dust emission
    biome(:,:),         &
    cult(:,:),          &
    veg (:,:,:),        & !leafe area index
    vmoist(:,:,:),      &
    vegmin2(:,:),        &
    soiltype(:,:),      & ! soil properties
    z0  (:,:),          & ! roughness length
    source  (:,:),      & ! preferential dust source
    alpha2   (:,:),      & ! ratio horiz/vertical flux
    feff_z0    (:,:),    & ! drag partition by roughness
    feff_veg   (:,:,:),    & ! drag partition by vegetation
    veff    (:,:,:),    &   ! effective vegetation
    mfac    (:,:),&          ! moisture factore
    ustar(:,:),&          ! (j,i,)
    srel_map(:,:,:),&          ! (j,i,nclass)
    mrel_map(:,:,:),&          ! (j,i,nclass)
    mrel_sum(:,:,:),&          ! (j,i,nclass)
    mrel_mx(:,:,:,:)          ! (j,i,nclass,nclass)
  END TYPE dust_subdomain
  TYPE (dust_subdomain), ALLOCATABLE, TARGET :: dust(:)

  !-- fixed dimensions for soil data
        INTEGER, PARAMETER :: SoilNumb = 5    & ! number of soil properties
  &                          ,SoilFrac = 1      ! number of soil fraction



  ! soil class properties
  INTEGER, PARAMETER :: &
    nats   = 45,        & ! amount of soil types
    nclass = 196          ! amount of particule classes
    ! !nclass = 392          ! amount of particule classes
    ! nclass = 14          ! amount of particule classes

  INTEGER :: &
    nmode          ! number of soil modes, depending on the soil data set
                   ! nmode = 3 for soilgrids data
                   ! nmode = 4 for the lookup table

  REAL, ALLOCATABLE, TARGET :: &
    median_dp(:)   ! median particle diameter of mode
                   ! IF (nmode = 3) median_dp = (707.0E-6,158.0E-6,15.0E-6,2.0E-6)
                   ! IF (nmode = 4) median_dp = (158.0E-6,15.0E-6,2.0E-6)

  ! dummy variable for input, allocate new when switch from 2d to 3d
  REAL(8), ALLOCATABLE ::  &
    read_input(:,:,:)       ! (j,i,time)

  ! REAL(8), ALLOCATABLE ::  &
  !   srel_map(:,:,:),&          ! (j,i,nclass)
  !   mrel_map(:,:,:),&          ! (j,i,nclass)
  !   mrel_sum(:,:,:),&          ! (j,i,nclass)
  !   mrel_mx(:,:,:,:)          ! (j,i,nclass,nclass)

  REAL (8)   :: &
    dp_meter(nclass) ! particle diameter [m]

  ! REAL(8), ALLOCATABLE ::  &
  !   ustar(:,:)          ! (j,i,nclass)

  REAL (8)   :: &
    Uth(Nclass)!,              & ! threshold friction velocity

  ! 2D Arrays
  REAL(8) ::                  &
    srel(nats,nclass),        & !
    srelV(nats,nclass),       & !
    su_srelV(nats,nclass)!,    & !

  ! internal switches
  LOGICAL     :: &
    lvegdaily,     &   ! daily or monthly input of lai
    laidaily,    & ! daily or monthly input of lai     MF
    lvegmin,     &        ! daily or monthly input of lai     MF
    lddebug       ! additional output


  !-- LM grid arrays
        REAL(8), POINTER ::   &
  &              SG_soilprop (:,:,:,:),     & ! soil properties
  &              SG_lai (:,:,:,:),          & ! leafe area  index
  &              SG_vegmin (:,:,:)             ! minimum of vegetation




  !---------------------------------------------------------------------
  ! !-- soil class properties
  !       INTEGER, PARAMETER :: nats   = 45     & ! amount of soil types, originally "nats = 12", "nats = 14"
  ! &                          ,nClass = 196      ! amount of particule classes
  !       REAL(8) :: srel(nats,Nclass)          & !
  ! &               ,srelV(nats,Nclass)         & !
  ! &               ,su_srelV(nats,nclass)      & !
  ! 	&               ,Uth(Nclass)                & ! threshold friction velocity
  ! 	&               ,Uth_bod(Nclass)              ! threshold friction velocity

  ! additional namelist variables for the offline model
#ifdef OFFLINE
  INTEGER ::   &
    ie_tot,    & ! number of lon points
    je_tot!,    & ! number of lat points

  REAL(8) ::    &
    startlon_tot, & ! lon of low left corner
    startlat_tot, & ! lat of low left corner
    timeinc,   &    ! increment of time in hours
    pollon,    &
    pollat,    &
    dlon,      &
    dlat,      &
    uconst,    &
    vconst,    &
    ustconst,  &
    dzconst


  CHARACTER(10) :: &
    ydate_ini, & ! Date when model Start YYYYMMDDHH
    ydate_end    ! Date when model END   YYYYMMDDHH

  CHARACTER(120) :: &
    windFile,       & ! wind input data
    youtdir,        & ! directory of the output file
    u_var_name,      &
    v_var_name,     &
    ust_var_name

  CHARACTER(3) :: &
    leveltype       ! 10m, or mlv

  LOGICAL :: &
    laccumulation
#endif

END MODULE mo_dust
