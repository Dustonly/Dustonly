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
    mineralmaptype, &   !0: none, 1:GMINER data
    threshold_scheme   ! 0 : Marticorena, 1 : Shao  !

  REAL :: &
    moistinc, &    ! time increment of soil moisture
    z0synop,  &
    z0const,  &
    SSR_min


  !-- Files with soil data
  CHARACTER(120) ::     &
    soiltypeFile,       & ! Filename of Soil Type Data
    mineraltypeFile,    & ! Filename of Mineralogical data
    psrcFile,           & ! Filename of preferential Dust Sources
    cultFile,           & ! Filename of Cultivation Class
    vegmonFile,         & ! Filename of monthly vegitation cover
    vegdayFile,         & ! Filename of daily vegetation cover
    vegminFile,         & ! Filename of min vegetation cover MF
    z0File,             & ! Filename of Roughness Length
    biomeFile,          & ! Filename of Vegetation Cover/Type Data
    moistFile


  ! description of dust particles for external use
  INTEGER, PARAMETER :: &
    DustBins = 5, & !8  ! number of dust particle fractions
    nmin = 18       !number of possible minerals (GMINER) plus first position for total dust in bins

  INTEGER :: bins

  INTEGER, DIMENSION(DustBins,nmin) :: DustInd          ! indices of dust particles, added nmin for the minerals
  DATA (DustInd(bins, 1), bins=1,DustBins)/ &
    1, 2, 3, 4, 5/
  DATA (DustInd(bins, 2), bins=1,DustBins)/ &
    6, 23, 40, 57, 74/
  DATA (DustInd(bins, 3), bins=1,DustBins)/ &
    7, 24, 41, 58, 75/
  DATA (DustInd(bins, 4), bins=1,DustBins)/ &
    8, 25, 42, 59, 76/
  DATA (DustInd(bins, 5), bins=1,DustBins)/ &
    9, 26, 43, 60, 77/
  DATA (DustInd(bins, 6), bins=1,DustBins)/ &
    10, 27, 44, 61, 78/
  DATA (DustInd(bins, 7), bins=1,DustBins)/ &
    11, 28, 45, 62, 79/
  DATA (DustInd(bins, 8), bins=1,DustBins)/ &
    12, 29, 46, 63, 80/
  DATA (DustInd(bins, 9), bins=1,DustBins)/ &
    13, 30, 47, 64, 81/
  DATA (DustInd(bins, 10), bins=1,DustBins)/ &
    14, 31, 48, 65, 82/
  DATA (DustInd(bins, 11), bins=1,DustBins)/ &
    15, 32, 49, 66, 83/
  DATA (DustInd(bins, 12), bins=1,DustBins)/ &
    16, 33, 50, 67, 84/
  DATA (DustInd(bins, 13), bins=1,DustBins)/ &
    17, 34, 51, 68, 85/
  DATA (DustInd(bins, 14), bins=1,DustBins)/ &
    18, 35, 52, 69, 86/
  DATA (DustInd(bins, 15), bins=1,DustBins)/ &
    19, 36, 53, 70, 87/
  DATA (DustInd(bins, 16), bins=1,DustBins)/ &
    20, 37, 54, 71, 88/
  DATA (DustInd(bins, 17), bins=1,DustBins)/ &
    21, 38, 55, 72, 89/
  DATA (DustInd(bins, 18), bins=1,DustBins)/ &
    22, 39, 56, 73, 90/

  REAL :: dustbin_top(DustBins)
  DATA dustbin_top(1) /1.E-6/,  &
       dustbin_top(2) /3.E-6/,  &
       dustbin_top(3) /9.E-6/,  &
       dustbin_top(4) /26.E-6/, &
       dustbin_top(5) /80.E-6/

  CHARACTER(20), DIMENSION(DustBins,nmin) :: DustName  !the position 1 is for non mineral bins
    DATA (DustName(1,minerals), minerals=1,nmin) / &
      'DP_01', 'DP_01_ili', 'DP_01_kao', 'DP_01_smc',      &
      'DP_01_cal', 'DP_01_qua', 'DP_01_hem', 'DP_01_fel',    &
      'DP_01_gyp', 'DP_01_ill','DP_01_kal','DP_01_sme',    &
      'DP_01_clc', 'DP_01_qur', 'DP_01_hmm','DP_01_fld',   &
      'DP_01_gps', 'DP_01_pho'/
    DATA (DustName(2,minerals), minerals=1,nmin) / &
      'DP_03', 'DP_03_ili', 'DP_03_kao', 'DP_03_smc',      &
      'DP_03_cal', 'DP_03_qua', 'DP_03_hem', 'DP_03_fel',    &
      'DP_03_gyp', 'DP_03_ill','DP_03_kal','DP_03_sme',    &
      'DP_03_clc', 'DP_03_qur', 'DP_03_hmm','DP_03_fld',   &
      'DP_03_gps', 'DP_03_pho'/
    DATA (DustName(3,minerals), minerals=1,nmin) / &
      'DP_09', 'DP_09_ili', 'DP_09_kao', 'DP_09_smc',      &
      'DP_09_cal', 'DP_09_qua', 'DP_09_hem', 'DP_09_fel',    &
      'DP_09_gyp', 'DP_09_ill','DP_09_kal','DP_09_sme',    &
      'DP_09_clc', 'DP_09_qur', 'DP_09_hmm','DP_09_fld',   &
      'DP_09_gps', 'DP_09_pho'/
    DATA (DustName(4,minerals), minerals=1,nmin) / &
      'DP_26', 'DP_26_ili', 'DP_26_kao', 'DP_26_smc',      &
      'DP_26_cal', 'DP_26_qua', 'DP_26_hem', 'DP_26_fel',    &
      'DP_26_gyp', 'DP_26_ill','DP_26_kal','DP_26_sme',    &
      'DP_26_clc', 'DP_26_qur', 'DP_26_hmm','DP_26_fld',   &
      'DP_26_gps', 'DP_26_pho'/
    DATA (DustName(5,minerals), minerals=1,nmin) / &
      'DP_80', 'DP_80_ili', 'DP_80_kao', 'DP_80_smc',      &
      'DP_80_cal', 'DP_80_qua', 'DP_80_hem', 'DP_80_fel',    &
      'DP_80_gyp', 'DP_80_ill','DP_80_kal','DP_80_sme',    &
      'DP_80_clc', 'DP_80_qur', 'DP_80_hmm','DP_80_fld',   &
      'DP_80_gps', 'DP_80_pho'/

  INTEGER, PARAMETER :: DustBins_m = 5  !8  ! number of dust particle fractions
  INTEGER :: DustInd_m(DustBins_m)            ! indices of dust particles

  REAL :: dustbin_top_m(DustBins_m)    !for the mineralogy DustBins
  DATA dustbin_top_m(1) /1.E-6/,  &
       dustbin_top_m(2) /3.E-6/,  &
       dustbin_top_m(3) /9.E-6/,  &
       dustbin_top_m(4) /26.E-6/, &
       dustbin_top_m(5) /80.E-6/

  CHARACTER(20) :: DustName_m(DustBins_m)
  DATA  DustName_m(1) /'DP_M_01'/,     &
        DustName_m(2) /'DP_M_03'/,     &
        DustName_m(3) /'DP_M_09'/,     &
        DustName_m(4) /'DP_M_26'/,     &
        DustName_m(5) /'DP_M_80'/


  TYPE dust_subdomain
  REAL(8), POINTER ::     &
    soilprop (:,:,:,:),   & !soil properties
    soilmap(:,:,:),       & ! sand, silt, clay map
    mineralmap(:,:,:),       & ! illite, kaolinite, smectite, feldpsar, calcite, hematite
    mineralmapmod(:,:,:),       & ! modified mineral map
    lai (:,:,:,:),        & !leafe area index
    vegmin (:,:,:),       & !minimum of vegetation
    alpha (:,:,:),        & !ratio horiz/vertical flux
    c_eff (:,:,:),        & !fraction efficace
    lai_eff (:,:,:,:),    & !effective surface for dust deflation from LAI condition
    w_str (:,:),        & !threshold soil moisture w' (Fecan, F. et al., 1999)
    umin2(:,:,:),         &
    d_emis(:,:,:),         & !dust emission
    d_emis_m(:,:,:,:),         & !dust mineralogical emission
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
    mrel_mx(:,:,:,:)        ! (j,i,nclass,nclass)
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

  ! dummy variable for input, allocate new wshen switch from 2d to 3d
  REAL(8), ALLOCATABLE ::  &
    read_input(:,:,:)       ! (j,i,time)

  REAL (8)   :: &
    dp_meter(nclass) ! particle diameter [m]

  REAL (8)   :: &
    Uth(Nclass)!,              & ! threshold friction velocity

  ! 2D Arrays
  REAL(8) ::                  &
    srel(nats,nclass),        & !   &
    srelV(nats,nclass),       & !
    su_srelV(nats,nclass)

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
