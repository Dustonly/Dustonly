!================================================================
! $Id: data_dust.f90,v 0 2018/11/08 Faust $
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
!---------------------------------------------------------------------
! Description:
! New module for dust emissions in muscat
! The aim is to simplify the dust Code
! create a simple structur: 1 Modul, 1 organizes routine, logical subroutines called in organize
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
! V0           2018-11-08 Start Coding
! V0.1         2018-11-08 add modules: dust_param_tegen
! V0.2         2018-11-16 add modules: dust_org, dust_tegen_data
!
!
! Code Description:
! Language: Fortran 90.
! Software Standards: COSMO Standards for Source Code Development
!---------------------------------------------------------------------
!
! MODULE dust_org
! !---------------------------------------------------------------------
! ! Description:
! ! This module contains general parameters for dust emisson
! !---------------------------------------------------------------------
! USE mo_dust, ONLY: &
!   DustBins  ! number of dust particle fractions
!
!
!   IMPLICIT NONE
!
!   ! fixed dimensions for soil data
!   INTEGER, PARAMETER :: &
!     SoilNumb = 5,       & ! number of soil properties
!     SoilFrac = 1          ! number of soil fraction
!
!   ! ! MUSCAT grid structure
!   ! TYPE dust_subdomain
!   !   REAL(8), POINTER ::     &
!   !     soilprop (:,:,:,:),   & !soil properties
!   !     lai (:,:,:,:),        & !fraction of vegetation
!   !     vegmin (:,:,:),       & !minimum of vegetation
!   !     alpha (:,:,:),        & !ratio horiz/vertical flux
!   !     c_eff (:,:,:),        & !fraction efficace
!   !     lai_eff (:,:,:,:),    & !effective surface for dust deflation from LAI condition
!   !     w_str (:,:,:),        & !threshold soil moisture w' (Fecan, F. et al., 1999)
!   !     umin2(:,:,:),         &
!   !     d_emis(:,:,:)          !dust emission
!   ! END TYPE dust_subdomain
!   ! TYPE (dust_subdomain), ALLOCATABLE, TARGET :: dust(:)
!
!   ! ! LM grid arrays
!   ! REAL(8), POINTER ::   &
!   !   SG_soilprop (:,:,:,:),     & ! soil properties
!   !   SG_veg (:,:,:,:),          & ! leafe area index
!   !   SG_vegmin (:,:,:)            ! minimum of vegetation
!
!   ! soil class properties
!   INTEGER, PARAMETER :: &
!     nats   = 45,        & ! amount of soil types
!     nclass = 196          ! amount of particule classes
!
!   ! dummy variable for input, allocate new when switch from 2d to 3d
!   REAL(8), ALLOCATABLE ::  &
!     read_input(:,:,:)       ! (j,i,time)
!
!   REAL (8)   :: &
!     Uth(Nclass)!,              & ! threshold friction velocity
!
!   ! 2D Arrays
!   REAL(8) ::                  &
!     srel(nats,nclass),        & !
!     srelV(nats,nclass),       & !
!     su_srelV(nats,nclass)!,    & !
!
!   ! new dust datatype with less in it.
!   TYPE dust_it
!   REAL(8), POINTER ::   &
!     biome(:,:),         &
!     cult(:,:),          &
!     veg (:,:,:),        & !leafe area index
!     vegmin(:,:)
!   END TYPE dust_it
!   TYPE (dust_it), ALLOCATABLE, TARGET :: dust_ini(:)
!
!
!
!   ! internal switches
!   LOGICAL     :: &
!     lvegdaily,     &   ! daily or monthly input of lai
!     laidaily,    & ! daily or monthly input of lai     MF
!     lvegmin        ! daily or monthly input of lai     MF
!
! END MODULE dust_org









MODULE dust_tegen_param
!---------------------------------------------------------------------
! Description:
! This module contains parameters for the Tegen dust emisson scheme
!---------------------------------------------------------------------
  USE mo_dust, ONLY: &
    DustBins  ! number of dust particle fractions

  IMPLICIT NONE

  INTEGER, PARAMETER ::      &
    combimax = 29,           &
    nbin = 24,               & !max number of bins
    nmode = 4,               &
    nspe = 14,               & !Nmode*3+2
    ntrace = DustBins          !max number of bins


  REAL(8), PARAMETER ::                       &
    a_rnolds    = 1331.647E0,                 & !Reynolds constant
    aeff        = 0.35E0,                     & !efficient fraction
    aw          = 1.121E0,                    & !soil moisture parameter
    abs_density = 2.6E0,                      & !absolute density of the soil (2.6 g/cm3)
    b_rnolds    = 0.38194E0,                  & !Reynolds constant
    bw          = 0.68E0,                     & !soil moisture parameter
    Cd          = 1.2507E-06,                 & !flux dimensioning parameter,1.00*roa/g !!BERND
    d_thrsld    = 0.00000231E0,               & !thresold value
    Dmin        = 0.00002E0,                  & !minimum particules diameter (cm)
    Dmax        = 0.155E0,                    & !maximum particules diameter (cm)
    Dstep       = 0.0460517018598807E0,       & !diameter increment (cm)
    e           = 2.7182818285E0,             & !
    g           = 981.E0,                     & !gravitational constant (cm.s-2)
    pi2         = 3.141592654E0,              & !pi=3.141592654
    porosity    = 0.4E0,                      & !porosity of the soil (40%)
    roa         = 0.001227E0,                 & !air density (g.cm-3)

    rop         = 2.65E0,                     & !particules density (g.cm-3) Sahara; 1.5g/cm3 for Diatomite
    rop_bod     = 2.1E0,                      & !particules density (g.cm-3) Bodele (Diatomite)

    u1fac       = 0.89E0,                     & !0.66; 1.0

    umin        = 13.75E0,                    & !minimum threshold friction windspeed (cm/s) 18.06/13.75/21.0
    Ustar_min   = 5.E0,                       & !minimum friction windspeed (cm/s)
    Ustar_max   = 150.E0,                     & !maximum friction windspeed (cm/s)
    Ustar_step  = 2.2674-02,                  & !friction windspeed increment
    VK          = 0.4E0,                      & !Von Karman constant: 0.4 (0.35 <--> 0.42 see Stull)
    w0          = 99.E0,                      & !soil layer water content threshold
    xnu         = 0.146E0,                    & !cinematic viscosity (cm2.s-1)
    x_rnolds    = 1.561228E0,                 & !Reynolds constant
    xeff        = 10.E0,                      & !efficient fraction
    ZZ          = 1000.E0                       !wind measurment height (cm)

  REAL(8), PARAMETER ::                       &
    veg_lim     = 0.6E0,                      & !lai limit for deflation !orgi 0.25 MF
    veg_lim2    = 0.5E0                         !lai limit to determine shrub or grass for tropical or temperate shrubland
  !  ----------------------------------------------------------
END MODULE dust_tegen_param









MODULE dust_tegen_data
!---------------------------------------------------------------------
! Description:
! This module contains parameters for the Tegen dust emisson scheme
!---------------------------------------------------------------------
  USE mo_dust, ONLY: &
    nats
  USE dust_tegen_param, ONLY: &
  combimax,nspe


  ! ------------
  !'mo_2000_veg_dust' contains some classification dependent biome and soil
  ! data needed for the dust emissions scheme.

  !----------VEGETATION PARAMETERS CLASSIFICATION DEPENDENT-----------

  !-   The 28 different biome types output by BIOME4
  !-------------------------------------------------------------------
  !      1       Tropical evergreen broadleaf forest
  !      2       Tropical semi-evergreen broadleaf forest
  !      3       Tropical deciduous broadleaf forest and woodland
  !      4       Temperate deciduous broadleaf forest
  !      5       Temperate evergreen needleleaf forest
  !      6       Warm-temperate evergreen broadleaf and mixed forest
  !      7       Cool mixed forest
  !      8       Cool evergreen needleleaf forest
  !      9       Cool-temperate evergreen needleleaf and mixed forest
  !      10      Cold evergreen needleleaf forest
  !      11      Cold deciduous forest
  !      12      Tropical savanna
  !      13      Tropical xerophytic shrubland
  !      14      Temperate xerophytic shrubland
  !      15      Temperate sclerophyll woodland and shrubland
  !      16      Temperate deciduous broadleaf savanna
  !      17      Temperate evergreen needleleaf open woodland
  !      18      Cold parkland
  !      19      Tropical grassland
  !      20      Temperate grassland
  !      21      Desert
  !      22      Graminoid and forb tundra
  !      23      Low and high shrub tundra
  !      24      Erect dwarf-shrub tundra
  !      25      Prostrate dwarf-shrub tundra
  !      26      Cushion-forb tundra
  !      27      Barren
  !      28      Ice
  !     (29)     Water (this is implied)
  !-------------------------------------------------------------------

  ! active vegetation types:
  !--------------------------

   IMPLICIT NONE

   INTEGER :: jspe

   INTEGER, DIMENSION(combimax) :: active

   data active(1)/1/
   data active(2)/1/
   data active(3)/1/
   data active(4)/1/
   data active(5)/1/
   data active(6)/1/
   data active(7)/1/
   data active(8)/1/
   data active(9)/1/
   data active(10)/1/
   data active(11)/1/
   data active(12)/1/
   data active(13)/1/
   data active(14)/1/
   data active(15)/1/
   data active(16)/1/
   data active(17)/1/
   data active(18)/1/
   data active(19)/1/
   data active(20)/1/
   data active(21)/1/
   data active(22)/1/
   data active(23)/1/
   data active(24)/1/
   data active(25)/1/
   data active(26)/1/
   data active(27)/1/
   data active(28)/0/
   data active(29)/0/

  ! Biomes including shrubs (active=1 only)
  !----------------------------------------

   REAL(8), DIMENSION(combimax) :: shrub

   data shrub(1)/1.E0/
   data shrub(2)/1.E0/
   data shrub(3)/1.E0/
   data shrub(4)/1.E0/
   data shrub(5)/1.E0/
   data shrub(6)/1.E0/
   data shrub(7)/1.E0/
   data shrub(8)/1.E0/
   data shrub(9)/1.E0/
   data shrub(10)/1.E0/
   data shrub(11)/1.E0/
   data shrub(12)/1.E0/
   data shrub(13)/1.E0/
   data shrub(14)/1.E0/
   data shrub(15)/1.E0/
   data shrub(16)/1.E0/
   data shrub(17)/1.E0/
   data shrub(18)/1.E0/
   data shrub(19)/0.E0/
   data shrub(20)/0.E0/
   data shrub(21)/0.E0/
   data shrub(22)/1.E0/
   data shrub(23)/1.E0/
   data shrub(24)/1.E0/
   data shrub(25)/1.E0/
   data shrub(26)/1.E0/
   data shrub(27)/0.E0/
   data shrub(28)/0.E0/
   data shrub(29)/0.E0/

  !  Vegetation albedo:
  !--------------------

   REAL(8), DIMENSION(combimax) :: vec_Av

   data vec_Av(1)/0.12E0/
   data vec_Av(2)/0.12E0/
   data vec_Av(3)/0.12E0/
   data vec_Av(4)/0.12E0/
   data vec_Av(5)/0.12E0/
   data vec_Av(6)/0.12E0/
   data vec_Av(7)/0.12E0/
   data vec_Av(8)/0.12E0/
   data vec_Av(9)/0.12E0/
   data vec_Av(10)/0.12E0/
   data vec_Av(11)/0.12E0/
   data vec_Av(12)/0.12E0/
   data vec_Av(13)/0.12E0/
   data vec_Av(14)/0.12E0/
   data vec_Av(15)/0.12E0/
   data vec_Av(16)/0.12E0/
   data vec_Av(17)/0.12E0/
   data vec_Av(18)/0.12E0/
   data vec_Av(19)/0.12E0/
   data vec_Av(20)/0.12E0/
   data vec_Av(21)/0.12E0/
   data vec_Av(22)/0.12E0/
   data vec_Av(23)/0.12E0/
   data vec_Av(24)/0.12E0/
   data vec_Av(25)/0.12E0/
   data vec_Av(26)/0.12E0/
   data vec_Av(27)/0.E0/
   data vec_Av(28)/0.E0/
   data vec_Av(29)/0.E0/

  ! Soil/liter Albedo:
  !--------------------

   REAL(8), DIMENSION(combimax):: vec_Ag

   data vec_Ag(1)/0.11E0/
   data vec_Ag(2)/0.11E0/
   data vec_Ag(3)/0.11E0/
   data vec_Ag(4)/0.11E0/
   data vec_Ag(5)/0.11E0/
   data vec_Ag(6)/0.11E0/
   data vec_Ag(7)/0.11E0/
   data vec_Ag(8)/0.11E0/
   data vec_Ag(9)/0.11E0/
   data vec_Ag(10)/0.11E0/
   data vec_Ag(11)/0.11E0/
   data vec_Ag(12)/0.15E0/
   data vec_Ag(13)/0.15E0/
   data vec_Ag(14)/0.15E0/
   data vec_Ag(15)/0.15E0/
   data vec_Ag(16)/0.15E0/
   data vec_Ag(17)/0.11E0/
   data vec_Ag(18)/0.11E0/
   data vec_Ag(19)/0.15E0/
   data vec_Ag(20)/0.15E0/
   data vec_Ag(21)/0.30E0/
   data vec_Ag(22)/0.11E0/
   data vec_Ag(23)/0.11E0/
   data vec_Ag(24)/0.11E0/
   data vec_Ag(25)/0.11E0/
   data vec_Ag(26)/0.11E0/
   data vec_Ag(27)/0.E0/
   data vec_Ag(28)/0.E0/
   data vec_Ag(29)/0.E0/

  ! Vegetation roughness length
  !-----------------------------

   REAL(8), DIMENSION(combimax):: vec_z0v

   data vec_z0v(1)/2.E0/
   data vec_z0v(2)/2.E0/
   data vec_z0v(3)/2.E0/
   data vec_z0v(4)/1.E0/
   data vec_z0v(5)/1.E0/
   data vec_z0v(6)/1.E0/
   data vec_z0v(7)/1.E0/
   data vec_z0v(8)/1.E0/
   data vec_z0v(9)/1.E0/
   data vec_z0v(10)/1.E0/
   data vec_z0v(11)/0.5E0/
   data vec_z0v(12)/0.25E0/
   data vec_z0v(13)/0.1E0/   !/0.25E0/
   data vec_z0v(14)/0.1E0/   !/0.25E0/
   data vec_z0v(15)/0.25E0/
   data vec_z0v(16)/0.25E0/
   data vec_z0v(17)/0.25E0/
   data vec_z0v(18)/0.25E0/
   data vec_z0v(19)/0.06E0/   !/0.10E0/
   data vec_z0v(20)/0.06E0/   !/0.10E0/
   data vec_z0v(21)/0.005E0/   !/0.05E0/
   data vec_z0v(22)/0.10E0/
   data vec_z0v(23)/0.10E0/
   data vec_z0v(24)/0.05E0/
   data vec_z0v(25)/0.05E0/
   data vec_z0v(26)/0.05E0/
   data vec_z0v(27)/0.E0/
   data vec_z0v(28)/0.E0/
   data vec_z0v(29)/0.E0/

  ! Vegetation LAI MAXIMUM & MINIMUM:
  !-----------------------------------
  ! MAXIMUM:

   REAL(8), DIMENSION(combimax) :: max_LAI

   data max_LAI(1)/0.12E0/
   data max_LAI(2)/0.12E0/
   data max_LAI(3)/0.12E0/
   data max_LAI(4)/0.12E0/
   data max_LAI(5)/0.12E0/
   data max_LAI(6)/0.12E0/
   data max_LAI(7)/0.12E0/
   data max_LAI(8)/0.12E0/
   data max_LAI(9)/0.12E0/
   data max_LAI(10)/0.12E0/
   data max_LAI(11)/0.12E0/
   data max_LAI(12)/0.12E0/
   data max_LAI(13)/0.12E0/
   data max_LAI(14)/0.12E0/
   data max_LAI(15)/0.12E0/
   data max_LAI(16)/0.12E0/
   data max_LAI(17)/0.12E0/
   data max_LAI(18)/0.12E0/
   data max_LAI(19)/0.12E0/
   data max_LAI(20)/0.12E0/
   data max_LAI(21)/0.12E0/
   data max_LAI(22)/0.12E0/
   data max_LAI(23)/0.12E0/
   data max_LAI(24)/0.12E0/
   data max_LAI(25)/0.12E0/
   data max_LAI(26)/0.12E0/
   data max_LAI(27)/0.E0/
   data max_LAI(28)/0.E0/
   data max_LAI(29)/0.E0/

  ! MINIMUM:

   REAL(8), DIMENSION(combimax) :: min_LAI

   data min_LAI(1)/0.12E0/
   data min_LAI(2)/0.12E0/
   data min_LAI(3)/0.12E0/
   data min_LAI(4)/0.12E0/
   data min_LAI(5)/0.12E0/
   data min_LAI(6)/0.12E0/
   data min_LAI(7)/0.12E0/
   data min_LAI(8)/0.12E0/
   data min_LAI(9)/0.12E0/
   data min_LAI(10)/0.12E0/
   data min_LAI(11)/0.12E0/
   data min_LAI(12)/0.12E0/
   data min_LAI(13)/0.12E0/
   data min_LAI(14)/0.12E0/
   data min_LAI(15)/0.12E0/
   data min_LAI(16)/0.12E0/
   data min_LAI(17)/0.12E0/
   data min_LAI(18)/0.12E0/
   data min_LAI(19)/0.12E0/
   data min_LAI(20)/0.12E0/
   data min_LAI(21)/0.12E0/
   data min_LAI(22)/0.12E0/
   data min_LAI(23)/0.12E0/
   data min_LAI(24)/0.12E0/
   data min_LAI(25)/0.12E0/
   data min_LAI(26)/0.12E0/
   data min_LAI(27)/0.E0/
   data min_LAI(28)/0.E0/
   data min_LAI(29)/0.E0/

  ! Vegetation fractional cover MAXIMUM & MINIMUM:
  ! MAXIMUM:

   REAL(8), DIMENSION(combimax) :: max_veg0

   data max_veg0(1)/0.12E0/
   data max_veg0(2)/0.12E0/
   data max_veg0(3)/0.12E0/
   data max_veg0(4)/0.12E0/
   data max_veg0(5)/0.12E0/
   data max_veg0(6)/0.12E0/
   data max_veg0(7)/0.12E0/
   data max_veg0(8)/0.12E0/
   data max_veg0(9)/0.12E0/
   data max_veg0(10)/0.12E0/
   data max_veg0(11)/0.12E0/
   data max_veg0(12)/0.12E0/
   data max_veg0(13)/0.12E0/
   data max_veg0(14)/0.12E0/
   data max_veg0(15)/0.12E0/
   data max_veg0(16)/0.12E0/
   data max_veg0(17)/0.12E0/
   data max_veg0(18)/0.12E0/
   data max_veg0(19)/0.12E0/
   data max_veg0(20)/0.12E0/
   data max_veg0(21)/0.12E0/
   data max_veg0(22)/0.12E0/
   data max_veg0(23)/0.12E0/
   data max_veg0(24)/0.12E0/
   data max_veg0(25)/0.12E0/
   data max_veg0(26)/0.12E0/
   data max_veg0(27)/0.E0/

  ! MINIMUM:

   REAL(8), DIMENSION(combimax) :: min_veg0

   data min_veg0(1)/0.12E0/
   data min_veg0(2)/0.12E0/
   data min_veg0(3)/0.12E0/
   data min_veg0(4)/0.12E0/
   data min_veg0(5)/0.12E0/
   data min_veg0(6)/0.12E0/
   data min_veg0(7)/0.12E0/
   data min_veg0(8)/0.12E0/
   data min_veg0(9)/0.12E0/
   data min_veg0(10)/0.12E0/
   data min_veg0(11)/0.12E0/
   data min_veg0(12)/0.12E0/
   data min_veg0(13)/0.12E0/
   data min_veg0(14)/0.12E0/
   data min_veg0(15)/0.12E0/
   data min_veg0(16)/0.12E0/
   data min_veg0(17)/0.12E0/
   data min_veg0(18)/0.12E0/
   data min_veg0(19)/0.12E0/
   data min_veg0(20)/0.12E0/
   data min_veg0(21)/0.12E0/
   data min_veg0(22)/0.12E0/
   data min_veg0(23)/0.12E0/
   data min_veg0(24)/0.12E0/
   data min_veg0(25)/0.12E0/
   data min_veg0(26)/0.12E0/
   data min_veg0(27)/0.E0/
   data min_veg0(28)/0.E0/
   data min_veg0(29)/0.E0/

  !----------------------------------------------------------------
  !solspe --> SOIL CARACTERISTICS:
  !
  !   ZOBLER texture classes
  !   SOLSPE: for 4 populations : values = 3*(Dmed sig p); ratio of fluxes; residual moisture
  !   Populations: Coarse sand, medium/fine sand, Silt, Clay
  !
  !
  ! MF
  ! change alpha values to more accurate ones
  ! Tegen Habil alpha
  ! Coarse Sand       -> 1.e-7
  ! Medium/Fine Sand  -> 1.e-6
  ! Silt              -> 1.e-5
  ! Clay (if < 45%)   -> 1.e-6
  ! clay (if > 45%)   -> 1.e-7
  !
  !----------------------------------------------------------------

   REAL(8), DIMENSION(nats,nspe) :: solspe

  !--     soil type 1 : Coarse
            data (solspe(1,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.43E0 ,          &
             0.0158E0, 2.E0, 0.4E0 ,           &
             0.0015E0, 2.E0, 0.17E0 ,          &
             0.0002E0, 2.E0, 0.E0 ,            &
             2.143E-06,  0.2E0/                   !MF 2.1E-06

  !--     soil type 2 : Medium
            data (solspe(2,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.E0 ,            &
             0.0158E0, 2.E0, 0.37E0 ,          &
             0.0015E0, 2.E0, 0.33E0 ,          &
             0.0002E0, 2.E0, 0.3E0 ,           &
             3.97E-06,   0.25E0/                   !MF 4.0E-06

  !--     soil type 3 : Fine
            data (solspe(3,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.E0 ,            &
             0.0158E0, 2.E0, 0.E0 ,            &
             0.0015E0, 2.E0, 0.33E0 ,          &
             0.0002E0, 2.E0, 0.67E0 ,          &
             3.367E-06,   0.5E0/                   !MF 1.0E-07

  !--     soil type 4 : Coarse Medium
            data (solspe(4,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.1E0 ,           &
             0.0158E0, 2.E0, 0.5E0 ,           &
             0.0015E0, 2.E0, 0.2E0 ,           &
             0.0002E0, 2.E0, 0.2E0 ,           &
             2.71E-06,  0.23E0/                   !MF 2.7E-06

  !--     soil type 5 : Coarse Fine
            data (solspe(5,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.E0 ,            &
             0.0158E0, 2.E0, 0.5E0 ,           &
             0.0015E0, 2.E0, 0.12E0 ,          &
             0.0002E0, 2.E0, 0.38E0 ,          &
             2.08E-06,  0.25E0/                   !MF 2.8E-06

  !--     soil type 6 : Medium Fine
            data (solspe(6,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.E0 ,            &
             0.0158E0, 2.E0, 0.27E0 ,          &
             0.0015E0, 2.E0, 0.25E0 ,          &
             0.0002E0, 2.E0, 0.48E0 ,          &
             2.818e-07,   0.36E0/                   !MF 1.0E-07

  !--     soil type 7 : Coarse, Medium, Fine
            data (solspe(7,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.23E0 ,          &
             0.0158E0, 2.E0, 0.23E0 ,          &
             0.0015E0, 2.E0, 0.19E0 ,          &
             0.0002E0, 2.E0, 0.35E0 ,          &
             2.503E-06,  0.25E0/                   !MF 2.5E-06

  !--     soil type 8 : Organic
            data (solspe(8,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.25E0 ,          &
             0.0158E0, 2.E0, 0.25E0 ,          &
             0.0015E0, 2.E0, 0.25E0 ,          &
             0.0002E0, 2.E0, 0.25E0 ,          &
             0.E0,     0.5E0/

  !--     soil type 9 : Ice
            data (solspe(9,jspe),jspe=1,nspe)/ &
             0.0707E0,  2.E0, 0.25E0 ,         &
             0.0158E0,  2.E0, 0.25E0 ,         &
             0.0015E0,  2.E0, 0.25E0 ,         &
             0.0002E0,  2.E0,  0.25E0 ,        &
             0.E0,      0.5E0/

  !--     soil type 10 : Potential Lakes (additional)
  !          GENERAL CASE
            data (solspe(10,jspe),jspe=1,nspe)/ &
             0.0707E0,  2.E0, 0.E0 ,            &
             0.0158E0,  2.E0, 0.E0,             &
             0.0015E0,  2.E0, 1.E0 ,            &
             0.0002E0, 2.E0,  0.E0 ,            &
             1.E-05,   0.25E0/

  !--     soil type 11 : Potential Lakes (clay)
  !          GENERAL CASE
            data (solspe(11,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.E0 ,             &
             0.0158E0, 2.E0, 0.E0 ,             &
             0.0015E0, 2.E0, 0.E0 ,             &
             0.0002E0, 2.E0, 1.E0 ,             &
             1.E-07,   0.25E0/                   !MF 1.0E-05


  !--     soil type 12 : Potential Lakes Australia
            data (solspe(12,jspe),jspe=1,nspe)/ &
             0.0707E0,  2.E0, 0.E0 ,            &
             0.0158E0,  2.E0, 0.E0 ,            &
             0.0027E0,  2E0,  1.E0 ,            &
             0.0002E0,  2E0,   0.E0 ,           &
             1.E-05,    0.25E0/

  !--     soil type 13 : Potential Lakes (additional)
  !       adapted to BODELE
           data (solspe(13,jspe),jspe=1,nspe)/ &
            0.02E0,    2.E0, 1.E0 ,            &
            0.0158E0,  2.E0, 0.E0 ,            &
            0.0015E0,  2.E0, 0.E0 ,            &
            0.0002E0,  2.E0, 0.E0 ,            &
            1.8E-06,   0.25E0/
  !       2.2E-06,  0.25/

  !--     soil type 14 : Potential Lakes (clay)
  !       adapted to BODELE
           data (solspe(14,jspe),jspe=1,nspe)/ &
            0.002E0,  2.E0, 0.7E0 ,            &
            0.0002E0, 2.E0, 0.3E0 ,            &
            0.0015E0, 2.E0, 0.E0  ,            &
            0.0002E0, 2.E0, 0.E0  ,            &
            1.E-05,   0.25E0/


  ! -----
  ! New soil types based on esdb
  ! 5 types + 26 mixed types

  !--     soil type 15 : Coarse
           data (solspe(15,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.4E0 ,           &
             0.0158E0, 2.E0, 0.375E0 ,         &
             0.0015E0, 2.E0, 0.135E0 ,         &
             0.0002E0, 2.E0, 0.09E0 ,          &
             1.86E-06,  0.2E0/

  !--     soil type 16 : Medium
           data (solspe(16,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.1E0 ,           &
             0.0158E0, 2.E0, 0.3E0 ,           &
             0.0015E0, 2.E0, 0.425E0 ,         &
             0.0002E0, 2.E0, 0.175E0 ,         &
             4.74E-06,  0.2E0/

  !--     soil type 17 : MedFine
           data (solspe(17,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.E0 ,            &
             0.0158E0, 2.E0, 0.075E0 ,         &
             0.0015E0, 2.E0, 0.75E0 ,          &
             0.0002E0, 2.E0, 0.175E0 ,         &
             7.75E-06,  0.2E0/

  !--     soil type 18 : Fine
           data (solspe(18,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.E0 ,            &
             0.0158E0, 2.E0, 0.2625E0 ,        &
             0.0015E0, 2.E0, 0.2625E0 ,        &
             0.0002E0, 2.E0, 0.475E0 ,         &
             3.36E-06,  0.2E0/

  !--     soil type 19 : VeryFine
           data (solspe(19,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.0E0 ,           &
             0.0158E0, 2.E0, 0.133E0 ,         &
             0.0015E0, 2.E0, 0.133E0 ,         &
             0.0002E0, 2.E0, 0.734E0 ,         &
             2.2E-06,  0.2E0/

  !--     soil type 20 : Coarse, Medium
           data (solspe(20,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.25E0 ,          &
             0.0158E0, 2.E0, 0.3375E0 ,        &
             0.0015E0, 2.E0, 0.28E0 ,          &
             0.0002E0, 2.E0, 0.1325E0 ,        &
             3.3E-06,  0.2E0/

  !--     soil type 21 : Coarse, MedFine
           data (solspe(21,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.2E0 ,           &
             0.0158E0, 2.E0, 0.225E0 ,         &
             0.0015E0, 2.E0, 0.4425E0 ,        &
             0.0002E0, 2.E0, 0.1325E0 ,        &
             4.8E-06,  0.2E0/

  !--     soil type 22 : Coarse, Fine
           data (solspe(22,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.2E0 ,           &
             0.0158E0, 2.E0, 0.31875E0 ,       &
             0.0015E0, 2.E0, 0.19875E0 ,       &
             0.0002E0, 2.E0, 0.2825E0 ,        &
             2.61E-06,  0.2E0/

  !--     soil type 23 : Coarse, VeryFine
           data (solspe(23,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.2E0 ,           &
             0.0158E0, 2.E0, 0.254E0 ,         &
             0.0015E0, 2.E0, 0.134E0 ,         &
             0.0002E0, 2.E0, 0.412E0 ,         &
             2.03E-06,  0.2E0/

  !--     soil type 24 : Medium, MedFine
           data (solspe(24,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.05E0 ,          &
             0.0158E0, 2.E0, 0.1875E0 ,        &
             0.0015E0, 2.E0, 0.5875E0 ,        &
             0.0002E0, 2.E0, 0.175E0 ,         &
             6.24E-06,  0.2E0/

  !--     soil type 25 : Medium, Fine
           data (solspe(25,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.05E0 ,          &
             0.0158E0, 2.E0, 0.28125E0 ,       &
             0.0015E0, 2.E0, 0.34375E0 ,       &
             0.0002E0, 2.E0, 0.325E0 ,         &
             4.05E-06,  0.2E0/

  !--     soil type 26 : Medium, VeryFine
           data (solspe(26,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.05E0 ,          &
             0.0158E0, 2.E0, 0.2165E0 ,        &
             0.0015E0, 2.E0, 0.279E0 ,         &
             0.0002E0, 2.E0, 0.4545E0 ,        &
             3.06E-06,  0.2E0/

  !--     soil type 27 : MedFine, Fine
           data (solspe(27,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.E0 ,            &
             0.0158E0, 2.E0, 0.16875E0 ,       &
             0.0015E0, 2.E0, 0.50625E0 ,       &
             0.0002E0, 2.E0, 0.325E0 ,         &
             5.56E-06,  0.2E0/

  !--     soil type 28 : MedFine, VeryFine
           data (solspe(28,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.E0 ,            &
             0.0158E0, 2.E0, 0.104E0 ,         &
             0.0015E0, 2.E0, 0.4415E0 ,        &
             0.0002E0, 2.E0, 0.4545E0 ,        &
             4.56E-06,  0.2E0/

  !--     soil type 29 : Fine, VeryFine
           data (solspe(29,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.E0 ,            &
             0.0158E0, 2.E0, 0.19775E0 ,       &
             0.0015E0, 2.E0, 0.19775E0 ,       &
             0.0002E0, 2.E0, 0.6045E0 ,        &
             2.24E-06,  0.2E0/

  !--     soil type 30 : Coarse, Medium, MedFine
           data (solspe(30,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.17E0 ,          &
             0.0158E0, 2.E0, 0.25E0 ,          &
             0.0015E0, 2.E0, 0.43E0 ,          &
             0.0002E0, 2.E0, 0.15E0 ,          &
             4.78E-06,  0.2E0/

  !--     soil type 31 : Coarse, Medium, Fine
           data (solspe(31,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.17E0 ,          &
             0.0158E0, 2.E0, 0.3125E0 ,        &
             0.0015E0, 2.E0, 0.2775E0 ,        &
             0.0002E0, 2.E0, 0.24E0 ,          &
             3.32E-06,  0.2E0/

  !--     soil type 32 : Coarse, Medium, VeryFine
           data (solspe(32,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.17E0 ,          &
             0.0158E0, 2.E0, 0.27E0 ,          &
             0.0015E0, 2.E0, 0.23E0 ,          &
             0.0002E0, 2.E0, 0.33E0 ,          &
             2.93E-06,  0.2E0/

  !--     soil type 33 : Coarse, MedFine, Fine
           data (solspe(33,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.13E0 ,          &
             0.0158E0, 2.E0, 0.2375E0 ,        &
             0.0015E0, 2.E0, 0.3825E0 ,        &
             0.0002E0, 2.E0, 0.25E0 ,          &
             4.32E-06,  0.2E0/

  !--     soil type 34 : Coarse, MedFine, VeryFine
           data (solspe(34,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.13E0 ,          &
             0.0158E0, 2.E0, 0.20E0 ,          &
             0.0015E0, 2.E0, 0.34E0 ,          &
             0.0002E0, 2.E0, 0.33E0 ,          &
             3.93E-06,  0.2E0/

  !--     soil type 35 : Coarse, Fine, VeryFine
           data (solspe(35,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.13E0 ,          &
             0.0158E0, 2.E0, 0.26E0 ,          &
             0.0015E0, 2.E0, 0.18E0 ,          &
             0.0002E0, 2.E0, 0.43E0 ,          &
             2.47E-06,  0.2E0/

  !--     soil type 36 : Medium, MedFine, Fine
           data (solspe(36,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.03E0 ,          &
             0.0158E0, 2.E0, 0.2125E0 ,        &
             0.0015E0, 2.E0, 0.485E0 ,         &
             0.0002E0, 2.E0, 0.275E0 ,         &
             5.28E-06,  0.2E0/

  !--     soil type 37 : Medium, MedFine, VeryFine
           data (solspe(37,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.03E0 ,          &
             0.0158E0, 2.E0, 0.17E0 ,          &
             0.0015E0, 2.E0, 0.44E0 ,          &
             0.0002E0, 2.E0, 0.36E0 ,          &
             4.89E-06,  0.2E0/

  !--     soil type 38 : Medium, Fine, VeryFine
           data (solspe(38,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.03E0 ,          &
             0.0158E0, 2.E0, 0.24E0 ,          &
             0.0015E0, 2.E0, 0.27E0 ,          &
             0.0002E0, 2.E0, 0.46E0 ,          &
             3.02E-06,  0.2E0/

  !--     soil type 39 : MedFine, Fine, VeryFine
           data (solspe(39,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.E0 ,            &
             0.0158E0, 2.E0, 0.16E0 ,          &
             0.0015E0, 2.E0, 0.38E0 ,          &
             0.0002E0, 2.E0, 0.46E0 ,          &
             4.02E-06,  0.2E0/

  !--     soil type 40 : Coarse, Medium, MedFine, Fine
           data (solspe(40,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.125E0 ,         &
             0.0158E0, 2.E0, 0.253125E0 ,      &
             0.0015E0, 2.E0, 0.393125E0 ,      &
             0.0002E0, 2.E0, 0.22875E0 ,       &
             4.43E-06,  0.2E0/

  !--     soil type 41 : Coarse, Medium, MedFine, VeryFine
           data (solspe(41,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.125E0 ,         &
             0.0158E0, 2.E0, 0.22075E0 ,       &
             0.0015E0, 2.E0, 0.36075E0 ,       &
             0.0002E0, 2.E0, 0.2935E0 ,        &
             4.13E-06,  0.2E0/

  !--     soil type 42 : Coarse, Medium, Fine, VeryFine
           data (solspe(42,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.125E0 ,         &
             0.0158E0, 2.E0, 0.267625E0 ,      &
             0.0015E0, 2.E0, 0.238875E0 ,      &
             0.0002E0, 2.E0, 0.3685E0 ,        &
             3.04E-06,  0.2E0/

  !--     soil type 43 : Coarse, MedFine, Fine, VeryFine
           data (solspe(43,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.1E0 ,           &
             0.0158E0, 2.E0, 0.211375E0 ,      &
             0.0015E0, 2.E0, 0.320125E0 ,      &
             0.0002E0, 2.E0, 0.3685E0 ,        &
             3.79E-06,  0.2E0/

  !--     soil type 44 : Medium, MedFine, Fine, VeryFine
           data (solspe(44,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.025E0 ,         &
             0.0158E0, 2.E0, 0.192625E0 ,      &
             0.0015E0, 2.E0, 0.392625E0 ,      &
             0.0002E0, 2.E0, 0.38975E0 ,       &
             4.51E-06,  0.2E0/

  !--     soil type 45 : Coarse, Medium, MedFine, Fine, VeryFine
           data (solspe(45,jspe),jspe=1,nspe)/ &
             0.0707E0, 2.E0, 0.1E0 ,           &
             0.0158E0, 2.E0, 0.2291E0 ,        &
             0.0015E0, 2.E0, 0.3411E0 ,        &
             0.0002E0, 2.E0, 0.3298E0 ,        &
             3.98E-06,  0.2E0/

END MODULE dust_tegen_data
