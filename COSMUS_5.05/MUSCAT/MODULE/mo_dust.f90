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
 USE data_modelconfig,   ONLY: DustMod ! Flag for Setting Soil Data
!======================================================================================
! Description:
! ------------
!  This module contains definitions of the provisional results required for the
!  dust-flux calculation in sub_dust_emis.f90.
!======================================================================================
!
!---------------------------------------------------------------------
! !-- control parameters
!       LOGICAL ::    &
! &       old_dust          ! if true use the original dust routines
!                           ! if false use new module src_dust
!
!       INTEGER ::    &
! &       nDust,      &    ! Flag for Dust Calculations
! &       psrcType         ! Flag for type of potential dust source   MF
!                          ! 0 : psrc, 1 : msgsrc, 2 : acDust
!
! !-- Files with soil data
!       CHARACTER(120) ::   &
! &           soiltypeFile       & ! Filename of Soil Type Data
! &          ,psrcFile           & ! Filename of preferential Dust Sources
! &          ,cultFile           & ! Filename of Cultivation Class
! &          ,vegmonFile         & ! Filename of monthly vegitation cover
! &          ,vegdayFile         & ! Filename of daily vegetation cover
! &          ,vegminFile         & ! Filename of min vegetation cover MF
! &          ,z0File             & ! Filename of Roughness Length
! &          ,biomeFile            ! Filename of Vegetation Cover/Type Data

!---------------------------------------------------------------------
! !-- description of dust particles
!       INTEGER, PARAMETER :: DustBins = 5  !8  ! number of dust particle fractions
!       INTEGER :: DustInd(DustBins)            ! indices of dust particles
!
!       CHARACTER(20) :: DustName(DustBins)
!       DATA  DustName(1) /'DP_01'/,     &
!             DustName(2) /'DP_03'/,     &
!             DustName(3) /'DP_09'/,     &
!             DustName(4) /'DP_26'/,     &
!             DustName(5) /'DP_80'/             !,     &
! !           DustName(6) /'DP_240'/,    &
! !           DustName(7) /'DP_720'/,    &
! !           DustName(8) /'DP_2200'/

!---------------------------------------------------------------------
!-- fixed dimensions for soil data
      INTEGER, PARAMETER :: SoilNumb = 5    & ! number of soil properties
&                          ,SoilFrac = 1      ! number of soil fraction

!-- internal switches
      LOGICAL     :: laidaily   ! daily or monthly input of lai     MF
      LOGICAL     :: lvegmin    ! daily or monthly input of lai     MF

!---------------------------------------------------------------------
!-- MUSCAT grid structure
      TYPE dust_subdomain
      REAL(8), POINTER ::   &
&              soilprop (:,:,:,:)   & !soil properties
&             ,lai (:,:,:,:)        & !leafe area index
&             ,vegmin (:,:,:)       & !minimum of vegetation
&             ,alpha (:,:,:)        & !ratio horiz/vertical flux
&             ,c_eff (:,:,:)        & !fraction efficace
&             ,lai_eff (:,:,:,:)    & !effective surface for dust deflation from LAI condition
&             ,w_str (:,:,:)        & !threshold soil moisture w' (Fecan, F. et al., 1999)
&             ,umin2(:,:,:)         &
&             ,d_emis(:,:,:)          !dust emission
      END TYPE dust_subdomain
      TYPE (dust_subdomain), ALLOCATABLE, TARGET :: dust(:)

!---------------------------------------------------------------------
!-- LM grid arrays
      REAL(8), POINTER ::   &
&              SG_soilprop (:,:,:,:),     & ! soil properties
&              SG_lai (:,:,:,:),          & ! leafe area  index
&              SG_vegmin (:,:,:)             ! minimum of vegetation

!---------------------------------------------------------------------
!-- soil class properties
      INTEGER, PARAMETER :: nats   = 45     & ! amount of soil types, originally "nats = 12", "nats = 14"
&                          ,nClass = 196      ! amount of particule classes
      REAL(8) :: srel(nats,Nclass)          & !
&               ,srelV(nats,Nclass)         & !
&               ,su_srelV(nats,nclass)      & !
	&               ,Uth(Nclass)                & ! threshold friction velocity
	&               ,Uth_bod(Nclass)              ! threshold friction velocity

!-- Vars for the Okin aprroach (Matthias Faust)
      REAL(8) :: &
        gap, gap1, gap2, &  ! gap size between plants
        gapheight           ! funktion of gap size / plant height


!-- New Vars for src_dust including one by one if necessary also in run_init
  ! INTEGER ::    &
  !   vegscheme   ! =0 no vegitation; =1 Okin scheme; =2 linear Tegen
  !
  !
  ! LOGICAL :: &
  !   lvegdaily,     &   ! daily or monthly input of lai     MF
  !   lwithz0,       &   ! =false without z0, =true with z0
  !   lwithbiom          ! =false without biomes, =true with biomes

  ! dummy variable for input, allocate new when switch from 2d to 3d
  REAL(8), ALLOCATABLE ::  &
    read_input(:,:,:)       ! (j,i,time)

  ! ! new dust datatype with less in it.
  ! TYPE dust_it
  ! REAL(8), POINTER ::   &
  !   soiltype(:,:),    & ! soil properties
  !   z0  (:,:),          & ! roughness length
  !   biom(:,:),        &
  !   veg (:,:,:),      &          !leafe area index
  !   vegmin(:,:)
  ! END TYPE dust_it
  ! TYPE (dust_it), ALLOCATABLE, TARGET :: dust_ini(:)
  !
  ! TYPE dust_fx
  ! REAL(8), POINTER ::   &
  !   source  (:,:),      & ! preferential dust source
  !   alpha   (:,:),      & ! ratio horiz/vertical flux
  !   feff    (:,:,:),    & ! drag partition
  !   d_emis  (:,:,:)       ! dust emission
  ! END TYPE dust_fx
  ! TYPE (dust_fx), ALLOCATABLE, TARGET :: dust_flux(:)

  INTEGER  ::  &
    dimveg       ! time dimension of vegetation  = 12 if monthly, = nSimuDays if daily

END MODULE mo_dust
