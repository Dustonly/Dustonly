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
 USE src_aerosol, ONLY: DustMod ! Flag for Setting Soil Data
!======================================================================================
! Description:
! ------------
!  This module contains definitions of the provisional results required for the
!  dust-flux calculation in sub_dust_emis.f90.
!======================================================================================
! 
!---------------------------------------------------------------------
!-- control parameters
      INTEGER :: nDust          ! Flag for Dust Calculations

!-- Files with soil data
      CHARACTER(120) ::   &
&           soiltypeFile       & ! Filename of Soil Type Data
&          ,psrcFile           & ! Filename of preferential Dust Sources
&          ,cultFile           & ! Filename of Cultivation Class
&          ,laimonFile         & ! Filename of monthly LAI
&          ,z0File             & ! Filename of Roughness Length
&          ,biomeFile            ! Filename of Vegetation Cover/Type Data

!---------------------------------------------------------------------
!-- description of dust particles
      INTEGER, PARAMETER :: DustBins = 5  !8  ! number of dust particle fractions
      INTEGER :: DustInd(DustBins)            ! indices of dust particles

      CHARACTER(20) :: DustName(DustBins)
      DATA  DustName(1) /'DP_01'/,     &
            DustName(2) /'DP_03'/,     &
            DustName(3) /'DP_09'/,     &
            DustName(4) /'DP_26'/,     &
            DustName(5) /'DP_80'/             !,     &
!           DustName(6) /'DP_240'/,    &
!           DustName(7) /'DP_720'/,    &
!           DustName(8) /'DP_2200'/

!---------------------------------------------------------------------
!-- fixed dimensions for soil data
      INTEGER, PARAMETER :: SoilNumb = 5    & ! number of soil properties
&                          ,SoilFrac = 1      ! number of soil fraction

!---------------------------------------------------------------------
!-- MUSCAT grid structure
      TYPE dust_subdomain
      REAL(8), POINTER ::   &
&              soilprop (:,:,:,:)   & !soil properties
&             ,lai (:,:,:,:)        & !leafe area index 
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
&              SG_soilprop (:,:,:,:)     & ! soil properties
&             ,SG_lai (:,:,:,:)            ! leafe area  index 

!---------------------------------------------------------------------
!-- soil class properties
      INTEGER, PARAMETER :: nats   = 14     & ! amount of soil types, originally "nats = 12"
&                          ,nClass = 196      ! amount of particule classes
      REAL(8) :: srel(nats,Nclass)          & ! 
&               ,srelV(nats,Nclass)         & ! 
&               ,su_srelV(nats,nclass)      & !
	&               ,Uth(Nclass)                & ! threshold friction velocity
	&               ,Uth_bod(Nclass)              ! threshold friction velocity

   END MODULE mo_dust
