PROGRAM dust_only

  USE mo_dust
  USE offline_org
  USE src_dust, ONLY: organize_dust

  ! USE   src_dust, ONLY: organize_dust  !new dust module !MF
  ! USE   mo_dust
  ! USE   sub_block
  ! USE   partition

  implicit none

  INTEGER        :: &
    ierr
  CHARACTER(120) :: &
    yerr

  ! INTEGER :: &
  !   ncdfID,        & ! id Var for the nc file
  !   timeID,      &
  !   rlonID,      &
  !   rlatID,      &
  !   lonID,       &
  !   latID,       &
  !   DE01ID,      &
  !   DE03ID,      &
  !   DE09ID,      &
  !   DE26ID,      &
  !   DE80ID,      &
  !   DETOTID,     &
  !   DEPM25ID,    &
  !   DEPM10ID,    &
  !   timeDim,     &
  !   rlonDim,     &
  !   rlatDim


  ierr = 0
  nDust = 1

  PRINT*, 'Start Dust emisson model'

  ! CHARACTER(20), ALLOCATABLE, TARGET :: species_name(:)

  ! Read in the INPUT FILE
  print*, 'Read NAMELIST'
  CALL read_namelist(ierr)
  IF (ierr /= 0) THEN
    print*, 'ERROR reading namelist'
    ierr =  100 + ierr
    print*, ierr
    STOP
  END IF

  print*, 'Define Grid'
  CALL def_grid(ierr)
  IF (ierr /= 0) THEN
    print*, 'ERROR def grid'
    ierr =  200 + ierr
    print*, ierr
    STOP
  END IF


  ! read in u wind
  print*, 'Read wind data'
  IF (uconst == 999.0) THEN
    CALL netcdf_in(TRIM(windFile),TRIM(u_var_name),u,lasttstep,.TRUE.,ierr,yerr)
    IF (ierr /= 0) THEN
      print*, 'ERROR netcdf in'
      ierr =  300 + ierr
      print*, ierr, yerr
      STOP
    END IF
  ENDIF

  ! read in v wind
  IF (vconst == 999.0) THEN
    CALL netcdf_in(TRIM(windFile),TRIM(u_var_name),v,lasttstep,.TRUE.,ierr,yerr)
    IF (ierr /= 0) THEN
      print*, 'ERROR netcdf in'
      ierr =  300 + ierr
      print*, ierr, yerr
      STOP
    END IF
  ENDIF
  ! CALL netcdf_in()

  ! ! simulation of muscat Species
  ! ALLOCATE(species_name(5))
  ! species_name = DustName
  ! tracer_name => species_name
  ! nt=5
  !

  ! - init routine
  CALL organize_dust('init',domain)
  IF (ierr /= 0) THEN
    print*, 'ERROR organize_dust:init'
    ierr =  300 + ierr
    print*, ierr
    STOP
  END IF

  ! init of the OUTPUT

  ! def Filename
  ! WRITE (ofilename,"(A6,I6.6,A2,I3.3,A3)") 'part_t',               &
  !        INT(nstart_part_gp*dt_part/60._wp,KIND=iintegers) ,'_p',k,'.nc'


  WRITE (ofilename,'(6A)') TRIM(youtdir),'/dust_emis_',ydate_ini,'-',ydate_end,'.nc'

  CALL netcdf_out('create',TRIM(ofilename),0,ierr)

  ! Time loop
  DO ntstep=firsttstep,lasttstep
    print'(A,I5,A,I5)', '    STEP  ',ntstep,'/',lasttstep-firsttstep

    IF (ntstep > firsttstep) THEN
      CALL organize_dust('calc',domain,flux=dust_flux)
    ENDIF
    IF (ierr /= 0) THEN
      print*, 'ERROR organize_dust:calc'
      ierr =  400 + ierr
      print*, ierr
      STOP
    END IF

    ! accumulate dust emisson g m-2 s-1 -> kg m-2
    IF (laccumulation) THEN
      dust_em_accum = dust_em_accum + dust(1)%d_emis*dt*1.E-3
    ELSE
      dust_em_accum = dust(1)%d_emis*dt*1.E-3
    END IF
    ! print*, maxval(dust(1)%d_emis*1.E-3)

    ! print*, dust(1)%d_emis(:,:,:)

    CALL netcdf_out('appand',TRIM(ofilename),ntstep,ierr)


  END DO

  print*,'done'


END PROGRAM



SUBROUTINE def_grid(ierr)
  USE offline_org
  USE mo_dust

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: &
    ierr

  REAL(8) :: &
    get_hstop


  ierr = 0




  ! init domain in MUSCAT style

  domain%igx0=0
  domain%igy0=0
  domain%igx1=ie_tot
  domain%igy1=je_tot
  domain%ix0=domain%igx0
  domain%iy0=domain%igy0
  domain%ix1=domain%igx1
  domain%iy1=domain%igy1
  domain%ntx=ie_tot
  domain%nty=je_tot
  domain%ngx=ie_tot
  domain%ngy=je_tot
  domain%ib=1
  domain%refine=0

  ALLOCATE(decomp(1))
  decomp = domain

  hstart=0.0
  hstop=get_hstop(ydate_ini,ydate_end)

  StartDate=ydate_ini

  dt = timeinc * 3600.0

  firsttstep = INT((hstart)/timeinc)
  lasttstep  = INT((hstop-hstart)/timeinc)


  ! ALLOCATE(geo(1))
  ! ! ALLOCATE(geo(1)%dxK(je_tot,ie_tot))
  ! ! ALLOCATE(geo(1)%dyK(je_tot,ie_tot))
  ! ALLOCATE(geo(1)%dz(1,je_tot,ie_tot))
  !
  ! ALLOCATE(meteo(1))
  ! ! ALLOCATE(meteo(1)%rho(1,je_tot,ie_tot,ScalCur))
  ! ! ALLOCATE(meteo(1)%QRSur(je_tot,ie_tot,ScalCur))
  ! ALLOCATE(meteo(1)%u(1,je_tot,ie_tot,lasttstep-firsttstep))
  ! ALLOCATE(meteo(1)%v(1,je_tot,ie_tot,lasttstep-firsttstep))

  ALLOCATE(u(je_tot,ie_tot,lasttstep-firsttstep))
  u=0.
  ALLOCATE(v(je_tot,ie_tot,lasttstep-firsttstep))
  v=0.
  ALLOCATE(dz(1,je_tot,ie_tot))
  dz=0


  ALLOCATE(dust_flux(ntz,je_tot,ie_tot,nt))
  dust_flux=0.

  AllOCATE(dust_em_accum(je_tot,ie_tot,nt))
  dust_em_accum=0


  IF (leveltype == '10m') dz=10

  ! IF (leveltype == 'mlv') THEN
  !
  ! ENDIF

  IF (uconst  /= 999.0) u  = uconst
  IF (vconst  /= 999.0) v  = vconst
  IF (dzconst /= 999.0) dz = dzconst

END SUBROUTINE def_grid


FUNCTION get_hstop(date_start,date_end)

  IMPLICIT NONE

  REAL(8) :: &
    get_hstop

  CHARACTER(*) :: &
    date_start,   &
    date_end

  INTEGER :: &
    syear,   &
    smon,    &
    sday,    &
    shour,   &
    eyear,   &
    emon,    &
    eday,    &
    ehour,   &
    days_this_month

  LOGICAL :: &
    lyear,   &
    lmon,    &
    lday,    &
    lhour

  ! read date into integers vars
  READ(date_start(1:4),*)  syear
  READ(date_start(5:6),*)  smon
  READ(date_start(7:8),*)  sday
  READ(date_start(9:10),*) shour

  READ(date_end(1:4),*)  eyear
  READ(date_end(5:6),*)  emon
  READ(date_end(7:8),*)  eday
  READ(date_end(9:10),*) ehour

  get_hstop = 0.0

  lyear=.FALSE.
  lmon =.FALSE.
  lday =.FALSE.
  lhour=.FALSE.

  IF (eyear == syear) lyear=.TRUE.
  IF (emon  == smon ) lmon =.TRUE.
  IF (eday  == sday ) lday =.TRUE.
  IF (ehour == shour) lhour=.TRUE.

  DO WHILE (.NOT. lyear .OR. .NOT. lmon .OR. .NOT. lday .OR. .NOT. lhour)

    ! add one hour
    shour = shour + 1
    get_hstop = get_hstop + 1.0

    ! add one day
    IF (shour >= 24) THEN
      shour = 0
      sday = sday + 1
    END IF

    ! how many days has this month
    IF (smon == 1 .or. smon == 3 .or. smon == 5 .or. smon == 7 .or. &
        smon == 8 .or. smon == 10 .or. smon == 12) THEN
      days_this_month = 31
    ELSE
      days_this_month = 30
    END IF

    ! only 28 days in february
    IF (smon == 2) THEN
      days_this_month = 28
      ! are we in a leap-year
      IF (MOD(syear,4) == 0) THEN
        days_this_month = 29
        ! special rule in the Gregorian calendar, sometimes the leap-year is skipped e.q. 1900
        IF (MOD(syear,100) == 0 .AND. MOD(syear,400) /= 0) days_this_month = 28
      END IF
    END IF

    ! add one month
    IF (sday > days_this_month) THEN
      sday = 1
      smon = smon + 1
    END IF

    ! add one year
    IF (smon > 12) THEN
      smon  = 1
      syear = syear + 1
    END IF

    ! ceck if end date is reached
    lyear=.FALSE.
    lmon =.FALSE.
    lday =.FALSE.
    lhour=.FALSE.
    IF (eyear == syear) lyear=.TRUE.
    IF (emon  == smon ) lmon =.TRUE.
    IF (eday  == sday ) lday =.TRUE.
    IF (ehour == shour) lhour=.TRUE.

  END DO

END FUNCTION get_hstop



SUBROUTINE read_namelist(ierr)

  USE   mo_dust

  implicit none

  INTEGER, INTENT(INOUT) :: ierr


  NAMELIST /SURFMODEL/ &
    ie_tot,            & ! number of lon points
    je_tot,            & ! number of lat points
    startlon_tot,           & ! lon of low left corner
    startlat_tot,           & ! lat of low left corner
    timeinc,           & ! increment of time in hours
    pollon,            &
    pollat,            &
    dlon,              &
    dlat,              &
    windFile,          &
    leveltype,         &
    ydate_ini,         & ! Date when model Start YYYYMMDDHH
    ydate_end,         & ! Date when model END   YYYYMMDDHH
    youtdir,           & ! directory of the output file
    lwithz0,           & ! =false without z0, =true with z0
    lwithbiom,         & ! =false without biomes, =true with biomes
    dust_scheme,       & ! 1=Tegen02
    veg_scheme,        & ! =0 no vegitation; =1 Okin scheme; =2 linear Tegen
    moist_scheme,      &
    psrcType,          & ! Flag for type of potential dust source ! 0 : psrc, 1 : msgsrc, 2 : acDust
    soiltypeFile,      & ! Filename of Soil Type Data
    psrcFile,          & ! Filename of preferential Dust Sources
    cultFile,          & ! Filename of Cultivation Class
    vegmonFile,        & ! Filename of monthly vegitation cover
    vegdayFile,        & ! Filename of daily vegetation cover
    vegminFile,        & ! Filename of min vegetation cover MF
    z0File,            & ! Filename of Roughness Length
    biomeFile,         & ! Filename of Vegetation Cover/Type Data
    moistFile,         &
    uconst,            &
    vconst,            &
    dzconst,           &
    u_var_name,        & ! Name of u variable in Input file
    v_var_name,        & ! Name of v variable in Input file
    laccumulation

  ! Defaults
  soiltypeFile  = 'without'   ! Soil Type Data
  psrcType      = 1           ! 0 : psrc, 1 : msgsrc, 2 : acDust  MF
  dust_scheme   = 1           ! 1=Tegen02
  veg_scheme    = 0           ! 0=no vegetation, 1=okin scheme, 2=linear (Tegen02)
  psrcFile      = 'without'   ! Potential Sources Data
  cultFile      = 'without'   ! Landuse Data
  z0File        = 'without'   ! Z0 Data
  biomeFile     = 'without'   ! Vegetation Classes
  moistFile     = 'without'
  vegmonFile    = 'without'   ! Leafe Area Index Data
  vegdayFile    = 'without'   ! Leafe Area Index Data   MF
  vegminFile    = 'without'   ! Leafe Area Index Data   MF
  lwithz0       = .FALSE.
  lwithbiom     = .FALSE.
  uconst        = 999.0
  vconst        = 999.0
  dzconst       = 999.0
  u_var_name    = 'u10'
  v_var_name    = 'v10'
  laccumulation = .FALSE.



  OPEN (UNIT=10, FILE='INPUT',FORM=  'FORMATTED', STATUS='OLD',IOSTAT=ierr)
  IF (ierr /= 0) THEN
    PRINT*, 'ERROR open INPUT file'
    ierr = 1
    RETURN
  ENDIF

  READ (10, SURFMODEL,IOSTAT=ierr)
  IF (ierr /= 0) THEN
    PRINT*, 'ERROR reading namelist'
    ierr = 2
    RETURN
  ENDIF

  CLOSE (10)




END SUBROUTINE read_namelist

SUBROUTINE netcdf_in(infile,varname,outvar,ntimes,timecheck,ierror,yerrmsg)

  ! Modules
  USE mo_dust
  USE offline_org
  USE netcdf
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
    ntimes         ! number of possible time steps in the model run
  LOGICAL,        INTENT(IN)    :: &
    timecheck     ! check if time fits model run

  REAL(8),        INTENT(OUT) :: &
    outvar(je_tot,ie_tot,ntimes)        ! Output var

  INTEGER,        INTENT(OUT)   :: &
    ierror

  CHARACTER(*),   INTENT(OUT)   :: &
    yerrmsg       ! error message
  ! CHARACTER(120) :: &
  !   yerrmsg       ! error message


  ! local Vars
  INTEGER  :: &
    i,j,t,    &  ! loop
    istart,   &  ! index of the start date in the var file
    istat,    &  ! local error code
    incID,     &  ! id of nc files
    itimeID,   &  ! id of time var
    ivarID,    &  ! id of the  var
    idimID,    &  ! id of a dimension
    idimlen       ! length of the above dimension

  REAL(8)  ::  &
    idate,     &  ! ydate read into float
    dayfrac,   &  ! fraction of a day
    var_scale, &    ! index of the start date in the var file
    var_offset

  REAL(8), ALLOCATABLE :: &
    times(:)

  REAL, ALLOCATABLE :: &
    var_read(:,:,:)



  ! start subroutine
  ! ---------------------------------------------------------

  ! Print short status
  PRINT*,'read_nc ', infile

  ! open the file
  istat = nf90_open(infile, nf90_nowrite, incid)
  IF (istat /= nf90_noerr) THEN
    ierror  = 10001
    yerrmsg = TRIM(nf90_strerror(istat))
    RETURN
  ENDIF

  ! get id of lat dimension
  istat = nf90_inq_dimid(incID, 'y', idimID)
  IF (istat /= nf90_noerr) THEN
    ierror  = 10002
    yerrmsg = TRIM(nf90_strerror(istat))
    RETURN
  ENDIF

  ! read the length of lat dimension
  istat = nf90_inquire_dimension(incID, idimID, len = idimlen)
  IF (istat /= nf90_noerr) THEN
    ierror  = 10003
    yerrmsg = TRIM(nf90_strerror(istat))
    RETURN
  ENDIF

  ! Check the Size
  IF (idimlen /= je_tot) THEN
    ierror = 10004
    yerrmsg = 'Error reading '//TRIM(varname)//' file: wrong lat dimension'
    RETURN
  END IF

   ! get id of lon dimension
   istat = nf90_inq_dimid(incID, 'x', idimID)
   IF (istat /= nf90_noerr) THEN
     ierror  = 10005
     yerrmsg = TRIM(nf90_strerror(istat))
     RETURN
   ENDIF

   ! read the length of lon dimension
   istat = nf90_inquire_dimension(incID, idimID, len = idimlen)
   IF (istat /= nf90_noerr) THEN
     ierror  = 10006
     yerrmsg = TRIM(nf90_strerror(istat))
     RETURN
   ENDIF

   ! Check the Size
   IF (idimlen /= ie_tot) THEN
     ierror = 10007
     yerrmsg = 'Error reading '//TRIM(varname)//' file: wrong lon dimension'
     RETURN
   END IF

  ! get id of time dimension
  istat = nf90_inq_dimid(incID, 'time', idimID)
  IF (istat /= nf90_noerr) THEN
    ierror  = 10008
    yerrmsg = TRIM(nf90_strerror(istat))
    RETURN
  ENDIF

  ! read the length of the time dimension
  istat = nf90_inquire_dimension(incID, idimID, len = idimlen)
  IF (istat /= nf90_noerr) THEN
    ierror  = 10009
    yerrmsg = TRIM(nf90_strerror(istat))
    RETURN
  ENDIF

  ! allocate the var ( times ) that hold the available dates in the nc file
  istat = 0
  ALLOCATE (times(idimlen), STAT=istat)
  IF (istat /= 0) THEN
    ierror = 10010
    yerrmsg = 'allocation of times failed'
    RETURN
  ENDIF

  ! get the id of the time var
  istat = nf90_inq_varid(incID, 'time', itimeID)
  IF (istat /= nf90_noerr) THEN
    ierror  = 10005
    yerrmsg = TRIM(nf90_strerror(istat))
    RETURN
  ENDIF

  ! get the dates in the nc file
  istat = nf90_get_var(incid, itimeID, times)
  IF (istat /= nf90_noerr) THEN
    ierror  = 10011
    yerrmsg = TRIM(nf90_strerror(istat))
    RETURN
  ENDIF

  ! transform date string to an integer
  READ(ydate_ini(1:8), *) idate
  READ(ydate_ini(9:10),*) dayfrac
  dayfrac = dayfrac/24.0
  idate=idate+dayfrac

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

  ! Allocate var to read
  istat = 0
  ALLOCATE (var_read(ie_tot,je_tot,ntimes), STAT=istat)
  IF (istat /= 0) THEN
    ierror = 10014
    yerrmsg = 'allocation of var_read failed'
    RETURN
  ENDIF


  ! get the id of the var
  istat = nf90_inq_varid(incID, varname , ivarID)
  IF (istat /= nf90_noerr) THEN
    ierror  = 10015
    yerrmsg = TRIM(nf90_strerror(istat))
    RETURN
  ENDIF

  ! get value of the scaling factor
  istat = nf90_get_att(incid, ivarID, 'scale_factor',var_scale)
  IF (istat /= nf90_noerr) THEN
    ! ierror  = 10016
    var_scale = 1
    ! yerrmsg = TRIM(nf90_strerror(istat))
    !RETURN
  ENDIF

  ! get value of the offset
  istat = nf90_get_att(incid, ivarID, 'add_offset',var_offset)
  IF (istat /= nf90_noerr) THEN
    ! ierror  = 10016
    var_offset = 0
    ! yerrmsg = TRIM(nf90_strerror(istat))
    !RETURN
  ENDIF

  ! get the var
  istat = nf90_get_var(incid, ivarID, var_read,start= (/1,1,istart/),count=(/ie_tot,je_tot,ntimes/))
  IF (istat /= nf90_noerr) THEN
    ierror  = 10017
    yerrmsg = TRIM(nf90_strerror(istat))
    RETURN
  ENDIF


  ! AllOCATE(outvar(je_tot,ie_tot,ntimes))

  ! wirte the var into outvar
  DO i=1,ie_tot
    DO j=1,je_tot
      DO t=1, ntimes
        ! digital values to physical values with the scaling factor
        ! u(j,i,t)=var_read(i,j,t)*var_scale+var_offset
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
END SUBROUTINE netcdf_in

SUBROUTINE netcdf_out(status,Filename,step,ierr)!,FileID,Var,ierr)

    USE mo_dust
    USE offline_org
    USE netcdf

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: &
      status                         ! action that will be done 'create', 'appand', 'close'

    CHARACTER(*), INTENT(IN) :: &    !OPTIONAL,
      Filename                       ! name of the nc File

    INTEGER, INTENT(IN) :: &
      step

    INTEGER, INTENT(OUT) :: &
      ierr

    INTEGER ::     &
      i,j,         &
      istat,       & ! netcdf status variable
      ! ncdfID,        & ! id Var for the nc file
      ! timeID,      &
      ! rlonID,      &
      ! rlatID,      &
      ! lonID,       &
      ! latID,       &
      ! DE01ID,      &
      ! DE03ID,      &
      ! DE09ID,      &
      ! DE26ID,      &
      ! DE80ID,      &
      ! DETOTID,      &
      ! DEPM25ID,      &
      ! DEPM10ID,      &
      ! timeDim,     &
      ! rlonDim,     &
      ! rlatDim,     &
      iztime_tmp!,  & ! tmp variable to hold time information

    REAL(8) :: &
      pi,      &
      geolon,  &
      geolat,  &
      rlon,    &
      rlat

    CHARACTER (LEN=40) :: ydate
    CHARACTER(120) :: yerrmsg



    ! Sec 1 create the nc file
    IF (status == 'create') THEN

      istat=nf90_create(Filename, NF90_SHARE, ncdfID)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10208
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF


      ! Define global attributes
      READ(ydate_ini(1:4),'(i4)') iztime_tmp
      istat = nf90_put_att(ncdfID, nf90_GLOBAL, 'ref_year', iztime_tmp)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10208
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF
      READ(ydate_ini(5:6),'(i2)') iztime_tmp
      istat = nf90_put_att(ncdfID, nf90_GLOBAL, 'ref_month', iztime_tmp)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10209
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF
      READ(ydate_ini(7:8),'(i2)') iztime_tmp
      istat = nf90_put_att(ncdfID, nf90_GLOBAL, 'ref_day', iztime_tmp)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10210
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF
      READ(ydate_ini(9:10),'(i2)') iztime_tmp
      istat = nf90_put_att(ncdfID, nf90_GLOBAL, 'ref_hour', iztime_tmp)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10211
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF
      ! READ(ydate_ini(11:12),'(i2)') iztime_tmp
      ! istat = nf90_put_att(ncdfID, nf90_GLOBAL, 'ref_min', iztime_tmp)
      ! IF (istat /= nf90_noerr) THEN
      !   ierr  = 10212
      !   yerrmsg = TRIM(nf90_strerror(istat))
      !   RETURN
      ! ENDIF
      ! READ(ydate_ini(13:14),'(i2)') iztime_tmp
      ! istat = nf90_put_att(ncdfID, nf90_GLOBAL, 'ref_sec', iztime_tmp)
      ! IF (istat /= nf90_noerr) THEN
      !   ierr  = 10213
      !   yerrmsg = TRIM(nf90_strerror(istat))
      !   RETURN
      ! ENDIF
      istat = nf90_put_att(ncdfID, nf90_GLOBAL, 'pollon', pollon)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10214
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF
      istat = nf90_put_att(ncdfID, nf90_GLOBAL, 'pollat', pollat)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10215
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF
      ! istat = nf90_put_att(ncdfID, nf90_GLOBAL, 'output_time_step_in_sec', &
      !                       NINT(hout_con*3600,iintegers))
      ! IF (istat /= nf90_noerr) THEN
      !   ierr  = 10216
      !   yerrmsg = TRIM(nf90_strerror(istat))
      !   RETURN
      ! ENDIF



      ! Define dimensions

      ! Time
      istat = nf90_def_dim(ncdfID, 'time', NF90_UNLIMITED, timeDim)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10217
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF

      ! rlon
      istat = nf90_def_dim(ncdfID, 'rlon', ie_tot, rlonDim)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10218
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF

      ! rlat
      istat = nf90_def_dim(ncdfID, 'rlat', je_tot, rlatDim)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10219
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF


      !----------------------------------------------------------------------
      !- Section 2.2: Define variables and their attributes
      !----------------------------------------------------------------------

      ! Time
      istat = nf90_def_var(ncdfID, 'time', nf90_FLOAT, (/ timeDim /), timeID)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10221
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, timeID, "standard_name", "time")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10222
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, timeID, "long_name", "time")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10223
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, timeID, "calendar", "gregorian")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10223
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      ydate = 'seconds since '// ydate_ini(1:4)//'-'//ydate_ini(5:6)//'-'  &
                 //ydate_ini(7:8)//' '//ydate_ini(9:10)//':00:00'!//':'//ydate_ini(11:12)//':'&
                 ! //ydate_ini(13:14)
      istat=nf90_put_att(ncdfID, timeID, "units", ydate)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10224
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF

      ! rlon
      istat = nf90_def_var(ncdfID, 'rlon', nf90_FLOAT, (/rlonDim/), rlonID )
      IF (istat /= nf90_noerr) THEN
        ierr  = 10225
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, rlonID, "standard_name", "grid_longitude")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10226
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, rlonID, "long_name", "rotated longitude")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10227
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, rlonID, "units", "degrees")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10228
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF

      ! rlat
      istat = nf90_def_var(ncdfID, 'rlat', nf90_FLOAT, (/rlatDim/), rlatID )
      IF (istat /= nf90_noerr) THEN
        ierr  = 10229
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, rlatID, "standard_name", "grid_latitude")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10230
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, rlatID, "long_name", "rotated latitude")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10231
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, rlatID, "units", "degrees")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10232
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF

      ! lon
      istat = nf90_def_var(ncdfID, "lon", nf90_FLOAT,(/rlonDim, rlatDim/), lonID)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10237
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, lonID, "standard_name", "longitude")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10238
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, lonID, "long_name", "longitude")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10239
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, lonID, "units", "degrees_east")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10240
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF


      ! lat
      istat = nf90_def_var(ncdfID, "lat", nf90_FLOAT,(/rlonDim, rlatDim/), latID)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10237
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, latID, "standard_name", "latitude")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10238
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, latID, "long_name", "latitude")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10239
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, latID, "units", "degrees_north")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10240
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF

      ! DE_01
      istat = nf90_def_var(ncdfID, "DE_01", nf90_FLOAT,(/rlonDim, rlatDim, timeDim/), DE01ID)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10237
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE01ID, "standard_name", "DE_01")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10238
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE01ID, "long_name", "Dust emisson < 1 µm")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10239
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE01ID, "units", "kg/m-2")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10240
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF

      ! DE_03
      istat = nf90_def_var(ncdfID, "DE_03", nf90_FLOAT,(/rlonDim, rlatDim, timeDim/), DE03ID)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10237
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE03ID, "standard_name", "DE_03")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10238
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE03ID, "long_name", "Dust emisson < 3 µm")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10239
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE03ID, "units", "kg/m-2")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10240
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF

      ! DE_09
      istat = nf90_def_var(ncdfID, "DE_09", nf90_FLOAT,(/rlonDim, rlatDim, timeDim/), DE09ID)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10237
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE09ID, "standard_name", "DE_09")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10238
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE09ID, "long_name", "Dust emisson < 9 µm")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10239
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE09ID, "units", "kg/m-2")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10240
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF

      ! DE_26
      istat = nf90_def_var(ncdfID, "DE_26", nf90_FLOAT,(/rlonDim, rlatDim, timeDim/), DE26ID)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10237
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE26ID, "standard_name", "DE_26")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10238
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE26ID, "long_name", "Dust emisson < 26 µm")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10239
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE26ID, "units", "kg/m-2")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10240
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF

      ! DE_80
      istat = nf90_def_var(ncdfID, "DE_80", nf90_FLOAT,(/rlonDim, rlatDim, timeDim/), DE80ID)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10237
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE80ID, "standard_name", "DE_80")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10238
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE80ID, "long_name", "Dust emisson < 80 µm")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10239
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DE80ID, "units", "kg/m-2")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10240
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF

      ! DE_TOT
      istat = nf90_def_var(ncdfID, "DE_TOT", nf90_FLOAT,(/rlonDim, rlatDim, timeDim/), DETOTID)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10237
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DETOTID, "standard_name", "DE_TOT")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10238
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DETOTID, "long_name", "Total dust emisson")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10239
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DETOTID, "units", "kg/m-2")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10240
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF

      ! DE_PM_2.5
      istat = nf90_def_var(ncdfID, "DE_PM_2.5", nf90_FLOAT,(/rlonDim, rlatDim, timeDim/), DEPM25ID)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10237
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DEPM25ID, "standard_name", "DE_PM_2.5")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10238
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DEPM25ID, "long_name", "Dust emisson < 2.5 µm")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10239
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DEPM25ID, "units", "kg/m-2")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10240
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF

      ! DE_PM_10
      istat = nf90_def_var(ncdfID, "DE_PM_10", nf90_FLOAT,(/rlonDim, rlatDim, timeDim/), DEPM10ID)
      IF (istat /= nf90_noerr) THEN
        ierr  = 10237
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DEPM10ID, "standard_name", "DE_PM_10")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10238
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DEPM10ID, "long_name", "Dust emisson < 10 µm")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10239
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF
      istat=nf90_put_att(ncdfID, DEPM10ID, "units", "kg/m-2")
      IF (istat /= nf90_noerr) THEN
        ierr  = 10240
        yerrmsg = TRIM(nf90_strerror(ierr))
        RETURN
      ENDIF



      istat=nf90_enddef(ncdfID)


      ! rlon
      DO i = 0, ie_tot-1
        rlon = (i*dlon)+startlon_tot
        istat = NF90_put_var(ncdfID,rlonID,rlon, start=(/i+1/))
      ENDDO
      ! rlat
      DO j = 0, je_tot-1
        rlat = (j*dlat)+startlat_tot
        istat = NF90_put_var(ncdfID,rlatID,rlat, start=(/j+1/))
      ENDDO
      ! lon/lat
      pi = 4. * atan(1.)
      DO i = 0, ie_tot-1
        DO j = 0, je_tot-1
          rlon = (i*dlon)+startlon_tot
          rlat = (j*dlat)+startlat_tot

          geolon=180./pi * atan((cos(pi/180.*rlat)*sin(pi/180.*rlon))/  &
                                (sin(pi/180.*pollat)*cos(pi/180.*rlat)* &
                                 cos(pi/180.*rlon)-sin(pi/180.*rlat)*   &
                                 cos(pi/180.*pollat))) + pollon + 180.

          geolat=180./pi * asin(sin(pi/180.*rlat)*sin(pi/180.*pollat) + &
                                cos(pi/180.*rlat)*cos(pi/180.*rlon)*cos(pi/180.*pollat))

          istat = NF90_PUT_VAR(ncdfID,lonID, geolon, start=(/i+1,j+1/) )
          IF (istat /= nf90_noerr) THEN
            ierr  = 10220
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          ! lat
          istat = NF90_PUT_VAR(ncdfID,latID, geolat, start=(/i+1,j+1/) )
          IF (istat /= nf90_noerr) THEN
            ierr  = 10220
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
        END DO
      END DO

      istat=nf90_close(ncdfID)

    ELSEIF (status == 'appand') THEN

      istat=nf90_open(Filename, NF90_WRITE, ncdfID)

      ! write time
      istat = NF90_put_var(ncdfID,timeID,step*dt,start=(/step+1/))
      IF (istat /= nf90_noerr) THEN
        ierr  = 10220
        yerrmsg = TRIM(nf90_strerror(istat))
        print*, yerrmsg
        RETURN
      ENDIF

      ! write DE_01
      istat = nf90_put_var(ncdfID,DE01ID, transpose(dust_em_accum(:,:,1)), start=(/1,1,step+1/) )
      IF (istat /= nf90_noerr) THEN
        ierr  = 10220
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF

      ! write DE_03
      istat = nf90_put_var(ncdfID,DE03ID, transpose(dust_em_accum(:,:,2)), start=(/1,1,step+1/) )
      IF (istat /= nf90_noerr) THEN
        ierr  = 10220
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF

      ! write DE_09
      istat = nf90_put_var(ncdfID,DE09ID, transpose(dust_em_accum(:,:,3)), start=(/1,1,step+1/) )
      IF (istat /= nf90_noerr) THEN
        ierr  = 10220
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF

      ! write DE_26
      istat = nf90_put_var(ncdfID,DE26ID, transpose(dust_em_accum(:,:,4)), start=(/1,1,step+1/) )
      IF (istat /= nf90_noerr) THEN
        ierr  = 10220
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF

      ! write DE_80
      istat = nf90_put_var(ncdfID,DE80ID, transpose(dust_em_accum(:,:,5)), start=(/1,1,step+1/) )
      IF (istat /= nf90_noerr) THEN
        ierr  = 10220
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF


      ! write DE_TOT
      istat = nf90_put_var(ncdfID,DETOTID, transpose(sum(dust_em_accum(:,:,:5), dim=3 )), start=(/1,1,step+1/) )
      IF (istat /= nf90_noerr) THEN
        ierr  = 10220
        yerrmsg = TRIM(nf90_strerror(istat))
        RETURN
      ENDIF


      istat=nf90_close(ncdfID)
    ENDIF





END SUBROUTINE netcdf_out
