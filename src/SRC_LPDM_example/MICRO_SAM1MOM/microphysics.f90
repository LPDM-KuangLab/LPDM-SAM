module microphysics

! module for original SAM bulk microphysics
! Marat Khairoutdinov, 2006

use grid, only: nx,ny,nzm,nz, dimx1_s,dimx2_s,dimy1_s,dimy2_s, rank ! subdomain grid information 
use params, only: doprecip, docloud
use micro_params
implicit none

!----------------------------------------------------------------------
!!! required definitions:

integer, parameter :: nmicro_fields = 2   ! total number of prognostic water vars

!!! microphysics prognostic variables are storred in this array:

real micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields)

integer, parameter :: flag_wmass(nmicro_fields) = (/1,1/)
integer, parameter :: index_water_vapor = 1 ! index for variable that has water vapor
integer, parameter :: index_cloud_ice = 1   ! index for cloud ice (sedimentation)
integer, parameter :: flag_precip(nmicro_fields) = (/0,1/)

! both variables correspond to mass, not number
integer, parameter :: flag_number(nmicro_fields) = (/0,0/)

! SAM1MOM 3D microphysical fields are output by default.
integer, parameter :: flag_micro3Dout(nmicro_fields) = (/0,0/)

real fluxbmk (nx, ny, 1:nmicro_fields) ! surface flux of tracers
real fluxtmk (nx, ny, 1:nmicro_fields) ! top boundary flux of tracers

!!! these arrays are needed for output statistics:

real mkwle(nz,1:nmicro_fields)  ! resolved vertical flux
real mkwsb(nz,1:nmicro_fields)  ! SGS vertical flux
real mkadv(nz,1:nmicro_fields)  ! tendency due to vertical advection
real mklsadv(nz,1:nmicro_fields)  ! tendency due to large-scale vertical advection
real mkdiff(nz,1:nmicro_fields)  ! tendency due to vertical diffusion
real mk0(nzm,1:nmicro_fields) !bloss: placeholder.  Could hold mean profiles
real mk_ref(nzm,1:nmicro_fields) !bloss: placeholder.  Could hold reference profiles

logical :: is_water_vapor(1:nmicro_fields)!bloss: placeholder

!======================================================================
! UW ADDITIONS

! number of fields output from micro_write_fields2D, micro_write_fields3D
integer :: nfields2D_micro=0 
integer :: nfields3D_micro=0 

!bloss: arrays with names/units for microphysical outputs in statistics.
character*3, dimension(nmicro_fields) :: mkname
character*80, dimension(nmicro_fields) :: mklongname
character*10, dimension(nmicro_fields) :: mkunits
real, dimension(nmicro_fields) :: mkoutputscale

!bloss: dummy arrays for effective radius and other variables
!    useful in computing cloud radiative properties, mainly
!    used by M2005 microphysics but included here so that
!    RRTM will compile against both schemes
real, allocatable, dimension(:,:,:) :: reffc, reffi, reffs, reffr, &
     CloudLiquidMassMixingRatio, CloudLiquidGammaExponent, CloudLiquidLambda, &
     CloudIceMassMixingRatio, SnowMassMixingRatio

! Flags related to M2005 cloud optics routines, included here
!   so that RRTM will compile against this microphysics as well.
logical :: dosnow_radiatively_active = .false.
logical :: dorrtm_cloud_optics_from_effrad_LegacyOption = .true.

! END UW ADDITIONS
!======================================================================

!------------------------------------------------------------------
! Optional (internal) definitions)

! make aliases for prognostic variables:
! note that the aliases should be local to microphysics

real q(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! total nonprecipitating water
real qp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! total precipitating water
equivalence (q(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,1))
equivalence (qp(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,2))

real qn(nx,ny,nzm)  ! cloud condensate (liquid + ice)

! variables accumulating precipitation and sedimentation tendencies for use in mse.f90
logical :: isallocatedMICRO = .false.
real, allocatable, dimension(:,:), public :: prec_accum, prec_ice_accum
real, allocatable, dimension(:,:,:), public :: qtot_sed, qice_sed

real qpsrc(nz)  ! source of precipitation microphysical processes
real qpfall(nz) ! source of precipitating water due to fall out in a given level
real qpevp(nz)  ! sink of precipitating water due to evaporation

real vrain, vsnow, vgrau, crain, csnow, cgrau  ! precomputed coefs for precip terminal velocity

CONTAINS

! required microphysics subroutines and function:
!----------------------------------------------------------------------
function micro_scheme_name()
  character(len=32) :: micro_scheme_name
  ! Return the scheme name, normally the same as the directory name with leading "MICRO_" removed  
  micro_scheme_name = "sam1mom" 
end function   micro_scheme_name
!----------------------------------------------------------------------
!!! Read microphysics options from prm file

subroutine micro_setparm()
  ! no user-definable options in SAM1MOM microphysics.
  integer :: ierr

  !======================================================================
  ! UW ADDITION
  if(.not.isallocatedMICRO) then
     ! allocate microphysical variables
     allocate(qtot_sed(nx,ny,nzm), qice_sed(nx,ny,nzm), & ! for budgets in mse.f90
          prec_accum(nx,ny), prec_ice_accum(nx,ny), & ! for budgets in mse.f90
          STAT=ierr)
     if(ierr.ne.0) then
        write(*,*) 'Failed to allocate microphysical arrays on proc ', rank
        call task_abort()
     else
        isallocatedMICRO = .true.
     end if
   end if
  ! END UW ADDITION
  !======================================================================

end subroutine micro_setparm

!----------------------------------------------------------------------
!!! Initialize microphysics:


subroutine micro_init()

  use grid, only: nrestart
  use vars, only: q0
  use params, only: dosmoke
  integer k

  a_bg = 1./(tbgmax-tbgmin)
  a_pr = 1./(tprmax-tprmin)
  a_gr = 1./(tgrmax-tgrmin)

  if(doprecip) call precip_init() 

  if(nrestart.eq.0) then

     micro_field = 0.
     do k=1,nzm
      q(:,:,k) = q0(k)
     end do
     qn = 0.
     fluxbmk = 0.
     fluxtmk = 0.

     if(docloud) then
       call cloud()
       call micro_diagnose()
     end if
     if(dosmoke) then
       call micro_diagnose()
     end if

  end if

  mkwle = 0.
  mkwsb = 0.
  mkadv = 0.
  mkdiff = 0.
  mklsadv = 0.

  qpsrc = 0.
  qpevp = 0.

  mkname(1) = 'QT'
  mklongname(1) = 'TOTAL WATER (VAPOR + CONDENSATE)'
  mkunits(1) = 'g/kg'
  mkoutputscale(1) = 1.e3

  mkname(2) = 'QP'
  mklongname(2) = 'PRECIPITATING WATER'
  mkunits(2) = 'g/kg'
  mkoutputscale(2) = 1.e3

  !bloss: added variables in upperbound.f90
  is_water_vapor = .false.
  is_water_vapor(1) = .true.
  mk0(:,:) = 0.

end subroutine micro_init

!----------------------------------------------------------------------
!!! fill-in surface and top boundary fluxes:
!
subroutine micro_flux()

  use vars, only: fluxbq, fluxtq

  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
  fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)

end subroutine micro_flux

!----------------------------------------------------------------------
!!! compute local microphysics processes (bayond advection and SGS diffusion):
!
subroutine micro_proc()

   use grid, only: nstep,dt,icycle
   use params, only: dosmoke

   call t_startf ('micro_proc')


   ! Update bulk coefficient
   if(doprecip.and.icycle.eq.1) call precip_init() 

   if(docloud) then
     call cloud()
     if(doprecip) call precip_proc()
     call micro_diagnose()
   end if
   if(dosmoke) then
     call micro_diagnose()
   end if

   call t_stopf ('micro_proc')

end subroutine micro_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine micro_diagnose()
 
   use vars

   real omn, omp
   integer i,j,k

   do k=1,nzm
    do j=1,ny
     do i=1,nx
       qv(i,j,k) = q(i,j,k) - qn(i,j,k)
       omn = max(0.,min(1.,(tabs(i,j,k)-tbgmin)*a_bg))
       qcl(i,j,k) = qn(i,j,k)*omn
       qci(i,j,k) = qn(i,j,k)*(1.-omn)
       omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
       qpl(i,j,k) = qp(i,j,k)*omp
       qpi(i,j,k) = qp(i,j,k)*(1.-omp)
     end do
    end do
   end do
       


end subroutine micro_diagnose

!----------------------------------------------------------------------
!!! function to compute terminal velocity for precipitating variables:
! In this particular case there is only one precipitating variable.

real function term_vel_qp(i,j,k,ind)
  
  use vars
  integer, intent(in) :: i,j,k,ind
  real wmax, omp, omg, qrr, qss, qgg

  term_vel_qp = 0.
  if(qp(i,j,k).gt.qp_threshold) then
    omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
    if(omp.eq.1.) then
       term_vel_qp = vrain*(rho(k)*qp(i,j,k))**crain
    elseif(omp.eq.0.) then
       omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
       qgg=omg*qp(i,j,k)
       qss=qp(i,j,k)-qgg
       term_vel_qp = (omg*vgrau*(rho(k)*qgg)**cgrau &
                                 +(1.-omg)*vsnow*(rho(k)*qss)**csnow)
    else
       omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
       qrr=omp*qp(i,j,k)
       qss=qp(i,j,k)-qrr
       qgg=omg*qss
       qss=qss-qgg
       term_vel_qp = (omp*vrain*(rho(k)*qrr)**crain &
                     +(1.-omp)*(omg*vgrau*(rho(k)*qgg)**cgrau &
                          +(1.-omg)*vsnow*(rho(k)*qss)**csnow))
    endif
  end if  
end function term_vel_qp

!----------------------------------------------------------------------
!!! compute sedimentation 
!
subroutine micro_precip_fall()
  
  use vars
  use params, only : pi

  real omega(nx,ny,nzm)
  real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real f0(nzm),df0(nzm)
  real dummy(1)
  integer ind
  integer i,j,k


  crain = b_rain / 4.
  csnow = b_snow / 4.
  cgrau = b_grau / 4.
  vrain = a_rain * gamr3 / 6. / (pi * rhor * nzeror) ** crain
  vsnow = a_snow * gams3 / 6. / (pi * rhos * nzeros) ** csnow
  vgrau = a_grau * gamg3 / 6. / (pi * rhog * nzerog) ** cgrau

! Initialize arrays that accumulate surface precipitation flux

 if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
   do j=1,ny
    do i=1,nx
     precsfc(i,j)=0.
    end do
   end do
   do k=1,nzm
    precflux(k) = 0.
   end do
 end if

 do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
    qpfall(k)=0.
    tlat(k) = 0.
 end do
   
 do k=1,nzm
  do j=1,ny
   do i=1,nx
       omega(i,j,k) = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
   end do
  end do
 end do

 if(dostatis) then
   do k=1,nzm
     do j=dimy1_s,dimy2_s
       do i=dimx1_s,dimx2_s
          df(i,j,k) = t(i,j,k)
       end do
     end do
   end do
 endif

 call precip_fall(qp, term_vel_qp, 2, omega, ind)

 if(dostatis) then
   call stat_varscalar(t,df,f0,df0,t2leprec)
   call setvalue(twleprec,nzm,0.)
   call stat_sw2(t,df,twleprec)
 endif

end subroutine micro_precip_fall

!----------------------------------------------------------------------
!!!! Collect microphysics history statistics (vertical profiles)
!
subroutine micro_statistics()
  
  use vars
  use hbuffer, only: hbuf_put
  use params, only : lcond

  real tmp(2), factor_xy 
  real qcz(nzm), qiz(nzm), qrz(nzm), qsz(nzm), qgz(nzm), omg
  integer i,j,k,n
  character(LEN=6) :: statname  !bloss: for conditional averages

  call t_startf ('micro_statistics')

  factor_xy = 1./float(nx*ny)

  do k=1,nzm
      tmp(1) = dz/rhow(k)
      tmp(2) = tmp(1) / dtn
      mkwsb(k,1) = mkwsb(k,1) * tmp(1) * rhow(k) * lcond
      mkwle(k,1) = mkwle(k,1)*tmp(2)*rhow(k)*lcond + mkwsb(k,1)
      if(docloud.and.doprecip) then
        mkwsb(k,2) = mkwsb(k,2) * tmp(1) * rhow(k) * lcond
        mkwle(k,2) = mkwle(k,2)*tmp(2)*rhow(k)*lcond + mkwsb(k,2)
      endif
  end do

  call hbuf_put('QTFLUX',mkwle(:,1),factor_xy)
  call hbuf_put('QTFLUXS',mkwsb(:,1),factor_xy)
  call hbuf_put('QPFLUX',mkwle(:,2),factor_xy)
  call hbuf_put('QPFLUXS',mkwsb(:,2),factor_xy)

  do k=1,nzm
    qcz(k) = 0.
    qiz(k) = 0.
    qrz(k) = 0.
    qsz(k) = 0.
    qgz(k) = 0.
    do j=1,ny
    do i=1,nx
      qcz(k)=qcz(k)+qcl(i,j,k)
      qiz(k)=qiz(k)+qci(i,j,k)
      qrz(k)=qrz(k)+qpl(i,j,k)
      omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
      qsz(k)=qsz(k)+qpi(i,j,k)*(1.-omg)
      qgz(k)=qgz(k)+qpi(i,j,k)*omg
    end do
    end do
  end do

  call hbuf_put('QC',qcz,1.e3*factor_xy)
  call hbuf_put('QI',qiz,1.e3*factor_xy)
  call hbuf_put('QR',qrz,1.e3*factor_xy)
  call hbuf_put('QS',qsz,1.e3*factor_xy)
  call hbuf_put('QG',qgz,1.e3*factor_xy)

  call hbuf_put('QTADV',mkadv(:,1)+qifall,factor_xy*86400000./dtn)
  call hbuf_put('QTDIFF',mkdiff(:,1),factor_xy*86400000./dtn)
  call hbuf_put('QTSINK',qpsrc,-factor_xy*86400000./dtn)
  call hbuf_put('QTSRC',qpevp,-factor_xy*86400000./dtn)
  call hbuf_put('QPADV',mkadv(:,2),factor_xy*86400000./dtn)
  call hbuf_put('QPDIFF',mkdiff(:,2),factor_xy*86400000./dtn)
  call hbuf_put('QPFALL',qpfall,factor_xy*86400000./dtn)
  call hbuf_put('QPSRC',qpsrc,factor_xy*86400000./dtn)
  call hbuf_put('QPEVP',qpevp,factor_xy*86400000./dtn)

  do n = 1,nmicro_fields
     call hbuf_put(trim(mkname(n))//'LSADV', &
          mklsadv(:,n),mkoutputscale(n)*factor_xy*86400.)
  end do

  do n = 1,ncondavg

     do k=1,nzm
        qcz(k) = 0.
        qiz(k) = 0.
        qrz(k) = 0.
        qsz(k) = 0.
        qgz(k) = 0.
        do j=1,ny
           do i=1,nx
              qcz(k)=qcz(k)+qcl(i,j,k)*condavg_mask(i,j,k,n)
              qiz(k)=qiz(k)+qci(i,j,k)*condavg_mask(i,j,k,n)
              qrz(k)=qrz(k)+qpl(i,j,k)*condavg_mask(i,j,k,n)
              omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
              qsz(k)=qsz(k)+qpi(i,j,k)*(1.-omg)*condavg_mask(i,j,k,n)
              qgz(k)=qgz(k)+qpi(i,j,k)*omg*condavg_mask(i,j,k,n)
           end do
        end do
     end do

     call hbuf_put('QC' // TRIM(condavgname(n)),qcz,1.e3)
     call hbuf_put('QI' // TRIM(condavgname(n)),qiz,1.e3)
     if(doprecip) then
        call hbuf_put('QR' // TRIM(condavgname(n)),qrz,1.e3)
        call hbuf_put('QS' // TRIM(condavgname(n)),qsz,1.e3)
        call hbuf_put('QG' // TRIM(condavgname(n)),qgz,1.e3)
     end if
  end do

  ncmn = 0.
  nrmn = 0.

  call t_stopf ('micro_statistics')

end subroutine micro_statistics

!----------------------------------------------------------------------
! called when stepout() called

subroutine micro_print()
end subroutine micro_print

!----------------------------------------------------------------------
!!! Initialize the list of microphysics statistics 
!
subroutine micro_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)

  use vars


   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count,trcount
   integer ntr, n


   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QTFLUX'
   deflist(count) = 'Nonprecipitating water flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QTFLUXS'
   deflist(count) = 'Nonprecipitating-water flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QPFLUX'
   deflist(count) = 'Precipitating-water turbulent flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QPFLUXS'
   deflist(count) = 'Precipitating-water turbulent flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0
   
   do n = 1,nmicro_fields
      count = count + 1
      trcount = trcount + 1
      namelist(count) = TRIM(mkname(n))//'LSADV'
      deflist(count) = 'Source of '//TRIM(mklongname(n))//' due to large-scale vertical advection'
      unitlist(count) = TRIM(mkunits(n))//'day'
      status(count) = 1    
      average_type(count) = 0
   end do

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QC'
   deflist(count) = 'Liquid cloud water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QI'
   deflist(count) = 'Icy cloud water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QR'
   deflist(count) = 'Rain water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QS'
   deflist(count) = 'Snow water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QG'
   deflist(count) = 'Graupel water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

  !bloss: setup to add an arbitrary number of conditional statistics
   do n = 1,ncondavg

      count = count + 1
      trcount = trcount + 1
      namelist(count) = 'QC' // TRIM(condavgname(n))
      deflist(count) = 'Mean Liquid cloud water in ' // TRIM(condavglongname(n))
      unitlist(count) = 'g/kg'
      status(count) = 1    
      average_type(count) = n

      count = count + 1
      trcount = trcount + 1
      namelist(count) = 'QI' // TRIM(condavgname(n))
      deflist(count) = 'Mean Icy cloud water in ' // TRIM(condavglongname(n))
      unitlist(count) = 'g/kg'
      status(count) = 1    
      average_type(count) = n

      if(doprecip) then
         count = count + 1
         trcount = trcount + 1
         namelist(count) = 'QR' // TRIM(condavgname(n))
         deflist(count) = 'Mean Rain water in ' // TRIM(condavglongname(n))
         unitlist(count) = 'g/kg'
         status(count) = 1    
         average_type(count) = n

         count = count + 1
         trcount = trcount + 1
         namelist(count) = 'QS' // TRIM(condavgname(n))
         deflist(count) = 'Mean Snow water in ' // TRIM(condavglongname(n))
         unitlist(count) = 'g/kg'
         status(count) = 1    
         average_type(count) = n

         count = count + 1
         trcount = trcount + 1
         namelist(count) = 'QG' // TRIM(condavgname(n))
         deflist(count) = 'Mean Graupel water in ' // TRIM(condavglongname(n))
         unitlist(count) = 'g/kg'
         status(count) = 1    
         average_type(count) = n
      end if

   end do

end subroutine micro_hbuf_init


!-----------------------------------------
subroutine micro_stat_2Dinit(ResetStorage)
  implicit none
  integer, intent(in) :: ResetStorage

  ! initialize microphysical 2D outputs as necessary

  if(ResetStorage.eq.1) then
    !bloss: If computing storage terms for individual hydrometeors,
    !         store initial profiles for computation of storage terms in budgets
    !bloss mkstor(:,:) = mk0(:,:)
  end if

end subroutine micro_stat_2Dinit
!-----------------------------------------
subroutine micro_write_fields2D(nfields1)
  implicit none
  integer, intent(inout) :: nfields1

  ! nothing here yet

end subroutine micro_write_fields2D

!-----------------------------------------
subroutine micro_write_fields3D(nfields1)
  implicit none
  integer, intent(inout) :: nfields1

  ! nothing here yet

end subroutine micro_write_fields3D

!-----------------------------------------------------------------------
! Supply function that computes total water in a domain:
!
real(8) function total_water()

  use vars, only : nstep,nprint,adz,dz,rho
  real(8) tmp
  integer i,j,k,m

  total_water = 0.
  do m=1,nmicro_fields
   if(flag_wmass(m).eq.1) then
    do k=1,nzm
      tmp = 0.
      do j=1,ny
        do i=1,nx
          tmp = tmp + micro_field(i,j,k,m)
        end do
      end do
      total_water = total_water + tmp*adz(k)*dz*rho(k)
    end do
   end if
  end do

end function total_water

! -------------------------------------------------------------------------------
! dummy effective radius functions:

logical function micro_provides_reffc()
  micro_provides_reffc = .false.
end function micro_provides_reffc

logical function micro_provides_reffi()
  micro_provides_reffi = .false.
end function micro_provides_reffi

function Get_reffc() ! liquid water
  real, pointer, dimension(:,:,:) :: Get_reffc
end function Get_reffc

function Get_reffi() ! ice
  real, pointer, dimension(:,:,:) :: Get_reffi
end function Get_reffi

function Get_nca() ! aerosol
  real, pointer, dimension(:,:,:) :: Get_nca
  Get_nca = 0.
end function Get_nca


end module microphysics



