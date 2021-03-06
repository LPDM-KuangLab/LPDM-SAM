module micro_params

   !bloss: added options that may be set in prm file namelist 
   !         -- initialized in micrphysics.f90
   logical, public :: &
        dototalwater = .true., &       ! use total water variable (vapor + cloud liquid)
        doicemicro = .true., &       ! use ice species (snow/cloud ice/graupel)
        dograupel = .true., &        ! use graupel
        dohail = .false., &          ! make graupel species have properties of hail
        dosb_warm_rain = .false., &  ! use Seifert & Beheng (2001) warm rain parameterization
        dopredictNc = .false., &     ! prediction of cloud droplet number
        dospecifyaerosol = .false., &! specify two modes of (sulfate) aerosol
        dosubgridw = .false., &      ! input estimate of subgrid w to microphysics
        doarcticicenucl = .false., & ! use arctic parameter values for ice nucleation
        docloudedgeactivation = .false. ! activate cloud droplets throughout the cloud

! Marat's scaling of autoconversion/activation based on grid spacing (dx).
   logical, public :: do_scale_dependence_of_autoconv = .true. 
   logical, public :: do_scale_dependence_of_activation = .true. 

   real, public :: Nc0 = 100. ! specified cloud droplet number conc (#/cm3)

   ! the following are used when dopredictNc==.true. and dospecifyaerosol==.false.
   ! Default is maritime value (/cm3), adapted from Rasmussen et al (2002) 
   !  by Hugh Morrison et al.  Values of 1000. and 0.5 suggested for continental
   real, public :: ccnconst = 120., ccnexpnt = 0.4
!bloss   real, public :: ccnconst = 1000., ccnexpnt = 0.5 ! continental suggestion from module_mp_graupel

   ! the following two modes of sulfate aerosol are used for cloud droplet activation
   !    when dopredictNc==.true. and dospecifyaerosol==.true.
   real, public :: &
        aer_rm1 = 0.052, & ! Mode 1: rm=geom mean radius (um)
        aer_n1 = 72.2, &   !         n=aer conc. (#/cm3)
        aer_sig1 = 2.04, & !         sig=geom standard deviation of aer size distn.
        aer_rm2 = 1.3, &   ! Mode 2: rm=geom mean radius (um)
        aer_n2 = 1.8, &    !         n=aer conc. (#/cm3)
        aer_sig2 = 2.5     !         sig=geom standard deviation of aer size distn.
   logical, public :: AerosolInputAsNumberPerMilligram = .false. ! self-descriptive

   ! option to fix value of pgam (exponent in cloud water gamma distn)
   logical, public :: dofix_pgam = .false.
   real, public ::    pgam_fixed = 10.3 ! Geoffroy et al (2010, doi:10.5194/acp-10-4835-2010)

   real, parameter :: rho_snow = 100., rho_water = 997., rho_cloud_ice = 500.

!bloss (Apr 09): Add option for output of cloud radar reflectivity.
!                Computed using quickbeam cloud radar simulator.
!                Will be output as histogram in statistics file 
!                (with nradar_bins bins between -40 and 20 dBZ) and
!                in 3D files as a full, instantaneous 3D field.
logical :: doreflectivity_cloudradar = .false.
integer :: binwidth_cloudradar = 5, & ! width of bins in dBZ for histogram output
     min_dBZbin_cloudradar = -40, max_dBZbin_cloudradar = 20 ! histogram edges
real*8 :: freq_cloudradar = 95., & ! Cloud radar frequency in GHz
     k2_cloudradar = -1., &  ! dielectric constand -- negative value lets this be computed within quickbeam
     missing_value_cloudradar = -9999. ! missing value for output.
integer :: surface_cloudradar = 0, & !
     usegasabs_cloudradar = 1, &
     doray_cloudradar = 0

! quickbeam is really expensive.  This parameter lets you call it less
!   frequently than the other statistics are computed
integer :: nskip_quickbeam = -1 

!bloss(24Apr2013): Add outputs of microphysical process rates
!   If icemicro==.true., this amount to an extra 68 outputs
!   in the statistics file.  Only 11 for warm cloud microphysics.
   logical, public :: do_output_micro_process_rates = .false. 
   integer :: nmicro_proc
   integer, parameter :: nmicro_process_rates = 74
   integer, parameter :: nmicro_process_rates_warm = 11
   character(len=8), dimension(nmicro_process_rates), parameter, public :: &
        micro_process_rate_names = (/ &
     'PCC     ', & ! COND/EVAP DROPLETS
     'PCCN    ', & ! CHANGE Q DROPLET ACTIVATION
     'NSUBC   ', & ! LOSS OF NC DURING EVAP
     'PRE     ', & ! EVAP OF RAIN
     'NSUBR   ', & ! LOSS OF NR DURING EVAP
     'PRA     ', & ! ACCRETION DROPLETS BY RAIN
     'NPRA    ', & ! CHANGE IN N DUE TO DROPLET ACC BY RAIN
     'PRC     ', & ! AUTOCONVERSION DROPLETS
     'NPRC    ', & ! CHANGE NC AUTOCONVERSION DROPLETS
     'NPRC1   ', & ! CHANGE NR AUTOCONVERSION DROPLETS
     'NRAGG   ', & ! SELF-COLLECTION OF RAIN
     'NSUBI   ', & ! LOSS OF NI DURING SUB.
     'NSUBS   ', & ! LOSS OF NS DURING SUB.
     'PRD     ', & ! DEP CLOUD ICE
     'PRDS    ', & ! DEP SNOW
     'NNUCCC  ', & ! CHANGE N DUE TO CONTACT FREEZ DROPLETS
     'MNUCCC  ', & ! CHANGE Q DUE TO CONTACT FREEZ DROPLETS
     'NNUCCD  ', & ! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
     'MNUCCD  ', & ! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
     'MNUCCR  ', & ! CHANGE Q DUE TO CONTACT FREEZ RAIN
     'NNUCCR  ', & ! CHANGE N DUE TO CONTACT FREEZ RAIN
     'NSAGG   ', & ! SELF-COLLECTION OF SNOW
     'PRAI    ', & ! CHANGE Q ACCRETION CLOUD ICE
     'PRCI    ', & ! CHANGE Q AUTOCONVERSION CLOUD ICE BY SNOW
     'PSACWS  ', & ! CHANGE Q DROPLET ACCRETION BY SNOW
     'NPSACWS ', & ! CHANGE N DROPLET ACCRETION BY SNOW
     'PSACWI  ', & ! CHANGE Q DROPLET ACCRETION BY CLOUD ICE
     'NPSACWI ', & ! CHANGE N DROPLET ACCRETION BY CLOUD ICE
     'NPRCI   ', & ! CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW
     'NPRAI   ', & ! CHANGE N ACCRETION CLOUD ICE
     'NMULTS  ', & ! ICE MULT DUE TO RIMING DROPLETS BY SNOW
     'NMULTR  ', & ! ICE MULT DUE TO RIMING RAIN BY SNOW
     'QMULTS  ', & ! CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
     'QMULTR  ', & ! CHANGE Q DUE TO ICE RAIN/SNOW
     'PRACS   ', & ! CHANGE Q RAIN-SNOW COLLECTION
     'NPRACS  ', & ! CHANGE N RAIN-SNOW COLLECTION
     'PSMLT   ', & ! CHANGE Q MELTING SNOW TO RAIN
     'EVPMS   ', & ! CHNAGE Q MELTING SNOW EVAPORATING
     'NSMLTS  ', & ! CHANGE N MELTING SNOW
     'NSMLTR  ', & ! CHANGE N MELTING SNOW TO RAIN
     'PIACR   ', & ! CHANGE QR, ICE-RAIN COLLECTION
     'NIACR   ', & ! CHANGE N, ICE-RAIN COLLECTION
     'PRACI   ', & ! CHANGE QI, ICE-RAIN COLLECTION
     'PIACRS  ', & ! CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
     'NIACRS  ', & ! CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW
     'PRACIS  ', & ! CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
     'EPRD    ', & ! SUBLIMATION CLOUD ICE
     'EPRDS   ', & ! SUBLIMATION SNOW
     'PRACG   ', & ! CHANGE IN Q COLLECTION RAIN BY GRAUPEL
     'PSACWG  ', & ! CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
     'PGSACW  ', & ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
     'PGRACS  ', & ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
     'PRDG    ', & ! DEP OF GRAUPEL
     'EPRDG   ', & ! SUB OF GRAUPEL
     'EVPMG   ', & ! CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
     'PGMLT   ', & ! CHANGE Q MELTING OF GRAUPEL
     'NPRACG  ', & ! CHANGE N COLLECTION RAIN BY GRAUPEL
     'NPSACWG ', & ! CHANGE N COLLECTION DROPLETS BY GRAUPEL
     'NSCNG   ', & ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
     'NGRACS  ', & ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
     'NGMLTG  ', & ! CHANGE N MELTING GRAUPEL
     'NGMLTR  ', & ! CHANGE N MELTING GRAUPEL TO RAIN
     'NSUBG   ', & ! CHANGE N SUB/DEP OF GRAUPEL
     'PSACR   ', & ! CONVERSION DUE TO COLL OF SNOW BY RAIN
     'NMULTG  ', & ! ICE MULT DUE TO ACC DROPLETS BY GRAUPEL
     'NMULTRG ', & ! ICE MULT DUE TO ACC RAIN BY GRAUPEL
     'QMULTG  ', & ! CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
     'QMULTRG ', & ! CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL
     'QHOMOC  ', & ! CHANGE Q DUE TO HOMOGENEOUS FREEZING OF CLOUD WATER
     'QHOMOR  ', & ! CHANGE Q DUE TO HOMOGENEOUS FREEZING OF RAIN
     'NHOMOC  ', & ! CHANGE N DUE TO HOMOGENEOUS FREEZING OF CLOUD WATER
     'NHOMOR  ', & ! CHANGE N DUE TO HOMOGENEOUS FREEZING OF RAIN
     'QMELTI  ', & ! CHANGE Q DUE TO MELTING OF CLOUD ICE
     'NMELTI  ' /) ! CHANGE N DUE TO MELTING OF CLOUD ICE

   character(len=80), dimension(nmicro_process_rates), parameter, public :: &
        micro_process_rate_longnames = (/ &
     'PCC     , COND/EVAP DROPLETS                                                    ', &
     'PCCN    , CHANGE Q DROPLET ACTIVATION                                           ', &
     'NSUBC   , LOSS OF NC DURING EVAP                                                ', &
     'PRE     , EVAP OF RAIN                                                          ', &
     'NSUBR   , LOSS OF NR DURING EVAP                                                ', &
     'PRA     , ACCRETION DROPLETS BY RAIN                                            ', &
     'NPRA    , CHANGE IN N DUE TO DROPLET ACC BY RAIN                                ', &
     'PRC     , AUTOCONVERSION DROPLETS                                               ', &
     'NPRC    , CHANGE NC AUTOCONVERSION DROPLETS                                     ', &
     'NPRC1   , CHANGE NR AUTOCONVERSION DROPLETS                                     ', &
     'NRAGG   , CHANGE IN NR DUE TO SELF-COLLECTION AND BREAKUP OF RAIN               ', &
     'NSUBI   , LOSS OF NI DURING SUB.                                                ', &
     'NSUBS   , LOSS OF NS DURING SUB.                                                ', &
     'PRD     , DEP CLOUD ICE                                                         ', &
     'PRDS    , DEP SNOW                                                              ', &
     'NNUCCC  , CHANGE N DUE TO CONTACT FREEZ DROPLETS                                ', &
     'MNUCCC  , CHANGE Q DUE TO CONTACT FREEZ DROPLETS                                ', &
     'NNUCCD  , CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)                       ', &
     'MNUCCD  , CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)                       ', &
     'MNUCCR  , CHANGE Q DUE TO CONTACT FREEZ RAIN                                    ', &
     'NNUCCR  , CHANGE N DUE TO CONTACT FREEZ RAIN                                    ', &
     'NSAGG   , SELF-COLLECTION OF SNOW                                               ', &
     'PRAI    , CHANGE Q ACCRETION CLOUD ICE                                          ', &
     'PRCI    , CHANGE Q AUTOCONVERSION CLOUD ICE BY SNOW                             ', &
     'PSACWS  , CHANGE Q DROPLET ACCRETION BY SNOW                                    ', &
     'NPSACWS , CHANGE N DROPLET ACCRETION BY SNOW                                    ', &
     'PSACWI  , CHANGE Q DROPLET ACCRETION BY CLOUD ICE                               ', &
     'NPSACWI , CHANGE N DROPLET ACCRETION BY CLOUD ICE                               ', &
     'NPRCI   , CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW                             ', &
     'NPRAI   , CHANGE N ACCRETION CLOUD ICE                                          ', &
     'NMULTS  , ICE MULT DUE TO RIMING DROPLETS BY SNOW                               ', &
     'NMULTR  , ICE MULT DUE TO RIMING RAIN BY SNOW                                   ', &
     'QMULTS  , CHANGE Q DUE TO ICE MULT DROPLETS/SNOW                                ', &
     'QMULTR  , CHANGE Q DUE TO ICE RAIN/SNOW                                         ', &
     'PRACS   , CHANGE Q RAIN-SNOW COLLECTION                                         ', &
     'NPRACS  , CHANGE N RAIN-SNOW COLLECTION                                         ', &
     'PSMLT   , CHANGE Q MELTING SNOW TO RAIN                                         ', &
     'EVPMS   , CHNAGE Q MELTING SNOW EVAPORATING                                     ', &
     'NSMLTS  , CHANGE N MELTING SNOW                                                 ', &
     'NSMLTR  , CHANGE N MELTING SNOW TO RAIN                                         ', &
     'PIACR   , CHANGE QR, ICE-RAIN COLLECTION                                        ', &
     'NIACR   , CHANGE N, ICE-RAIN COLLECTION                                         ', &
     'PRACI   , CHANGE QI, ICE-RAIN COLLECTION                                        ', &
     'PIACRS  , CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW                          ', &
     'NIACRS  , CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW                           ', &
     'PRACIS  , CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW                          ', &
     'EPRD    , SUBLIMATION CLOUD ICE                                                 ', &
     'EPRDS   , SUBLIMATION SNOW                                                      ', &
     'PRACG   , CHANGE IN Q COLLECTION RAIN BY GRAUPEL                                ', &
     'PSACWG  , CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL                            ', &
     'PGSACW  , CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW            ', &
     'PGRACS  , CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW                ', &
     'PRDG    , DEP OF GRAUPEL                                                        ', &
     'EPRDG   , SUB OF GRAUPEL                                                        ', &
     'EVPMG   , CHANGE Q MELTING OF GRAUPEL AND EVAPORATION                           ', &
     'PGMLT   , CHANGE Q MELTING OF GRAUPEL                                           ', &
     'NPRACG  , CHANGE N COLLECTION RAIN BY GRAUPEL                                   ', &
     'NPSACWG , CHANGE N COLLECTION DROPLETS BY GRAUPEL                               ', &
     'NSCNG   , CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW     ', &
     'NGRACS  , CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW         ', &
     'NGMLTG  , CHANGE N MELTING GRAUPEL                                              ', &
     'NGMLTR  , CHANGE N MELTING GRAUPEL TO RAIN                                      ', &
     'NSUBG   , CHANGE N SUB/DEP OF GRAUPEL                                           ', &
     'PSACR   , CONVERSION DUE TO COLL OF SNOW BY RAIN                                ', &
     'NMULTG  , ICE MULT DUE TO ACC DROPLETS BY GRAUPEL                               ', &
     'NMULTRG , ICE MULT DUE TO ACC RAIN BY GRAUPEL                                   ', &
     'QMULTG  , CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL, SINK OF QC, SOURCE OF QI   ', &
     'QMULTRG , CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL, SINK OF QR, SOURCE OF QI       ', &
     'QHOMOC  , CHANGE Q DUE TO HOMOGENEOUS FREEZING OF CLOUD WATER TO FORM CLOUD ICE ', &
     'QHOMOR  , CHANGE Q DUE TO HOMOGENEOUS FREEZING OF RAIN TO FORM GRAUPEL          ', &
     'NHOMOC  , CHANGE N DUE TO HOMOGENEOUS FREEZING OF CLOUD WATER TO FORM CLOUD ICE ', &
     'NHOMOR  , CHANGE N DUE TO HOMOGENEOUS FREEZING OF RAIN TO FORM GRAUPEL          ', &
     'QMELTI  , CHANGE Q DUE TO MELTING OF CLOUD ICE TO FORM RAIN                     ', &
     'NMELTI  , CHANGE N DUE TO MELTING OF CLOUD ICE TO FORM RAIN                     ' /)




 end module micro_params
