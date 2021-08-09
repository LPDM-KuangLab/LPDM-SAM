# The Lagrangian parcel dispersion model

- [Introduction of the model](#introduction)

- [Model description and validation](#description-validation)

- [Implementing the source code](#implementing-src)

- [How to run LPDM](#run-lpdm)

- [Post-processing LPDM output](#post-processing)

- [LPDM related publications](#publication)


## Introduction of the model  <a name="introduction"></a>

This Lagrangian parcel dispersion model (LPDM) is a Fortran module that can be embedded in a Large-Eddy Simulations (LES) models or Cloud-Resolving Models (CRMs). The LPDM can release Lagrangian particles in the LES/CRMs, track their locations, and output their trajectories. This model has been widely used in a variety of topics in atmospheric convection (see [LPDM related publications](#publication)).

## Model description and validation  <a name="description-validation"></a>

Could we share the link of the manuscript (A Lagrangian Particle Model to study Atmospheric Dynamics: the LPDMX) here?

## Implementing the source code <a name="implementing-src"></a>

An example version of SAM source code with LPDM can be found at /LPDM_tutorial/src/SRC_LPDM_example. This SRC_LPDM_example is based on SAM 6.10.6. It is very likely that users want to use a different version of SAM source code. Therefore, we documented here how to implement LPDM model into your customized SAM code.

#### 1. Add lpdm_mod.f90
copy lpdm_mod.f90 (src/SRC_LPDM_example/lpdm_mod.f90) into your SRC directory.

#### 2. Modify grid.f90

add the following codes into grid.f90 (see src/SRC_LPDM_example/grid.f90 line 254-258)

``` fortran
! parameters for lpdm
logical:: lpdm_do =.false.       ! run LPDM jinie
integer:: lpdm_ndt =5      ! lpdm output time interval (dt)
integer:: lpdm_sub =1     ! seperate output on subdomain on x
integer:: lpdm_start =0    ! time step when to start the LPDM
```

#### 3. Modify main.f90

add the following line (src/SRC_LPDM_example/main.f90 line 11)

``` fortran
use lpdm_mod
```

add the following two lines (src/SRC_LPDM_example/main.f90 line 80-81) after call setforcing()

``` fortran
if(lpdm_do)  call lpdm_ini
if(lpdm_do)  call lpdm_output
```

add the following lines (src/SRC_LPDM_example/main.f90 line 350-357) after call stepout(nstatsteps)

``` fortran
! before running lpdm, update the boundaries informations

 if (lpdm_do .and. nstep.ge.lpdm_start) then
   call boundaries(0)
   call boundaries(2)
   call lpdm_update
   call lpdm_output
 endif

```

#### 4. Modify printout.f90

add the following line (src/SRC_LPDM_example/printout.f90 line 6)

``` fortran
use lpdm_mod, only: lpdm_top,lpdm_num
```

add the following lines (src/SRC_LPDM_example/printout.f90 line 253-260)

``` fortran
print*,'run LPDM = ', lpdm_do
if (lpdm_do) then
print*,'lpdm_top = ',lpdm_top
print*,'lpdm_num = ',lpdm_num
print*,'lpdm_ndt = ',lpdm_ndt
print*,'lpdm_sub = ',lpdm_sub
print*,'lpdm_start = ', lpdm_start
endif
```

#### 5. Modify setparm.f90

add the following line (src/SRC_LPDM_example/setparm.f90 line 49) to add the following variables into the "NAMELIST /PARAMETERS/".

``` fortran
lpdm_do, lpdm_ndt, lpdm_sub, lpdm_start ! for LPDM
```

## How to run LPDM  <a name="run-lpdm"></a>

A full description of how to run SAM is documented in [Harvard Climate Modeling](https://wiki.harvard.edu/confluence/display/climatemodeling/SAM), including how to setup turn on LPDM and properly set LPDM related parameters. The following content here are basically copied from [Harvard Climate Modeling](https://wiki.harvard.edu/confluence/display/climatemodeling/SAM):

#### During compiling
During compiling the model (before building the model), to set the Lagrangian particle dispersion model, go to SRC/lpdm_mod.f90 and set lpdm_top and lpdm_num. Each Lagrangian particle represents an air parcel with a given mass of m (m depends on lpdm_top and lpdm_num and is the same for all the Lagrangian particles).
``` fortran
integer, parameter :: lpdm_top  = 70  ! top level below which particles are initially released, mirror reflection b.c.
integer, parameter :: lpdm_num   = 700  ! the number of particles initially in one column
```

To support lpdm, nx_gl and ny_gl has the upper limit of 999, and  nsubdomain_x *nsubdomain_y has the upper limit of 214 (or edit source code)

#### Before submitting a job
To run lpdm model, you need to specify parameters in the prm file under your specific case directory

```
lpdm_do =  .true ,  ! true or false, set it to be true if you want to include particle model
lpdm_ndt = 4,   ! frequency to output particle location
lpdm_sub=64,  ! total number of processors used, should be consistent with the total number of cores specified
```


For the current LPDM model, it does not have the restart function, I will add this function in the future. Therefore if your run stops
before it finishes, you will need to start from the beginning. REMEMBER: Remove the existing lpdm output first!!


## Post-processing LPDM output  <a name="post-processing"></a>

LPDM outputs are stored in [SAM output directory]/OUT_LPDM

Here we provide an example script (utils/lpdm_loadxyz.m) which can convert the .xyz output into .mat files. Please follow the comments in lpdm_loadxyz.m to adjust paths and parameters properly to your own case.

## LPDM related publications  <a name="publication"></a>

Tian, Y. and Kuang, Z.: Dependence of entrainment in shallow cumulus convection on vertical velocity and distance to cloud edge, Geophys- ical Research Letters, 43, 4056–4065, 2016.

Tian, Y. and Kuang, Z.: Why Does Deep Convection Have Different Sensitivities to Temperature Perturbations in the Lower versus Upper Troposphere?, Journal Of The Atmospheric Sciences, 76, 27–40, 2019.

Tian, Y., Z, K., Singh, M. S., and Nie, J.: The Vertical Momentum Budget of Shallow Cumulus Convection: Insights From a Lagrangian Perspective, Journal of Advances in Modeling Earth Systems, 2019.

Torri, G. and Kuang, Z.: A Lagrangian Study of Precipitation-Driven Downdrafts, Journal of the Atmospheric Sciences, 73, 839–853, 2016.

Torri, G., Kuang, Z., and Tian, Y.: Mechanisms for convection triggering by cold pools, Geophysical Research Letters, 42, 1943–1950, 2015.
