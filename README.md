This is my fork of FVCOM-ICM, used for the Salish Sea Model (SSM) in
conjunction with a modified version of the FVCOM hydrodynamic model. I have
not been keeping up with changes since 2022; my goals have been to make
small cleanups, optimizations, compatibility/usability improvements, and
perform some code verficiation. The only modeling enhancement made so far
has been to allow the air/sea exchange method 6 (hybrid with Wanninkhof et
al, 2013) to work when the model is running for more than one calendar
year.

I am only hosting the model source code here, not any input files it
needs.

# Building

## METIS

Build the bundled version of METIS. In METIS_source, create a file named
local.mk to override any compiler options set in the makefile. Then type
`make`.

## The model

First, review the configurable options in the Makefile, as which are turned
on affect the necessary dependencies. The netcdf fortran libraries and an
MPI library are typical. I have also added some unit tests to the code that
require [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit) to
build. Note that I haven't been keeping up-to-date on this, and version 4.4
is known to work.

Similar to METIS, create a file named local.mk to set and preserve custom
compiler options and flags. Then type `make`.

For many uses, two separate builds of the model are needed. The default
build assumes initial alkalinity values are specified in an input file, but
for cold starts of the SSM this data is absent. An alternate build of the
    model can be made by enabling `FLAG_5`, which conventionally also
    involves changing the name of the compiled executable. The typical
    method is to save the first compiled executable, run `make clean`,
    adjust compiler options, and rebuild to get the second version. This is
    really one of those usability issues I should eventually fix.

# Original Readme

Here is the download package (source code and other utilities) for the FVCOM-ICMv4.0 model, which 
was originally developed at University of Massachusetts with support from U.S. Army Corps of Engineers USACE. 
Pacific Northwest National Laboratory (PNNL) completed subsequent developments and improvements as part of an effort 
to establish a robust unstructured biogeochemical model for Puget Sound, Georgia Basin, and the Strait of Juan 
de Fuca (Salish Sea) for U.S. EPA and Washington State Department of Ecology (Ecology). 
The unstructured biogeochemical model FVCOM-ICM has been through numerous phases of testing, and model 
performance has been well documented via peer-reviewed journal publications and several study reports. 
Initial testing and Phase 1 and 2 with model applications to the Puget Sound domain were completed in 2012. 
The most recent updated version, 4.0 includes explicit calculations for kinetics of zooplankton and submerged 
aquatic vegetations. The model is going through continuous developments and improvements. Thus, users should 
use the recommended procedures together with their own best judgments in model applications.
