#  Minor changes to turn on -DMULTIPROCESSOR
#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/bash
            DEF_FLAGS     = -P -C -traditional
#             EXEC          = FVCOM_ICM_v4_pH_TAinitialFromS
            EXEC           = FVCOM_ICM_v4_pH_TAinitialFromInput

#==========================================================================
#  BEGIN USER DEFINITION SECTION
#==========================================================================
#        SELECT MODEL OPTIONS
#          SELECT FROM THE FOLLOWING OPTIONS BEFORE COMPILING CODE
#          SELECT/UNSELECT BY COMMENTING/UNCOMMENTING LINE (#)
#          CODE MUST BE CLEANED (with "make clean") AND
#          RECOMPILED IF NEW SET OF OPTIONS IS DESIRED
#--------------------------------------------------------------------------

          FLAG_0 = -DDOUBLE_PRECISION

#--------------------------------------------------------------------------
#        MULTI_PROCESSOR    INCLUDES PARALLELIZATION WITH MPI
#                           REQUIRES LINKING MPI LIBRARIES OR COMPILING
#                           WITH A PRELINKED SCRIPT (mpif90/mpf90/etc)
#                           DEFAULT: NO PARALLEL CAPABILITY
#                           UNCOMMENT TO INCLUDE MPI PARALLEL CAPABILITY
#--------------------------------------------------------------------------

           FLAG_1 = -DMULTIPROCESSOR
           PARLIB = -L../METIS_source -lmetis

#--------------------------------------------------------------------------
#        NETCDF OUTPUT      DUMP OUTPUT INTO NETCDF FILES (yes/no)
#                           REQUIRES SYSTEM DEPENDENT NETCDF LIBRARIES
#                           COMPILED WITH SAME F90 COMPILER
#                           SET PATH TO LIBRARIES WITH IOLIBS
#                           SET PATH TO INCLUDE FILES (netcdf.mod) WITH IOINCS
#                           DEFAULT: NO NETCDF OUTPUT
#                           UNCOMMENT TO INCLUDE NETCDF OUTPUT CAPABILITY
#--------------------------------------------------------------------------
           FLAG_2       =  -DNETCDF_IO
           IOLIBS       =  -L../netcdf_install/lib -lnetcdf # /home/BIAN431/mypsdiag/netcdf_install/lib  -lnetcdf
           IOINCS       =  -I../netcdf_install/include      # /home/BIAN431/mypsdiag/netcdf_install/include

#
#Flag to turn on the various outputs for program checking
#
#
            FLAG_3       = -DKILLTREES

#--------------------------------------------------------------------------
#        TYPE OF DIC SEDIMENT FLUX   to choose one of  two options for
#                                    DIC flux from sediments to water column
#           options are
#              CSOD_CALC:  uses QUAL2K approach (converts carbonaceous
#                          sediment O2 demand into DIC flux)
#              MASS_CALC:  mass balance approach (DIC flux = organic C
#                          reaching sedims - CH4 flux from seds - burial)
#--------------------------------------------------------------------------
            FLAG_4       = -DCSOD_CALC

#--------------------------------------------------------------------------
#        INITIAL TALK CONDITIONS AS A FUNCTION OF INITAL S and INITIAL RHO
#--------------------------------------------------------------------------
#              FLAG_5       = -DTA_INI_FROM_S

#--------------------------------------------------------------------------
#        SELECT COMPILER/PLATFORM SPECIFIC DEFINITIONS
#          SELECT FROM THE FOLLOWING PLATFORMS OR USE "OTHER" TO DEFINE
#          THE FOLLOWING VARIABLES:
#          CPP:  PATH TO C PREPROCESSOR
#           FC:  PATH TO FORTRAN COMPILER (OR MPI COMPILE SCRIPT)
#          OPT:  COMPILER OPTIONS
#       MPILIB:  PATH TO MPI LIBRARIES (IF NOT LINKED THROUGH COMPILE SCRIPT)
#--------------------------------------------------------------------------

#           FLAG_6 = -DWET_DRY	#flag to use wetting and drying

#           FLAG_7 = -DSED_TEMP	#flag for calculating sediment temperature
				#using diffusion from water column

           FLAG_8 = -DINTERP_HYDRO

           FLAG_9 = -DOUTPUT_NETCDF  #flag to turn on output of netcdf results

#           FLAG_10 = -DNEWSIGMA  #option to use spatially varying sigma levels

#           FLAG_11 = -DNO_BIOLOGY -DONLY_AIRSEA

#--------------------------------------------------------------------------
#  Intel Compiler Definitions
#--------------------------------------------------------------------------
#        CPP      = /usr/bin/cpp
#        CPPFLAGS = $(DEF_FLAGS) -DINTEL
#        FC       = /share/apps/intel/Compiler/11.0/069/bin/intel64/ifort
#        DEBFLGS  = #-check all
#        OPT      = #-O3 -xN -axN -tpp7
#        CLIB     = #-static-libcxa
#--------------------------------------------------------------------------
#  Intel Compiler Definitions
#--------------------------------------------------------------------------
          CPP      = /usr/bin/cpp
          CPPFLAGS = $(DEF_FLAGS) -DINTEL
          FC       = mpif90
          #DEBFLGS  = -r8 -warn all -i4 -g -nozero -O0  -check all -g -debug -ftrapuv -debug-parameters -fp-stack-check -fpe0 -traceback -automatic

          #this is a less restrictive debugging option
#          DEBFLGS = -g  -i4 -r8   -fp-stack-check -debug -traceback -warn all -warn unused -check uninit -debug-parameters #-nozero

          OPT      = -r8 -i4 -O2 # -ftrapuv -fpe0 -nozero
#          OPT      = #-heap-arrays  -i4 -O2 -vec-report0 -r8 -xHost -no-prec-div  -no-prec-sqrt -assume nominus0 -assume noprotect_parens -assume norealloc_lhs  -shared-intel -warn unused -check uninit #  -p -g  -fno-inline-functions
          CLIB     = #-static-libcxa
#--------------------------------------------------------------------------
#   Linux/Portland Group Definitions
#--------------------------------------------------------------------------
#         CPP      = /usr/bin/cpp
#         CPPFLAGS = $(DEF_FLAGS)
#         FC       = pgf90
#         DEBFLGS  = -Mbounds -g -Mprof=func
#         OPT      = #-fast  -Mvect=assoc,cachesize:512000,sse
#--------------------------------------------------------------------------
#  gfortran defs
#--------------------------------------------------------------------------
#         CPP      = /usr/bin/cpp
#         CPPFLAGS = $(DEF_FLAGS)  -DGFORTRAN
#         FC       = gfortran  -O3
#         DEBFLGS  =
#         OPT      =
#         CLIB     =
#==========================================================================
#  END USER DEFINITION SECTION
#==========================================================================

         FFLAGS = $(DEBFLGS) $(OPT)
         MDEPFLAGS = --cpp --fext=f90 --file=-
         RANLIB = ranlib

#--------------------------------------------------------------------------
#  CAT Preprocessing Flags
#--------------------------------------------------------------------------

           CPPARGS = $(CPPFLAGS) $(DEF_FLAGS) $(FLAG_0) $(FLAG_1) $(FLAG_2) $(FLAG_3) 		\
                     $(FLAG_4) $(FLAG_5) $(FLAG_6) $(FLAG_7) $(FLAG_8) $(FLAG_9) $(FLAG_10) $(FLAG_11)

#--------------------------------------------------------------------------
#  Libraries
#--------------------------------------------------------------------------

            LIBS  = $(CLIB)  $(PARLIB) $(IOLIBS) $(MPILIB)
#            LIBS  = -L../SPARSKIT2 -lskit $(PARLIB)
#             LIBS  = $(PARLIB) $(PETSC_LIB)
#           INCS  = $(IOINCS)$(PETSC_FC_INCLUDES)$(IOINCS)
           INCS  = $(IOINCS)

#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90 .F .F90

.F.o:
	$(CPP) $(CPPARGS) $(INCS) $*.F > $*.f90
	$(FC)  -c $(FFLAGS) $(INCS) $*.f90
#	\rm $*.f90

#--------------------------------------------------------------------------
#  FVCOM Source Code.
#--------------------------------------------------------------------------

MAIN  = mod_prec.F	mod_types.F	mod_utils.F	mod_input.F	mod_control.F	\
	mod_lims.F	mod_bcs.F	mod_sizes.F	mod_hydrovars.F	\
	mod_buffers.F 	mod_fileinfo.F					\
	mod_wqm.F	mod_zoop.F	mod_wqminit.F	dens2.F		\
	utilities.F	mod_par.F	mod_wd.F	mod_sed_sav_exchange_vars.F	\
	bracket.F	mod_sed_sf_exchange_vars.F 	mod_sf.F	\
	mod_sed_df_exchange_vars.F	mod_df.F	mod_owq.F	\
	mod_sav.F	mod_ba.F	mod_sed.F	mod_algal.F	\
	mod_CO2SYSconstants.F mod_pH_pCO2.F mod_DIC_TALK.F  		\
	mod_bcmap.F	mod_obcs.F	mod_ncd.F	ncdio.F		\
	mod_tge.F	cell_area.F	pdomdec.F			\
	domdec.F	genmap.F	bcs_force.F	wqm_inputs.F	\
	mod_kin.F	mod_hydro.F		tvds.F			\
	adv_wqm.F	vdif_wqm.F	vertvl.F	viscofh.F	\
	bcond_wqm.F	fct_nut.F	mod_filenames.F	wqm_main.F

 SRCS = $(MAIN)

 OBJS = $(SRCS:.F=.o)

#--------------------------------------------------------------------------
#  Linking Directives
#--------------------------------------------------------------------------

$(EXEC):	$(OBJS)
		$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXEC) $(OBJS) $(LIBS)

#--------------------------------------------------------------------------
#  Target to create dependecies.
#--------------------------------------------------------------------------

depend makedepends:
		/pic/projects/psgb_model/makedepf90_working/bin/makedepf90   $(SRCS) > makedepends


#--------------------------------------------------------------------------
#  Tar Up Code
#--------------------------------------------------------------------------

tarfile:
	tar cvf fvcom-icm.tar *.F  makefile makedepends

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod *.f90

clobber:	clean
		/bin/rm -f *.f90 *.o $(EXEC)

#--------------------------------------------------------------------------
#  Common rules for all Makefiles - do not edit.
#--------------------------------------------------------------------------

emptyrule::

#--------------------------------------------------------------------------
#  Empty rules for directories that do not have SUBDIRS - do not edit.
#--------------------------------------------------------------------------

install::
	@echo "install in $(CURRENT_DIR) done"

install.man::
	@echo "install.man in $(CURRENT_DIR) done"

Makefiles::

includes::
include makedepends
# DO NOT DELETE
