!mod_ncd.F
!************************************************************************
!**                                                                    **
!**                           FVCOM-ICM_4.0                            **
!**                                                                    **
!**               A Finite Volume Based Integrated Compartment         **
!**                         Water Quality Model                        **      
!**        The original unstructured-grid ICM code was developed by    ** 
!**    the FVCOM development team at the University of Massachusetts   ** 
!**         through a contract with U.S. Army Corps of Engineers       ** 
!**         [Dr. Changsheng Chen (PI), Dr. Jianhua Qi and              ** 
!**                      Dr. Geoffrey W. Cowles]                       **
!**                                                                    **
!**                Subsequent Development and Maintenance by           ** 
!**                   PNNL/UW Salish Sea Modeling Center               **
!**                                                                    **
!**                 Tarang Khangaonkar    :  PNNL (2008 - Present)     **
!**                 Lakshitha Premathilake:  PNNL (2019 - Present)     **
!**                 Adi Nugraha           :  PNNL/UW (2018 - Present)  **
!**                 Kurt Glaesmann        :  PNNL (2008 - Present)     **
!**                 Laura Bianucci        :  PNNL/DFO(2015 - Present)  **
!**                 Wen Long              :  PNNL (2012-2016)          **
!**                 Taeyum Kim            :  PNNL (2008-2011)          **
!**                 Rochelle G Labiosa    :  PNNL (2009-2010)          **
!**                                                                    **
!**                                                                    **
!**                     Adopted from CE-QUAL-ICM  Model                **
!**                           Developed by:                            **
!**                                                                    **
!**             Carl F. Cerco      : Water quality scheme              **
!**             Raymond S. Chapman : Numerical solution scheme         **
!**             Thomas M. Cole     : Computer algorithms & coding      **
!**             Hydroqual          : Sediment compartment              **
!**                                                                    **
!**                    Water Quality Modeling Group                    **
!**                    U.S. Army Corps of Engineers                    **
!**                    Waterways Experiment Station                    **
!**                    Vicksburg, Mississippi 39180                    **
!**                                                                    **
!************************************************************************
!
!==============================================================================|
!   NETCDF UTILITIES                                                              !
!==============================================================================|
!
!tools for read and out netcdf results
!

Module MOD_NCD
  !
      Use NETCDF
      Use MOD_PREC, Only: SP, CDF_PREC !

	  Implicit None
	  Save
#if defined (OUTPUT_NETCDF)
!--Control Variables----------------------------------------------!
   logical :: cdf_out           	!!true to activate netcdf input/output
   integer :: nout_vars !=107        	!!number of variables to output 
   integer :: cdf_int           	!!output every cdf_int iterations
   integer :: cdf_stk           	!!cdf_stk outputs are put in each file
									!!CDF_STK=0: ALL OUTPUTS IN SINGLE FILE
   integer :: stck_cnt          	!!counts number of outputs in each file
   integer :: out_cnt           	!!counts number of outputs
   character(len=1024) :: cdfname 	!!netcdf file name

   character(len=256), allocatable, dimension(:) :: cdf_vdp 




!--NetCDF IDs----------------------------------------------------!

   !--NetCDF File
   integer :: nc_ofid

   !--Dimensions
   integer :: node_did,nele_did
   integer :: scl_did,siglay_did,siglev_did

   integer :: three_did,four_did,obc_did,obc2_did, nine_did, sav_did
   integer :: time_did
   integer :: nsav_did 

   !--Grid Variables
   integer :: nprocs_vid,partition_vid
   integer :: idens_vid
   integer :: x_vid,y_vid!,lat_vid,lon_vid

   integer :: nv_vid,nbe_vid
   !integer :: aw0_vid,awx_vid,awy_vid
   !integer :: a1u_vid,a2u_vid
   integer :: siglay_vid,siglev_vid!,siglay_shift_vid

#if defined (NEWSIGMA)
   integer :: siglev3dn_vid
#endif

   !--Flow Variables
   integer :: time_vid
   integer :: iint_vid
   integer :: u_vid
   integer :: v_vid
   !integer :: w_vid
   integer :: wd_vid
   !integer :: ww_vid
   integer :: wts_vid

   integer :: s1_vid
   integer :: t1_vid
   integer :: el_vid
   integer :: h_vid
!   integer :: km_vid
   integer :: kh_vid
!   integer :: viscofh_vid
   !integer :: ua_vid
   !integer :: va_vid
   integer :: d_vid
   !integer :: dtfa_vid
!   integer :: xflux_obc_vid
!   integer :: uard_obcn_vid

	Integer ::	 T_vid 			&!(1:MGL, 1:KBM1)			! Correction to assign values for all layers
				,SALT_vid		&! (1:MGL, 1:KBM1)		! Correction to assign values for all layers
				!,D_vid 			&!(1:MGL)	&! Added D_GL, EL_GL, H_GL
				!,EL_vid			&! (1:MGL)	&! 
				!,H_vid 			&!(1:MGL)	&! 
				,CCHL1_vid 		&!
				,total_netPP_vid&!
				,DOXG_vid		&! C2 (1:MGL, 1:KBM1, 27)	!
				,LDOC_vid		&! C2 (1:MGL, 1:KBM1, 9)  		!LDOC
				,RDOC_vid		&!C2 (1:MGL, 1:KBM1, 10)		!RDOC
				,LPOC_vid 		&!C2 (1:MGL, 1:KBM1, 11)		!LPOC
				,RPOC_vid		&!C2 (1:MGL, 1:KBM1, 12)		!RPOC
				,B1_vid			&!C2 (1:MGL, 1:KBM1, 4)		!B1
				,B2_vid			&!C2 (1:MGL, 1:KBM1, 5)		!B2
				,NH4_vid		&!C2 (1:MGL, 1:KBM1, 13)		!NH4
				,NO3_vid		&!C2 (1:MGL, 1:KBM1, 14)		!NO3
				,LDON_vid		&!C2 (1:MGL, 1:KBM1, 16) 		!LDON
				,RDON_vid		&!C2 (1:MGL, 1:KBM1, 17)		!RDON
				,LPON_vid		&!C2 (1:MGL, 1:KBM1, 18)		!LPON
				,RPON_vid		&!C2 (1:MGL, 1:KBM1, 19)		!RPON
				,PO4_vid		&!C2 (1:MGL, 1:KBM1, 20)		!PO4
				,TDIC_vid		&!C2 (1:MGL, 1:KBM1, 33)		!TDIC	
				,TALK_vid		&!C2 (1:MGL, 1:KBM1, 34)		!TALK 	
				,pH_vid			&!(1:MGL, 1:KBM1)			!pH		
				,pCO2_vid		&!(1:MGL, 1:KBM1)		        !pCO2
				,P1_vid			&
				,P2_vid			&
				,BM1_vid		&
				,BM2_vid		&
				,NL1_vid		&
				,NL2_vid		&
				,PL1_vid		&
				,PL2_vid		&
				,B1SZ_vid		&
				,B2SZ_vid		&
				,B1LZ_vid		&
				,B2LZ_vid		&
				,PR1_vid		&
				,PR2_vid		&
            !If (LIGHT_EXTINCTION) Then
				,IAVG_vid		&
				,FI1_vid		&
				,FI2_vid		&
            !End If
				,IK1_vid		&
				,IK2_vid		&
				,GPP_vid 		&!
				,NPP_vid 		&!
				,DDOC_vid		&
				,DCOD_vid		&
				,NITRIF_vid		&
				,DOSZ_vid		&
				,DOLZ_vid		&
				,DOPR_vid		&
				,DORALG_vid		&
				,RESP_vid		&
				,FTCOD_vid		&
				,COD_vid		&
				,REAERDO_vid	&
			!Added the following for debuguing benthic fluxes
				,BFLUX_vid		&
			!Added the following to save TDIC and TALK fluxes as output
				,DICUPT_vid		&
				,DICBMP_vid		&
				,DICPRD_vid		&
				,DICMNL_vid		&
				,DICDEN_vid		&
				!,DICNIT_GL = DICNIT	&!
				,DICGAS_vid		&
				,DICSED_vid		&
				,ALKNH4_vid		&
				,ALKNO3_vid		&
				,ALKNIT_vid		&
				,ALKDEN_vid		&
				,ALKREM_vid		&
				,ALKNH4SED_vid	&
				,ALKNO3SED_vid	&
			!Added the following to save TDIC and TALK advection + horiz diff AND vert diff fluxes as output
!				,ADVFLUX_vid	&
!				,VDIFFLUX_vid	&
			!If (SEDIMENT_CALC) Then
				,JPOC_vid		&
				,JPON_vid		&
				,JPOP_vid		&
				,JPOS_vid		&
				,CTEMP_vid		&
				!,CH40_GL = 0.d0 &!temporarily set to zero for no CH4 in water column model !0.d0 instead of 0.0 to make it equal to serial run-
				,SODTM1S_vid	&
				,JNH4TM1S_vid	&
				,JNO3TM1S_vid	&
				,BENDEN_vid		&
				,JCH4TM1S_vid	&
				,JCH4GTM1S_vid	&
				,JHSTM1S_vid	&
				,JPO4TM1S_vid	&
				,JSITM1S_vid	&
				,NH41TM1S_vid	&
				,NH42TM1S_vid	&
				,NO31TM1S_vid	&
				,NO32TM1S_vid	&
				,PO41TM1S_vid	&
				,PO42TM1S_vid	&
				,SI1TM1S_vid	&
				,SI2TM1S_vid	&
				,CH41TM1S_vid	&
				,CH42TM1S_vid	&
				,HS1TM1S_vid	&
				,HS2TM1S_vid	&
				,CPOP_vid		&
				,CPON_vid		&
				,CPOC_vid		&
				,CPOS_vid		&
				,HSED1_vid		&
				,BENSTRTM1S_vid &
			! If (SOLIDS_CALC) Then
				,TSS_vid		&
				,SSI_vid		&
				,TURB_vid       &
				,KESS_vid       &
				,PARAD_vid      &
			! end if 
			! If (SAV_CALC) Then
				,LEAF_vid		&
				,STEM_vid		&
				,ROOT_vid       &
				,TUBER_vid      &
				,EP_vid         &
			! end if !AN
			! If (ZOO_CALC) Then
				,SZ_vid		    &
				,LZ_vid
			! end if 


		!ENDIF
  !T,    & !1 water column temperature (degC)
  !SALT, & !2 salinity in water column (psu)
  !SSI,  & !3 inorganic solids (g/m^3)
  !B1,   & !4 algal group 1 biomass concentration (gC/m^3)
  !B2,   & !5 algal group 2 biomass concentration (gC/m^3)
  !B3,   & !6 algal group 3 biomass concentration (gC/m^3)
  !SZ,   & !7 micro(small) zooplankton concentration (gC/m^3)    !
  !LZ,   & !8 meso(large) zooplankton concentration (gC/m^3)     !
  !LDOC, & !9 labile dissolved organic carbon concentration (gC/m^3)
  !RDOC, & !10 refractory dissolved organic carbon concentration (gC/m^3)
  !LPOC, & !11 labile particulate organic carbon concentration (gC/m^3)
  !RPOC, & !12 refractory particulate organic carbon concentration (gC/m^3)
  !NH4,  & !13 ammonia concentration in water column (gN/m^3)
  !NO3,  & !14 nitrate concentration in water column (gN/m^3)
  !UREA, & !15 urea concentration (gN/m^3)    !
  !LDON, & !16 labile dissolved organic nitrogen concentration (gN/m^3)
  !RDON, & !17 labile dissolved organic nitrogen concentration (gN/m^3)
  !LPON, & !18 labile particulate organic nitrogen concentration (gN/m^3)
  !RPON, & !19 refractory particulate organic nitrogen concentration (gN/m^3)
  !PO4,  & !20 phosphate concentation in water column (gP/m^3)
  !LDOP, & !21 labile dissolved organic phosphorus concentration (gP/m^3)
  !RDOP, & !22 refractory dissolved organic phosphorus concentration (gP/m^3)
  !LPOP, & !23 labile particulate organic phosphorus concentration (gP/m^3)
  !RPOP, & !24 refractory particulate organic phosphorus concentration (gP/m^3)
  !PIP,  & !25 particulate inorganic phosphorus (adsorbed to SSI) (gP/m^3)
  !COD,  & !26 chemical oxygen demand in water column (gO2-equiv/m^3)
  !DOXG, & !27 dissolved oxygen in water column (gO2/m^3)
  !SIUPB,& !28 particulate biogenic silica (unavaiable ) concentration in water column (gSi/m^3) 
  !SIAT, & !29 total silica (available) concentration in water column (gSi/m^3) 
  !PIB1, & !Phosphorus in algae 1 (gP/m^3)
  !PIB2, & !Phosphorus in algae 2 (gP/m^3)
  !PIB3    !Phosphorus in algae 3 (gP/m^3)
  !TDIC,  & !33 Dissolved inorg carbon (mmolC/m^3)    				                !
  !TALK, & !34 Total alkalinity (mmol/m^3) 		   		  			!

! Added following for debugging algae
!--Info Variables
   character(len=120) :: institution='DOE/PNNL/MSL/BSRC'
   character(len=120) :: netcdf_timestring = 'seconds after 00:00:00'

#endif /* end of OUTPUT_NETCDF */

   Integer :: NC_FID

Contains


  !functions:
  ! function 	GETDIM()			:   read dimensions
  !
  !subroutines:
  !
  ! subroutine	GETSVAR()			:	read static variable
  ! subroutine	GETDVAR()			:   read time varying dynamic variable
  ! subroutine	PUTSVAR()			:   write static variable
  ! subroutine	PUTDVAR()			:   write time varying dynamic variable

#if defined (OUTPUT_NETCDF)
  ! subroutine  set_ncd_io()        :   read assimilation parameters from input
  ! subroutine  write_netcdf_setup():   set up dimensioning and write grid
  ! subroutine  out_netcdf_data()        :   write time dependent data
  ! subroutine  putvar()            :   collect variable to global and dump
  ! subroutine  handle_ncerr()      :   deal with netcdf error

#endif
    !==============================================================================!
      Integer Function GETDIM (FID, SSIZE, DIMNAME)
    !==============================================================================!
    !  Read dimensions
    !==============================================================================!
    !
         Implicit None
         Integer, Intent (In) :: FID
         Integer, Intent (In) :: SSIZE
         Character (Len=SSIZE), Intent (In) :: DIMNAME
         Integer :: LENGTH
         Integer :: IERR
         Integer :: DIMID
         Character (Len=NF90_MAX_NAME) :: TEMPNAME
    !
         IERR = NF90_INQ_DIMID (FID, TRIM(DIMNAME), DIMID)
         If (IERR /= NF90_NOERR) Then
            Write (*,*) 'Error getting dimension id: ', TRIM (DIMNAME)
            Write (*,*) TRIM (NF90_STRERROR(IERR))
            Stop
         End If
    !
         IERR = NF90_INQUIRE_DIMENSION (FID, DIMID, TEMPNAME, LENGTH)
         If (IERR /= NF90_NOERR) Then
            Write (*,*) 'Error getting dimension: ', TRIM (DIMNAME)
            Write (*,*) TRIM (NF90_STRERROR(IERR))
            Stop
         End If
    !
         GETDIM = LENGTH
    !
      End Function GETDIM
  !==============================================================================!
  !
  !==============================================================================!
  !
      Subroutine GETSVAR (FID, NLEN, VARNAME, I1, I2, TEMP)
    !============================================================================!
    !  Read Static variables
    !==============================================================================!
    !
         Implicit None
         Integer, Intent (In) :: FID
         Integer, Intent (In) :: NLEN
         Character (Len=NLEN), Intent (In) :: VARNAME
         Integer, Intent (In) :: I1, I2
         Real (SP), Intent (Out) :: TEMP (I1, I2)
         Integer :: IERR
         Integer :: VARID
         Integer, Allocatable :: DIMS (:)
    !
         If (I2 == 1) Then
            Allocate (DIMS(1))
            DIMS (1) = 1
         Else
            Allocate (DIMS(2))
            DIMS (1) = 1
            DIMS (2) = 1
         End If
    !
         IERR = NF90_INQ_VARID (FID, TRIM(VARNAME), VARID)
         If (IERR /= NF90_NOERR) Then
            Write (*,*) 'error getting variable id: ', TRIM (VARNAME)
            Write (*,*) TRIM (NF90_STRERROR(IERR))
            Stop
         End If
    !
         IERR = NF90_GET_VAR (FID, VARID, TEMP, DIMS)
         If (IERR /= NF90_NOERR) Then
            Write (*,*) 'error getting variable: ', TRIM (VARNAME)
            Write (*,*) TRIM (NF90_STRERROR(IERR))
            Stop
         End If
    ! Add explicit deallocate
         Deallocate (DIMS)
    !
         Return
      End Subroutine GETSVAR
  !==============================================================================!
  !
  !==============================================================================!
      Subroutine GETDVAR (FID, NLEN, VARNAME, I1, I2, TEMP, NT)
    !============================================================================!
    !  Read time dynamic variables
    !==============================================================================!
    !
         Implicit None
         Integer, Intent (In) :: FID
         Integer, Intent (In) :: NLEN
         Character (Len=NLEN), Intent (In) :: VARNAME
         Integer, Intent (In) :: I1, I2
         Real (CDF_PREC), Intent (Out) :: TEMP (I1, I2)
         Integer :: IERR
         Integer :: VARID
         Integer :: NT
         Integer, Allocatable :: DIMS (:)
    !
         If (I2 == 1) Then
            Allocate (DIMS(2))
            DIMS (1) = 1
            DIMS (2) = NT
         Else
            Allocate (DIMS(3))
            DIMS (1) = 1
            DIMS (2) = 1
            DIMS (3) = NT
         End If
    !
         IERR = NF90_INQ_VARID (FID, TRIM(VARNAME), VARID)
         If (IERR /= NF90_NOERR) Then
            Write (*,*) 'error getting variable id: ', TRIM (VARNAME)
            Write (*,*) TRIM (NF90_STRERROR(IERR))
            Stop
         End If
    !
         IERR = NF90_GET_VAR (FID, VARID, TEMP, DIMS)
         If (IERR /= NF90_NOERR) Then
            Write (*,*) 'error getting variable: ', TRIM (VARNAME)
            Write (*,*) TRIM (NF90_STRERROR(IERR))
            Stop
         End If
		! Add explicit deallocate
         Deallocate (DIMS)
		!
         Return
      End Subroutine GETDVAR
  !==============================================================================!
  !
  !==============================================================================!
      Subroutine PUTDVAR (FID, NLEN, VARNAME, I1, TEMP, NT)
    !============================================================================!
    !  Write dynamic time variables
    !==============================================================================!
    !
         Implicit None
         Integer, Intent (In) :: FID
         Integer, Intent (In) :: NLEN
         Character (Len=NLEN), Intent (In) :: VARNAME
         Integer, Intent (In) :: I1
         Real (SP), Intent (In) :: TEMP (I1)
         Integer :: IERR
         Integer :: VARID
         Integer, Allocatable :: DIMS (:)
         Integer :: NT
    !
         If (I1 == 1) Then
            Allocate (DIMS(1))
            DIMS (1) = NT
         Else
            Allocate (DIMS(2))
            DIMS (1) = 1
            DIMS (2) = NT
         End If
    !
         IERR = NF90_INQ_VARID (FID, TRIM(VARNAME), VARID)
         If (IERR /= NF90_NOERR) Then
            Write (*,*) 'error getting variable id: ', TRIM (VARNAME)
            Write (*,*) TRIM (NF90_STRERROR(IERR))
            Stop
         End If
    !
         IERR = NF90_PUT_VAR (FID, VARID, TEMP, DIMS)
         If (IERR /= NF90_NOERR) Then
            Write (*,*) 'error getting variable: ', TRIM (VARNAME)
            Write (*,*) TRIM (NF90_STRERROR(IERR))
            Stop
         End If
		! Add explicit deallocate
         Deallocate (DIMS)
    !
         Return
      End Subroutine PUTDVAR
  !==============================================================================!
  !
  !==============================================================================!
      Subroutine PUTSVAR (FID, NLEN, VARNAME, I1, TEMP)
    !============================================================================!
    !  Write static variables
    !==============================================================================!
    !
         Implicit None
         Integer, Intent (In) :: FID
         Integer, Intent (In) :: NLEN
         Character (Len=NLEN), Intent (In) :: VARNAME
         Integer, Intent (In) :: I1
         Real (SP), Intent (In) :: TEMP (I1)
         Integer :: IERR
         Integer :: VARID
         Integer, Dimension (1) :: DIMS
    !
         DIMS (1) = 1
    !
         IERR = NF90_INQ_VARID (FID, TRIM(VARNAME), VARID)
         If (IERR /= NF90_NOERR) Then
            Write (*,*) 'error getting variable id: ', TRIM (VARNAME)
            Write (*,*) TRIM (NF90_STRERROR(IERR))
            Stop
         End If
    !
         IERR = NF90_PUT_VAR (FID, VARID, TEMP, DIMS)
         If (IERR /= NF90_NOERR) Then
            Write (*,*) 'error getting variable: ', TRIM (VARNAME)
            Write (*,*) TRIM (NF90_STRERROR(IERR))
            Stop
         End If
    !
         Return
      End Subroutine PUTSVAR

#if defined (OUTPUT_NETCDF)
!==============================================================================|
!  Collect Data to Global Array and Write to Netcdf File                       |
!==============================================================================|

   subroutine putvar(i1,i2,n1,n1gl,kt,k1,map_type,var,nc_fid,vid,myid,nprocs,ipt)

!------------------------------------------------------------------------------|

#  if defined (MULTIPROCESSOR)
   use mod_par, only:emap,nmap ,gather !,bmap
#  endif
   use mod_types
   Use MOD_CONTROL, Only: MSR
   implicit none
   integer, intent(in) :: i1,i2,n1,n1gl,kt,k1,nc_fid,vid,myid,nprocs,ipt
   character(len=*),intent(in)   :: map_type
   real(sp), dimension(i1:i2,kt) :: var

   real(sp), allocatable, dimension(:,:) :: temp,gtemp
   integer :: ierr,k1m1
   integer, allocatable :: dims(:)


   k1m1 = k1
   if(k1m1 == 1)then
     allocate(dims(2))
     dims(1) = 1
     dims(2) = stck_cnt
   else
     allocate(dims(3))
     dims(1) = 1
     dims(2) = 1
     dims(3) = stck_cnt
   end if


   if(map_type(1:1) /= "e" .and. map_type(1:1) /= "n")then 
     write(ipt,*)'map_type input to putvar should be "e" OR "n"'
     call pstop
   end if

   if(nprocs==1)then
     allocate(temp(n1,k1m1))  ; temp(1:n1,1:k1m1) = var(1:n1,1:k1m1)
   end if

#  if defined (MULTIPROCESSOR)

   if(nprocs > 1)then
     allocate(gtemp(n1gl,kt))
     if(map_type(1:1) == "e")then
       call gather(i1,i2,n1,n1gl,kt,myid,nprocs,emap,var,gtemp)
     else if (map_type(1:1) == "n")then
       call gather(i1,i2,n1,n1gl,kt,myid,nprocs,nmap,var,gtemp)
!     else  !type is "b"
!       call gather(i1,i2,n1,n1gl,kt,myid,nprocs,bmap,var,gtemp)
     end if
     allocate(temp(n1gl,k1m1)) ; temp(1:n1gl,1:k1m1) = gtemp(1:n1gl,1:k1m1)
     deallocate(gtemp)
   end if

#  endif

   if(myid /= 1) return

   ierr = nf90_put_var(nc_fid,vid,temp,START=dims)

   if(ierr /= nf90_noerr)then
     IF(msr) THEN
       WRITE(*,*)'Enter here....',vid 
     END IF
     call handle_ncerr(ierr,"error writing variable to netcdf",ipt)
   end if
   deallocate(dims)

   return
   end subroutine putvar


!------------------------------------------------------------------------------|
!  CHECK NETCDF ERROR AND WRITE MESSAGE TO FILE UNIT IPT                       |
!------------------------------------------------------------------------------|
   SUBROUTINE handle_ncerr(status,errmsge,ipt)
   integer, intent(in) :: status,ipt
   character(len=*)    :: errmsge
   if(status /=nf90_noerr)then
     write(ipt,*)trim(errmsge)
     write(ipt,*)trim(nf90_strerror(status))
     call pstop
   end if
   END SUBROUTINE handle_ncerr

!==============================================================================|


!==============================================================================|
!  READ IN PARAMETERS CONTROLLING INPUT/OUTPUT FROM RUNTIME PARAM FILE         |
!==============================================================================|
   SUBROUTINE set_ncd_io(ipt)
   use mod_prec, only:
   use mod_input, only: scan_file
   use netcdf
   use mod_control, only: casename, MSR
   use mod_fileinfo, only: inctr
   implicit none
   integer, intent(in) ::ipt
!--Local Vars----------------------------|
   real(sp)           :: realvec(150)
   integer            :: intvec(150)
   integer            :: iscan
   character(len=120) :: fname
   character(len=256), dimension(400) :: charvec
   integer            :: i

!----------------------------------------|


   out_cnt = 0

   fname = "./"//trim(casename)//"_run.dat"

!------------------------------------------------------------------------------|
!   cdf_out: netcdf activation flag
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_OUT",LVAL = CDF_OUT, IUNIT=INCTR)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_out: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if

!------------------------------------------------------------------------------|
!  cdf_int: output is performed every cdf_int iterations
!------------------------------------------------------------------------------|

   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_INT",ISCAL = CDF_INT, IUNIT=INCTR)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_int: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if

!------------------------------------------------------------------------------|
!  cdf_stk: # dumps / file
!------------------------------------------------------------------------------|

   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_STK",ISCAL = CDF_STK, IUNIT=INCTR)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_stk: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if


!------------------------------------------------------------------------------|
!     cdf_vdp: list of variables to write to output file
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_VDP",CVEC = CHARVEC,NSZE = NOUT_VARS, IUNIT=INCTR)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_vdp: ',iscan
     call pstop
   end if
   if(nout_vars <= 0)then
     write(ipt,*)'incorrect number of netcdf cdf_vdp variables specified'
     write(ipt,*)'in input file',nout_vars
     call pstop
   end if

 ! CDF_VDP  = u v wts kh el salt temp d CCHL1 total_netPP DOXG LDOC RDOC LPOC RPOC  B1 B2 NH4\\
 !NO3 LDON RDON LPON RPON PO4 TDIC TALK pH pCO2 P1 P2 BM1 BM2 NL1 NL2 PL1 PL2 B1SZ B2SZ B1LZ\\
 !B2LZ PR1 PR2 IAVG FI1 FI2 IK1 IK2 GPP NPP DDOC DCOD NITRIF DOSZ DOLZ DOPR DORALG RESP FTCOD\\
 !COD REAERDO BFLUX DICUPT DICBMP DICPRD DICMNL DICDEN DICGAS DICSED ALKNH4 ALKNO3 ALKNIT\\
 !ALKDEN ALKREM ALKNH4SED ALKNO3SED JPOC JPON JPOP JPOS CTEMP SODTM1S JNH4TM1S JNO3TM1S BENDEN\\
 !JCH4TM1S JCH4GTM1S JHSTM1S JPO4TM1S JSITM1S NH41TM1S NH42TM1S NO31TM1S NO32TM1S PO41TM1S\\
 !PO42TM1S SI1TM1S SI2TM1S CH41TM1S CH42TM1S HS1TM1S HS2TM1S CPOP CPON CPOC CPOS HSED1\\
 !BENSTRTM1S

  allocate(cdf_vdp(nout_vars))
  cdf_vdp(1:nout_vars)= charvec(1:nout_vars)


!------------------------------------------------------------------------------|
!            SCREEN REPORT OF SET VARIABLES                                    !
!------------------------------------------------------------------------------|
   if(msr)then
     write(ipt,*)''
     write(ipt,*)'!        netcdf parameters                  '
     if(cdf_out)then
       write(ipt,*)'!  output netcdf         :  active'
       write(ipt,*)'!  output every # its    : ',cdf_int
       write(ipt,*)'!  # dumps / file        : ',cdf_stk
       write(ipt,*)'!  # variables to write  : ',nout_vars
       do i=1,nout_vars
         write(ipt,999)i,trim(cdf_vdp(i))
       end do
     else
       write(ipt,*)'!  # output netcdf       :  disabled'
     end if
   end if


   return
   999 format(' !  variable #',i4,'        :',a13)
   END SUBROUTINE set_ncd_io
!===================================================================


!==============================================================================|
!  Write NetCDF Header and Static Variables                                    |
!==============================================================================|
   SUBROUTINE write_netcdf_setup(filecnt,ipt)

#  if defined (MULTIPROCESSOR)
   use mod_par , only:el_pid,ngid
#  endif
   use netcdf


   !use mod_types, only:
   use mod_utils, only:get_timestamp
   use mod_hydrovars, only: his_outdir,stn_outdir,vx,vy,vxmin, vymin,z,zz,uu,vv,wts,h, d, Z2D,SIGVAR
   use mod_control, only: casetitle, MSR, PAR, CASENAME, SERIAL
   use mod_sizes,only: MGL, NGL, NSAVM
   use mod_lims, only: KB, KBM1, MLOC, NLOC, myid,nprocs, NTLOC,MTLOC,IINT
   use mod_tge, only: nv

   use mod_wqm, only: SEDIMENT_CALC, SOLIDS_CALC, SAV_CALC

   USE MOD_WQMINIT, ONLY: CARBONATE_CALC, COD_CALC, ZOO_CALC
   implicit none
   integer, intent(in)   :: filecnt,ipt
   integer, dimension(3) :: dynm3de_lev,dynm3de_lay,dynm3dobc_lay
   integer, dimension(3) :: dynm3dn_lev,dynm3dn_lay,dynm3dobc_lay2,dynm3dn_3,dynm3dn_9, &
  & dynm3dn_nsav
   integer, dimension(2) :: stat3de_lev,stat3de_lay
   integer, dimension(2) :: stat3dn_lev,stat3dn_lay
   integer, dimension(2) :: specdim
   integer, dimension(2) :: dynm2de,dynm2dn,dynm2dobc
   integer, dimension(1) :: stat2de,stat2dn
   integer, dimension(1) :: stat_lev,stat_lay,dynmtime ,stat_scl
   character(len=100)    :: netcdf_convention
   character(len=100)    :: timestamp ,temp
   integer               :: i,j,ierr,i1,i2
   integer               :: maxnode,maxnodep,maxelem,maxelemp,itmp
   real(sp), allocatable :: tmp(:,:),tvec(:)
   character(len=4)      :: nchar

!==============================================================================|

!==============================================================================|
!  Set up Constants and Initialize Counters                                    |
!==============================================================================|

!--Initialize Stack Count
   stck_cnt = 1

!--NetCDF Convention String
   netcdf_convention = 'CF-1.0'

!--Time Stamp for History
   call get_timestamp(temp)
   timestamp = 'model started at: '//trim(temp)

!==============================================================================|
!  OPEN FILE AND DEFINE VARIABLES                                              |
!==============================================================================|
   if(msr)then

!--Define NetCDF Output Filename
   write(nchar,'(I4)')filecnt
   if(filecnt < 10)then
     cdfname = trim(his_outdir)//trim(casename)//'_FVCOMICM_0000'//trim(adjustl(nchar))//'.nc'
   elseif(filecnt < 100)then
     cdfname = trim(his_outdir)//trim(casename)//'_FVCOMICM_000'//trim(adjustl(nchar))//'.nc'
   elseif(filecnt < 1000)then
     cdfname = trim(his_outdir)//trim(casename)//'_FVCOMICM_00'//trim(adjustl(nchar))//'.nc'
   elseif(filecnt < 10000)then
     cdfname = trim(his_outdir)//trim(casename)//'_FVCOMICM_0'//trim(adjustl(nchar))//'.nc'
    elseif(filecnt < 100000)then
     cdfname = trim(his_outdir)//trim(casename)//'_FVCOMICM_'//trim(adjustl(nchar))//'.nc'
   else
     write(ipt,*)'error in netcdf module'
     write(ipt,*)'# history files > 100000'
     stop
   endif

!--Create File
   ierr = nf90_create(path=trim(cdfname),cmode=nf90_clobber,ncid=nc_ofid)
   if(ierr /= nf90_eexist)then
     call handle_ncerr(ierr,"file creation error",ipt)
   else
     write(ipt,*)'file :',trim(cdfname),' already exists'
     write(ipt,*)'exiting'
     stop
   end if

!--Description of File Contents
   ierr = nf90_put_att(nc_ofid,nf90_global,"title"      ,trim(casetitle))
   ierr = nf90_put_att(nc_ofid,nf90_global,"institution",trim(institution))
   ierr = nf90_put_att(nc_ofid,nf90_global,"source"     ,trim('FVCOM-ICM 3.0'))
   ierr = nf90_put_att(nc_ofid,nf90_global,"history"    ,trim(timestamp))
   !ierr = nf90_put_att(nc_ofid,nf90_global,"references" ,trim(fvcom_website))
   ierr = nf90_put_att(nc_ofid,nf90_global,"Conventions",trim(netcdf_convention))

!--Define Fixed Model Dimensions
   ierr = nf90_def_dim(nc_ofid,"scalar" ,1      ,scl_did    )
   ierr = nf90_def_dim(nc_ofid,"node"   ,mgl    ,node_did   )
   ierr = nf90_def_dim(nc_ofid,"nele"   ,ngl    ,nele_did   )
   ierr = nf90_def_dim(nc_ofid,"siglay" ,kbm1   ,siglay_did )
   ierr = nf90_def_dim(nc_ofid,"siglev" ,kb     ,siglev_did )
   ierr = nf90_def_dim(nc_ofid,"three"  ,3      ,three_did  )
   ierr = nf90_def_dim(nc_ofid,"four"   ,4      ,four_did   )
   ierr = nf90_def_dim(nc_ofid,"nine"   ,9      ,nine_did   )    !
!   ierr = nf90_def_dim(nc_ofid,"obc"    ,IOBCN_RGL,obc_did )
!   ierr = nf90_def_dim(nc_ofid,"obc2"   ,IOBCN_RGL,obc2_did  )
!   Print*, 'IOBCN_GL', IOBCN_RGL
!--Define Unlimited Model Dimension
   ierr = nf90_def_dim(nc_ofid,"time"   ,nf90_unlimited,time_did)
   ierr = nf90_def_dim(nc_ofid,"nsav"   ,nsavm      ,nsav_did   ) ! AN

!--Set Up Data Dimensioning - Static Vars
   stat_scl     = (/scl_did/)             !!scalar variable
   stat_lay     = (/siglay_did/)          !!vertical variables at layers
   stat_lev     = (/siglev_did/)          !!vertical variables at levels
   stat2de      = (/nele_did/)            !!2d element vars
   stat2dn      = (/node_did/)            !!2d nodal vars
   stat3de_lay  = (/nele_did,siglay_did/) !!3d element vars at layers
   stat3de_lev  = (/nele_did,siglev_did/) !!3d element vars at levels
   stat3dn_lay  = (/node_did,siglay_did/) !!3d node    vars at layers
   stat3dn_lev  = (/node_did,siglev_did/) !!3d node    vars at levels

!--Set Up Data Dimensioning - Dynamic Vars
   dynm2de      = (/nele_did,time_did/)            !!2d element vars
   dynm2dn      = (/node_did,time_did/)            !!2d nodal vars
   dynm2dobc    = (/obc_did,time_did/)             !!2d obc vars
   dynm3de_lay  = (/nele_did,siglay_did,time_did/) !!3d elem vars at layers
   dynm3de_lev  = (/nele_did,siglev_did,time_did/) !!3d elem vars at levels
   dynm3dn_lay  = (/node_did,siglay_did,time_did/) !!3d node vars at layers
   dynm3dn_lev  = (/node_did,siglev_did,time_did/) !!3d node vars at levels
   dynm3dobc_lay  = (/obc_did,siglay_did,time_did/) !!3d node vars at layers, bound IOBCN_GL+1
   dynm3dobc_lay2 = (/obc2_did,siglay_did,time_did/) !!3d node vars at layers, bound IOBCN_GL
   dynmtime     = (/time_did/)

   dynm3dn_3  = (/node_did,three_did,time_did/) !!3d node vars at layers
   dynm3dn_9  = (/node_did,nine_did,time_did/) !!3d node vars at layers
   dynm3dn_nsav= (/node_did,nsav_did,time_did/) !!3d node vars at layers

!--Define Coordinate Variables and Attributes

   !!====NPROCS: Number of Processors=======================!
   ierr = nf90_def_var(nc_ofid,"nprocs",nf90_int,stat_scl,nprocs_vid)
   ierr = nf90_put_att(nc_ofid,nprocs_vid,"long_name","number of processors")

   !!====PARTITION: Partion Number of Element===============!
   ierr = nf90_def_var(nc_ofid,"partition",nf90_int,stat2de,partition_vid)
   ierr = nf90_put_att(nc_ofid,partition_vid,"long_name","partition")

   !!====Initial Density (Used for Constructing 3D Domain)==!
   ierr = nf90_def_var(nc_ofid,"Initial_Density",nf90_float,stat3dn_lay,idens_vid)
   ierr = nf90_put_att(nc_ofid,idens_vid,"long_name","Initial Density")

   !!====X Grid Coordinate at Nodes (VX) (Meters)===========!
   ierr = nf90_def_var(nc_ofid,"x",nf90_float,stat2dn,x_vid)
   ierr = nf90_put_att(nc_ofid,x_vid,"long_name","nodal x-coordinate")
   ierr = nf90_put_att(nc_ofid,x_vid,"units","meters")

   !!====Y Grid Coordinate at Nodes (VY) (Meters)===========!
   ierr = nf90_def_var(nc_ofid,"y",nf90_float,stat2dn,y_vid)
   ierr = nf90_put_att(nc_ofid,y_vid,"long_name","nodal y-coordinate")
   ierr = nf90_put_att(nc_ofid,y_vid,"units","meters")


   !!!====Longitudinal Coordinate at Nodes (LON) (degrees)===!
   !ierr = nf90_def_var(nc_ofid,"lon",nf90_float,stat2dn,lon_vid)
   !ierr = nf90_put_att(nc_ofid,lon_vid,"long_name","Longitude")
   !ierr = nf90_put_att(nc_ofid,lon_vid,"standard_name","longitude")
   !ierr = nf90_put_att(nc_ofid,lon_vid,"units","degrees_east")

   !!!====Latitudinal  Coordinate at Nodes (LAT) (degrees)===!
   !ierr = nf90_def_var(nc_ofid,"lat",nf90_float,stat2dn,lat_vid)
   !ierr = nf90_put_att(nc_ofid,lat_vid,"long_name","Latitude")
   !ierr = nf90_put_att(nc_ofid,lat_vid,"standard_name","latitude")
   !ierr = nf90_put_att(nc_ofid,lat_vid,"units","degrees_north")

   !!====Sigma Coordinate for Sigma Layers (ZZ)  (-)========!
   ierr = nf90_def_var(nc_ofid,"siglay",nf90_float,stat_lay,siglay_vid)
   ierr = nf90_put_att(nc_ofid,siglay_vid,"long_name","Sigma Layers")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"standard_name","ocean_sigma_coordinate")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"positive","up")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"valid_min","-1")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"valid_max","0")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"formula_terms","siglay:siglay eta:zeta depth:depth")

   !!====Shifted Sigma Layer Coordinate for Viz ============!
!   ierr = nf90_def_var(nc_ofid,"siglay_shift",nf90_float,stat_lay,siglay_shift_vid)
!   ierr = nf90_put_att(nc_ofid,siglay_shift_vid,"long_name","Shifted Sigma Layers")

   !!====Sigma Coordinate for Sigma Levels (Z)   (-)========!
   ierr = nf90_def_var(nc_ofid,"siglev",nf90_float,stat_lev,siglev_vid)
   ierr = nf90_put_att(nc_ofid,siglev_vid,"long_name","Sigma Levels")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"standard_name","ocean_sigma_coordinate")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"positive","up")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"valid_min","-1")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"valid_max","0")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"formula_terms","siglev:siglev eta:zeta depth:depth")

#if defined (NEWSIGMA)

   ierr = nf90_def_var(nc_ofid,"siglev3dn",nf90_float,stat3dn_lev,siglev3dn_vid)
   ierr = nf90_put_att(nc_ofid,siglev3dn_vid,"long_name","Sigma Levels")
   ierr = nf90_put_att(nc_ofid,siglev3dn_vid,"standard_name","ocean_sigma_coordinate")
   ierr = nf90_put_att(nc_ofid,siglev3dn_vid,"positive","up")
   ierr = nf90_put_att(nc_ofid,siglev3dn_vid,"valid_min","-1")
   ierr = nf90_put_att(nc_ofid,siglev3dn_vid,"valid_max","0")
   ierr = nf90_put_att(nc_ofid,siglev3dn_vid,"formula_terms","siglev:siglev eta:zeta depth:depth")

#endif

!--Define Mesh Relevant Variables and Attributes

   !!====Bathymetry at Nodes (H) (meters)===================!
   ierr = nf90_def_var(nc_ofid,"h",nf90_float,stat2dn,h_vid)
   ierr = nf90_put_att(nc_ofid,h_vid,"long_name","Bathymetry")
   ierr = nf90_put_att(nc_ofid,h_vid,"units","meters")
   ierr = nf90_put_att(nc_ofid,h_vid,"positive","down")
   ierr = nf90_put_att(nc_ofid,h_vid,"standard_name","depth")
   ierr = nf90_put_att(nc_ofid,h_vid,"grid","fvcom_grid")

   !!====Nodes surrounding each Element (NV)================!
   specdim = (/nele_did,three_did/)
   ierr = nf90_def_var(nc_ofid,"nv",nf90_float,specdim,nv_vid)
   ierr = nf90_put_att(nc_ofid,nv_vid,"long_name","nodes surrounding element")

   !!!====Momentum Stencil Interpolation Coefficients========!
   !specdim = (/nele_did,four_did/)
   !ierr = nf90_def_var(nc_ofid,"a1u",nf90_float,specdim,a1u_vid)
   !ierr = nf90_put_att(nc_ofid,a1u_vid,"long_name","a1u")
   !ierr = nf90_def_var(nc_ofid,"a2u",nf90_float,specdim,a2u_vid)
   !ierr = nf90_put_att(nc_ofid,a2u_vid,"long_name","a2u")

   !!!====Element Based Interpolation Coefficients===========!
   !specdim = (/nele_did,three_did/)
   !ierr = nf90_def_var(nc_ofid,"aw0",nf90_float,specdim,aw0_vid)
   !ierr = nf90_put_att(nc_ofid,aw0_vid,"long_name","aw0")
   !ierr = nf90_def_var(nc_ofid,"awx",nf90_float,specdim,awx_vid)
   !ierr = nf90_put_att(nc_ofid,awx_vid,"long_name","awx")
   !ierr = nf90_def_var(nc_ofid,"awy",nf90_float,specdim,awy_vid)
   !ierr = nf90_put_att(nc_ofid,awy_vid,"long_name","awy")

!--Define Model Time Variables and Attributes
   ierr = nf90_def_var(nc_ofid,"time",nf90_float,dynmtime,time_vid)
   ierr = nf90_put_att(nc_ofid,time_vid,"long_name","Time")
   ierr = nf90_put_att(nc_ofid,time_vid,"units",trim(netcdf_timestring))
   ierr = nf90_put_att(nc_ofid,time_vid,"calendar","none")
   ierr = nf90_def_var(nc_ofid,"iint",nf90_int,dynmtime,iint_vid)
   ierr = nf90_put_att(nc_ofid,iint_vid,"long_name","internal mode iteration number")

!--Define Time Dependent Flow Variables (selected by user from input file)
   do i=1,nout_vars

     select case(trim(cdf_vdp(i)))

     case("u")  !!===============u=======================================!
     ierr = nf90_def_var(nc_ofid,"u",nf90_float,dynm3de_lay,u_vid)
     ierr = nf90_put_att(nc_ofid,u_vid,"long_name","Eastward Water Velocity")
     ierr = nf90_put_att(nc_ofid,u_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,u_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,u_vid,"type","data")

     case("v")  !!===============v=======================================!
     ierr = nf90_def_var(nc_ofid,"v",nf90_float,dynm3de_lay,v_vid)
     ierr = nf90_put_att(nc_ofid,v_vid,"long_name","Northward Water Velocity")
     ierr = nf90_put_att(nc_ofid,v_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,v_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,v_vid,"type","data")

     !case("w") !!===============w=========================================!
     !ierr = nf90_def_var(nc_ofid,"w",nf90_float,dynm3de_lev,w_vid)
     !ierr = nf90_put_att(nc_ofid,w_vid,"long_name","Pseudo Vertical Velocity")
     !ierr = nf90_put_att(nc_ofid,w_vid,"units","meters s-1")
     !ierr = nf90_put_att(nc_ofid,w_vid,"grid","fvcom_grid")
     !ierr = nf90_put_att(nc_ofid,w_vid,"type","data")

     !case("ww") !!===============ww======================================!
     !ierr = nf90_def_var(nc_ofid,"ww",nf90_float,dynm3de_lay,ww_vid)
     !ierr = nf90_put_att(nc_ofid,ww_vid,"long_name","Upward Water Velocity")
     !ierr = nf90_put_att(nc_ofid,ww_vid,"units","meters s-1")
     !ierr = nf90_put_att(nc_ofid,ww_vid,"grid","fvcom_grid")
     !ierr = nf90_put_att(nc_ofid,ww_vid,"type","data")

     case("wts") !!===============wts=====================================!
     ierr = nf90_def_var(nc_ofid,"wts",nf90_float,dynm3dn_lev,wts_vid)
     ierr = nf90_put_att(nc_ofid,wts_vid,"long_name","Upward Water Velocity at node")
     ierr = nf90_put_att(nc_ofid,wts_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,wts_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,wts_vid,"type","data")


     !case("dtfa") !!===============dtfa===================================!
     !ierr = nf90_def_var(nc_ofid,"dtfa",nf90_float,dynm2dn,dtfa_vid)
     !ierr = nf90_put_att(nc_ofid,dtfa_vid,"long_name","FSH + Water Height")
     !ierr = nf90_put_att(nc_ofid,dtfa_vid,"units","meters")
     !ierr = nf90_put_att(nc_ofid,dtfa_vid,"positive","down")
     !ierr = nf90_put_att(nc_ofid,dtfa_vid,"grid","fvcom_grid")
     !ierr = nf90_put_att(nc_ofid,dtfa_vid,"type","data")

!     case("xflux_obc") !!===============xfluxobc==========================!
!     ierr = nf90_def_var(nc_ofid,"xflux_obc",nf90_float,dynm3dobc_lay2,xflux_obc_vid)
!     ierr = nf90_put_att(nc_ofid,xflux_obc_vid,"long_name","Xflux at OBC")
!     ierr = nf90_put_att(nc_ofid,xflux_obc_vid,"units","tracer unit * meters+3 s-1 ")
!     ierr = nf90_put_att(nc_ofid,xflux_obc_vid,"grid","fvcom_grid")
!     ierr = nf90_put_att(nc_ofid,xflux_obc_vid,"type","data")

!!
!     case("uard_obcn") !!===============uardobcn============================!
!     ierr = nf90_def_var(nc_ofid,"uard_obcn",nf90_float,dynm2dobc,uard_obcn_vid)
!     ierr = nf90_put_att(nc_ofid,uard_obcn_vid,"long_name","UARD at OBC")
!     ierr = nf90_put_att(nc_ofid,uard_obcn_vid,"units","meters s-1")
!     ierr = nf90_put_att(nc_ofid,uard_obcn_vid,"grid","fvcom_grid")
!     ierr = nf90_put_att(nc_ofid,uard_obcn_vid,"type","data")

     !case("km") !!===============km======================================!
     !ierr = nf90_def_var(nc_ofid,"km",nf90_float,dynm3dn_lev,km_vid)
     !ierr = nf90_put_att(nc_ofid,km_vid,"long_name","Turbulent Eddy Viscosity")
     !ierr = nf90_put_att(nc_ofid,km_vid,"units","meters2 s-1")
     !ierr = nf90_put_att(nc_ofid,km_vid,"grid","fvcom_grid")
     !ierr = nf90_put_att(nc_ofid,km_vid,"type","data")

     case("kh") !!===============kh======================================!
     ierr = nf90_def_var(nc_ofid,"kh",nf90_float,dynm3dn_lev,kh_vid)
     ierr = nf90_put_att(nc_ofid,kh_vid,"long_name","Turbulent Eddy Diffusivity")
     ierr = nf90_put_att(nc_ofid,kh_vid,"units","meters2 s-1")
     ierr = nf90_put_att(nc_ofid,kh_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,kh_vid,"type","data")

     !case("viscofh") !!===============viscofh======================================!
     !ierr = nf90_def_var(nc_ofid,"viscofh",nf90_float,dynm3dn_lay,viscofh_vid)
     !ierr = nf90_put_att(nc_ofid,viscofh_vid,"long_name","Giruzibtak Eddy Diffusivity")
     !ierr = nf90_put_att(nc_ofid,viscofh_vid,"units","meters2 s-1")
     !ierr = nf90_put_att(nc_ofid,viscofh_vid,"grid","fvcom_grid")
     !ierr = nf90_put_att(nc_ofid,viscofh_vid,"type","data")


     case("temp") !!===============t1======================================!
     ierr = nf90_def_var(nc_ofid,"temp",nf90_float,dynm3dn_lay,t1_vid)
     ierr = nf90_put_att(nc_ofid,t1_vid,"long_name","temperature")
     ierr = nf90_put_att(nc_ofid,t1_vid,"standard_name","sea_water_temperature")
     ierr = nf90_put_att(nc_ofid,t1_vid,"units","degrees_C")
     ierr = nf90_put_att(nc_ofid,t1_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,t1_vid,"type","data")

     case("salt") !!===============s1======================================!
     ierr = nf90_def_var(nc_ofid,"salinity",nf90_float,dynm3dn_lay,s1_vid)
     ierr = nf90_put_att(nc_ofid,s1_vid,"long_name","salinity")
     ierr = nf90_put_att(nc_ofid,s1_vid,"standard_name","sea_water_salinity")
     ierr = nf90_put_att(nc_ofid,s1_vid,"units","1e-3")
     ierr = nf90_put_att(nc_ofid,s1_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,s1_vid,"type","data")

     case("el") !!===============el======================================!
     ierr = nf90_def_var(nc_ofid,"zeta",nf90_float,dynm2dn,el_vid)
     ierr = nf90_put_att(nc_ofid,el_vid,"long_name","Water Surface Elevation")
     ierr = nf90_put_att(nc_ofid,el_vid,"units","meters")
     ierr = nf90_put_att(nc_ofid,el_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,el_vid,"standard_name","sea_surface_elevation")
     ierr = nf90_put_att(nc_ofid,el_vid,"type","data")

     case("d") !!===============d=======================================!
     ierr = nf90_def_var(nc_ofid,"depth",nf90_float,dynm2dn,d_vid)
     ierr = nf90_put_att(nc_ofid,d_vid,"long_name","Water Depth")
     ierr = nf90_put_att(nc_ofid,d_vid,"units","meters")
     ierr = nf90_put_att(nc_ofid,d_vid,"positive","down")
     ierr = nf90_put_att(nc_ofid,d_vid,"type","data")

     !case("ua") !!===============ua======================================!
     !ierr = nf90_def_var(nc_ofid,"ua",nf90_float,dynm2de,ua_vid)
     !ierr = nf90_put_att(nc_ofid,ua_vid,"long_name","Vertically Averaged x-velocity")
     !ierr = nf90_put_att(nc_ofid,ua_vid,"units","meters s-1")
     !ierr = nf90_put_att(nc_ofid,ua_vid,"type","data")

     !case("va") !!===============va======================================!
     !ierr = nf90_def_var(nc_ofid,"va",nf90_float,dynm2de,va_vid)
     !ierr = nf90_put_att(nc_ofid,va_vid,"long_name","Vertically Averaged y-velocity")
     !ierr = nf90_put_att(nc_ofid,va_vid,"units","meters s-1")
     !ierr = nf90_put_att(nc_ofid,va_vid,"type","data")

#if defined (WET_DRY)
     case("wd") !!===============WET DRY FLAG============================!
     ierr = nf90_def_var(nc_ofid,"wd",nf90_float,dynm2dn,wd_vid)
     ierr = nf90_put_att(nc_ofid,wd_vid,"long_name","Wet Dry Flag")
	 ierr = nf90_put_att(nc_ofid,wd_vid,"standard_name","Wet_Dry_Flag")
     ierr = nf90_put_att(nc_ofid,wd_vid,"units","-")
     ierr = nf90_put_att(nc_ofid,wd_vid,"type","data")
#endif

!example  case("var") !!===============var====================================!
!example  ierr = nf90_def_var(nc_ofid,"truevar",nf90_float,dimensions,var_vid)
!example  ierr = nf90_put_att(nc_ofid,var_vid,"long_name","A Good Descriptive Name")
!example  ierr = nf90_put_att(nc_ofid,var_vid,"units","UDUNITS compatible units")
!example  ierr = nf90_put_att(nc_ofid,var_vid,"standard_name","CF-convention standard name")
!example  ierr = nf90_put_att(nc_ofid,var_vid,"type","data")

     !    new variable definition
     !1.) add new definition above here by copying example above and modifying
     !2.) copy dimensions from variable which has same dimensions as var
     !3.) change variable name if necessary to something more descriptive
     !   e.g. model name for temperature is t1, use temp instead
     !4.) give the variable a reasonable "long_name"
     !5.) look up the variables standard_name from the cf-convention standard_name list
     !   http://www.cgd.ucar.edu/cms/eaton/cf-metadata/standard_name.html
     !   if it does not exist, do not provide a standard name attribute
     !6.) set variable units conforming to udunits standard
     !   http://my.unidata.ucar.edu/content/software/udunits/index.html


		case ("CCHL1")

			ierr = nf90_def_var(nc_ofid,"CCHL1",nf90_float,dynm3dn_lay,CCHL1_vid)
			ierr = nf90_put_att(nc_ofid,CCHL1_vid,"long_name","C to CHL ratio")
			ierr = nf90_put_att(nc_ofid,CCHL1_vid,"standard_name","Carbn_to_CHLA_ratio")
			ierr = nf90_put_att(nc_ofid,CCHL1_vid,"units","gC/gCHLA")
			ierr = nf90_put_att(nc_ofid,CCHL1_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,CCHL1_vid,"type","data")

		case ("total_netPP")

			ierr = nf90_def_var(nc_ofid,"total_netPP",nf90_float,dynm2dn,total_netPP_vid)
			ierr = nf90_put_att(nc_ofid,total_netPP_vid,"long_name","total_netPP")
			ierr = nf90_put_att(nc_ofid,total_netPP_vid,"standard_name","total_netPP")
			ierr = nf90_put_att(nc_ofid,total_netPP_vid,"units","gC/m^2/day")
			ierr = nf90_put_att(nc_ofid,total_netPP_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,total_netPP_vid,"type","data")

		case ( "DOXG")

			ierr = nf90_def_var(nc_ofid,"DOXG",nf90_float,dynm3dn_lay,DOXG_vid)
			ierr = nf90_put_att(nc_ofid,DOXG_vid,"long_name","Dissolved Oxygen")
			ierr = nf90_put_att(nc_ofid,DOXG_vid,"standard_name","DOXG")
			ierr = nf90_put_att(nc_ofid,DOXG_vid,"units","MG/L")
			ierr = nf90_put_att(nc_ofid,DOXG_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DOXG_vid,"type","data")

		case (  "LDOC")
			ierr = nf90_def_var(nc_ofid,"LDOC",nf90_float,dynm3dn_lay,LDOC_vid)
			ierr = nf90_put_att(nc_ofid,LDOC_vid,"long_name","LDOC")
			ierr = nf90_put_att(nc_ofid,LDOC_vid,"standard_name","LDOC")
			ierr = nf90_put_att(nc_ofid,LDOC_vid,"units","gC meters-3")
			ierr = nf90_put_att(nc_ofid,LDOC_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,LDOC_vid,"type","data")


		case ( "RDOC")
			ierr = nf90_def_var(nc_ofid,"RDOC",nf90_float,dynm3dn_lay,RDOC_vid)
			ierr = nf90_put_att(nc_ofid,RDOC_vid,"long_name","RDOC")
			ierr = nf90_put_att(nc_ofid,RDOC_vid,"standard_name","RDOC")
			ierr = nf90_put_att(nc_ofid,RDOC_vid,"units","gC meters-3")
			ierr = nf90_put_att(nc_ofid,RDOC_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,RDOC_vid,"type","data")


		case ( "LPOC")
			ierr = nf90_def_var(nc_ofid,"LPOC",nf90_float,dynm3dn_lay,LPOC_vid)
			ierr = nf90_put_att(nc_ofid,LPOC_vid,"long_name","LPOC")
			ierr = nf90_put_att(nc_ofid,LPOC_vid,"standard_name","LPOC")
			ierr = nf90_put_att(nc_ofid,LPOC_vid,"units","MG-C/L")
			ierr = nf90_put_att(nc_ofid,LPOC_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,LPOC_vid,"type","data")

		case ( "RPOC")
			ierr = nf90_def_var(nc_ofid,"RPOC",nf90_float,dynm3dn_lay,RPOC_vid)
			ierr = nf90_put_att(nc_ofid,RPOC_vid,"long_name","RPOC")
			ierr = nf90_put_att(nc_ofid,RPOC_vid,"standard_name","RPOC")
			ierr = nf90_put_att(nc_ofid,RPOC_vid,"units","gC meters-3")
			ierr = nf90_put_att(nc_ofid,RPOC_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,RPOC_vid,"type","data")

		case ( "B1")
			ierr = nf90_def_var(nc_ofid,"B1",nf90_float,dynm3dn_lay,B1_vid)
			ierr = nf90_put_att(nc_ofid,B1_vid,"long_name","B1")
			ierr = nf90_put_att(nc_ofid,B1_vid,"standard_name","B1")
			ierr = nf90_put_att(nc_ofid,B1_vid,"units","gC meters-3")
			ierr = nf90_put_att(nc_ofid,B1_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,B1_vid,"type","data")

		case ( "B2")
			ierr = nf90_def_var(nc_ofid,"B2",nf90_float,dynm3dn_lay,B2_vid)
			ierr = nf90_put_att(nc_ofid,B2_vid,"long_name","B2")
			ierr = nf90_put_att(nc_ofid,B2_vid,"standard_name","B2")
			ierr = nf90_put_att(nc_ofid,B2_vid,"units","gC meters-3")
			ierr = nf90_put_att(nc_ofid,B2_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,B2_vid,"type","data")

		case ( "NH4")
			ierr = nf90_def_var(nc_ofid,"NH4",nf90_float,dynm3dn_lay,NH4_vid)
			ierr = nf90_put_att(nc_ofid,NH4_vid,"long_name","NH4")
			ierr = nf90_put_att(nc_ofid,NH4_vid,"standard_name","NH4")
			ierr = nf90_put_att(nc_ofid,NH4_vid,"units","gN meters-3")
			ierr = nf90_put_att(nc_ofid,NH4_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,NH4_vid,"type","data")

		case ( "NO3")
			ierr = nf90_def_var(nc_ofid,"NO3",nf90_float,dynm3dn_lay,NO3_vid)
			ierr = nf90_put_att(nc_ofid,NO3_vid,"long_name","NO3")
			ierr = nf90_put_att(nc_ofid,NO3_vid,"standard_name","NO3")
			ierr = nf90_put_att(nc_ofid,NO3_vid,"units","gN meters-3")
			ierr = nf90_put_att(nc_ofid,NO3_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,NO3_vid,"type","data")

		case ( "LDON")
			ierr = nf90_def_var(nc_ofid,"LDON",nf90_float,dynm3dn_lay,LDON_vid)
			ierr = nf90_put_att(nc_ofid,LDON_vid,"long_name","LDON")
			ierr = nf90_put_att(nc_ofid,LDON_vid,"standard_name","LDON")
			ierr = nf90_put_att(nc_ofid,LDON_vid,"units","gN meters-3")
			ierr = nf90_put_att(nc_ofid,LDON_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,LDON_vid,"type","data")


		case ( "RDON")
			ierr = nf90_def_var(nc_ofid,"RDON",nf90_float,dynm3dn_lay,RDON_vid)
			ierr = nf90_put_att(nc_ofid,RDON_vid,"long_name","RDON")
			ierr = nf90_put_att(nc_ofid,RDON_vid,"standard_name","RDON")
			ierr = nf90_put_att(nc_ofid,RDON_vid,"units","gN meters-3")
			ierr = nf90_put_att(nc_ofid,RDON_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,RDON_vid,"type","data")


		case ( "LPON")
			ierr = nf90_def_var(nc_ofid,"LPON",nf90_float,dynm3dn_lay,LPON_vid)
			ierr = nf90_put_att(nc_ofid,LPON_vid,"long_name","LPON")
			ierr = nf90_put_att(nc_ofid,LPON_vid,"standard_name","LPON")
			ierr = nf90_put_att(nc_ofid,LPON_vid,"units","gN meters-3")
			ierr = nf90_put_att(nc_ofid,LPON_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,LPON_vid,"type","data")

		case ( "RPON")
			ierr = nf90_def_var(nc_ofid,"RPON",nf90_float,dynm3dn_lay,RPON_vid)
			ierr = nf90_put_att(nc_ofid,RPON_vid,"long_name","RPON")
			ierr = nf90_put_att(nc_ofid,RPON_vid,"standard_name","RPON")
			ierr = nf90_put_att(nc_ofid,RPON_vid,"units","gN meters-3")
			ierr = nf90_put_att(nc_ofid,RPON_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,RPON_vid,"type","data")

		case ( "PO4")
			ierr = nf90_def_var(nc_ofid,"PO4",nf90_float,dynm3dn_lay,PO4_vid)
			ierr = nf90_put_att(nc_ofid,PO4_vid,"long_name","PO4")
			ierr = nf90_put_att(nc_ofid,PO4_vid,"standard_name","PO4")
			ierr = nf90_put_att(nc_ofid,PO4_vid,"units","gP meters-3")
			ierr = nf90_put_att(nc_ofid,PO4_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,PO4_vid,"type","data")

		case ( "P1")
			ierr = nf90_def_var(nc_ofid,"P1",nf90_float,dynm3dn_lay,P1_vid)
			ierr = nf90_put_att(nc_ofid,P1_vid,"long_name","P1")
			ierr = nf90_put_att(nc_ofid,P1_vid,"standard_name","P1")
			ierr = nf90_put_att(nc_ofid,P1_vid,"units","day-1")
			ierr = nf90_put_att(nc_ofid,P1_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,P1_vid,"type","data")
		case ( "P2")
			ierr = nf90_def_var(nc_ofid,"P2",nf90_float,dynm3dn_lay,P2_vid)
			ierr = nf90_put_att(nc_ofid,P2_vid,"long_name","P2")
			ierr = nf90_put_att(nc_ofid,P2_vid,"standard_name","P2")
			ierr = nf90_put_att(nc_ofid,P2_vid,"units","day-1")
			ierr = nf90_put_att(nc_ofid,P2_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,P2_vid,"type","data")
		case ( "BM1")
			ierr = nf90_def_var(nc_ofid,"BM1",nf90_float,dynm3dn_lay,BM1_vid)
			ierr = nf90_put_att(nc_ofid,BM1_vid,"long_name","BM1")
			ierr = nf90_put_att(nc_ofid,BM1_vid,"standard_name","BM1")
			ierr = nf90_put_att(nc_ofid,BM1_vid,"units","day-1")
			ierr = nf90_put_att(nc_ofid,BM1_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,BM1_vid,"type","data")
		case ( "BM2")
			ierr = nf90_def_var(nc_ofid,"BM2",nf90_float,dynm3dn_lay,BM2_vid)
			ierr = nf90_put_att(nc_ofid,BM2_vid,"long_name","BM2")
			ierr = nf90_put_att(nc_ofid,BM2_vid,"standard_name","BM2")
			ierr = nf90_put_att(nc_ofid,BM2_vid,"units","day-1")
			ierr = nf90_put_att(nc_ofid,BM2_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,BM2_vid,"type","data")
		case ( "NL1")
			ierr = nf90_def_var(nc_ofid,"NL1",nf90_float,dynm3dn_lay,NL1_vid)
			ierr = nf90_put_att(nc_ofid,NL1_vid,"long_name","NL1")
			ierr = nf90_put_att(nc_ofid,NL1_vid,"standard_name","NL1")
			ierr = nf90_put_att(nc_ofid,NL1_vid,"units","-")
			ierr = nf90_put_att(nc_ofid,NL1_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,NL1_vid,"type","data")
		case ( "NL2")
			ierr = nf90_def_var(nc_ofid,"NL2",nf90_float,dynm3dn_lay,NL2_vid)
			ierr = nf90_put_att(nc_ofid,NL2_vid,"long_name","NL2")
			ierr = nf90_put_att(nc_ofid,NL2_vid,"standard_name","NL2")
			ierr = nf90_put_att(nc_ofid,NL2_vid,"units","-")
			ierr = nf90_put_att(nc_ofid,NL2_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,NL2_vid,"type","data")
		case ( "PL1")
			ierr = nf90_def_var(nc_ofid,"PL1",nf90_float,dynm3dn_lay,PL1_vid)
			ierr = nf90_put_att(nc_ofid,PL1_vid,"long_name","PL1")
			ierr = nf90_put_att(nc_ofid,PL1_vid,"standard_name","PL1")
			ierr = nf90_put_att(nc_ofid,PL1_vid,"units","-")
			ierr = nf90_put_att(nc_ofid,PL1_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,PL1_vid,"type","data")
		case ( "PL2")
			ierr = nf90_def_var(nc_ofid,"PL2",nf90_float,dynm3dn_lay,PL2_vid)
			ierr = nf90_put_att(nc_ofid,PL2_vid,"long_name","PL2")
			ierr = nf90_put_att(nc_ofid,PL2_vid,"standard_name","PL2")
			ierr = nf90_put_att(nc_ofid,PL2_vid,"units","-")
			ierr = nf90_put_att(nc_ofid,PL2_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,PL2_vid,"type","data")
!: the following are related to zoop
		case ( "B1SZ")
			ierr = nf90_def_var(nc_ofid,"B1SZ",nf90_float,dynm3dn_lay,B1SZ_vid)
			ierr = nf90_put_att(nc_ofid,B1SZ_vid,"long_name","B1SZ")
			ierr = nf90_put_att(nc_ofid,B1SZ_vid,"standard_name","B1SZ")
			ierr = nf90_put_att(nc_ofid,B1SZ_vid,"units","gC meter-3 day-1")
			ierr = nf90_put_att(nc_ofid,B1SZ_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,B1SZ_vid,"type","data")
		case ( "B2SZ")
			ierr = nf90_def_var(nc_ofid,"B2SZ",nf90_float,dynm3dn_lay,B2SZ_vid)
			ierr = nf90_put_att(nc_ofid,B2SZ_vid,"long_name","B2SZ")
			ierr = nf90_put_att(nc_ofid,B2SZ_vid,"standard_name","B2SZ")
			ierr = nf90_put_att(nc_ofid,B2SZ_vid,"units","gC meter-3 day-1")
			ierr = nf90_put_att(nc_ofid,B2SZ_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,B2SZ_vid,"type","data")
		case ( "B1LZ")
			ierr = nf90_def_var(nc_ofid,"B1LZ",nf90_float,dynm3dn_lay,B1LZ_vid)
			ierr = nf90_put_att(nc_ofid,B1LZ_vid,"long_name","B1LZ")
			ierr = nf90_put_att(nc_ofid,B1LZ_vid,"standard_name","B1LZ")
			ierr = nf90_put_att(nc_ofid,B1LZ_vid,"units","gC meter-3 day-1")
			ierr = nf90_put_att(nc_ofid,B1LZ_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,B1LZ_vid,"type","data")
		case ( "B2LZ")
			ierr = nf90_def_var(nc_ofid,"B2LZ",nf90_float,dynm3dn_lay,B2LZ_vid)
			ierr = nf90_put_att(nc_ofid,B2LZ_vid,"long_name","B2LZ")
			ierr = nf90_put_att(nc_ofid,B2LZ_vid,"standard_name","B2LZ")
			ierr = nf90_put_att(nc_ofid,B2LZ_vid,"units","gC meter-3 day-1")
			ierr = nf90_put_att(nc_ofid,B2LZ_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,B2LZ_vid,"type","data")
!
		case ( "PR1")
			ierr = nf90_def_var(nc_ofid,"PR1",nf90_float,dynm3dn_lay,PR1_vid)
			ierr = nf90_put_att(nc_ofid,PR1_vid,"long_name","PR1")
			ierr = nf90_put_att(nc_ofid,PR1_vid,"standard_name","PR1")
			ierr = nf90_put_att(nc_ofid,PR1_vid,"units","gC meter-3 day-1")
			ierr = nf90_put_att(nc_ofid,PR1_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,PR1_vid,"type","data")
		case ( "PR2")
			ierr = nf90_def_var(nc_ofid,"PR2",nf90_float,dynm3dn_lay,PR2_vid)
			ierr = nf90_put_att(nc_ofid,PR2_vid,"long_name","PR2")
			ierr = nf90_put_att(nc_ofid,PR2_vid,"standard_name","PR2")
			ierr = nf90_put_att(nc_ofid,PR2_vid,"units","gC meter-3 day-1")
			ierr = nf90_put_att(nc_ofid,PR2_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,PR2_vid,"type","data")
		case ( "IAVG")
			ierr = nf90_def_var(nc_ofid,"IAVG",nf90_float,dynm3dn_lay,IAVG_vid)
			ierr = nf90_put_att(nc_ofid,IAVG_vid,"long_name","IAVG")
			ierr = nf90_put_att(nc_ofid,IAVG_vid,"standard_name","IAVG")
			ierr = nf90_put_att(nc_ofid,IAVG_vid,"units","Einstein meters-2 day-1")
			ierr = nf90_put_att(nc_ofid,IAVG_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,IAVG_vid,"type","data")
		case ( "FI1")
			ierr = nf90_def_var(nc_ofid,"FI1",nf90_float,dynm3dn_lay,FI1_vid)
			ierr = nf90_put_att(nc_ofid,FI1_vid,"long_name","FI1")
			ierr = nf90_put_att(nc_ofid,FI1_vid,"standard_name","FI1")
			ierr = nf90_put_att(nc_ofid,FI1_vid,"units","-")
			ierr = nf90_put_att(nc_ofid,FI1_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,FI1_vid,"type","data")
		case ( "FI2")
			ierr = nf90_def_var(nc_ofid,"FI2",nf90_float,dynm3dn_lay,FI2_vid)
			ierr = nf90_put_att(nc_ofid,FI2_vid,"long_name","FI2")
			ierr = nf90_put_att(nc_ofid,FI2_vid,"standard_name","FI2")
			ierr = nf90_put_att(nc_ofid,FI2_vid,"units","-")
			ierr = nf90_put_att(nc_ofid,FI2_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,FI2_vid,"type","data")
		case ( "IK1")
			ierr = nf90_def_var(nc_ofid,"IK1",nf90_float,dynm3dn_lay,IK1_vid)
			ierr = nf90_put_att(nc_ofid,IK1_vid,"long_name","IK1")
			ierr = nf90_put_att(nc_ofid,IK1_vid,"standard_name","IK1")
			ierr = nf90_put_att(nc_ofid,IK1_vid,"units","Einstein meters-2 day-1")
			ierr = nf90_put_att(nc_ofid,IK1_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,IK1_vid,"type","data")
		case ( "IK2")
			ierr = nf90_def_var(nc_ofid,"IK2",nf90_float,dynm3dn_lay,IK2_vid)
			ierr = nf90_put_att(nc_ofid,IK2_vid,"long_name","IK2")
			ierr = nf90_put_att(nc_ofid,IK2_vid,"standard_name","IK2")
			ierr = nf90_put_att(nc_ofid,IK2_vid,"units","Einstein meters-2 day-2")
			ierr = nf90_put_att(nc_ofid,IK2_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,IK2_vid,"type","data")
		case ( "GPP")
			ierr = nf90_def_var(nc_ofid,"GPP",nf90_float,dynm3dn_lay,GPP_vid)
			ierr = nf90_put_att(nc_ofid,GPP_vid,"long_name","GPP")
			ierr = nf90_put_att(nc_ofid,GPP_vid,"standard_name","GPP")
			ierr = nf90_put_att(nc_ofid,GPP_vid,"units","gC/m^2/day")
			ierr = nf90_put_att(nc_ofid,GPP_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,GPP_vid,"type","data")
		case ( "NPP")
			ierr = nf90_def_var(nc_ofid,"NPP",nf90_float,dynm3dn_lay,NPP_vid)
			ierr = nf90_put_att(nc_ofid,NPP_vid,"long_name","NPP")
			ierr = nf90_put_att(nc_ofid,NPP_vid,"standard_name","NPP")
			ierr = nf90_put_att(nc_ofid,NPP_vid,"units","gC/m^2/day")
			ierr = nf90_put_att(nc_ofid,NPP_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,NPP_vid,"type","data")
		case ( "DDOC")
			ierr = nf90_def_var(nc_ofid,"DDOC",nf90_float,dynm3dn_lay,DDOC_vid)
			ierr = nf90_put_att(nc_ofid,DDOC_vid,"long_name","DDOC")
			ierr = nf90_put_att(nc_ofid,DDOC_vid,"standard_name","DDOC")
			ierr = nf90_put_att(nc_ofid,DDOC_vid,"units","gO2/m^3/day")
			ierr = nf90_put_att(nc_ofid,DDOC_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DDOC_vid,"type","data")
		case ( "DCOD")
			ierr = nf90_def_var(nc_ofid,"DCOD",nf90_float,dynm3dn_lay,DCOD_vid)
			ierr = nf90_put_att(nc_ofid,DCOD_vid,"long_name","DCOD")
			ierr = nf90_put_att(nc_ofid,DCOD_vid,"standard_name","DCOD")
			ierr = nf90_put_att(nc_ofid,DCOD_vid,"units","gO2-equiv/m^3/day")
			ierr = nf90_put_att(nc_ofid,DCOD_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DCOD_vid,"type","data")
		case ( "NITRIF")
			ierr = nf90_def_var(nc_ofid,"NITRIF",nf90_float,dynm3dn_lay,NITRIF_vid)
			ierr = nf90_put_att(nc_ofid,NITRIF_vid,"long_name","NITRIF")
			ierr = nf90_put_att(nc_ofid,NITRIF_vid,"standard_name","NITRIF")
			ierr = nf90_put_att(nc_ofid,NITRIF_vid,"units","gO2-equiv/m^3/day")
			ierr = nf90_put_att(nc_ofid,NITRIF_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,NITRIF_vid,"type","data")
		case ( "DOSZ")
			ierr = nf90_def_var(nc_ofid,"DOSZ",nf90_float,dynm3dn_lay,DOSZ_vid)
			ierr = nf90_put_att(nc_ofid,DOSZ_vid,"long_name","DOSZ")
			ierr = nf90_put_att(nc_ofid,DOSZ_vid,"standard_name","DOSZ")
			ierr = nf90_put_att(nc_ofid,DOSZ_vid,"units","gO2/m^3/day")
			ierr = nf90_put_att(nc_ofid,DOSZ_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DOSZ_vid,"type","data")
		case ( "DOLZ")
			ierr = nf90_def_var(nc_ofid,"DOLZ",nf90_float,dynm3dn_lay,DOLZ_vid)
			ierr = nf90_put_att(nc_ofid,DOLZ_vid,"long_name","DOLZ")
			ierr = nf90_put_att(nc_ofid,DOLZ_vid,"standard_name","DOLZ")
			ierr = nf90_put_att(nc_ofid,DOLZ_vid,"units","gO2/m^3/day")
			ierr = nf90_put_att(nc_ofid,DOLZ_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DOLZ_vid,"type","data")
		case ( "DOPR")
			ierr = nf90_def_var(nc_ofid,"DOPR",nf90_float,dynm3dn_lay,DOPR_vid)
			ierr = nf90_put_att(nc_ofid,DOPR_vid,"long_name","DOPR")
			ierr = nf90_put_att(nc_ofid,DOPR_vid,"standard_name","DOPR")
			ierr = nf90_put_att(nc_ofid,DOPR_vid,"units","gO2/m^3/day")
			ierr = nf90_put_att(nc_ofid,DOPR_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DOPR_vid,"type","data")
		case ( "DORALG")
			ierr = nf90_def_var(nc_ofid,"DORALG",nf90_float,dynm3dn_lay,DORALG_vid)
			ierr = nf90_put_att(nc_ofid,DORALG_vid,"long_name","DORALG")
			ierr = nf90_put_att(nc_ofid,DORALG_vid,"standard_name","DORALG")
			ierr = nf90_put_att(nc_ofid,DORALG_vid,"units","gO2/m^3/day")
			ierr = nf90_put_att(nc_ofid,DORALG_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DORALG_vid,"type","data")
		case ( "RESP")
			ierr = nf90_def_var(nc_ofid,"RESP",nf90_float,dynm3dn_lay,RESP_vid)
			ierr = nf90_put_att(nc_ofid,RESP_vid,"long_name","RESP")
			ierr = nf90_put_att(nc_ofid,RESP_vid,"standard_name","RESP")
			ierr = nf90_put_att(nc_ofid,RESP_vid,"units","gO2/m^3/day")
			ierr = nf90_put_att(nc_ofid,RESP_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,RESP_vid,"type","data")

		case ( "FTCOD")
			IF(COD_CALC)THEN
			ierr = nf90_def_var(nc_ofid,"FTCOD",nf90_float,dynm3dn_lay,FTCOD_vid)
			ierr = nf90_put_att(nc_ofid,FTCOD_vid,"long_name","FTCOD")
			ierr = nf90_put_att(nc_ofid,FTCOD_vid,"standard_name","FTCOD")
			ierr = nf90_put_att(nc_ofid,FTCOD_vid,"units","gO2-equiv/m^3/day")
			ierr = nf90_put_att(nc_ofid,FTCOD_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,FTCOD_vid,"type","data")
			ENDIF
		case ( "COD")
			IF(COD_CALC)THEN
			ierr = nf90_def_var(nc_ofid,"COD",nf90_float,dynm3dn_lay,COD_vid)
			ierr = nf90_put_att(nc_ofid,COD_vid,"long_name","COD")
			ierr = nf90_put_att(nc_ofid,COD_vid,"standard_name","COD")
			ierr = nf90_put_att(nc_ofid,COD_vid,"units","gO2-equiv/m^3/")
			ierr = nf90_put_att(nc_ofid,COD_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,COD_vid,"type","data")
			ENDIF
		case ( "REAERDO")
			ierr = nf90_def_var(nc_ofid,"REAERDO",nf90_float,dynm3dn_lay,REAERDO_vid)
			ierr = nf90_put_att(nc_ofid,REAERDO_vid,"long_name","REAERDO")
			ierr = nf90_put_att(nc_ofid,REAERDO_vid,"standard_name","REAERDO")
			ierr = nf90_put_att(nc_ofid,REAERDO_vid,"units","gO2/m^3/day")
			ierr = nf90_put_att(nc_ofid,REAERDO_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,REAERDO_vid,"type","data")
		case ( "BFLUX")
			ierr = nf90_def_var(nc_ofid,"BFLUX",nf90_float,dynm3dn_9,BFLUX_vid)  !should be 9 layers only
			ierr = nf90_put_att(nc_ofid,BFLUX_vid,"long_name","BFLUX")
			ierr = nf90_put_att(nc_ofid,BFLUX_vid,"standard_name","BFLUX")
			ierr = nf90_put_att(nc_ofid,BFLUX_vid,"units","???")
			ierr = nf90_put_att(nc_ofid,BFLUX_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,BFLUX_vid,"type","data")
!: the following are related to CARBONATE_CALC
		case ( "TDIC")
		case ( "TALK")
		case ( "pH")
		case ( "pCO2")
		case ( "DICUPT")
		case ( "DICBMP")
		case ( "DICPRD")
		case ( "DICMNL")
		case ( "DICDEN")
		case ( "DICGAS")
		case ( "DICSED")
		case ( "ALKNH4")
		case ( "ALKNO3")
		case ( "ALKNIT")
		case ( "ALKDEN")
		case ( "ALKREM")
		case ( "ALKNH4SED")
		case ( "ALKNO3SED")
		!case ( "ADVFLUX")
		!	ierr = nf90_def_var(nc_ofid,"ADVFLUX",nf90_float,dynm3dn_lay,ADVFLUX_vid)
		!	ierr = nf90_put_att(nc_ofid,ADVFLUX_vid,"long_name","ADVFLUX")
		!	ierr = nf90_put_att(nc_ofid,ADVFLUX_vid,"standard_name","ADVFLUX")
		!	ierr = nf90_put_att(nc_ofid,ADVFLUX_vid,"units","????")
		!	ierr = nf90_put_att(nc_ofid,ADVFLUX_vid,"grid","fvcom_grid")
		!	ierr = nf90_put_att(nc_ofid,ADVFLUX_vid,"type","data")
		!case ( "VDIFFLUX")
		!	ierr = nf90_def_var(nc_ofid,"VDIFFLUX",nf90_float,dynm3dn_lay,VDIFFLUX_vid)
		!	ierr = nf90_put_att(nc_ofid,VDIFFLUX_vid,"long_name","VDIFFLUX")
		!	ierr = nf90_put_att(nc_ofid,VDIFFLUX_vid,"standard_name","VDIFFLUX")
		!	ierr = nf90_put_att(nc_ofid,VDIFFLUX_vid,"units","???")
		!	ierr = nf90_put_att(nc_ofid,VDIFFLUX_vid,"grid","fvcom_grid")
		!	ierr = nf90_put_att(nc_ofid,VDIFFLUX_vid,"type","data")
!: the following are related to SEDIMENT_CALC
		case ( "JPOC")
		case ( "JPON")
		case ( "JPOP")
		case ( "JPOS")
		case ( "CTEMP")
		case ( "SODTM1S")
		case ( "JNH4TM1S")
		case ( "JNO3TM1S")
		case ( "BENDEN")
		case ( "JCH4TM1S")
		case ( "JCH4GTM1S")
		case ( "JHSTM1S")
		case ( "JPO4TM1S")
		case ( "JSITM1S")
		case ( "NH41TM1S")
		case ( "NH42TM1S")
		case ( "NO31TM1S")
		case ( "NO32TM1S")
		case ( "PO41TM1S")
		case ( "PO42TM1S")
		case ( "SI1TM1S")
		case ( "SI2TM1S")
		case ( "CH41TM1S")
		case ( "CH42TM1S")
		case ( "HS1TM1S")
		case ( "HS2TM1S")
		case ( "CPOP")
		case ( "CPON")
		case ( "CPOC")
		case ( "CPOS")
		case ( "HSED1")
		case ( "BENSTRTM1S")

		case ( "TSS")
		case ( "SSI")
		case ( "KESS")
		case ( "TURB")
		case ( "PARAD")

		case ( "LEAF")
		case ( "STEM")
		case ( "ROOT")
		case ( "TUBER")
		case ( "EP")

		case ( "SZ")
		case ( "LZ")
     case default
       write(ipt,*)'variable',cdf_vdp(i),' not set up for netcdf output'
       write(ipt,*)'modify module mod_ncdio.f'
       call pstop
end select

 IF(CARBONATE_CALC)THEN
		select case (trim(cdf_vdp(i)))
		case ( "TDIC")
			ierr = nf90_def_var(nc_ofid,"TDIC",nf90_float,dynm3dn_lay,TDIC_vid)
			ierr = nf90_put_att(nc_ofid,TDIC_vid,"long_name","TDIC")
			ierr = nf90_put_att(nc_ofid,TDIC_vid,"standard_name","TDIC")
			ierr = nf90_put_att(nc_ofid,TDIC_vid,"units","mmolC meters-3")
			ierr = nf90_put_att(nc_ofid,TDIC_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,TDIC_vid,"type","data")
		case ( "TALK")
			ierr = nf90_def_var(nc_ofid,"TALK",nf90_float,dynm3dn_lay,TALK_vid)
			ierr = nf90_put_att(nc_ofid,TALK_vid,"long_name","TALK")
			ierr = nf90_put_att(nc_ofid,TALK_vid,"standard_name","TALK")
			ierr = nf90_put_att(nc_ofid,TALK_vid,"units","mmol meters-3")
			ierr = nf90_put_att(nc_ofid,TALK_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,TALK_vid,"type","data")
		case ( "pH")
			ierr = nf90_def_var(nc_ofid,"pH",nf90_float,dynm3dn_lay,pH_vid)
			ierr = nf90_put_att(nc_ofid,pH_vid,"long_name","pH")
			ierr = nf90_put_att(nc_ofid,pH_vid,"standard_name","pH")
			ierr = nf90_put_att(nc_ofid,pH_vid,"units","SeaWater Scale")
			ierr = nf90_put_att(nc_ofid,pH_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,pH_vid,"type","data")
		case ( "pCO2")
			ierr = nf90_def_var(nc_ofid,"pCO2",nf90_float,dynm3dn_lay,pCO2_vid)
			ierr = nf90_put_att(nc_ofid,pCO2_vid,"long_name","pCO2")
			ierr = nf90_put_att(nc_ofid,pCO2_vid,"standard_name","pCO2")
			ierr = nf90_put_att(nc_ofid,pCO2_vid,"units","")
			ierr = nf90_put_att(nc_ofid,pCO2_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,pCO2_vid,"type","data")
		case ( "DICUPT")
			ierr = nf90_def_var(nc_ofid,"DICUPT",nf90_float,dynm3dn_lay,DICUPT_vid)
			ierr = nf90_put_att(nc_ofid,DICUPT_vid,"long_name","DICUPT")
			ierr = nf90_put_att(nc_ofid,DICUPT_vid,"standard_name","DICUPT")
			ierr = nf90_put_att(nc_ofid,DICUPT_vid,"units","mmolC/m^3/day")
			ierr = nf90_put_att(nc_ofid,DICUPT_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DICUPT_vid,"type","data")
		case ( "DICBMP")
			ierr = nf90_def_var(nc_ofid,"DICBMP",nf90_float,dynm3dn_lay,DICBMP_vid)
			ierr = nf90_put_att(nc_ofid,DICBMP_vid,"long_name","DICBMPC")
			ierr = nf90_put_att(nc_ofid,DICBMP_vid,"standard_name","DICBMP")
			ierr = nf90_put_att(nc_ofid,DICBMP_vid,"units","mmolC/m3/day")
			ierr = nf90_put_att(nc_ofid,DICBMP_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DICBMP_vid,"type","data")
		case ( "DICPRD")
			ierr = nf90_def_var(nc_ofid,"DICPRD",nf90_float,dynm3dn_lay,DICPRD_vid)
			ierr = nf90_put_att(nc_ofid,DICPRD_vid,"long_name","DICPRD")
			ierr = nf90_put_att(nc_ofid,DICPRD_vid,"standard_name","DICPRD")
			ierr = nf90_put_att(nc_ofid,DICPRD_vid,"units"," mmolC/m3/day")
			ierr = nf90_put_att(nc_ofid,DICPRD_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DICPRD_vid,"type","data")
		case ( "DICMNL")
			ierr = nf90_def_var(nc_ofid,"DICMNL",nf90_float,dynm3dn_lay,DICMNL_vid)
			ierr = nf90_put_att(nc_ofid,DICMNL_vid,"long_name","DICMNL")
			ierr = nf90_put_att(nc_ofid,DICMNL_vid,"standard_name","DICMNL")
			ierr = nf90_put_att(nc_ofid,DICMNL_vid,"units","mmolC/m3/day")
			ierr = nf90_put_att(nc_ofid,DICMNL_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DICMNL_vid,"type","data")
		case ( "DICDEN")
			ierr = nf90_def_var(nc_ofid,"DICDEN",nf90_float,dynm3dn_lay,DICDEN_vid)
			ierr = nf90_put_att(nc_ofid,DICDEN_vid,"long_name","DICDEN")
			ierr = nf90_put_att(nc_ofid,DICDEN_vid,"standard_name","DICDEN")
			ierr = nf90_put_att(nc_ofid,DICDEN_vid,"units","mmolC/m3/day")
			ierr = nf90_put_att(nc_ofid,DICDEN_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DICDEN_vid,"type","data")
		case ( "DICGAS")
			ierr = nf90_def_var(nc_ofid,"DICGAS",nf90_float,dynm3dn_lay,DICGAS_vid)
			ierr = nf90_put_att(nc_ofid,DICGAS_vid,"long_name","DICGAS")
			ierr = nf90_put_att(nc_ofid,DICGAS_vid,"standard_name","DICGAS")
			ierr = nf90_put_att(nc_ofid,DICGAS_vid,"units","mmolC/m3/day")
			ierr = nf90_put_att(nc_ofid,DICGAS_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DICGAS_vid,"type","data")
		case ( "DICSED")
			ierr = nf90_def_var(nc_ofid,"DICSED",nf90_float,dynm3dn_lay,DICSED_vid)
			ierr = nf90_put_att(nc_ofid,DICSED_vid,"long_name","DICSED")
			ierr = nf90_put_att(nc_ofid,DICSED_vid,"standard_name","DICSED")
			ierr = nf90_put_att(nc_ofid,DICSED_vid,"units","mmolC/m3/day")
			ierr = nf90_put_att(nc_ofid,DICSED_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,DICSED_vid,"type","data")
		case ( "ALKNH4")
			ierr = nf90_def_var(nc_ofid,"ALKNH4",nf90_float,dynm3dn_lay,ALKNH4_vid)
			ierr = nf90_put_att(nc_ofid,ALKNH4_vid,"long_name","ALKNH4")
			ierr = nf90_put_att(nc_ofid,ALKNH4_vid,"standard_name","ALKNH4")
			ierr = nf90_put_att(nc_ofid,ALKNH4_vid,"units","mmol/m^3/day")
			ierr = nf90_put_att(nc_ofid,ALKNH4_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,ALKNH4_vid,"type","data")
		case ( "ALKNO3")
			ierr = nf90_def_var(nc_ofid,"ALKNO3",nf90_float,dynm3dn_lay,ALKNO3_vid)
			ierr = nf90_put_att(nc_ofid,ALKNO3_vid,"long_name","ALKNO3")
			ierr = nf90_put_att(nc_ofid,ALKNO3_vid,"standard_name","ALKNO3")
			ierr = nf90_put_att(nc_ofid,ALKNO3_vid,"units","mmol/m^3/day")
			ierr = nf90_put_att(nc_ofid,ALKNO3_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,ALKNO3_vid,"type","data")
		case ( "ALKNIT")
			ierr = nf90_def_var(nc_ofid,"ALKNIT",nf90_float,dynm3dn_lay,ALKNIT_vid)
			ierr = nf90_put_att(nc_ofid,ALKNIT_vid,"long_name","ALKNIT")
			ierr = nf90_put_att(nc_ofid,ALKNIT_vid,"standard_name","ALKNIT")
			ierr = nf90_put_att(nc_ofid,ALKNIT_vid,"units","mmol/m^3/day")
			ierr = nf90_put_att(nc_ofid,ALKNIT_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,ALKNIT_vid,"type","data")
		case ( "ALKDEN")
			ierr = nf90_def_var(nc_ofid,"ALKDEN",nf90_float,dynm3dn_lay,ALKDEN_vid)
			ierr = nf90_put_att(nc_ofid,ALKDEN_vid,"long_name","ALKDEN")
			ierr = nf90_put_att(nc_ofid,ALKDEN_vid,"standard_name","ALKDEN")
			ierr = nf90_put_att(nc_ofid,ALKDEN_vid,"units","mmol/m^3/day")
			ierr = nf90_put_att(nc_ofid,ALKDEN_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,ALKDEN_vid,"type","data")
		case ( "ALKREM")
			ierr = nf90_def_var(nc_ofid,"ALKREM",nf90_float,dynm3dn_lay,ALKREM_vid)
			ierr = nf90_put_att(nc_ofid,ALKREM_vid,"long_name","ALKREM")
			ierr = nf90_put_att(nc_ofid,ALKREM_vid,"standard_name","ALKREM")
			ierr = nf90_put_att(nc_ofid,ALKREM_vid,"units","mmol/m^3/day")
			ierr = nf90_put_att(nc_ofid,ALKREM_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,ALKREM_vid,"type","data")
		case ( "ALKNH4SED")
			ierr = nf90_def_var(nc_ofid,"ALKNH4SED",nf90_float,dynm3dn_lay,ALKNH4SED_vid)
			ierr = nf90_put_att(nc_ofid,ALKNH4SED_vid,"long_name","ALKNH4SED")
			ierr = nf90_put_att(nc_ofid,ALKNH4SED_vid,"standard_name","ALKNH4SED")
			ierr = nf90_put_att(nc_ofid,ALKNH4SED_vid,"units","mmol/m^3/day")
			ierr = nf90_put_att(nc_ofid,ALKNH4SED_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,ALKNH4SED_vid,"type","data")
		case ( "ALKNO3SED")
			ierr = nf90_def_var(nc_ofid,"ALKNO3SED",nf90_float,dynm3dn_lay,ALKNO3SED_vid)
			ierr = nf90_put_att(nc_ofid,ALKNO3SED_vid,"long_name","ALKNO3SED")
			ierr = nf90_put_att(nc_ofid,ALKNO3SED_vid,"standard_name","ALKNO3SED")
			ierr = nf90_put_att(nc_ofid,ALKNO3SED_vid,"units","mmol/m^3/day")
			ierr = nf90_put_att(nc_ofid,ALKNO3SED_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,ALKNO3SED_vid,"type","data")
end select
ENDIF

	IF (SEDIMENT_CALC) THEN
		select case(trim(cdf_vdp(i)))
		case ( "JPOC")
			ierr = nf90_def_var(nc_ofid,"JPOC",nf90_float,dynm3dn_3,JPOC_vid)  !3 layers
			ierr = nf90_put_att(nc_ofid,JPOC_vid,"long_name","JPOC")
			ierr = nf90_put_att(nc_ofid,JPOC_vid,"standard_name","JPOC")
			ierr = nf90_put_att(nc_ofid,JPOC_vid,"units","mgC/m^2/day")
			ierr = nf90_put_att(nc_ofid,JPOC_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,JPOC_vid,"type","data")
		case ( "JPON")
			ierr = nf90_def_var(nc_ofid,"JPON",nf90_float,dynm3dn_3,JPON_vid)
			ierr = nf90_put_att(nc_ofid,JPON_vid,"long_name","JPON")
			ierr = nf90_put_att(nc_ofid,JPON_vid,"standard_name","JPON")
			ierr = nf90_put_att(nc_ofid,JPON_vid,"units","mgN/m^2/day")
			ierr = nf90_put_att(nc_ofid,JPON_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,JPON_vid,"type","data")
		case ( "JPOP")
			ierr = nf90_def_var(nc_ofid,"JPOP",nf90_float,dynm3dn_3,JPOP_vid)
			ierr = nf90_put_att(nc_ofid,JPOP_vid,"long_name","JPOP")
			ierr = nf90_put_att(nc_ofid,JPOP_vid,"standard_name","JPOP")
			ierr = nf90_put_att(nc_ofid,JPOP_vid,"units","mgP/m^2/day")
			ierr = nf90_put_att(nc_ofid,JPOP_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,JPOP_vid,"type","data")
		case ( "JPOS")
			ierr = nf90_def_var(nc_ofid,"JPOS",nf90_float,dynm2dn,JPOS_vid)
			ierr = nf90_put_att(nc_ofid,JPOS_vid,"long_name","JPOS")
			ierr = nf90_put_att(nc_ofid,JPOS_vid,"standard_name","JPOS")
			ierr = nf90_put_att(nc_ofid,JPOS_vid,"units","gO2-equiv/m^2/day")
			ierr = nf90_put_att(nc_ofid,JPOS_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,JPOS_vid,"type","data")
		case ( "CTEMP")
			ierr = nf90_def_var(nc_ofid,"CTEMP",nf90_float,dynm2dn,CTEMP_vid)
			ierr = nf90_put_att(nc_ofid,CTEMP_vid,"long_name","CTEMP")
			ierr = nf90_put_att(nc_ofid,CTEMP_vid,"standard_name","CTEMP")
			ierr = nf90_put_att(nc_ofid,CTEMP_vid,"units","degC")
			ierr = nf90_put_att(nc_ofid,CTEMP_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,CTEMP_vid,"type","data")
		case ( "SODTM1S")
			ierr = nf90_def_var(nc_ofid,"SODTM1S",nf90_float,dynm2dn,SODTM1S_vid)
			ierr = nf90_put_att(nc_ofid,SODTM1S_vid,"long_name","SODTM1S")
			ierr = nf90_put_att(nc_ofid,SODTM1S_vid,"standard_name","SODTM1S")
			ierr = nf90_put_att(nc_ofid,SODTM1S_vid,"units","gO2/m^2/day")
			ierr = nf90_put_att(nc_ofid,SODTM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,SODTM1S_vid,"type","data")
		case ( "JNH4TM1S")
			ierr = nf90_def_var(nc_ofid,"JNH4TM1S",nf90_float,dynm2dn,JNH4TM1S_vid)
			ierr = nf90_put_att(nc_ofid,JNH4TM1S_vid,"long_name","JNH4TM1S")
			ierr = nf90_put_att(nc_ofid,JNH4TM1S_vid,"standard_name","JNH4TM1S")
			ierr = nf90_put_att(nc_ofid,JNH4TM1S_vid,"units","mgN/m^2/day")
			ierr = nf90_put_att(nc_ofid,JNH4TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,JNH4TM1S_vid,"type","data")
		case ( "JNO3TM1S")
			ierr = nf90_def_var(nc_ofid,"JNO3TM1S",nf90_float,dynm2dn,JNO3TM1S_vid)
			ierr = nf90_put_att(nc_ofid,JNO3TM1S_vid,"long_name","JNO3TM1S")
			ierr = nf90_put_att(nc_ofid,JNO3TM1S_vid,"standard_name","JNO3TM1S")
			ierr = nf90_put_att(nc_ofid,JNO3TM1S_vid,"units","mgN/m^2/day")
			ierr = nf90_put_att(nc_ofid,JNO3TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,JNO3TM1S_vid,"type","data")
		case ( "BENDEN")
			ierr = nf90_def_var(nc_ofid,"BENDEN",nf90_float,dynm2dn,BENDEN_vid)
			ierr = nf90_put_att(nc_ofid,BENDEN_vid,"long_name","BENDEN")
			ierr = nf90_put_att(nc_ofid,BENDEN_vid,"standard_name","BENDEN")
			ierr = nf90_put_att(nc_ofid,BENDEN_vid,"units","mgN/m^2/day")
			ierr = nf90_put_att(nc_ofid,BENDEN_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,BENDEN_vid,"type","data")
		case ( "JCH4TM1S")
			ierr = nf90_def_var(nc_ofid,"JCH4TM1S",nf90_float,dynm2dn,JCH4TM1S_vid)
			ierr = nf90_put_att(nc_ofid,JCH4TM1S_vid,"long_name","JCH4TM1S")
			ierr = nf90_put_att(nc_ofid,JCH4TM1S_vid,"standard_name","JCH4TM1S")
			ierr = nf90_put_att(nc_ofid,JCH4TM1S_vid,"units","mgO2/m^2/day")
			ierr = nf90_put_att(nc_ofid,JCH4TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,JCH4TM1S_vid,"type","data")
		case ( "JCH4GTM1S")
			ierr = nf90_def_var(nc_ofid,"JCH4GTM1S",nf90_float,dynm2dn,JCH4GTM1S_vid)
			ierr = nf90_put_att(nc_ofid,JCH4GTM1S_vid,"long_name","JCH4GTM1S")
			ierr = nf90_put_att(nc_ofid,JCH4GTM1S_vid,"standard_name","JCH4GTM1S")
			ierr = nf90_put_att(nc_ofid,JCH4GTM1S_vid,"units","mgO2/m^2/day")
			ierr = nf90_put_att(nc_ofid,JCH4GTM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,JCH4GTM1S_vid,"type","data")
		case ( "JHSTM1S")
			ierr = nf90_def_var(nc_ofid,"JHSTM1S",nf90_float,dynm2dn,JHSTM1S_vid)
			ierr = nf90_put_att(nc_ofid,JHSTM1S_vid,"long_name","JHSTM1S")
			ierr = nf90_put_att(nc_ofid,JHSTM1S_vid,"standard_name","JHSTM1S")
			ierr = nf90_put_att(nc_ofid,JHSTM1S_vid,"units","mgO2/m^2/day")
			ierr = nf90_put_att(nc_ofid,JHSTM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,JHSTM1S_vid,"type","data")
		case ( "JPO4TM1S")
			ierr = nf90_def_var(nc_ofid,"JPO4TM1S",nf90_float,dynm2dn,JPO4TM1S_vid)
			ierr = nf90_put_att(nc_ofid,JPO4TM1S_vid,"long_name","JPO4TM1S")
			ierr = nf90_put_att(nc_ofid,JPO4TM1S_vid,"standard_name","JPO4TM1S")
			ierr = nf90_put_att(nc_ofid,JPO4TM1S_vid,"units","mgP/m^2/day")
			ierr = nf90_put_att(nc_ofid,JPO4TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,JPO4TM1S_vid,"type","data")
		case ( "JSITM1S")
			ierr = nf90_def_var(nc_ofid,"JSITM1S",nf90_float,dynm2dn,JSITM1S_vid)
			ierr = nf90_put_att(nc_ofid,JSITM1S_vid,"long_name","JSITM1S")
			ierr = nf90_put_att(nc_ofid,JSITM1S_vid,"standard_name","JSITM1S")
			ierr = nf90_put_att(nc_ofid,JSITM1S_vid,"units","mgO2/m^2/day")
			ierr = nf90_put_att(nc_ofid,JSITM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,JSITM1S_vid,"type","data")
		case ( "NH41TM1S")
			ierr = nf90_def_var(nc_ofid,"NH41TM1S",nf90_float,dynm2dn,NH41TM1S_vid)
			ierr = nf90_put_att(nc_ofid,NH41TM1S_vid,"long_name","NH41TM1S")
			ierr = nf90_put_att(nc_ofid,NH41TM1S_vid,"standard_name","NH41TM1S")
			ierr = nf90_put_att(nc_ofid,NH41TM1S_vid,"units","mgN/m^3")
			ierr = nf90_put_att(nc_ofid,NH41TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,NH41TM1S_vid,"type","data")
		case ( "NH42TM1S")
			ierr = nf90_def_var(nc_ofid,"NH42TM1S",nf90_float,dynm2dn,NH42TM1S_vid)
			ierr = nf90_put_att(nc_ofid,NH42TM1S_vid,"long_name","NH42TM1S")
			ierr = nf90_put_att(nc_ofid,NH42TM1S_vid,"standard_name","NH42TM1S")
			ierr = nf90_put_att(nc_ofid,NH42TM1S_vid,"units","mgN/m^3")
			ierr = nf90_put_att(nc_ofid,NH42TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,NH42TM1S_vid,"type","data")
		case ( "NO31TM1S")
			ierr = nf90_def_var(nc_ofid,"NO31TM1S",nf90_float,dynm2dn,NO31TM1S_vid)
			ierr = nf90_put_att(nc_ofid,NO31TM1S_vid,"long_name","NO31TM1S")
			ierr = nf90_put_att(nc_ofid,NO31TM1S_vid,"standard_name","NO31TM1S")
			ierr = nf90_put_att(nc_ofid,NO31TM1S_vid,"units","mgN/m^3")
			ierr = nf90_put_att(nc_ofid,NO31TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,NO31TM1S_vid,"type","data")
		case ( "NO32TM1S")
			ierr = nf90_def_var(nc_ofid,"NO32TM1S",nf90_float,dynm2dn,NO32TM1S_vid)
			ierr = nf90_put_att(nc_ofid,NO32TM1S_vid,"long_name","NO32TM1S")
			ierr = nf90_put_att(nc_ofid,NO32TM1S_vid,"standard_name","NO32TM1S")
			ierr = nf90_put_att(nc_ofid,NO32TM1S_vid,"units","mgN/m^3")
			ierr = nf90_put_att(nc_ofid,NO32TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,NO32TM1S_vid,"type","data")
		case ( "PO41TM1S")
			ierr = nf90_def_var(nc_ofid,"PO41TM1S",nf90_float,dynm2dn,PO41TM1S_vid)
			ierr = nf90_put_att(nc_ofid,PO41TM1S_vid,"long_name","PO41TM1S")
			ierr = nf90_put_att(nc_ofid,PO41TM1S_vid,"standard_name","PO41TM1S")
			ierr = nf90_put_att(nc_ofid,PO41TM1S_vid,"units","mgP/m^3")
			ierr = nf90_put_att(nc_ofid,PO41TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,PO41TM1S_vid,"type","data")
		case ( "PO42TM1S")
			ierr = nf90_def_var(nc_ofid,"PO42TM1S",nf90_float,dynm2dn,PO42TM1S_vid)
			ierr = nf90_put_att(nc_ofid,PO42TM1S_vid,"long_name","PO42TM1S")
			ierr = nf90_put_att(nc_ofid,PO42TM1S_vid,"standard_name","PO42TM1S")
			ierr = nf90_put_att(nc_ofid,PO42TM1S_vid,"units","mgP/m^3")
			ierr = nf90_put_att(nc_ofid,PO42TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,PO42TM1S_vid,"type","data")
		case ( "SI1TM1S")
			ierr = nf90_def_var(nc_ofid,"SI1TM1S",nf90_float,dynm2dn,SI1TM1S_vid)
			ierr = nf90_put_att(nc_ofid,SI1TM1S_vid,"long_name","SI1TM1S")
			ierr = nf90_put_att(nc_ofid,SI1TM1S_vid,"standard_name","SI1TM1S")
			ierr = nf90_put_att(nc_ofid,SI1TM1S_vid,"units","mgSi/m^3")
			ierr = nf90_put_att(nc_ofid,SI1TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,SI1TM1S_vid,"type","data")
		case ( "SI2TM1S")
			ierr = nf90_def_var(nc_ofid,"SI2TM1S",nf90_float,dynm2dn,SI2TM1S_vid)
			ierr = nf90_put_att(nc_ofid,SI2TM1S_vid,"long_name","SI2TM1S")
			ierr = nf90_put_att(nc_ofid,SI2TM1S_vid,"standard_name","SI2TM1S")
			ierr = nf90_put_att(nc_ofid,SI2TM1S_vid,"units","mgSi/m^3")
			ierr = nf90_put_att(nc_ofid,SI2TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,SI2TM1S_vid,"type","data")
		case ( "CH41TM1S")
			ierr = nf90_def_var(nc_ofid,"CH41TM1S",nf90_float,dynm2dn,CH41TM1S_vid)
			ierr = nf90_put_att(nc_ofid,CH41TM1S_vid,"long_name","CH41TM1S")
			ierr = nf90_put_att(nc_ofid,CH41TM1S_vid,"standard_name","CH41TM1S")
			ierr = nf90_put_att(nc_ofid,CH41TM1S_vid,"units","mgO2/m^3")
			ierr = nf90_put_att(nc_ofid,CH41TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,CH41TM1S_vid,"type","data")
		case ( "CH42TM1S")
			ierr = nf90_def_var(nc_ofid,"CH42TM1S",nf90_float,dynm2dn,CH42TM1S_vid)
			ierr = nf90_put_att(nc_ofid,CH42TM1S_vid,"long_name","CH42TM1S")
			ierr = nf90_put_att(nc_ofid,CH42TM1S_vid,"standard_name","CH42TM1S")
			ierr = nf90_put_att(nc_ofid,CH42TM1S_vid,"units","????")
			ierr = nf90_put_att(nc_ofid,CH42TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,CH42TM1S_vid,"type","data")
		case ( "HS1TM1S")
			ierr = nf90_def_var(nc_ofid,"HS1TM1S",nf90_float,dynm2dn,HS1TM1S_vid)
			ierr = nf90_put_att(nc_ofid,HS1TM1S_vid,"long_name","HS1TM1S")
			ierr = nf90_put_att(nc_ofid,HS1TM1S_vid,"standard_name","HS1TM1S")
			ierr = nf90_put_att(nc_ofid,HS1TM1S_vid,"units","mgO2/m^3")
			ierr = nf90_put_att(nc_ofid,HS1TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,HS1TM1S_vid,"type","data")
		case ( "HS2TM1S")
			ierr = nf90_def_var(nc_ofid,"HS2TM1S",nf90_float,dynm2dn,HS2TM1S_vid)
			ierr = nf90_put_att(nc_ofid,HS2TM1S_vid,"long_name","HS2TM1S")
			ierr = nf90_put_att(nc_ofid,HS2TM1S_vid,"standard_name","HS2TM1S")
			ierr = nf90_put_att(nc_ofid,HS2TM1S_vid,"units","mgO2/m^3")
			ierr = nf90_put_att(nc_ofid,HS2TM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,HS2TM1S_vid,"type","data")
		case ( "CPOP")
			ierr = nf90_def_var(nc_ofid,"CPOP",nf90_float,dynm3dn_3,CPOP_vid)		!3 layers
			ierr = nf90_put_att(nc_ofid,CPOP_vid,"long_name","CPOP")
			ierr = nf90_put_att(nc_ofid,CPOP_vid,"standard_name","CPOP")
			ierr = nf90_put_att(nc_ofid,CPOP_vid,"units","mgP/m^3")
			ierr = nf90_put_att(nc_ofid,CPOP_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,CPOP_vid,"type","data")
		case ( "CPON")
			ierr = nf90_def_var(nc_ofid,"CPON",nf90_float,dynm3dn_3,CPON_vid)		!3 layers
			ierr = nf90_put_att(nc_ofid,CPON_vid,"long_name","CPON")
			ierr = nf90_put_att(nc_ofid,CPON_vid,"standard_name","CPON")
			ierr = nf90_put_att(nc_ofid,CPON_vid,"units","mgN/m^3")
			ierr = nf90_put_att(nc_ofid,CPON_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,CPON_vid,"type","data")
		case ( "CPOC")
			ierr = nf90_def_var(nc_ofid,"CPOC",nf90_float,dynm3dn_3,CPOC_vid)		!3 layers
			ierr = nf90_put_att(nc_ofid,CPOC_vid,"long_name","CPOC")
			ierr = nf90_put_att(nc_ofid,CPOC_vid,"standard_name","CPOC")
			ierr = nf90_put_att(nc_ofid,CPOC_vid,"units","mgC/m^3")
			ierr = nf90_put_att(nc_ofid,CPOC_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,CPOC_vid,"type","data")
		case ( "CPOS")
			ierr = nf90_def_var(nc_ofid,"CPOS",nf90_float,dynm2dn,CPOS_vid)
			ierr = nf90_put_att(nc_ofid,CPOS_vid,"long_name","CPOS")
			ierr = nf90_put_att(nc_ofid,CPOS_vid,"standard_name","CPOS")
			ierr = nf90_put_att(nc_ofid,CPOS_vid,"units","mgSi/m^3")
			ierr = nf90_put_att(nc_ofid,CPOS_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,CPOS_vid,"type","data")
		case ( "HSED1")
			ierr = nf90_def_var(nc_ofid,"HSED1",nf90_float,dynm2dn,HSED1_vid)
			ierr = nf90_put_att(nc_ofid,HSED1_vid,"long_name","HSED1")
			ierr = nf90_put_att(nc_ofid,HSED1_vid,"standard_name","HSED1")
			ierr = nf90_put_att(nc_ofid,HSED1_vid,"units","m")
			ierr = nf90_put_att(nc_ofid,HSED1_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,HSED1_vid,"type","data")
		case ( "BENSTRTM1S")
			ierr = nf90_def_var(nc_ofid,"BENSTRTM1S",nf90_float,dynm2dn,BENSTRTM1S_vid)
			ierr = nf90_put_att(nc_ofid,BENSTRTM1S_vid,"long_name","HSED1")
			ierr = nf90_put_att(nc_ofid,BENSTRTM1S_vid,"standard_name","HSED1")
			ierr = nf90_put_att(nc_ofid,BENSTRTM1S_vid,"units","-")
			ierr = nf90_put_att(nc_ofid,BENSTRTM1S_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,BENSTRTM1S_vid,"type","data")

		case default
			!do nothing
end select
ENDIF

	IF (ZOO_CALC) THEN !Outputted ZOO
		select case(trim(cdf_vdp(i)))
		case ( "SZ")
			ierr = nf90_def_var(nc_ofid,"SZ",nf90_float,dynm3dn_lay,SZ_vid)
			ierr = nf90_put_att(nc_ofid,SZ_vid,"long_name","SZ")
			ierr = nf90_put_att(nc_ofid,SZ_vid,"standard_name","SZ")
			ierr = nf90_put_att(nc_ofid,SZ_vid,"units","gC/m^3")
			ierr = nf90_put_att(nc_ofid,SZ_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,SZ_vid,"type","data")
		case ( "LZ")
			ierr = nf90_def_var(nc_ofid,"LZ",nf90_float,dynm3dn_lay,LZ_vid)
			ierr = nf90_put_att(nc_ofid,LZ_vid,"long_name","LZ")
			ierr = nf90_put_att(nc_ofid,LZ_vid,"standard_name","LZ")
			ierr = nf90_put_att(nc_ofid,LZ_vid,"units","gC/m^3")
			ierr = nf90_put_att(nc_ofid,LZ_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,LZ_vid,"type","data")
		case default
			!do nothing
end select
ENDIF

	IF (SOLIDS_CALC) THEN 
		select case(trim(cdf_vdp(i)))
		case ("TSS")
			ierr = nf90_def_var(nc_ofid,"TSS",nf90_float,dynm3dn_lay,TSS_vid)
			ierr = nf90_put_att(nc_ofid,TSS_vid,"long_name","TSS")
			ierr = nf90_put_att(nc_ofid,TSS_vid,"standard_name","TSS")
			ierr = nf90_put_att(nc_ofid,TSS_vid,"units","g meters-3")
			ierr = nf90_put_att(nc_ofid,TSS_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,TSS_vid,"type","data")
		case ("KESS")
			ierr = nf90_def_var(nc_ofid,"KESS",nf90_float,dynm3dn_lay,KESS_vid)
			ierr = nf90_put_att(nc_ofid,KESS_vid,"long_name","KESS")
			ierr = nf90_put_att(nc_ofid,KESS_vid,"standard_name","KESS")
			ierr = nf90_put_att(nc_ofid,KESS_vid,"units","meter-1")
			ierr = nf90_put_att(nc_ofid,KESS_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,KESS_vid,"type","data")
		case ( "TURB")
			ierr = nf90_def_var(nc_ofid,"TURB",nf90_float,dynm3dn_lay,TURB_vid)
			ierr = nf90_put_att(nc_ofid,TURB_vid,"long_name","TURB")
			ierr = nf90_put_att(nc_ofid,TURB_vid,"standard_name","TURB")
			ierr = nf90_put_att(nc_ofid,TURB_vid,"units","NTU")
			ierr = nf90_put_att(nc_ofid,TURB_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,TURB_vid,"type","data")
		case ( "SSI")
			ierr = nf90_def_var(nc_ofid,"SSI",nf90_float,dynm3dn_lay,SSI_vid)
			ierr = nf90_put_att(nc_ofid,SSI_vid,"long_name","SSI")
			ierr = nf90_put_att(nc_ofid,SSI_vid,"standard_name","SSI")
			ierr = nf90_put_att(nc_ofid,SSI_vid,"units","g meters-3")
			ierr = nf90_put_att(nc_ofid,SSI_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,SSI_vid,"type","data")
	    case ( "PARAD")
			ierr = nf90_def_var(nc_ofid,"PARAD",nf90_float,dynm3dn_lay,PARAD_vid)
			ierr = nf90_put_att(nc_ofid,PARAD_vid,"long_name","PARAD")
			ierr = nf90_put_att(nc_ofid,PARAD_vid,"standard_name","PARAD")
			ierr = nf90_put_att(nc_ofid,PARAD_vid,"units","Einstein meters-2 day-1")
			ierr = nf90_put_att(nc_ofid,PARAD_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,PARAD_vid,"type","data")
		case default
			!do nothing
end select
ENDIF

	 IF (SAV_CALC) THEN !Outputted SAV
		select case(trim(cdf_vdp(i)))
		case ( "LEAF")
			ierr = nf90_def_var(nc_ofid,"LEAF",nf90_float,dynm3dn_nsav,LEAF_vid)
			ierr = nf90_put_att(nc_ofid,LEAF_vid,"long_name","LEAF")
			ierr = nf90_put_att(nc_ofid,LEAF_vid,"standard_name","LEAF")
			ierr = nf90_put_att(nc_ofid,LEAF_vid,"units","gC/m^2")
			ierr = nf90_put_att(nc_ofid,LEAF_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,LEAF_vid,"type","data")
		case ( "STEM")
			ierr = nf90_def_var(nc_ofid,"STEM",nf90_float,dynm3dn_nsav,STEM_vid)
			ierr = nf90_put_att(nc_ofid,STEM_vid,"long_name","STEM")
			ierr = nf90_put_att(nc_ofid,STEM_vid,"standard_name","STEM")
			ierr = nf90_put_att(nc_ofid,STEM_vid,"units","gC/m^2")
			ierr = nf90_put_att(nc_ofid,STEM_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,STEM_vid,"type","data")
		case ( "ROOT")
			ierr = nf90_def_var(nc_ofid,"ROOT",nf90_float,dynm3dn_nsav,ROOT_vid)
			ierr = nf90_put_att(nc_ofid,ROOT_vid,"long_name","ROOT")
			ierr = nf90_put_att(nc_ofid,ROOT_vid,"standard_name","ROOT")
			ierr = nf90_put_att(nc_ofid,ROOT_vid,"units","gC/m^2")
			ierr = nf90_put_att(nc_ofid,ROOT_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,ROOT_vid,"type","data")
		case ( "TUBER")
			ierr = nf90_def_var(nc_ofid,"TUBER",nf90_float,dynm3dn_nsav,TUBER_vid)
			ierr = nf90_put_att(nc_ofid,TUBER_vid,"long_name","TUBER")
			ierr = nf90_put_att(nc_ofid,TUBER_vid,"standard_name","TUBER")
			ierr = nf90_put_att(nc_ofid,TUBER_vid,"units","gC/m^2")
			ierr = nf90_put_att(nc_ofid,TUBER_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,TUBER_vid,"type","data")
	    case ( "EP")
			ierr = nf90_def_var(nc_ofid,"EP",nf90_float,dynm3dn_nsav,EP_vid)
			ierr = nf90_put_att(nc_ofid,EP_vid,"long_name","EP")
			ierr = nf90_put_att(nc_ofid,EP_vid,"standard_name","EP")
			ierr = nf90_put_att(nc_ofid,EP_vid,"units","gC-epiphytes/gm^2-leaf")
			ierr = nf90_put_att(nc_ofid,EP_vid,"grid","fvcom_grid")
			ierr = nf90_put_att(nc_ofid,EP_vid,"type","data")
		case default
			!do nothing
end select
ENDIF
end do


!--Exit Define Mode
   ierr = nf90_enddef(nc_ofid)
   ierr = nf90_close(nc_ofid)

   end if !(msr)

!==============================================================================|
!  WRITE VARIABLES TO FILE                                                     |
!==============================================================================|
   if(msr)then
     ierr = nf90_open(trim(cdfname),nf90_write,nc_ofid)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"file open error",ipt)
     end if
   end if

   !!!====Longitude at Nodes (LON) ==========================!
   !i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   !call putvar(i1,i2,MLOC,mgl,1,1,"n",vx+vxmin,nc_ofid,lon_vid,myid,nprocs,ipt)

   !!!====Latitude  at Nodes (LAT) ==========================!
   !i1 = lbound(vy,1) ; i2 = ubound(vy,1)
   !call putvar(i1,i2,MLOC,mgl,1,1,"n",vy+vymin,nc_ofid,lat_vid,myid,nprocs,ipt)

   !!====Number of Processors (NPROCS) =====================!
   if(msr)then
   ierr = nf90_put_var(nc_ofid,nprocs_vid,nprocs)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing nprocs variable to netcdf",ipt)
   end if
#  if defined (MULTIPROCESSOR)
   ierr = nf90_put_var(nc_ofid,partition_vid,el_pid)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing el_pid variable to netcdf",ipt)
   end if
#  endif
   end if

   !!!====Initial Density Field==============================!
   !i1 = lbound(rho1,1) ; i2 = ubound(rho1,1)
   !call putvar(i1,i2,MLOC,mgl,kb,kb-1,"n",rho1,nc_ofid,idens_vid,myid,nprocs,ipt)


   !!====X Grid Coordinate at Nodes (VX)====================!

   i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   call putvar(i1,i2,MLOC,mgl,1,1,"n",vx+vxmin,nc_ofid,x_vid,myid,nprocs,ipt)

   !!====Y Grid Coordinate at Nodes (VY)====================!
   i1 = lbound(vy,1) ; i2 = ubound(vy,1)
   call putvar(i1,i2,MLOC,mgl,1,1,"n",vy+vymin,nc_ofid,y_vid,myid,nprocs,ipt)

   !!====Bathymetry at Nodes (H)============================!
   i1 = lbound(h,1) ; i2 = ubound(h,1)
   call putvar(i1,i2,MLOC,mgl,1,1,"n",h,nc_ofid,h_vid,myid,nprocs,ipt)

   !!====Nodes surrounding each Element (NV)================!
   allocate(tmp(0:NTLOC,3))
   if(serial)then
     tmp(0:NTLOC,1:3) = real(nv(0:NTLOC,1:3),sp)
   end if
#  if defined (MULTIPROCESSOR)
   if(par)then
   do j=1,3
   do i=1,NLOC
     tmp(i,j) = real(ngid(nv(i,j)),sp)
   end do
   end do
   end if
#  endif
   i1 = lbound(tmp,1) ; i2 = ubound(tmp,1)
   call putvar(i1,i2,NLOC,ngl,3,3,"e",tmp,nc_ofid,nv_vid,myid,nprocs,ipt)
   deallocate(tmp)


   !!!====Momentum Stencil Interpolation Coefficients========!
   !i1 = lbound(a1u,1) ; i2 = ubound(a1u,1)
   !call putvar(i1,i2,NLOC,ngl,4,4,"e",a1u,nc_ofid,a1u_vid,myid,nprocs,ipt)
   !i1 = lbound(a2u,1) ; i2 = ubound(a2u,1)
   !call putvar(i1,i2,NLOC,ngl,4,4,"e",a2u,nc_ofid,a2u_vid,myid,nprocs,ipt)

   !!!====Element Based Interpolation Coefficients===========!
   !i1 = lbound(aw0,1) ; i2 = ubound(aw0,1)
   !call putvar(i1,i2,NLOC,ngl,3,3,"e",aw0,nc_ofid,aw0_vid,myid,nprocs,ipt)
   !i1 = lbound(awx,1) ; i2 = ubound(awx,1)
   !call putvar(i1,i2,NLOC,ngl,3,3,"e",awx,nc_ofid,awx_vid,myid,nprocs,ipt)
   !i1 = lbound(awy,1) ; i2 = ubound(awy,1)
   !call putvar(i1,i2,NLOC,ngl,3,3,"e",awy,nc_ofid,awy_vid,myid,nprocs,ipt)

   !!====Sigma Layers (zz)==================================!
   if(msr)then
   allocate(tvec(kbm1))
   tvec(1:kbm1) = zz(1:kbm1)
   !tvec(1:kbm1) = zz(1,1:kbm1)   !T.W., changed back to v2.3, for WQ linkage purpose
   ierr = nf90_put_var(nc_ofid,siglay_vid,tvec)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing variable to netcdf",ipt)
   end if
   deallocate(tvec)

   !allocate(tvec(kbm1))
   !tvec(1:kbm1) = z(2:kb)
   !!tvec(1:kbm1) = z(1,2:kb)  !T.W., changed back v2.3
   !ierr = nf90_put_var(nc_ofid,siglay_shift_vid,tvec)
   !if(ierr /= nf90_noerr)then
   !  call handle_ncerr(ierr,"error writing variable to netcdf",ipt)
   !end if
   !deallocate(tvec)

   allocate(tvec(kb))
   tvec(1:kb) = z(1:kb)
!   tvec(1:kb) = z(1,1:kb)  !T.W., changed back to v2.3
   ierr = nf90_put_var(nc_ofid,siglev_vid,tvec)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing variable to netcdf",ipt)
   end if
   deallocate(tvec)
   endif

#if defined (NEWSIGMA)

		allocate(tmp(0:MTLOC,1:KB))
		if(serial)then
			tmp(0:MTLOC,1:KB) = real(Z2D(0:MTLOC,1:KB),sp)
		end if
#  	if defined (MULTIPROCESSOR)
		if(par)then
			do j=1,KB
			do i=1,MTLOC
				tmp(i,j) =real(Z2D(i,j),sp)
			end do
			end do
		end if
#	endif
		i1 = lbound(tmp,1) ; i2 = ubound(tmp,1)
		call putvar(i1,i2,MLOC,MGL,KB,KB,"n",tmp,nc_ofid,siglev3dn_vid,myid,nprocs,ipt)
		deallocate(tmp)

#endif

!==============================================================================|
!  close the file                                                              |
!==============================================================================|

   if(msr) ierr = nf90_close(nc_ofid)

   return
   end subroutine write_netcdf_setup
!==============================================================================|


   subroutine out_netcdf_data(ipt)
!==============================================================================|
!   Write Time Dependent NetCDF Data to File                                   |
!==============================================================================|
   use netcdf
#  if defined (WET_DRY)
   use mod_wd, only:iswetn
#  endif
   use mod_lims, only: NLOC, MLOC
   USE MOD_CONTROL, ONLY: MSR
   USE MOD_HYDROVARS, ONLY: UU,VV,WTS,thour, KH, EL,h,d!,Z,ZZ
   USE MOD_WQM, only: t,salt,CCHL1,total_netPP,DOXG, LDOC, RDOC, 						&
					  LPOC, RPOC, B1,B2, NH4, NO3, LDON, RDON, LPON, RPON, 				&
					  PO4,TDIC, TALK, pH, pCO2, P1, P2, BM1, BM2, NL1, NL2,PL1, PL2 	&
					  ,PR1, PR2, FI1, FI2, IK1, IK2, GPP, NPP, DDOC, 					&
					  DCOD, NITRIF, DOPR, DORALG, RESP, FTCOD, COD, 					&
					  REAERDO, BFLUX, DDOC, DCOD, NITRIF, DORALG, RESP,FTCOD, 			&
					  SEDIMENT_CALC, SOLIDS_CALC, SSI, SAV_CALC !,&
				!	,ADVFLUX	&
				!   ,VDIFFLUX

	USE MOD_SED, ONLY: JPOC, JPON, JPOP, JPOS, CTEMP,SODTM1S, JNH4TM1S,JNO3TM1S, BENDEN, JCH4TM1S 	&
				,JCH4GTM1S,JHSTM1S,JPO4TM1S,JSITM1S	,NH41TM1S,NH42TM1S,NO31TM1S ,NO32TM1S 			&
				,PO41TM1S,PO42TM1S	,SI1TM1S ,SI2TM1S	,CH41TM1S	,CH42TM1S	,HS1TM1S			&
				,HS2TM1S,CPOP,CPON,CPOC	,CPOS,HSED1	,BENSTRTM1S

	USE MOD_CO2SYS, ONLY: 	DICUPT, DICBMP, DICPRD, DICMNL, DICDEN, DICGAS, &
							DICSED, ALKNH4, ALKNO3, ALKNIT,ALKDEN, ALKREM,	&
							ALKNH4SED,ALKNO3SED

	USE MOD_ZOOP, ONLY: B1SZ, B2SZ, B1LZ, B2LZ, DOSZ, DOLZ, SZ, LZ

	USE MOD_OWQ, ONLY: IAVG, TSS, KESS, TURB, PARAD !
	!USE MOD_OWQ, ONLY: IAVG, KESS, TURB, PARAD
	USE MOD_WQMINIT, ONLY: CARBONATE_CALC,COD_CALC, ZOO_CALC

	USE MOD_SED, ONLY: SODTM1S
	USE MOD_SIZES, only: MGL, NGL, NSAVM!
	USE MOD_LIMS, only: MYID,NPROCS, MTLOC,IINT,KB!,

	USE MOD_SAV, ONLY: LEAF, STEM, ROOT, TUBER, EP

   implicit none
   integer,intent(in) :: ipt
   integer :: i,ierr,i1,i2,k,icheck
   integer :: dims(1)
!  real*4, allocatable :: ftemp(:)
   real(sp), allocatable :: ftemp(:)

!==============================================================================|

!--Update Counter
   out_cnt = out_cnt + 1
   stck_cnt = stck_cnt + 1

!--Write Header information if first output of file
   if(cdf_stk == 0)then
     if(out_cnt == 1) call write_netcdf_setup(1,IPT=IPT)
   else
     icheck = mod(out_cnt-1,cdf_stk)
     if(icheck ==0 .or. out_cnt==1)call write_netcdf_setup((out_cnt-1)/cdf_stk+1,IPT=IPT)
   endif

!--Open File
   if(msr)then
     ierr = nf90_open(trim(cdfname),nf90_write,nc_ofid)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"file open error",ipt=IPT)
     end if
   end if

!--Dump Time/IINT to File
   dims(1) = stck_cnt
!
   IF(msr) then
     ierr    = nf90_put_var(nc_ofid,iint_vid,iint,START=dims)
     ierr    = nf90_put_var(nc_ofid,time_vid,thour*3600.,START=dims)
   ENDIF

!--Write Variables to File
   if(msr) write(ipt,*)'dumping to netcdf file: ',trim(cdfname),stck_cnt

   do i=1,nout_vars

	 !if(msr) write(ipt,*) 'variable ', i, '=', trim(cdf_vdp(i)), 'nout_vars=', nout_vars

     select case(trim(cdf_vdp(i)))

		case("u")  !!===============U=======================================!
			i1 = lbound(uu,1) ; i2 = ubound(uu,1)
			call putvar(i1,i2,NLOC,ngl,kb,kb-1,"e",uu,nc_ofid,u_vid,myid,nprocs,ipt)
		case("v")  !!===============V=======================================!
			i1 = lbound(vv,1) ; i2 = ubound(vv,1)
			call putvar(i1,i2,NLOC,ngl,kb,kb-1,"e",vv,nc_ofid,v_vid,myid,nprocs,ipt)

		case("wts") !!===============WTS=====================================!
			i1 = lbound(wts,1) ; i2 = ubound(wts,1)
			call putvar(i1,i2,MLOC,mgl,kb,kb,"n",wts,nc_ofid,wts_vid,myid,nprocs,ipt)

		!case("dtfa") !!===============DTFA===================================!
			!  i1 = lbound(dtfa,1) ; i2 = ubound(dtfa,1)
			!  call putvar(i1,i2,MLOC,mgl,1,1,"n",dtfa,nc_ofid,dtfa_vid,myid,nprocs,ipt)
		!case("xflux_obc") !!===============XFLUX_OBC============================!
		!	 !  Change to "b" and other stuff, like changed IOBCN_RGL to IOBCN+1
		!       i1 = lbound(xflux_obc,1) ; i2 = ubound(xflux_obc,1)
		!       call putvar(i1,i2,IOBCN,IOBCN_GL,kb,kb-1,"b",XFLUX_OBC,nc_ofid,&
		!        xflux_obc_vid,myid,nprocs,ipt)

		!     case("uard_obcn") !!===============UARD_OBCN===========================!
		!       i1 = lbound(uard_obcn,1) ; i2 = ubound(uard_obcn,1)
		!       call putvar(i1,i2,IOBCN,IOBCN_GL,1,1,"b",UARD_OBCN,nc_ofid,&
		!        uard_obcn_vid,myid,nprocs,ipt)

		!case("km") !!===============KM======================================!
			!  i1 = lbound(km,1) ; i2 = ubound(km,1)
			!  call putvar(i1,i2,MLOC,mgl,kb,kb,"n",km,nc_ofid,km_vid,myid,nprocs,ipt)

		case("kh") !!===============KH======================================!
			i1 = lbound(kh,1) ; i2 = ubound(kh,1)
			call putvar(i1,i2,MLOC,mgl,kb,kb,"n",kh,nc_ofid,kh_vid,myid,nprocs,ipt)

		!case("viscofh")!!===========VISCOFH=================================!
			!  i1 = lbound(viscofh,1); i2 = ubound(viscofh,1)
			!  call putvar(i1,i2,MLOC,mgl,kb,kb-1,"n",viscofh,nc_ofid,viscofh_vid,myid,nprocs,ipt)

		case("el") !!===============EL======================================!
			i1 = lbound(el,1) ; i2 = ubound(el,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",el,nc_ofid,el_vid,myid,nprocs,ipt)
		case("d") !!===============D=======================================!
			i1 = lbound(d,1) ; i2 = ubound(d,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",d,nc_ofid,d_vid,myid,nprocs,ipt)
		case("temp") !!===============T1======================================!
			i1 = lbound(t,1) ; i2 = ubound(t,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",t,nc_ofid,t1_vid,myid,nprocs,ipt)
		case("salt") !!===============S1======================================!
			i1 = lbound(salt,1) ; i2 = ubound(salt,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",salt,nc_ofid,s1_vid,myid,nprocs,ipt)

		case("wd") !!===============WETDRY==================================!
# if defined (WET_DRY)
			allocate(ftemp(0:MTLOC)) ; ftemp = iswetn
			i1 = lbound(ftemp,1) ; i2 = ubound(ftemp,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",ftemp,nc_ofid,wd_vid,myid,nprocs,ipt)
			deallocate(ftemp)
# endif

!example     case("s1") !!===============S1======================================!
!example       i1 = lbound(s1,1) ; i2 = ubound(s1,1)
!example       call putvar(i1,i2,MLOC,mgl,kb,kb-1,"n",s1,nc_ofid,s1_vid,myid,nprocs,ipt)

     !new variable output - add a new variable (e.g. 'var') to output
     !1.) copy example section above
     !2.) modify case for your variable 'case("var")'
     !3.) modify bounds for your variable
     !4.) modify putvar for your variable by finding a putvar for a variable
     !    with same dimensions and type ("e" or "n")
     !5.) modify variable vid with your variables vid (e.g. "var_vid")

		case ("CCHL1")
			i1 = lbound(CCHL1,1) ; i2 = ubound(CCHL1,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",CCHL1,nc_ofid,CCHL1_vid,myid,nprocs,ipt)


		case ("total_netPP")
			i1 = lbound(total_netPP,1) ; i2 = ubound(total_netPP,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",total_netPP,nc_ofid,total_netPP_vid,myid,nprocs,ipt)


		case ( "DOXG")
			i1 = lbound(DOXG,1) ; i2 = ubound(DOXG,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DOXG,nc_ofid,DOXG_vid,myid,nprocs,ipt)

		case (  "LDOC")
			i1 = lbound(LDOC,1) ; i2 = ubound(LDOC,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",LDOC,nc_ofid,LDOC_vid,myid,nprocs,ipt)

		case ( "RDOC")
			i1 = lbound(RDOC,1) ; i2 = ubound(RDOC,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",RDOC,nc_ofid,RDOC_vid,myid,nprocs,ipt)


		case ( "LPOC")
			i1 = lbound(LPOC,1) ; i2 = ubound(LPOC,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",LPOC,nc_ofid,LPOC_vid,myid,nprocs,ipt)


		case ( "RPOC")
			i1 = lbound(RPOC,1) ; i2 = ubound(RPOC,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",RPOC,nc_ofid,RPOC_vid,myid,nprocs,ipt)

		case ( "B1")
			i1 = lbound(B1,1) ; i2 = ubound(B1,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",B1,nc_ofid,B1_vid,myid,nprocs,ipt)


		case ( "B2")
			i1 = lbound(B2,1) ; i2 = ubound(B1,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",B2,nc_ofid,B2_vid,myid,nprocs,ipt)

		case ( "NH4")
			i1 = lbound(NH4,1) ; i2 = ubound(NH4,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",NH4,nc_ofid,NH4_vid,myid,nprocs,ipt)


		case ( "NO3")
			i1 = lbound(NO3,1) ; i2 = ubound(NO3,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",NO3,nc_ofid,NO3_vid,myid,nprocs,ipt)


		case ( "LDON")
			i1 = lbound(LDON,1) ; i2 = ubound(LDON,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",LDON,nc_ofid,LDON_vid,myid,nprocs,ipt)

		case ( "RDON")
			i1 = lbound(RDON,1) ; i2 = ubound(RDON,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",RDON,nc_ofid,RDON_vid,myid,nprocs,ipt)

		case ( "LPON")
			i1 = lbound(LPON,1) ; i2 = ubound(LPON,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",LPON,nc_ofid,LPON_vid,myid,nprocs,ipt)


		case ( "RPON")
			i1 = lbound(RPON,1) ; i2 = ubound(RPON,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",RPON,nc_ofid,RPON_vid,myid,nprocs,ipt)


		case ( "PO4")
			i1 = lbound(PO4,1) ; i2 = ubound(PO4,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",PO4,nc_ofid,PO4_vid,myid,nprocs,ipt)


		case ( "P1")
			i1 = lbound(P1,1) ; i2 = ubound(P1,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",P1,nc_ofid,P1_vid,myid,nprocs,ipt)

		case ( "P2")
			i1 = lbound(P2,1) ; i2 = ubound(P2,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",P2,nc_ofid,P2_vid,myid,nprocs,ipt)

		case ( "BM1")
			i1 = lbound(BM1,1) ; i2 = ubound(BM1,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",BM1,nc_ofid,BM1_vid,myid,nprocs,ipt)

		case ( "BM2")
			i1 = lbound(BM2,1) ; i2 = ubound(BM2,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",BM2,nc_ofid,BM2_vid,myid,nprocs,ipt)

		case ( "NL1")
			i1 = lbound(NL1,1) ; i2 = ubound(NL1,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",NL1,nc_ofid,NL1_vid,myid,nprocs,ipt)

		case ( "NL2")
			i1 = lbound(NL2,1) ; i2 = ubound(NL2,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",NL2,nc_ofid,NL2_vid,myid,nprocs,ipt)

		case ( "PL1")
			i1 = lbound(PL1,1) ; i2 = ubound(PL1,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",PL1,nc_ofid,PL1_vid,myid,nprocs,ipt)

		case ( "PL2")
			i1 = lbound(PL2,1) ; i2 = ubound(PL2,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",PL2,nc_ofid,PL2_vid,myid,nprocs,ipt)
!The following are related to Zoop
		case ( "B1SZ")
			i1 = lbound(B1SZ,1) ; i2 = ubound(B1SZ,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",B1SZ,nc_ofid,B1SZ_vid,myid,nprocs,ipt)

		case ( "B2SZ")
			i1 = lbound(B2SZ,1) ; i2 = ubound(B2SZ,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",B2SZ,nc_ofid,B2SZ_vid,myid,nprocs,ipt)

		case ( "B1LZ")
			i1 = lbound(B1LZ,1) ; i2 = ubound(B1LZ,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",B1LZ,nc_ofid,B1LZ_vid,myid,nprocs,ipt)

		case ( "B2LZ")
			i1 = lbound(B2LZ,1) ; i2 = ubound(B2LZ,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",B2LZ,nc_ofid,B2LZ_vid,myid,nprocs,ipt)
!--
		case ( "PR1")
			i1 = lbound(PR1,1) ; i2 = ubound(PR1,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",PR1,nc_ofid,PR1_vid,myid,nprocs,ipt)

		case ( "PR2")
			i1 = lbound(PR2,1) ; i2 = ubound(PR2,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",PR2,nc_ofid,PR2_vid,myid,nprocs,ipt)

		case ( "IAVG")
			i1 = lbound(IAVG,1) ; i2 = ubound(IAVG,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",IAVG,nc_ofid,IAVG_vid,myid,nprocs,ipt)

		case ( "FI1")
			i1 = lbound(FI1,1) ; i2 = ubound(FI1,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",FI1,nc_ofid,FI1_vid,myid,nprocs,ipt)

		case ( "FI2")
			i1 = lbound(FI2,1) ; i2 = ubound(FI2,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",FI2,nc_ofid,FI2_vid,myid,nprocs,ipt)

		case ( "IK1")
			i1 = lbound(IK1,1) ; i2 = ubound(IK1,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",IK1,nc_ofid,IK1_vid,myid,nprocs,ipt)

		case ( "IK2")
			i1 = lbound(IK2,1) ; i2 = ubound(IK2,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",IK2,nc_ofid,IK2_vid,myid,nprocs,ipt)

		case ( "GPP")
			i1 = lbound(GPP,1) ; i2 = ubound(GPP,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",GPP,nc_ofid,GPP_vid,myid,nprocs,ipt)

		case ( "NPP")
			i1 = lbound(NPP,1) ; i2 = ubound(NPP,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",NPP,nc_ofid,NPP_vid,myid,nprocs,ipt)

		case ( "DDOC")
			i1 = lbound(DDOC,1) ; i2 = ubound(DDOC,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DDOC,nc_ofid,DDOC_vid,myid,nprocs,ipt)

		case ( "DCOD")
			i1 = lbound(DCOD,1) ; i2 = ubound(DCOD,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DCOD,nc_ofid,DCOD_vid,myid,nprocs,ipt)

		case ( "NITRIF")
			i1 = lbound(NITRIF,1) ; i2 = ubound(NITRIF,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",NITRIF,nc_ofid,NITRIF_vid,myid,nprocs,ipt)

		case ( "DOSZ")
			i1 = lbound(DOSZ,1) ; i2 = ubound(DOSZ,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DOSZ,nc_ofid,DOSZ_vid,myid,nprocs,ipt)

		case ( "DOLZ")
			i1 = lbound(DOLZ,1) ; i2 = ubound(DOLZ,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DOLZ,nc_ofid,DOLZ_vid,myid,nprocs,ipt)

		case ( "DOPR")
			i1 = lbound(DOPR,1) ; i2 = ubound(DOPR,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DOPR,nc_ofid,DOPR_vid,myid,nprocs,ipt)

		case ( "DORALG")
			i1 = lbound(DORALG,1) ; i2 = ubound(DORALG,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DORALG,nc_ofid,DORALG_vid,myid,nprocs,ipt)

		case ( "RESP")
			i1 = lbound(RESP,1) ; i2 = ubound(RESP,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",RESP,nc_ofid,RESP_vid,myid,nprocs,ipt)

		case ( "FTCOD")
			IF(COD_CALC)THEN
			i1 = lbound(FTCOD,1) ; i2 = ubound(FTCOD,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",FTCOD,nc_ofid,FTCOD_vid,myid,nprocs,ipt)
			ENDIF
		case ( "COD")
			IF(COD_CALC)THEN
			i1 = lbound(COD,1) ; i2 = ubound(COD,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",COD,nc_ofid,COD_vid,myid,nprocs,ipt)
			ENDIF
		case ( "REAERDO")
			i1 = lbound(REAERDO,1) ; i2 = ubound(REAERDO,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",REAERDO,nc_ofid,REAERDO_vid,myid,nprocs,ipt)

		case ( "BFLUX")
			i1 = lbound(BFLUX,1) ; i2 = ubound(BFLUX,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,9,"n",BFLUX,nc_ofid,BFLUX_vid,myid,nprocs,ipt)	  !9 variables for BFLUX

!carbonate_calc
!The following are related to pH and carbonate_calc
		case ( "TDIC" )
		case ( "TALK" )
		case ( "pH" )
		case ( "pCO2" )
		case ( "DICUPT" )
		case ( "DICBMP" )
		case ( "DICPRD" )
		case ( "DICMNL" )
		case ( "DICDEN" )
		case ( "DICGAS" )
		case ( "DICSED" )
		case ( "ALKNH4" )
		case ( "ALKNO3" )
		case ( "ALKNIT" )
		case ( "ALKDEN" )
		case ( "ALKREM" )
		case ( "ALKNH4SED" )
		case ( "ALKNO3SED" )
!
		!case ( "ADVFLUX")
		!	i1 = lbound(ADVFLUX,1) ; i2 = ubound(ADVFLUX,1)
		!	call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",ADVFLUX,nc_ofid,ADVFLUX_vid,myid,nprocs,ipt)

		!case ( "VDIFFLUX")
		!	i1 = lbound(VDIFFLUX,1) ; i2 = ubound(VDIFFLUX,1)
		!	call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",VDIFFLUX,nc_ofid,VDIFFLUX_vid,myid,nprocs,ipt)
!The following are part of sediment module mod_sed
		case ( "JPOC" )
		case ( "JPON" )
		case ( "JPOP" )
		case ( "JPOS" )
		case ( "CTEMP" )
		case ( "SODTM1S" )
		case ( "JNH4TM1S" )
		case ( "JNO3TM1S" )
		case ( "BENDEN" )
		case ( "JCH4TM1S" )
		case ( "JCH4GTM1S" )
		case ( "JHSTM1S" )
		case ( "JPO4TM1S" )
		case ( "JSITM1S" )
		case ( "NH41TM1S" )
		case ( "NH42TM1S" )
		case ( "NO31TM1S" )
		case ( "NO32TM1S" )
		case ( "PO41TM1S" )
		case ( "PO42TM1S" )
		case ( "SI1TM1S" )
		case ( "SI2TM1S" )
		case ( "CH41TM1S" )
		case ( "CH42TM1S" )
		case ( "HS1TM1S" )
		case ( "HS2TM1S" )
		case ( "CPOP" )
		case ( "CPON" )
		case ( "CPOC" )
		case ( "CPOS" )
		case ( "HSED1" )
		case ( "BENSTRTM1S" )

		case ( "TSS" )
		case ( "SSI" )
		case ( "KESS" )
		case ( "TURB" )
		case ( "PARAD" )

		case ( "LEAF")
		case ( "STEM")
		case ( "ROOT")
		case ( "TUBER")
		case ( "EP")

		case ( "SZ")
		case ( "LZ")

		case default
			if(msr)then
				write(ipt,*)'variable',cdf_vdp(i),' not set up for netcdf output'
				write(ipt,*)'modify module MOD_NCDIO.f'
				call pstop
			end if
     end select

	 IF(CARBONATE_CALC)THEN
		select case(trim(cdf_vdp(i)))

		case ( "TDIC")
			i1 = lbound(TDIC,1) ; i2 = ubound(TDIC,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",TDIC,nc_ofid,TDIC_vid,myid,nprocs,ipt)


		case ( "TALK")
			i1 = lbound(TALK,1) ; i2 = ubound(TALK,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",TALK,nc_ofid,TALK_vid,myid,nprocs,ipt)

		case ( "pH")
			i1 = lbound(pH,1) ; i2 = ubound(pH,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",pH,nc_ofid,pH_vid,myid,nprocs,ipt)

		case ( "pCO2")

			i1 = lbound(pCO2,1) ; i2 = ubound(pCO2,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",pCO2,nc_ofid,pCO2_vid,myid,nprocs,ipt)

		case ( "DICUPT")
			i1 = lbound(DICUPT,1) ; i2 = ubound(DICUPT,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DICUPT,nc_ofid,DICUPT_vid,myid,nprocs,ipt)

		case ( "DICBMP")
			i1 = lbound(DICBMP,1) ; i2 = ubound(DICBMP,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DICBMP,nc_ofid,DICBMP_vid,myid,nprocs,ipt)

		case ( "DICPRD")
			i1 = lbound(DICPRD,1) ; i2 = ubound(DICPRD,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DICPRD,nc_ofid,DICPRD_vid,myid,nprocs,ipt)

		case ( "DICMNL")
			i1 = lbound(DICMNL,1) ; i2 = ubound(DICMNL,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DICMNL,nc_ofid,DICMNL_vid,myid,nprocs,ipt)

		case ( "DICDEN")
			i1 = lbound(DICDEN,1) ; i2 = ubound(DICDEN,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DICDEN,nc_ofid,DICDEN_vid,myid,nprocs,ipt)

		case ( "DICGAS")

			i1 = lbound(DICGAS,1) ; i2 = ubound(DICGAS,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DICGAS,nc_ofid,DICGAS_vid,myid,nprocs,ipt)	!

		case ( "DICSED")
			i1 = lbound(DICSED,1) ; i2 = ubound(DICSED,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",DICSED,nc_ofid,DICSED_vid,myid,nprocs,ipt)	!

		case ( "ALKNH4")
			i1 = lbound(ALKNH4,1) ; i2 = ubound(ALKNH4,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",ALKNH4,nc_ofid,ALKNH4_vid,myid,nprocs,ipt)

		case ( "ALKNO3")
			i1 = lbound(ALKNO3,1) ; i2 = ubound(ALKNO3,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",ALKNO3,nc_ofid,ALKNO3_vid,myid,nprocs,ipt)

		case ( "ALKNIT")
			i1 = lbound(ALKNIT,1) ; i2 = ubound(ALKNIT,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",ALKNIT,nc_ofid,ALKNIT_vid,myid,nprocs,ipt)

		case ( "ALKDEN")
			i1 = lbound(ALKDEN,1) ; i2 = ubound(ALKDEN,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",ALKDEN,nc_ofid,ALKDEN_vid,myid,nprocs,ipt)

		case ( "ALKREM")
			i1 = lbound(ALKREM,1) ; i2 = ubound(ALKREM,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",ALKREM,nc_ofid,ALKREM_vid,myid,nprocs,ipt)

		case ( "ALKNH4SED")
			i1 = lbound(ALKNH4SED,1) ; i2 = ubound(ALKNH4SED,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",ALKNH4SED,nc_ofid,ALKNH4SED_vid,myid,nprocs,ipt)

		case ( "ALKNO3SED")
			i1 = lbound(ALKNO3SED,1) ; i2 = ubound(ALKNO3SED,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",ALKNO3SED,nc_ofid,ALKNO3SED_vid,myid,nprocs,ipt)

		case default
			!do nothing
		end select

	 ENDIF
	 !write(*,*)'i=',i,'cdf_vdp(i)=',trim(cdf_vdp(i))

	 IF(SEDIMENT_CALC)THEN
		select case(trim(cdf_vdp(i)))
	 	case ("JPOC")
			i1 = lbound(JPOC,1) ; i2 = ubound(JPOC,1)
			call putvar(i1,i2,MLOC,mgl,3,3,"n",JPOC,nc_ofid,JPOC_vid,myid,nprocs,ipt)

		case ("JPON")
			i1 = lbound(JPON,1) ; i2 = ubound(JPON,1)
			call putvar(i1,i2,MLOC,mgl,3,3,"n",JPON,nc_ofid,JPON_vid,myid,nprocs,ipt)

		case ( "JPOP")
			i1 = lbound(JPOP,1) ; i2 = ubound(JPOP,1)
			call putvar(i1,i2,MLOC,mgl,3,3,"n",JPOP,nc_ofid,JPOP_vid,myid,nprocs,ipt)

		case ( "JPOS")
			i1 = lbound(JPOS,1) ; i2 = ubound(JPOS,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",JPOS,nc_ofid,JPOS_vid,myid,nprocs,ipt)

		case ( "CTEMP")
			i1 = lbound(CTEMP,1) ; i2 = ubound(CTEMP,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",CTEMP,nc_ofid,CTEMP_vid,myid,nprocs,ipt)

		case ( "SODTM1S")
			allocate(ftemp(1:MTLOC)) ; ftemp = SODTM1S
			i1 = lbound(ftemp,1) ; i2 = ubound(ftemp,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",ftemp,nc_ofid,SODTM1S_vid,myid,nprocs,ipt)
			deallocate(ftemp)

		case ( "JNH4TM1S")
			i1 = lbound(JNH4TM1S,1) ; i2 = ubound(JNH4TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",JNH4TM1S,nc_ofid,JNH4TM1S_vid,myid,nprocs,ipt)

		case ( "JNO3TM1S")
			i1 = lbound(JNO3TM1S,1) ; i2 = ubound(JNO3TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",JNO3TM1S,nc_ofid,JNO3TM1S_vid,myid,nprocs,ipt)

		case ( "BENDEN")
			i1 = lbound(BENDEN,1) ; i2 = ubound(BENDEN,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",BENDEN,nc_ofid,BENDEN_vid,myid,nprocs,ipt)

		case ( "JCH4TM1S")
			i1 = lbound(JCH4TM1S,1) ; i2 = ubound(JCH4TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",JCH4TM1S,nc_ofid,JCH4TM1S_vid,myid,nprocs,ipt)

		case ( "JCH4GTM1S")
			i1 = lbound(JCH4GTM1S,1) ; i2 = ubound(JCH4GTM1S,1)
			!call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",JCH4GTM1S,nc_ofid,JCH4GTM1S_vid,myid,nprocs,ipt)
      call putvar(i1,i2,MLOC,mgl,1,1,"n",JCH4GTM1S,nc_ofid,JCH4GTM1S_vid,myid,nprocs,ipt)  

		case ( "JHSTM1S")
			i1 = lbound(JHSTM1S,1) ; i2 = ubound(JHSTM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",JHSTM1S,nc_ofid,JHSTM1S_vid,myid,nprocs,ipt)

		case ( "JPO4TM1S")
			i1 = lbound(JPO4TM1S,1) ; i2 = ubound(JPO4TM1S,1)
			!call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",JPO4TM1S,nc_ofid,JPO4TM1S_vid,myid,nprocs,ipt)
      call putvar(i1,i2,MLOC,mgl,1,1,"n",JPO4TM1S,nc_ofid,JPO4TM1S_vid,myid,nprocs,ipt)  

		case ( "JSITM1S")
			i1 = lbound(JSITM1S,1) ; i2 = ubound(JSITM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",JSITM1S,nc_ofid,JSITM1S_vid,myid,nprocs,ipt)

		case ( "NH41TM1S")
			i1 = lbound(NH41TM1S,1) ; i2 = ubound(NH41TM1S,1)
			!call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",NH41TM1S,nc_ofid,NH41TM1S_vid,myid,nprocs,ipt)
      call putvar(i1,i2,MLOC,mgl,1,1,"n",NH41TM1S,nc_ofid,NH41TM1S_vid,myid,nprocs,ipt) 

		case ( "NH42TM1S")
			i1 = lbound(NH42TM1S,1) ; i2 = ubound(NH42TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",NH42TM1S,nc_ofid,NH42TM1S_vid,myid,nprocs,ipt)

		case ( "NO31TM1S")
			i1 = lbound(NO31TM1S,1) ; i2 = ubound(NO31TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",NO31TM1S,nc_ofid,NO31TM1S_vid,myid,nprocs,ipt)

		case ( "NO32TM1S")
			i1 = lbound(NO32TM1S,1) ; i2 = ubound(NO32TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",NO32TM1S,nc_ofid,NO32TM1S_vid,myid,nprocs,ipt)

		case ( "PO41TM1S")
			i1 = lbound(PO41TM1S,1) ; i2 = ubound(PO41TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",PO41TM1S,nc_ofid,PO41TM1S_vid,myid,nprocs,ipt)

		case ( "PO42TM1S")
			i1 = lbound(PO42TM1S,1) ; i2 = ubound(PO42TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",PO42TM1S,nc_ofid,PO42TM1S_vid,myid,nprocs,ipt)

		case ( "SI1TM1S")
			i1 = lbound(SI1TM1S,1) ; i2 = ubound(SI1TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",SI1TM1S,nc_ofid,SI1TM1S_vid,myid,nprocs,ipt)

		case ( "SI2TM1S")
			i1 = lbound(SI2TM1S,1) ; i2 = ubound(SI2TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",SI2TM1S,nc_ofid,SI2TM1S_vid,myid,nprocs,ipt)

		case ( "CH41TM1S")
			i1 = lbound(CH41TM1S,1) ; i2 = ubound(CH41TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",CH41TM1S,nc_ofid,CH41TM1S_vid,myid,nprocs,ipt)

		case ( "CH42TM1S")
			i1 = lbound(CH42TM1S,1) ; i2 = ubound(CH42TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",CH42TM1S,nc_ofid,CH42TM1S_vid,myid,nprocs,ipt)

		case ( "HS1TM1S")
			i1 = lbound(HS1TM1S,1) ; i2 = ubound(HS1TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",HS1TM1S,nc_ofid,HS1TM1S_vid,myid,nprocs,ipt)

		case ( "HS2TM1S")
			i1 = lbound(HS2TM1S,1) ; i2 = ubound(HS2TM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",HS2TM1S,nc_ofid,HS2TM1S_vid,myid,nprocs,ipt)

		case ( "CPOP")
			i1 = lbound(CPOP,1) ; i2 = ubound(CPOP,1)
			call putvar(i1,i2,MLOC,mgl,3,3,"n",CPOP,nc_ofid,CPOP_vid,myid,nprocs,ipt)
		case ( "CPON")
			i1 = lbound(CPON,1) ; i2 = ubound(CPON,1)
			call putvar(i1,i2,MLOC,mgl,3,3,"n",CPON,nc_ofid,CPON_vid,myid,nprocs,ipt)
		case ( "CPOC")
			i1 = lbound(CPOC,1) ; i2 = ubound(CPOC,1)
			call putvar(i1,i2,MLOC,mgl,3,3,"n",CPOC,nc_ofid,CPOC_vid,myid,nprocs,ipt)
		case ( "CPOS")
			i1 = lbound(CPOS,1) ; i2 = ubound(CPOS,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",CPOS,nc_ofid,CPOS_vid,myid,nprocs,ipt)
		case ( "HSED1")
			i1 = lbound(HSED1,1) ; i2 = ubound(HSED1,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",HSED1,nc_ofid,HSED1_vid,myid,nprocs,ipt)
		case ( "BENSTRTM1S")
			i1 = lbound(BENSTRTM1S,1) ; i2 = ubound(BENSTRTM1S,1)
			call putvar(i1,i2,MLOC,mgl,1,1,"n",BENSTRTM1S,nc_ofid,BENSTRTM1S_vid,myid,nprocs,ipt)

		case default
			!do nothing
		end select
	 ENDIF
     IF (ZOO_CALC) THEN 
		select case(trim(cdf_vdp(i)))
		case ( "SZ")
			i1 = lbound(SZ,1) ; i2 = ubound(SZ,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",SZ,nc_ofid,SZ_vid,myid,nprocs,ipt)
		case ( "LZ")
			i1 = lbound(LZ,1) ; i2 = ubound(LZ,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",LZ,nc_ofid,LZ_vid,myid,nprocs,ipt)
		case default
			!do nothing
		end select
	 ENDIF

	 !IF (SOLIDS_CALC) THEN 
		select case(trim(cdf_vdp(i)))
		case ( "KESS")
			i1 = lbound(KESS,1) ; i2 = ubound(KESS,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",KESS,nc_ofid,KESS_vid,myid,nprocs,ipt)
		case ( "TURB")
			i1 = lbound(TURB,1) ; i2 = ubound(TURB,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",TURB,nc_ofid,TURB_vid,myid,nprocs,ipt)
		case ( "SSI")
			i1 = lbound(SSI,1) ; i2 = ubound(SSI,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",SSI,nc_ofid,SSI_vid,myid,nprocs,ipt)
		case ( "TSS")
			i1 = lbound(TSS,1) ; i2 = ubound(TSS,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",TSS,nc_ofid,TSS_vid,myid,nprocs,ipt)
		case ( "PARAD")
			i1 = lbound(PARAD,1) ; i2 = ubound(PARAD,1)
			call putvar(i1,i2,MLOC,mgl,kb-1,kb-1,"n",PARAD,nc_ofid,PARAD_vid,myid,nprocs,ipt)
		case default
			!do nothing
		end select
	 !ENDIF

	 IF (SAV_CALC) THEN 
		select case(trim(cdf_vdp(i)))
		case ( "LEAF")
			i1 = lbound(LEAF,1) ; i2 = ubound(LEAF,1)
			call putvar(i1,i2,MLOC,mgl,nsavm,nsavm,"n",LEAF,nc_ofid,LEAF_vid,myid,nprocs,ipt)
		case ( "STEM")
			i1 = lbound(STEM,1) ; i2 = ubound(STEM,1)
			call putvar(i1,i2,MLOC,mgl,nsavm,nsavm,"n",STEM,nc_ofid,STEM_vid,myid,nprocs,ipt)
		case ( "ROOT")
			i1 = lbound(ROOT,1) ; i2 = ubound(ROOT,1)
			call putvar(i1,i2,MLOC,mgl,nsavm,nsavm,"n",ROOT,nc_ofid,ROOT_vid,myid,nprocs,ipt)
		case ( "TUBER")
			i1 = lbound(TUBER,1) ; i2 = ubound(TUBER,1)
			call putvar(i1,i2,MLOC,mgl,nsavm,nsavm,"n",TUBER,nc_ofid,TUBER_vid,myid,nprocs,ipt)
		case ( "EP")
			i1 = lbound(EP,1) ; i2 = ubound(EP,1)
			call putvar(i1,i2,MLOC,mgl,nsavm,nsavm,"n",EP,nc_ofid,EP_vid,myid,nprocs,ipt)
		case default
			!do nothing
		end select
	 ENDIF
   ENDDO

!==============================================================================|
!  CONSTANT OUTPUTS                                                            |
!==============================================================================|


!==============================================================================|
!  CLOSE THE FILE                                                              |
!==============================================================================|

   if(msr) ierr = nf90_close(nc_ofid)

   return
   end subroutine out_netcdf_data

#endif /* end of OUTPUT_NETCDF */

End Module MOD_NCD
