!mod_hydrovars.F
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
!  VARIABLES                                                                   |
!==============================================================================|
!
!
Module MOD_HYDROVARS !
  !
      Use MOD_PREC, Only: SP
      Use MOD_LIMS, Only: NSTATIONMAX, KB, KBM1, NTLOC, MTLOC
      Use MOD_TYPES, Only: BC
  !
      Use MOD_SIZES, Only: MGL, NGL, NOBTY, NCP !!

      Implicit None
      Save
  !
  !--Constants-------------------------------------------------------------------!
      Real (SP), Parameter :: GRAV = 9.81
      Real (SP), Parameter :: PI = 3.141592653
      Real (SP), Parameter :: PI2 = 6.283185307
      Real (SP), Parameter :: ZERO = 0.0
      Real (SP), Parameter :: ONE_THIRD = 1.0 / 3.0
  !
  !--------------------------Temporary Array------------------------------------------!
  !
      Integer, Allocatable :: NVG (:, :)
  !
  !--------------------------Global Grid Variables------------------------------------!
  !
      Real (SP), Allocatable :: XG (:)!!GLOBAL X-COORD AT NODE
      Real (SP), Allocatable :: YG (:)!!GLOBAL X-COORD AT NODE
      Real (SP), Allocatable :: HG (:)!!GLOBAL DEPTH AT NODE
      Real (SP), Allocatable :: XCG (:)!!GLOBAL X-COORD AT FACE CENTER
      Real (SP), Allocatable :: YCG (:)!!GLOBAL X-COORD AT FACE CENTER
  !
  !--------------------------Grid Metrics---------------------------------------------!
  !
      Real (SP) :: VXMIN, VYMIN, VXMAX, VYMAX
      Real (SP), Allocatable :: XC (:)!!X-COORD AT FACE CENTER
      Real (SP), Allocatable :: YC (:)!!Y-COORD AT FACE CENTER
      Real (SP), Allocatable :: VX (:)!!X-COORD AT GRID POINT
      Real (SP), Allocatable :: VY (:)!!Y-COORD AT GRID POINT
      Real (SP), Allocatable :: ART (:)!!AREA OF ELEMENT
      Real (SP), Allocatable :: ART1 (:)!!AREA OF NODE-BASE CONTROl VOLUME
      Real (SP), Allocatable :: ART2 (:)!!AREA OF ELEMENTS AROUND NODE
!
  !----------------1-d arrays for the sigma coordinate -------------------------------!
  !
      Real (SP), Allocatable :: Z (:)		!!SIGMA COORDINATE VALUE
      Real (SP), Allocatable :: ZZ (:)		!!INTRA LEVEL SIGMA VALUE
      Real (SP), Allocatable :: DZ (:)		!!DELTA-SIGMA VALUE
      Real (SP), Allocatable :: DZZ (:)		!!DELTA OF INTRA LEVEL SIGMA

!#if defined (NEWSIGMA)
	  Real (SP), Allocatable :: Z2DG (:,:)	!!SIGMA COORDINATE VALUE global, spatially varying
	  Real (SP), Allocatable :: Z2D (:,:)	!!SIGMA COORDINATE VALUE, spatially varying
      Real (SP), Allocatable :: ZZ2D (:,:)	!!INTRA LEVEL SIGMA VALUE, spatially varying
      Real (SP), Allocatable :: DZ2D (:,:)	!!DELTA-SIGMA VALUE, spatially varying
      Real (SP), Allocatable :: DZZ2D (:,:)	!!DELTA OF INTRA LEVEL SIGMA, spatially varying
!#endif

  !
  !---------------2-d flow variable arrays at elements-------------------------------!
  !
      Real (SP), Allocatable :: H1 (:)!!BATHYMETRIC DEPTH
  !
  !---------------2-d flow variable arrays at nodes----------------------------------!
  !
      Real (SP), Allocatable :: H (:)!!BATHYMETRIC DEPTH
      Real (SP), Allocatable :: D (:)!!CURRENT DEPTH
      Real (SP), Allocatable :: DT (:)!!DEPTH AT PREVIOUS TIME STEP
      Real (SP), Allocatable :: DT1 (:)!!DEPTH AT PREVIOUS TIME STEP
      Real (SP), Allocatable :: EL (:)!!CURRENT SURFACE ELEVATION
      Real (SP), Allocatable :: ET (:)!!SURFACE ELEVATION AT PREVIOUS TIME STEP
      Real (SP), Allocatable :: DTFA (:)!!ADJUSTED DEPTH FOR MASS CONSERVATION
  !
  !---------------- internal mode   arrays-(element based)----------------------------!
  !
      Real (SP), Allocatable :: UU (:, :)!X-VELOCITY
      Real (SP), Allocatable :: VV (:, :)!Y-VELOCITY
      Real (SP), Allocatable :: UUT (:, :)!X-VELOCITY FROM PREVIOUS TIMESTEP
      Real (SP), Allocatable :: VVT (:, :)!Y-VELOCITY FROM PREVIOUS TIMESTEP
  !         REAL(SP), ALLOCATABLE :: WWT(:,:)        !Z-VELOCITY FROM PREVIOUS TIMESTEP !Never used
      Real (SP), Allocatable :: WTST (:, :)!Vertical velocity in sigma from PREVIOUS TIMESTEP
      Real (SP), Allocatable :: UARD_OBCNT (:)
      Real (SP), Allocatable :: XFLUX_OBCT (:, :)
      Real (SP), Allocatable :: DTFAT (:)!tykim
  !         REAL(SP), ALLOCATABLE :: TT_T(:,:)       !never used
  !         REAL(SP), ALLOCATABLE :: SALTT(:,:)      !never used
  !
  !-----------------------3d variable arrays-(node based)-----------------------------!
  !
      Real (SP), Allocatable :: WTS (:, :)!!VERTICAL VELOCITY IN SIGMA SYSTEM
      Real (SP), Allocatable :: UARD_OBCN (:)
      Real (SP), Allocatable :: XFLUX_OBC (:, :)
  !         REAL(SP), ALLOCATABLE :: WTTS(:,:)       !!VERTICAL VELOCITY IN SIGMA SYSTEM     !never used
      Real (SP), Allocatable :: KH (:, :)!!TURBULENT DIFFUSIVITY
  !
      Real (SP), Allocatable :: XFLUX_OBC_WQM (:, :, :)!
  !
  !--------------------hydrodynamics-------------------------------------------------
  !
      Real (SP), Allocatable :: VISCOFH (:, :)
  !
      Real (SP), Allocatable :: UNC1 (:, :)
      Real (SP), Allocatable :: VNC1 (:, :)
      Real (SP), Allocatable :: WNC1 (:, :)
      Real (SP), Allocatable :: WTSNC1 (:, :)
      Real (SP), Allocatable :: UARD_OBCNNC1 (:)
      Real (SP), Allocatable :: XFLUX_OBCNC1 (:, :)
      Real (SP), Allocatable :: DTFANC1 (:)
      Real (SP), Allocatable :: KHNC1 (:, :)
      Real (SP), Allocatable :: TNC1 (:, :)
      Real (SP), Allocatable :: SNC1 (:, :)
      Real (SP), Allocatable :: ELNC1 (:)
  !
      Real (SP), Allocatable :: UNC2 (:, :)
      Real (SP), Allocatable :: VNC2 (:, :)
      Real (SP), Allocatable :: WNC2 (:, :)
      Real (SP), Allocatable :: WTSNC2 (:, :)
      Real (SP), Allocatable :: UARD_OBCNNC2 (:)
      Real (SP), Allocatable :: XFLUX_OBCNC2 (:, :)
      Real (SP), Allocatable :: DTFANC2 (:)
      Real (SP), Allocatable :: KHNC2 (:, :)
      Real (SP), Allocatable :: TNC2 (:, :)
      Real (SP), Allocatable :: SNC2 (:, :)
      Real (SP), Allocatable :: ELNC2 (:)
!
      Real (SP), Allocatable :: RHO (:, :)!
  !
      Integer (4) :: num_hyd_ints !number of records in each hydrodynamics netcdf file
  !
  !
  !
      Type (BC) :: TIME_MAP
      Real (SP) :: THOUR1 !!SIMULATION TIME AT END OF CURRENT EXTERNAL STEP (IEXT) IN HOURS
      Real (SP) :: THOUR
  !
  !Added following for hydrodynamics input
      Character (Len=1024) :: NCFILE_DIR, NCFILE_PREFIX, NCFILE_SUFFIX, &
     & NCFILE_NUMBER
      Character (Len=1024) :: FORMAT_STR
      Character (Len=1024) :: hydro_dir, hydro_prefix, hydro_suffix ! suffix of filename, e.g. '.nc'
      Integer (4) :: hydro_filenumwidth, hydro_filenumstart, hydro_Nrec ! number of records in each of hydrodynamics file
      Real (SP) :: hydro_dlt ! time step in hydrodynamics file (in seconds), e.g. 100 for 100sec
  !
  !The following are for outputs should really go into mod_wqm.F or mod_output.F
  !
      Real (SP) :: t_his_start, t_his_end, t_his_dlt !starting time, ending time, and interval of history outputs (days)
  !
      Integer (4) :: Nstation, NstationNum_GL (NSTATIONMAX)!maximum number of station is NstationMax!
      Real (SP) :: t_stn_start, t_stn_end, t_stn_dlt !starting time, ending time, and interval of station outputs (days)
  !
      Character (Len=1024) :: STNFN !file name for station output
      Character (Len=1024) :: HISFN !file name for history output

	  Character (Len=1024) :: HIS_OUTDIR !output dir of the history output
	  Character (Len=1024) :: STN_OUTDIR !output dir of the station output
      Character (Len=1024) :: HISFN_PREFIX !prefix of history output file
      Character (Len=1024) :: HISFN_EXT !extention name of history output file
      Character (Len=1024) :: HISFN_FINAL
      Logical :: HISFN_SPLIT_BYLEVEL = .False. !Ture or False for splitting history output in files level by level (default is .FALSE.)
  !
      Namelist / hydro_netcdf / hydro_dir, hydro_prefix, hydro_suffix, &
     & hydro_filenumwidth, hydro_filenumstart, hydro_Nrec, hydro_dlt
      Namelist / wqm_history /HIS_OUTDIR, HISFN, t_his_start, t_his_end, &
     & t_his_dlt, HISFN_SPLIT_BYLEVEL
      Namelist / wqm_stations /STN_OUTDIR, STNFN, Nstation, NstationNum_GL, &
     & t_stn_start, t_stn_end, t_stn_dlt

      Integer (4) :: IFNC !file number index for hydrodynamics netcdf files, set to hydro_filenumstart initially for cold start, set otherwise
  !for restart
      Integer (4) :: NTRECNC !time record index for a particular hydrodynamics netcdf file, reset to 1 upon opening new file.
      Integer (4) :: NTHYDRO !overall time record index for all netcdf files, increment by 1 each time a hydrodynamics record is read

	  CHARACTER(Len=3):: SIGVARC !' ON' or 'OFF'  !where the model uses variable sigma

	  Logical :: SIGVAR		!SIGVASR=.TRUE. if SIGVAR=' ON'

Contains
  !
  !Subroutines:
  !	SUBROUTINE HYDRO_GEOM_ALLOC()
  !	SUBROUTINE HYDRO_GEOM_DEALLOC()
  ! 	SUBROUTINE HYDRO_ALLOC()
  ! 	SUBROUTINE HYDRO_DEALLOC()
  !
  !Subroutine to allocate geometry related variables from hydrodyanmics (netcdf)
      Subroutine HYDRO_GEOM_ALLOC
         Allocate (XG(0:MGL))
         XG = 0.0_SP
         Allocate (YG(0:MGL))
         YG = 0.0_SP
         Allocate (NVG(0:NGL, 4))
         NVG = 0
         Allocate (HG(0:MGL))
         HG = 0.0_SP
         Allocate (Z(KB))
         Z = 0.0_SP !!SIGMA COORDINATE VALUE
         Allocate (ZZ(KB))
         ZZ = 0.0_SP !!INTRA LEVEL SIGMA VALUE
         Allocate (DZ(KB))
         DZ = 0.0_SP !!DELTA-SIGMA VALUE
         Allocate (DZZ(KB))
         DZZ = 0.0_SP !!DELTA OF INTRA LEVEL SIGMA

!#if defined (NEWSIGMA)

		 Allocate (Z2DG(0:MGL,KB))
		 Z2DG=0.0_SP
		 Allocate (Z2D(0:MTLOC,KB))
         Z2D = 0.0_SP !!SIGMA COORDINATE VALUE
         Allocate (ZZ2D(0:MTLOC,KB))
         ZZ2D = 0.0_SP !!INTRA LEVEL SIGMA VALUE
         Allocate (DZ2D(0:MTLOC,KB))
         DZ2D = 0.0_SP !!DELTA-SIGMA VALUE
         Allocate (DZZ2D(0:MTLOC,KB))
         DZZ2D = 0.0_SP !!DELTA OF INTRA LEVEL SIGMA

!#endif
    !
         Allocate (XCG(0:NGL))
         XCG = 0.0_SP
         Allocate (YCG(0:NGL))
         YCG = 0.0_SP
    !
      End Subroutine HYDRO_GEOM_ALLOC
  !
  !Subroutine to allocate geometry related variables from hydrodyanmics (netcdf)
      Subroutine HYDRO_GEOM_DEALLOC
    !
         If (ALLOCATED(XG)) 	DEALLOCATE (XG)
         If (ALLOCATED(YG)) 	DEALLOCATE (YG)
         If (ALLOCATED(NVG)) 	DEALLOCATE (NVG)
         If (ALLOCATED(HG)) 	DEALLOCATE (HG)
         If (ALLOCATED(Z)) 		DEALLOCATE (Z)
         If (ALLOCATED(ZZ)) 	DEALLOCATE (ZZ)
         If (ALLOCATED(DZ)) 	DEALLOCATE (DZ)
         If (ALLOCATED(DZZ)) 	DEALLOCATE (DZZ)

!#if defined (NEWSIGMA)
		 IF(Allocated (Z2DG))	DEALLOCATE(Z2DG)
		 IF(Allocated (Z2D))	DEALLOCATE(Z2D)
         IF(Allocated (ZZ2D))	DEALLOCATE(ZZ2D)
         IF(Allocated (DZ2D))	DEALLOCATE(DZ2D)
         IF(Allocated (DZZ2D))	DEALLOCATE(DZZ2D)
!#endif
         If (ALLOCATED(XCG)) 	DEALLOCATE (XCG)
         If (ALLOCATED(YCG)) 	DEALLOCATE (YCG)
    !
      End Subroutine HYDRO_GEOM_DEALLOC
  !
      Subroutine HYDRO_ALLOC
    !
    !local geometry
         Allocate (VX(0:MTLOC))
         VX = 0.0_SP !!X-COORD AT GRID POINT
         Allocate (VY(0:MTLOC))
         VY = 0.0_SP !!X-COORD AT GRID POINT
         Allocate (XC(0:NTLOC))
         XC = 0.0_SP !!X-COORD AT FACE CENTER
         Allocate (YC(0:NTLOC))
         YC = 0.0_SP !!X-COORD AT FACE CENTER
         Allocate (H(0:MTLOC))
         H = 0.0_SP !!BATHYMETRIC DEPTH
         Allocate (H1(0:NTLOC))
         H1 = 0.0_SP !!BATHYMETRIC DEPTH
    !
    !Moved from cell_area.F to here
         Allocate (ART(0:NTLOC))
         ART = 0.0_SP !!AREA OF ELEMENT
         Allocate (ART1(0:MTLOC))
         ART1 = 0.0_SP !!AREA OF NODE-BASE CONTROl VOLUME
         Allocate (ART2(0:MTLOC))
         ART2 = 0.0_SP !!AREA OF ELEMENTSAROUND NODE
    !
    !
    !
         Allocate (UNC1(0:NTLOC, KB))
         UNC1 = 0.0
         Allocate (VNC1(0:NTLOC, KB))
         VNC1 = 0.0
         Allocate (WNC1(0:MTLOC, KB))
         WNC1 = 0.0
         Allocate (WTSNC1(0:MTLOC, KB))
         WTSNC1 = 0.0
         Allocate (UARD_OBCNNC1(0:NOBTY+1))
         UARD_OBCNNC1 = 0.0
         Allocate (XFLUX_OBCNC1(0:NOBTY, KBM1))
         XFLUX_OBCNC1 = 0.0
         Allocate (DTFANC1(0:MTLOC))
         DTFANC1 = 0.0
         Allocate (KHNC1(0:MTLOC, KB))
         KHNC1 = 0.0
         Allocate (TNC1(0:MTLOC, KBM1))
         TNC1 = 0.0
         Allocate (SNC1(0:MTLOC, KBM1))
         SNC1 = 0.0
         Allocate (ELNC1(0:MTLOC))
         ELNC1 = 0.0
    !
         Allocate (UNC2(0:NTLOC, KB))
         UNC2 = 0.0
         Allocate (VNC2(0:NTLOC, KB))
         VNC2 = 0.0
         Allocate (WNC2(0:MTLOC, KB))
         WNC2 = 0.0
         Allocate (WTSNC2(0:MTLOC, KB))
         WTSNC2 = 0.0
         Allocate (UARD_OBCNNC2(0:NOBTY+1))
         UARD_OBCNNC2 = 0.0
         Allocate (XFLUX_OBCNC2(0:NOBTY, KBM1))
         XFLUX_OBCNC2 = 0.0
         Allocate (DTFANC2(0:MTLOC))
         DTFANC2 = 0.0
         Allocate (KHNC2(0:MTLOC, KB))
         KHNC2 = 0.0
         Allocate (TNC2(0:MTLOC, KBM1))
         TNC2 = 0.0
         Allocate (SNC2(0:MTLOC, KBM1))
         SNC2 = 0.0
         Allocate (ELNC2(0:MTLOC))
         ELNC2 = 0.0
    !
         Allocate (RHO(0:MTLOC, KBM1))!
         RHO = 1000.0 !
    !
         Allocate (UU(0:NTLOC, KB))
         UU = 0.0
         Allocate (VV(0:NTLOC, KB))
         VV = 0.0
         Allocate (UUT(0:NTLOC, KB))
         UUT = 0.0
         Allocate (VVT(0:NTLOC, KB))
         VVT = 0.0
         Allocate (WTST(0:MTLOC, KB))
         WTST = 0.0
         Allocate (UARD_OBCNT(0:NOBTY+1))
         UARD_OBCNT = 0.0
         Allocate (XFLUX_OBCT(0:NOBTY, KBM1))
         XFLUX_OBCT = 0.0
         Allocate (DTFAT(0:MTLOC))
         DTFAT = 0.0
    !
         Allocate (WTS(0:MTLOC, KB))
         WTS = 0.0
         Allocate (UARD_OBCN(0:NOBTY+1))
         UARD_OBCN = 0.0
         Allocate (XFLUX_OBC(0:NOBTY, KBM1))
         XFLUX_OBC = 0.0
         Allocate (KH(0:MTLOC, KB))
         KH = 0.0
    !
         Allocate (XFLUX_OBC_WQM(0:NOBTY, KBM1, NCP))!Calculation occurs in adv_wqm.F and is then used in bcond_wqm.F
         XFLUX_OBC_WQM = 0.0 !                       !as it was saved and used in FVCOM
    !
         Allocate (EL(0:MTLOC))
         EL = 0.0
         Allocate (ET(0:MTLOC))
         ET = 0.0
         Allocate (D(0:MTLOC))
         D = 0.0
         Allocate (DT(0:MTLOC))
         DT = 0.0
         Allocate (DT1(0:NTLOC))
         DT1 = 0.0
         Allocate (DTFA(0:MTLOC))
         DTFA = 0.0
    !
    !Allocate (VISCOFH(0:NTLOC, KB))  !
         Allocate (VISCOFH(0:MTLOC, KBM1))
         VISCOFH = 0.0
    !
      End Subroutine HYDRO_ALLOC
  !
      Subroutine HYDRO_DEALLOC
    !
    !local geometry
         If (ALLOCATED(VX)) DEALLOCATE (VX)!1
         If (ALLOCATED(VY)) DEALLOCATE (VY)!2
         If (ALLOCATED(XC)) DEALLOCATE (XC)!3
         If (ALLOCATED(YC)) DEALLOCATE (YC)!4
         If (ALLOCATED(H)) DEALLOCATE (H)!5
         If (ALLOCATED(H1)) DEALLOCATE (H1)!6
    !
         If (ALLOCATED(ART)) DEALLOCATE (ART)!7
         If (ALLOCATED(ART1)) DEALLOCATE (ART1)!8
         If (ALLOCATED(ART2)) DEALLOCATE (ART2)!9
    !
         If (ALLOCATED(UNC1)) DEALLOCATE (UNC1)!10
         If (ALLOCATED(VNC1)) DEALLOCATE (VNC1)!11
         If (ALLOCATED(WNC1)) DEALLOCATE (WNC1)!12
         If (ALLOCATED(WTSNC1)) DEALLOCATE (WTSNC1)!13
         If (ALLOCATED(UARD_OBCNNC1)) DEALLOCATE (UARD_OBCNNC1)!14
         If (ALLOCATED(XFLUX_OBCNC1)) DEALLOCATE (XFLUX_OBCNC1)!15
         If (ALLOCATED(DTFANC1)) DEALLOCATE (DTFANC1)!16
         If (ALLOCATED(KHNC1)) DEALLOCATE (KHNC1)!17
         If (ALLOCATED(TNC1)) DEALLOCATE (TNC1)!18
         If (ALLOCATED(SNC1)) DEALLOCATE (SNC1)!19
         If (ALLOCATED(ELNC1)) DEALLOCATE (ELNC1)!20
    !
         If (ALLOCATED(UNC2)) DEALLOCATE (UNC2)!21
         If (ALLOCATED(VNC2)) DEALLOCATE (VNC2)!22
         If (ALLOCATED(WNC2)) DEALLOCATE (WNC2)!23
         If (ALLOCATED(WTSNC2)) DEALLOCATE (WTSNC2)!24
         If (ALLOCATED(UARD_OBCNNC2)) DEALLOCATE (UARD_OBCNNC2)!25
         If (ALLOCATED(XFLUX_OBCNC2)) DEALLOCATE (XFLUX_OBCNC2)!26
         If (ALLOCATED(DTFANC2)) DEALLOCATE (DTFANC2)!27
         If (ALLOCATED(KHNC2)) DEALLOCATE (KHNC2)!28
         If (ALLOCATED(TNC2)) DEALLOCATE (TNC2)!29
         If (ALLOCATED(SNC2)) DEALLOCATE (SNC2)!30
         If (ALLOCATED(ELNC2)) DEALLOCATE (ELNC2)!31
    !
         If (ALLOCATED(RHO)) DEALLOCATE (RHO)!
    !
         If (ALLOCATED(UU)) DEALLOCATE (UU)!32
         If (ALLOCATED(VV)) DEALLOCATE (VV)!33
    !
         If (ALLOCATED(UUT)) DEALLOCATE (UUT)!34
         If (ALLOCATED(VVT)) DEALLOCATE (VVT)!35
         If (ALLOCATED(WTST)) DEALLOCATE (WTST)!36
         If (ALLOCATED(UARD_OBCNT)) DEALLOCATE (UARD_OBCNT)!37
         If (ALLOCATED(XFLUX_OBCT)) DEALLOCATE (XFLUX_OBCT)!38
         If (ALLOCATED(DTFAT)) DEALLOCATE (DTFAT)!39
    !
         If (ALLOCATED(WTS)) DEALLOCATE (WTS)!40
         If (ALLOCATED(UARD_OBCN)) DEALLOCATE (UARD_OBCN)!41
         If (ALLOCATED(XFLUX_OBC)) DEALLOCATE (XFLUX_OBC)!42
         If (ALLOCATED(KH)) DEALLOCATE (KH)!43
    !
         If (ALLOCATED(XFLUX_OBC_WQM)) DEALLOCATE (XFLUX_OBC_WQM)!
    !
         If (ALLOCATED(EL)) DEALLOCATE (EL)!44
         If (ALLOCATED(ET)) DEALLOCATE (ET)!45
    !
         If (ALLOCATED(D)) DEALLOCATE (D)!46
         If (ALLOCATED(DT)) DEALLOCATE (DT)!47
         If (ALLOCATED(DT1)) DEALLOCATE (DT1)!48
         If (ALLOCATED(DTFA)) DEALLOCATE (DTFA)!49
    !
         If (ALLOCATED(VISCOFH)) DEALLOCATE (VISCOFH)!50
    !
      End Subroutine HYDRO_DEALLOC
  !
  !
End Module MOD_HYDROVARS
