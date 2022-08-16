!wqm_inputs.F
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
!Subroutine INPUTS()
!
Subroutine INPUTS ()
  !
      Use MOD_FILEINFO, Only: DIA, CBC, S1, S2, S3, BFI, BAI, MET, BFO, &
     & KEI, ATM, STL, AGR, SVI, CON, SNP, PLT, APL, TFL, OPL, SFI, SFO, &
     & ICI, MRL, MBL, RSI, UNIT_LINKAGE, CNAME, CONFN, KFL, DFO
!
      Use MOD_LIMS, Only: MTLOC, NTLOC, NPROCS, MYID, NSTATIONMAX, &
     & MLOC, KB, KBM1
      Use MOD_SIZES, Only: MGL, NGL, NOBTY, NCP, NSAVM
      Use MOD_PREC, Only: SP
  !
      Use MOD_TGE, Only: ISBCE, ISONB, NV, TRIANGLE_GRID_EDGE, &
     & TGE_ALLOC
  !
      Use MOD_BCMAP, Only: BCMAP
  !
      Use MOD_HYDROVARS, Only: ART, ART1, DZ,DZ2D, H, D, DT, DT1, EL, ET, &
     & DTFA, UU, VV, UUT, VVT, WTST, UARD_OBCNT, XFLUX_OBCT, DTFAT, &
     & WTS, UARD_OBCN, XFLUX_OBC, KH, VISCOFH, UNC1, VNC1, WNC1, &
     & WTSNC1, UARD_OBCNNC1, XFLUX_OBCNC1, DTFANC1, KHNC1, TNC1, SNC1, &
     & ELNC1, UNC2, VNC2, WNC2, WTSNC2, UARD_OBCNNC2, XFLUX_OBCNC2, &
     & DTFANC2, KHNC2, TNC2, SNC2, ELNC2, NCFILE_PREFIX, NCFILE_SUFFIX, &
     & NCFILE_NUMBER, FORMAT_STR, hydro_dir, hydro_prefix, &
     & hydro_suffix, hydro_filenumwidth, hydro_filenumstart, &
     & hydro_Nrec, t_his_dlt, Nstation, NstationNum_GL, t_stn_dlt, &
     & hydro_netcdf, wqm_history, wqm_stations, IFNC, NTRECNC, NTHYDRO, &
     & HYDRO_GEOM_ALLOC, HYDRO_ALLOC, RHO, HIS_OUTDIR,STN_OUTDIR,SIGVAR,SIGVARC ! 
  !
  !Took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: SERIAL, MSR, PAR, DTI
!
      Use MOD_CO2SYS, Only: CO2SYSCONST_ALLOC
      Use MOD_PH_PCO2, Only: CALC_PH_PCO2
!
# if defined (KILLTREES)
      Use MOD_WQMINIT, Only: HYDC
# endif
  Use MOD_WQM, Only: AHMDLT, FILGTH, JDAY, NAPL, NDLT, NSNP, TMEND, &
       & TMSTRT, NDIA, NKFL, NMBL, NOINT, NOPL, NPLT, NRSO, NTFL, REDS1C, &
       & REDS1N, REDS1P, REDS2C, TH, XYDF, ZDFBCK, ZDFMUL, ATMFN, KEIFN, &
       & METFN, NHYDF, NTVDF, REDCBC, REDCBN, REDCBP, REDS2N, REDS2P, &
       & REDS3C, REDS3N, REDS3P, S1FN, S2FN, S3FN, SVIFN, ZOOFN, ALOFN, &
       & ATMOS_LOADS, BAIFN, BAOFN, BENTHIC_FLUXES, BFIFN, BFOFN, &
       & BOUNDARY_CONC, CBCFN, DFIFN, DFOFN, SAV_CALC, SAV_LOADS, &
       & SETTLING, SOURCE_ONE, SOURCE_THR, SOURCE_TWO, SVOFN, ZFOFN, AC, &
       & ATMPTR, BAOPTR, BFIPTR, CBCPTR, HYDPTR, IJDAY, KEIPTR, &
       & KINETIC_FLUXES, LIGHT_EXTINCTION, METPTR, NAC, S1PTR, S2PTR, &
       & S3PTR, SAVPTR, SEDIMENT_CALC, BBMI, C1, C2, CTEMP, DFEEDI, &
       & AANOX, ANDC, AOCR, AONT, AREAR, BREAR, CREAR, KADPO4, KADSA, &
       & KHCOAG, KHNDN, KHNNT, KHOCOD, KHODOC, KHONT, KLDC, KTCOD, KTHDR, &
       & KTMNL, KTNT1, KTNT2, KTSUA, TMNT, TRCOD, TRHDR, TRMNL, TRSUA, &
       & KLDN, KLPC, KRDC, KRDN, KRPC, KLDP, KLPN, KLPP, KRDP, KRPN, &
       & KCOD, KDCALG, KLCALG, KRPP, KSUA, KDNALG, KDPALG, KLNALG, &
       & KLPALG, KRCOAG, NCB, NTM, S1LB, S1LN, S2LB, S2LN, S3LB, S3LN, &
       & WS1, WS2, WSL, WSR, WSS, FLOW, KHSO, KSDOC, KSNH4, KSNO3, KSO, &
       & KSPO4, KSSA, MTCNO3, SEDNO3, TRSDOC, TRSNH4, TRSNO3, TRSO, &
       & TRSPO4, TRSSA, WS3, WSU, AVERAGE_PLOTS, CONSERVE_MASS, &
       & DIAGNOSTICS, JDAYMBL, KFLDP, MASS_BALANCE, MBLDP, NIT, PIP_CALC, &
       & QUALITY_DIAG, SAV_PLOTS, SEDIMENT_DIAG, SOLIDS_CALC, &
       & STEP_BOUNDARY, TRANSPORT_FLUXES, XY_DIFFUSION, Z_DIFFUSION, AC1, &
       & ACCHL1, ACCHL2, ACCHL3, AFI1, AFI2, AFI3, AKE, ANL1, ANL2, APL1, &
       & APL2, ASL1, ASL2, DLT, ELTMS, ELTMSKFL, ELTMSPLT, ELTMSTFL, &
       & FNDLT, NHMR, NWQMR, NXATM, NXBFI, NXCBC, NXMET, NXS1, NXS2, &
       & NXS3, NXSAV, AASRAT, ABBM, ABENCH4A, ABENCH4G, ABENCOD, ABENDO, &
       & ABENDOC, ABENNH4, ABENNO3, ABENPO4, ABENSA, ABLITE, ACFIX, &
       & ACPIP, ACPOS, AFIB, AGPP, ANL3, ANLB, ANPP, ANPPB, APCFWS, APL3, &
       & APLB, APNFWS, APPFWS, APSFWS, ARESP, ASL3, ASSFWS, ABM1, ABM2, &
       & ACPOC, ACPON, ACPOP, ADFEED, AFLXCSF, AFLXNSF, AFLXPSF, AJNSF, &
       & AJPSF, AP1, AP2, APR1, APR2, ARPOCSF, ARPONSF, ARPOPSF, ASASF, &
       & ASFCFEC, ASFCPSF, ASFEED, ASFGCIN, ASODSF, ASSIPSF, ASSISASF, &
       & ASSISF, ASSISUSF, ASUSF, A_T, AALGDOC, AALGDON, AALGDOP, &
       & AALGNH4, AALGNO3, AALGPO4, AALGPOC, AALGPON, AALGPOP, AALGRES, &
       & AALGUP, ABM3, ADCOD, ADENIT, ADENNO3, ADO, ADOPR, ADORALG, &
       & AHDRPOC, AHDRPON, AHDRPOP, AMNLDOC, AMNLDON, AMNLDOP, ANFIX, &
       & ANT, AP3, APR3, APSD, ASAP, ABADO, ABADOC, ABANH4, ABANO3, &
       & ABAPO4, ABAPOC, ABAPON, ABAPOP, ABMB, ADDOC, ANC1, ANC2, ANC3, &
       & ANITRIF, APB, APRB, ASC1, ASC2, ASC3, CMASS, DTC, INFLOW, IWCMC, &
	   & LEAF_GL, ROOT_GL, STEM_GL, TUBER_GL, EP_GL, &
# if defined (KILLTREES)
       & HYDFN,NB,NSB,NQF,NL,KE,KEISS,KEVSS,KEDOC,&
# endif
       & IWCMN, IWCMP, IWCMS, JCS1MAX, JCS2MAX, JCS3MAX, NHQF, WQM_ALLOC, &
       & T, SALT, TDIC, TALK ,SEDTEMP_CALC ! 
  !   and also get TDIC and TALK to be used in initialization of CO2star
!
      Use MOD_OWQ, Only: INITKE, INTKE, NXKEI, LAT, LONG, GSTAR440, &
     & PHISTAR676, BNSTARNTU, BNSTARCHL, PSTARINT, PSTAR440, &
     & PSTARCOLOR, SCDOM, EZEROINT, TCHL2CHL, DOFFSET, DEPTHKE, NWAVEL, &
     & COLOR, TURB, ALAMB, GLAMB, PLAMB, PHILAMB, BLAMB, WAVEL, EZERO, &
     & OWQ_ALLOC, OWQ_READ
      Use MOD_SF, Only: SFEEDER, NSPECIES, HYPOXFX_SF, SFEED, SFEED_GL, &
     & SF_ALLOC, SF_READ, SF_INIT
  !
      Use MOD_DF, Only: DFEEDER, DFEEDER_OUTPUT, HYPOXFX_DF, &
     & DFEEDM1S_GL, DFEEDM1S, DF_READ, DF_INIT, DF_ALLOC
  !
      Use MOD_WQMINIT, Only: RSODP, DLTDP, SNPDP, TFLDP, PLTDP, APLDP, &
     & OPLDP, DIADP, F, S1LNMAX, S2LNMAX, S3LNMAX, NXSNP, NXPLT, NXAPL, &
     & NXTFL, NXKFL, NXTVD, NXOPL, NXMBL, NXDIA,  &
     & MXDLT, MINSTEP, SNPC, RSIC, RSOC, RSOD, BCC, S1C, S2C, &
     & S3C, MDC, PLTC, FLC, MBLC, BFC, VBC, QPLTC, XYDFC, XYDFU, ZDFC, &
     & ICOC, ATMC, SAVLC, SEDC,SEDTC, AUTOC, SPLTC, TFLC, DIAC, STLC, APLTC, &
     & KFLC, OPLC, BFOC, BAOC, DFOC, SAVMC, SAVPLTC, DFLC, SFLC, &
     & SFLOXC, DFLOXC, SLC, BNDTC, CONSC, ICIC, SPVARM, PRINTM, TITLE, &
     & OLDTITLE, MAPFN, GEOFN, ICIFN, RSIFN, AGRFN, STLFN, MRLFN, &
     & KFLFN, ICOFN, SNPFN, PLTFN, APLFN, DIAFN, TFLFN, RSOFN, OPLFN, &
     & MBLFN, SFIFN, SFOFN, RESTART_OUT, SNAPSHOTS, END_RUN, &
     & MODIFY_ICONC, VOLUME_BALANCE, QUICKEST, UPWIND, ICOND_OUT, &
     & UNI_ICON_IN, UNI_ICON_IN_SED_VAR, BIN_ICON_IN, RES_ICON_IN, &
     & AUTO_STEPPING, PLOTS, OXYGEN_PLOTS, ZOO_CALC, RESTART_IN, &
     & TEMPERATURE_CALC, ALGAE_CALC, CARBON_CALC, NITROGEN_CALC, &
     & PHOSPHORUS_CALC, COD_CALC, OXYGEN_CALC, SILICA_CALC, LEFT_FLOWB, &
     & RIGHT_FLOWB, IFLOWP, C1MIN, C1MAX, DOVDAYS, OINT, ACC, CIC, &
     & SFEEDI, DLTD, SNPD, SNPF, DLTVAL, DLTMAX, DLTFTN, PLTD, PLTF, &
     & APLTD, APLF, TFLD, TFLF, KFLD, KFLF, OPLD, OPLF, DIAD, DIAF, &
     & MBLD, MBLF, WQMINIT_ALLOC, CARBONATE_CALC
  !
      Use MOD_ALGAL, Only: ALG_READ
  !
      Use MOD_BA, Only: BBM, BBM_GL, BALGAE_CALC, BA_OUTPUT, BALC
  !
      Use MOD_SED, Only: BENTHIC_OUTPUT, STEADY_STATE_SED_IC, &
     & BENSTRTM1S, SODTM1S, NH41TM1S, NH4T2TM1S, NO31TM1S, NO3T2TM1S, &
     & HS1TM1S, HST2TM1S, SI1TM1S, SIT2TM1S, PO41TM1S, PO4T2TM1S, &
     & CH4T2TM1S, CH41TM1S, SO4T2TM1S, SODTM1S, PO41TM1S, NH41TM1S, &
     & NO31TM1S, HS1TM1S, SI1TM1S, CTEMP_GL, CPOP_GL, CPON_GL, CPOC_GL, &
     & CPOS_GL, PO4T2TM1S_GL, NH4T2TM1S_GL, NO3T2TM1S_GL, HST2TM1S_GL, &
     & CH4T2TM1S_GL, CH41TM1S_GL, SO4T2TM1S_GL, SIT2TM1S_GL, &
     & BENSTRTM1S_GL, SODTM1S_GL, PO41TM1S_GL, NH41TM1S_GL, &
     & NO31TM1S_GL, HS1TM1S_GL, SI1TM1S_GL, CPOC, CPON, CPOP, CPOS, &
     & CPOPI, CPONI, CPOCI, CPOSI, PO41TI, PO4T2I, NO31TI, NO3T2I, &
     & NH41TI, NH4T2I, CH4T2I, CH41TI, SO41TI, SO4T2I, HS1TI, HST2I, &
     & SI1TI, SIT2I, CTEMPI, BENSTI, SODI, SED_INIT, SED_INIT2, &
     & SED_READ, SED_CALC
  !

  !
	   Use MOD_SAV, Only: ALEAF, AROOT, ASTEM, ATUBER, AEP, &
     & APLEAF, ABMLEAF, ABMTUBER, APEP, ABMEP, APREP, ASLSH, ANLSAV, &
     & APLSAV, ANLEPI, APLEPI, AFNSED, AFPSED, AFHS, AEPATN, AWATATN, &
     & AFISH, AFIEP, ANPPSAV, ANPPEPI, ADOSAVW, ADOCSAVW, APOCSAVW, &
     & ANH4SAVW, ANO3SAVW, ADONSAVW, APONSAVW, APO4SAVW, ADOPSAVW, &
     & APOPSAVW, ADOEPIW, ADOCEPIW, APOCEPIW, ANH4EPIW, ANO3EPIW, &
     & ADONEPIW, APONEPIW, APO4EPIW, ADOPEPIW, APOPEPIW, ASEDDOSAV, &
     & ASEDPOCSAV, ASEDPONSAV, ASEDPOPSAV, ASEDNH4SAV, ASEDPO4SAV, &
     & SAV_ALLOC, SAV_LOADS_ALLOC, SAV_READ, SAV_INIT, SAV_INIT_ICI, &
	 & LEAF, STEM, ROOT, TUBER, EP, &
	 & LEAFI, STEMI, ROOTI, TUBERI, EPI,NSAVCELL_GL,SAVNDS_GL
  !
#  if defined (MULTIPROCESSOR)
      Use MOD_PAR, Only: EL_PID, NGID, NHN, HN_LST, EC, NC, EXCHANGE
!
#  endif

#  if defined (WET_DRY)
	  USE MOD_WD, ONLY: WD_SET_PARAM , WD_ALLOC, WD_SET_DATA
#  endif
  !
      Use MOD_ZOOP, Only: ANCSZ, ANCLZ, APCSZ, APCLZ, ACLSZ, ACLLZ, &
     & ARSZ, ARLZ, ABMSZ, ABMLZ, AMSZ, AMLZ, APRSZLZ, AGSZ, AGLZ, &
     & ADOCSZ, APOCSZ, ADOCLZ, APOCLZ, ANH4SZ, ADONSZ, APONSZ, ANH4LZ, &
     & ADONLZ, APONLZ, APO4SZ, ADOPSZ, APOPSZ, APO4LZ, ADOPLZ, APOPLZ, &
     & APRSZ, APRLZ, APISZ, APILZ, AB1SZ, AB2SZ, AB3SZ, AB1LZ, AB2LZ, &
     & AB3LZ, ADOSZ, ADOLZ, ASASZ, ASUSZ, ASALZ, ASULZ, ZOOP_ALLOC, &
     & ZOOP_READ
  !
      Use MOD_DIC_TALK, Only: PCO2AIR_READ_INITIAL !
  !
      Implicit None
      Integer :: J, JF, JG, JCON, K, L, N, NC_ID2
      Integer :: I, II, ISF, JP, JNP
      Integer :: JS1, JS2, JS3
      Integer :: YEARSTRT
      Character (Len=1024) :: NCFILE
      Character (Len=4) :: YEAR
      Real (SP), Allocatable :: FTEMP (:)!,FTEMP2(:)
      Real (SP), Allocatable :: RTMP (:, :, :)
      Real (SP), Allocatable :: RTMP2 (:, :), RTMP21 (:, :), RTMP22 (:, &
     & :), RTMP23 (:, :), RTMP24 (:, :), RTMP25 (:, :), RTMP26 (:, :)
      Real (SP) :: TAasFsalt, auxS, auxT ! 
      REAL(SP), DIMENSION(:),ALLOCATABLE :: INI_LEAF, INI_STEM, INI_ROOT, INI_TUBER, INI_EPY
      INTEGER :: SVID
      Character (Len=120) :: ERROR_MSG

	  Character (LEN=1024) :: MKHISOUTDIR
	  Character (LEN=1024) :: MKSTNOUTDIR
  !
  !************************************************************************
  !**                              Inputs                                **
  !************************************************************************
  !
  !: added UNIT_LINKAGE and wqm_linage.in for input control of linkage
  !       to FVCOM
  !***** Open Linkage control file (Linkage to FVCOM)
  !
      Open (UNIT_LINKAGE, File='wqm_linkage.in', Status='OLD')
  !
      NstationNum_GL = 0
      Nstation = 0
  !
      Read (UNIT_LINKAGE, Nml=hydro_netcdf, Err=10)!read $hydro namelist
      Read (UNIT_LINKAGE, Nml=wqm_history, Err=20)!read history  output namelist
      Read (UNIT_LINKAGE, Nml=wqm_stations, Err=30)!read stations output namelist
  !
  !---Below is what wqm_linkage.in should look like-----
  !$hydro_netcdf
  !     hydro_dir='/home/long075/fvcom/outputs/'                 !directory name where the fvcom
  !     hydro_prefix='psm_'
  !     hydro_suffix='.nc'                                       !extentionname
  !     hydro_filenumstart= 1                                    !number in first netcdf input file
  !     hydro_Nrec   = 72                                        !number of records in each netcdf file
  !     hydro_filenumwidth = 4                                   !number of digits in file numbe following hydro_prefix
  !$end
  !$wqm_history
  !    HIS_OUTDIR='./output/'
  !    HISFN = 'output/wqm_history.out'
  !    t_his_start = 1
  !    t_his_end  =  365
  !    t_his_dlt  =  0.25
  !    HISFN_SPLIT_BYLEVEL = .FALSE.
  !
  !$end
  !$wqm_stations
  !    STN_OUTDIR='./output/'
  !    STNFN       = 'output/wqm_station.out'                    !output file name for stations
  !    Nstation    = 12                                          !total number of stations
  !    NStationLoc = 1 2 3 4 5 6 7 8 9 10 11 12                  !node numbers in fvcom model grid for outputing station timeseries
  !    t_stn_start = 1
  !    t_stn_end   = 365
  !    t_stn_dlt   = 0.25
  !$end
  !-------------------------------------------------------------
  !
      Go To 91!if nothing wrong, goto 91 directly
  !
10    ERROR_MSG = 'Error when reading hydro, pls check wqm_linkage.in'
      Go To 90
  !
20    ERROR_MSG = 'ERROR when reading wqm_history, pls check wqm_linkag&
     &e.in'
      Go To 90
  !
30    ERROR_MSG = 'ERROR when reading wqm_stations, pls check wqm_linka&
     &ge.in'
      Go To 90
  !
90    Close (UNIT_LINKAGE)
      Call ERROR_CHECK (ERROR_MSG,-1)
  !
91    Close (UNIT_LINKAGE)
      If (MSR) WRITE (*,*) 'Finished reading wqm_linkage.in ...'
  !
  !find out the FVCOM netcdf output file name (first file)
  !
      NCFILE_PREFIX = TRIM (hydro_prefix)
      NCFILE_SUFFIX = TRIM (hydro_suffix)
      Write (FORMAT_STR, '(A4,I1,A1)') '(I0.', hydro_filenumwidth, ')'
      NCFILE_NUMBER = '' !empty string for NCFILE_NUMBER
  !Populate it with the number hydro_filenumstart
      Write (NCFILE_NUMBER(1:hydro_filenumwidth), TRIM(FORMAT_STR)) &
     & hydro_filenumstart
      NCFILE = TRIM (hydro_dir) // TRIM (NCFILE_PREFIX) // TRIM &
     & (NCFILE_NUMBER) // TRIM (NCFILE_SUFFIX)
  !
  !***** Open control file
  !
      Open (CON, File=CONFN, Status='OLD')
  !
  !***** Title cards
  !
      Read (CON, 1000)
      Read (CON, 1010) (TITLE(J), J=1, 6)
  !

  !
      Read (CON, 1020) MGL, NGL, KB, NOBTY
  !

  !

	  SIGVARC=' OFF'  !bydefault set it to off

      Read (CON, 1031) TMSTRT, TMEND, YEARSTRT, SIGVARC

#if defined (NEWSIGMA)
	 SIGVAR=SIGVARC==' ON'  !spatially variable sigma coordinates rather than uniform
#else
	 SIGVAR=.FALSE.
#endif

      Read (CON, 1040) NDLT
      Read (CON, 1030) (DLTD(J), J=1, NDLT)
      Read (CON, 1030) (DLTVAL(J), J=1, NDLT)
	  DTI=DLTVAL(1)  !
      Read (CON, 1030) (DLTMAX(J), J=1, NDLT)
      Read (CON, 1030) (DLTFTN(J), J=1, NDLT)
  !
  ! AHMDLT IS FVCOM output interval
      Read (CON, 1030) AHMDLT, FILGTH

      JDAY = TMSTRT

  !--SELECT NETCDF FILE TO READ FOR GRID CONSTRUCTION
      If (MSR) WRITE (YEAR, '(I4.4)') YEARSTRT

  !--DETERMINE NUMBER OF ELEMENTS AND NODES IN THE MODEL
  !
  !read global dimensions
      Call NCD_READ_GRID (NCFILE)
  !
  !allocate global variables for geometry related hydrodynamics
      Call HYDRO_GEOM_ALLOC
  !
      If (MSR) WRITE (*,*) '!  # OF NODES            :', MGL
      If (MSR) WRITE (*,*) '!  # OF ELEMENTS         :', NGL
      If (MSR) WRITE (*,*) '!  # OF SIGMA LEVELS     :', KB !KBM1
      If (MSR) WRITE (*,*)
      If (MSR) WRITE (*,*) '!      MESH READING      :    FINISHED'
  !
  !Read global grid coordiantes, sigma levels etc
      Call NCD_READ_SHAPE (NCFILE)
  !
  !
  !--DECOMPOSE DOMAIN BY ELEMENTS USING METIS
  !
#  if defined (MULTIPROCESSOR)
      Allocate (EL_PID(NGL))
      EL_PID = 1
      If (PAR) Call DOMDEC (NGL, NPROCS, EL_PID, MSR)
#  endif
  !
  !
  !--GENERATE GLOBAL<==>LOCAL ELEMENT/NODE MAPPING
  !
#  if defined (MULTIPROCESSOR)
      Call GENMAP !removed "IF(PAR)" 
#  endif
  !
  !Allocate variables related to local hydrodynamics
      Call HYDRO_ALLOC
  !
  !Allocate water quality related variables for local domain
  !
  ! CALL WQM_ALLOC !Moved after CON file is read
  !
      Call ZOOP_ALLOC
  !
      Call WQMINIT_ALLOC
!
  !: allocate all CO2SYS parameters
      Call CO2SYSCONST_ALLOC
!
  !
  !--MAP OPEN BOUNDARY CONDITION NODES TO LOCAL DOMAIN
  !
      Call BCMAP
  !
  !
  !--INPUT AND SETUP BOUNDARY FORCING (HEAT/RIVERS/WIND/etc)
  !
      Call BCS_FORCE
  !
  !
  !--SHIFT GRID/CORIOLIS/BATHYMETRY TO LOCAL DOMAIN
  !
  !
      Call TGE_ALLOC !allocate NV (element to node connectivity)
      Call PDOMDEC !calculate local bathymetry, grid locations etc (using NV)
  !
  !
  !--SET UP GRID METRICS (FLUX EDGES/CONTROL VOLUMES/ETC)
  !
      Call TRIANGLE_GRID_EDGE !Set up fluxes and control Volumes
  !
      Call CELL_AREA !Calculate Element and Control Volume Areas
  !
  !--SET ISBCE AND ISONB CORRECTLY IN HALO CELLS/NODES
#  if defined (MULTIPROCESSOR)
      Allocate (FTEMP(0:NTLOC))
      FTEMP = ISBCE
      If (PAR) Call EXCHANGE (EC, NTLOC, 1, MYID, NPROCS, FTEMP)
      ISBCE = FTEMP
      Deallocate (FTEMP)
      Allocate (FTEMP(0:MTLOC))
      FTEMP = ISONB
      If (PAR) Call EXCHANGE (NC, MTLOC, 1, MYID, NPROCS, FTEMP)
      ISONB = FTEMP
      Deallocate (FTEMP)
#  endif
  !
  !
  !--EXCHANGE SHAPE FACTOR INFORMATION
  !
#  if defined (MULTIPROCESSOR)
  !   IF(PAR)CALL EXCHANGE(EC,NTLOC,4,MYID,NPROCS,A1U,A2U)
  !   IF(PAR)CALL EXCHANGE(EC,NTLOC,3,MYID,NPROCS,AWX,AWY,AW0)
  !   IF(PAR)CALL EXCHANGE(EC,NTLOC,1,MYID,NPROCS,ALPHA)
      If (PAR) Call EXCHANGE (EC, NTLOC, 1, MYID, NPROCS, ART)
#  endif
  !
  !***** Output control
  !
  !
      Read (CON, 1050) SNPC, NSNP
      Read (CON, 1030) (SNPD(J), J=1, NSNP)
      Read (CON, 1030) (SNPF(J), J=1, NSNP)
      Read (CON, 1055) PLTC, QPLTC, SPLTC, SAVPLTC, NPLT
      Read (CON, 1030) (PLTD(J), J=1, NPLT)
      Read (CON, 1030) (PLTF(J), J=1, NPLT)
      Read (CON, 1050) APLTC, NAPL
      Read (CON, 1030) (APLTD(J), J=1, NAPL)
      Read (CON, 1030) (APLF(J), J=1, NAPL)
      Read (CON, 1050) TFLC, NTFL
      Read (CON, 1030) (TFLD(J), J=1, NTFL)
      Read (CON, 1030) (TFLF(J), J=1, NTFL)
      Read (CON, 1052) KFLC, NKFL
      Read (CON, 1030) (KFLD(J), J=1, NKFL)
      Read (CON, 1030) (KFLF(J), J=1, NKFL)
      Read (CON, 1052) OPLC, NOPL, NOINT
      Read (CON, 1030) (OINT(J), J=1, NOINT)
      Read (CON, 1030) (OPLD(J), J=1, NOPL)
      Read (CON, 1030) (OPLF(J), J=1, NOPL)
      Read (CON, 1050) MBLC, NMBL
      Read (CON, 1030) (MBLD(J), J=1, NMBL)
      Read (CON, 1030) (MBLF(J), J=1, NMBL)
      Read (CON, 1050) DIAC, NDIA
      Read (CON, 1030) (DIAD(J), J=1, NDIA)
      Read (CON, 1030) (DIAF(J), J=1, NDIA)
      Read (CON, 1050) RSOC, NRSO, RSIC
      Read (CON, 1030) (RSOD(J), J=1, NRSO)
  !
  !
  !***** Hydrodynamic solution scheme
  !

      Read (CON, 1060) SLC, CONSC, TH, MINSTEP

  !
  !***** Controls
  !
      Read (CON, 1071) SEDC, AUTOC, VBC, BFOC, STLC, ICIC, ICOC, SAVMC, &
     & BAOC, SEDTC !

#if !defined (SED_TEMP)
	 SEDTC ='OFF'  !Set it to OFF if  SED_TEMP not defined
#endif
  !
  ! suspension and deposit feeders and benthic algae
  !
      Read (CON, 1072) SFLC, DFLC, SFLOXC, DFLOXC, BALC, DFOC !: added DFLOXC, and DFOC
  !
  !
  !***** Dead sea case
  !
      Read (CON, 1070) FLC, XYDFC, ZDFC
  !
  !
  !***** Dispersion
  !
      Read (CON, 1030) XYDF, ZDFMUL, ZDFBCK
  !
  !
  !
  !***** Dispersion type (closure (ON) or constant (OFF))
  !
      Read (CON, 1070) XYDFU
  !
  !
  !***** Constituent control cards
  !
      Read (CON, 1070) BCC, S1C, S2C, S3C, MDC, BFC, ATMC, SAVLC
  !
  !
  !
      Read (CON, 1080) REDS1C, REDS1N, REDS1P, REDS2C, REDS2N, REDS2P, &
     & REDS3C, REDS3N, REDS3P
  !
  !
      Read (CON, 1080) REDCBC, REDCBN, REDCBP
  !
  !
  !
      Read (CON, 1065) BNDTC
  !
  !
      Read (CON, 1070) (ACC(JCON), JCON=1, NCP)
  !
  !
  !
  !***** Input filenames
  !
      Read (CON, 1020) NHYDF, NTVDF
  !
  !
      Read (CON, 1090) MAPFN
  !
  !
      Read (CON, 1090) GEOFN
  !
  !
  !
      Read (CON, 1090) ICIFN
  !
  !
  !
      Read (CON, 1090) RSIFN
  !
  !
      Read (CON, 1090) AGRFN
  !
  !
  !
      Read (CON, 1090) ZOOFN
  !
  !
      Read (CON, 1090) SFIFN
  !
  !
      Read (CON, 1090) STLFN
  !
  !
  !
      Read (CON, 1090) MRLFN 
  !
  !
      Read (CON, 1090) KEIFN
  !
  !
      Read (CON, 1090) (METFN(J), J=1, NTVDF)
      If (MSR) Then
         Write (*,*) (METFN(J), J=1, NTVDF)
      End If
      Read (CON, 1090) (S1FN(J), J=1, NTVDF)
      Read (CON, 1090) (S2FN(J), J=1, NTVDF)
      Read (CON, 1090) (S3FN(J), J=1, NTVDF)
      Read (CON, 1090) (ATMFN(J), J=1, NTVDF)
      Read (CON, 1090) (SVIFN(J), J=1, NTVDF)
      Read (CON, 1090) (BFIFN(J), J=1, NTVDF)
      Read (CON, 1090) (CBCFN(J), J=1, NTVDF)
      Read (CON, 1090) BAIFN !Benthic algae input file
      Read (CON, 1090) DFIFN !Deposition feeder input file name
  !***** Output filenames
  !
      Read (CON, 1090) ICOFN !MNOEL  2-5-93
      Read (CON, 1090) SNPFN
      Read (CON, 1090) RSOFN
      Read (CON, 1090) PLTFN
      Read (CON, 1090) APLFN
      Read (CON, 1090) DIAFN
      Read (CON, 1090) TFLFN
      Read (CON, 1090) KFLFN
      Read (CON, 1090) OPLFN
      Read (CON, 1090) MBLFN
      Read (CON, 1090) ALOFN
      Read (CON, 1090) ZFOFN
      Read (CON, 1090) BFOFN
      Read (CON, 1090) SVOFN
      Read (CON, 1090) SFOFN
      Read (CON, 1090) BAOFN !  benthic algae output file name
      Read (CON, 1090) DFOFN !  Deposition feeder output file name
      Close (CON)
  !
  !***** Initialize I/O logical control variables
  !
      BOUNDARY_CONC = BCC == ' ON'
      MODIFY_ICONC = MDC == ' ON'
      SOURCE_ONE = S1C == ' ON' !
      BENTHIC_FLUXES = BFC == ' ON'
      SOURCE_TWO = S2C == ' ON' !
      SOURCE_THR = S3C == ' ON' !
      ATMOS_LOADS = ATMC == ' ON'
      SAV_LOADS = SAVLC == ' ON'
      SAV_CALC = SAVMC == ' ON'
  !
  !
  !Make sure SAV_LOADS and SAV_CALC are not on at the same time
  !we can have only one of them or neither of them
      If (SAV_LOADS .And. SAV_CALC) Then
         ERROR_MSG = 'Oops, SAV_LOADS and SAV_CALC can''t be on at same&
        & time!!'
         ERROR_MSG = TRIM (ERROR_MSG) // 'Please check SAVMC and SAVLC &
        &in input file'
         Call ERROR_CHECK (ERROR_MSG,-1)
      End If
  !
      SETTLING = STLC == ' ON'
      ICOND_OUT = ICOC == ' ON' !
      SEDIMENT_CALC = SEDC == ' ON'
#if defined (SED_TEMP)
	  SEDTEMP_CALC = SEDTC == ' ON'  !Added flag to turn on /off sediment temperature calculation
#else
	  SEDTEMP_CALC = .FALSE.  !Do not do sediment temperature calculation (use bottom water temperature)
	                          !or air temperature (for dry land)
#endif
									 !
  !
  !PRINT error if BOTH SEDC and BFC are ' ON'
  !
      If (MSR) Then
         If (BENTHIC_FLUXES .And. SEDIMENT_CALC) Then !fix the format of this error message later
            Write (*,*) 'Oops, should not have BFC and SEDC ON at same &
           &time'
            Call PSTOP
         End If
      End If
!
  !IF TDIC sedim flux calculated from CSOD, we must have SEDC=ON and BFC = OFF.
#if defined(CSOD_CALC)
      If (BENTHIC_FLUXES) Then !fix the format of this error message later
     !WRITE(*,*)'Oops, should have BFC=OFF and SEDC=ON to run with CSOD_CALC'
     !CALL PSTOP
!
         Write (*,*) 'NOTE: you have BFC=ON and SEDC=OFF (no sediment c&
        &alc, even if CSOD_CALC is defined)'
      End If
#endif
!
      BALGAE_CALC = BALC == ' ON' ! benthic algae calculation
      KINETIC_FLUXES = KFLC == ' ON'
      LIGHT_EXTINCTION = ((ACC(4) == ' ON') .Or. (ACC(5) == ' ON') .Or. &
     & (ACC(6) == ' ON'))
      ZOO_CALC = ((ACC(7) == ' ON') .Or. (ACC(8) == ' ON'))
  !   UNI_ICON_IN      = ICIC == ' UNIFORM'              !MNOEL 2-5-93
      UNI_ICON_IN = ((ICIC == ' UNIFORM') .Or. (ICIC == ' WU_SEDV'))!: to read spatially variable sedims ini cond
      UNI_ICON_IN_SED_VAR = ICIC == ' WU_SEDV' !: to read spatially variable sedims ini cond
      BIN_ICON_IN = ICIC == '  BINARY' !
      RES_ICON_IN = ICIC == ' RESTART'
  ! deposit and suspension feeder control
  !
      DFEEDER = DFLC == ' ON'
      SFEEDER = SFLC == ' ON'
      HYPOXFX_SF = SFLOXC == ' ON'
      HYPOXFX_DF = DFLOXC == ' ON'
  !
  !***** Active Constituents
  !
      NAC = 0
      Do JCON = 1, NCP
         If (ACC(JCON) == ' ON') Then
            NAC = NAC + 1
            AC (NAC) = JCON
         End If
      End Do
  !
  !-------
  ! sanity check: if only one of TALK/TDIC are defined, define the other one
      If (ACC(33) == ' ON' .And. ACC(34) == 'OFF') Then
         NAC = NAC + 1
         AC (NAC) = 34
      End If
      If (ACC(33) == 'OFF' .And. ACC(34) == ' ON') Then
         AC (NAC) = 33 !Keep TDIC at position 33
         NAC = NAC + 1
         AC (NAC) = 34 !Keep TALK at position 34
      End If
  !-------
  !
      VOLUME_BALANCE = VBC == ' ON'
      FLOW = FLC == ' ON'
      PLOTS = PLTC == ' ON'
      MASS_BALANCE = MBLC == ' ON'
      OXYGEN_PLOTS = OPLC == ' ON'
      SNAPSHOTS = SNPC == ' ON'
      BENTHIC_OUTPUT = BFOC == ' ON'
      DFEEDER_OUTPUT = DFOC == ' ON' ! deposition feeder output
      BA_OUTPUT = BAOC == ' ON' ! benthic algae output
      TRANSPORT_FLUXES = TFLC == ' ON' .Or. MBLC == ' ON'
      RESTART_OUT = RSOC == ' ON'
      RESTART_IN = RSIC == ' ON'
      DIAGNOSTICS = DIAC == ' ON'
      Z_DIFFUSION = ZDFC == ' ON'
      AVERAGE_PLOTS = APLTC == ' ON'
      QUALITY_DIAG = QPLTC == ' ON'
      AUTO_STEPPING = AUTOC == ' ON'
      XY_DIFFUSION = XYDFC == ' ON'
      UPWIND = SLC == '  UPWIND'
      QUICKEST = SLC == 'QUICKEST'
      STEP_BOUNDARY = BNDTC == '    STEP'
      CONSERVE_MASS = CONSC == '    MASS'
      SEDIMENT_DIAG = SPLTC == ' ON' .And. (SEDIMENT_CALC .Or. &
     & BENTHIC_FLUXES)
      SAV_PLOTS = SAVPLTC == ' ON' .And. (SAV_LOADS .Or. SAV_CALC)
      END_RUN = .False.
  !
      If (ACC(1) == ' ON') TEMPERATURE_CALC = .True.
      If (ACC(3) == ' ON') SOLIDS_CALC = .True.
      If ((ACC(4) == ' ON') .Or. (ACC(5) == ' ON') .Or. (ACC(6) == ' ON&
     &')) ALGAE_CALC = .True.
      If (ACC(9) == ' ON') CARBON_CALC = .True.
      If (ACC(13) == ' ON') NITROGEN_CALC = .True.
      If (ACC(20) == ' ON') PHOSPHORUS_CALC = .True.
      If (ACC(26) == ' ON') COD_CALC = .True.
      If (ACC(27) == ' ON') OXYGEN_CALC = .True.
      If (ACC(29) == ' ON') SILICA_CALC = .True.
      If (ACC(25) == ' ON') PIP_CALC = .True.
!
      If (ACC(33) == ' ON') CARBONATE_CALC = .True. !
!
	  ! Get NASVM for matrix allocation 
	  If (SAV_CALC) Then
     !read SVIFN file and find maximum number of SAV species - NSAVM
     !
         Open (SVI, File=SVIFN(1), Status='OLD')
     !
         Read (SVI, '(///)')!jump 4th lines (first skip 3 lines and then read nothing from 4th line)
         Read (SVI, '(A72)') (TITLE(J), J=1, 6)!read six lines
         Read (SVI, '(//8X,I8)') NSAVM !go down two lines and read the third line
	     If (MSR) write (*,*) '  --> Max number of SAV =', NSAVM
         Close (SVI)
      End If
!
  !Have to allocate after LIGHT_EXTINCTION is set to True or False
      Call WQM_ALLOC
  !
  !***** Open remaining input files
  !
      HYDPTR = 1
      METPTR = 1
      CBCPTR = 1
      S1PTR = 1
      S2PTR = 1
      S3PTR = 1
      BFIPTR = 1
      BAOPTR = 1
      KEIPTR = 1
      ATMPTR = 1
      SAVPTR = 1
  !
      Open (MET, File=METFN(METPTR), Status='OLD')
      If (SETTLING) OPEN (STL, FILE=STLFN, STATUS='OLD')
  !   IF (LIGHT_EXTINCTION) OPEN (KEI,FILE=KEIFN,STATUS='OLD')    
      If (ALGAE_CALC) OPEN (AGR, FILE=AGRFN, STATUS='OLD')!LIGHT_EXTINCTION -> ALGAE_CALC
      If (SEDIMENT_CALC .Or. BENTHIC_FLUXES) OPEN (BFI, &
     & FILE=BFIFN(BFIPTR), STATUS='OLD')
  !
      If (BALGAE_CALC) OPEN (BAI, FILE=BAIFN, STATUS='OLD')
  !
      If (BOUNDARY_CONC) OPEN (CBC, FILE=CBCFN(CBCPTR), STATUS='OLD')
      If (SOURCE_ONE) OPEN (S1, FILE=S1FN(S1PTR), STATUS='OLD')
      If (SOURCE_TWO) OPEN (S2, FILE=S2FN(S2PTR), STATUS='OLD')
      If (SOURCE_THR) OPEN (S3, FILE=S3FN(S3PTR), STATUS='OLD')
      If (ATMOS_LOADS) OPEN (ATM, FILE=ATMFN(ATMPTR), STATUS='OLD')
  !

  !
      If (Nstation > NSTATIONMAX) Then
         ERROR_MSG = ''
         Write (ERROR_MSG, '(I5)') NSTATIONMAX
         ERROR_MSG = 'Nstation exceeding maximum value allowed:' // &
        & TRIM (ERROR_MSG)
         Call ERROR_CHECK (ERROR_MSG,-1)
      End If
  !
      If (Nstation <= 0 .Or. MAXVAL(NstationNum_GL(1:Nstation), 1) > &
     & MGL .Or. MINVAL(NstationNum_GL(1:Nstation), 1) < 1) Then
         ERROR_MSG = 'Error in wqm_stations section, please check wqm_l&
        &inkage.in'
         Call ERROR_CHECK (ERROR_MSG,-1)
      End If
  !
      If (t_stn_dlt < 0.0) Then
         ERROR_MSG = 'Error wqm_stations section, t_stn_dlt must be > 0&
        &. Pls check wqm_linkage.in'
         Call ERROR_CHECK (ERROR_MSG,-1)
      End If
  !
      If (t_his_dlt <= 0.0) Then
         ERROR_MSG = 'Error in wqm_history section, t_his_dlt must be >&
        & 0. Pls check wqm_linkage.in'
         Call ERROR_CHECK (ERROR_MSG,-1)
      End If

	  !check if output dir exists

#  if !defined (DOS)
			MKHISOUTDIR = "mkdir -p "//TRIM(his_outdir)
			MKSTNOUTDIR = "mkdir -p "//TRIM(stn_outdir)

			IF(MSR)THEN
#  		if !defined (CRAY)
				CALL SYSTEM( TRIM(MKHISOUTDIR) )
				CALL SYSTEM( TRIM(MKSTNOUTDIR) )
#		else
				CALL CRAY_SYSTEM_CALL(TRIM(MKHISOUTDIR))
				CALL CRAY_SYSTEM_CALL(TRIM(MKSTNOUTDIR))
#  		endif
			END IF
#  endif


  !
  ! :
  ! Read two time steps in FVCOM netcdf output files. Note that this may have to be moved
  ! to after the control file reading is finished. For model restart, one would
  ! need to find out what time record to read based on restarting time, as opposed to
  ! reading the first two records of all files combined.

  ! coarsegrid = read every timestep, open new file every 2 simulation hours
  !  DAY = DAY + 1
      IJDAY = AINT (JDAY)
      If (MSR) WRITE (*,*) 'IJDAY, JDAY', IJDAY, JDAY
  !
  !
  !  R_HOUR = R_HOUR + 1
  !  R_HOUR_INT = R_HOUR_INT + 1
  !  R_HOUR_first = R_HOUR
  !`  WRITE(FILENUMBER,'(I4.4)') R_HOUR_first
  !
      IFNC = 0 !initialize indices for file
      NTRECNC = 0 !indices for time record in each file
      NTHYDRO = 0 !indices for overall hydrodyanmics
  !
  ! The first record
  !
      IFNC = 1 !
      NTRECNC = 1 !read from the first record for the current file
      NTHYDRO = 1 !the first record over all


  !
      NCFILE_NUMBER = '' !empty string for NCFILE_NUMBER
  !Populate it with the number with hydro_filenumstart
  !
      Write (NCFILE_NUMBER(1:hydro_filenumwidth), TRIM(FORMAT_STR)) &
     & (IFNC+hydro_filenumstart-1)
      NCFILE = TRIM (hydro_dir) // TRIM (NCFILE_PREFIX) // TRIM &
     & (NCFILE_NUMBER) // TRIM (NCFILE_SUFFIX)
  !
      Call NCD_READ_OPEN (NCFILE, UNC1, VNC1, WTSNC1, UARD_OBCNNC1, &
     & XFLUX_OBCNC1, DTFANC1, KHNC1, ELNC1, TNC1, SNC1, NTRECNC)!reading NTRECNC'th record in NCFILE

	  !Commented these variables for the previous timestep

	  UU = UNC1
      VV = VNC1
      WTS = WTSNC1
      UARD_OBCN = UARD_OBCNNC1
      XFLUX_OBC = XFLUX_OBCNC1
      DTFA = DTFANC1
      KH = KHNC1
      T = TNC1
      SALT = SNC1
      EL = ELNC1
      D = H + EL
	  ET = ELNC1   !
      DT = H + ET

#if defined (WET_DRY)
	    CALL WD_SET_PARAM
		CALL WD_ALLOC
		CALL WD_SET_DATA
#endif

  !
  !The seocnd record
  !
      NTRECNC = NTRECNC + 1 !point to next record
      NTHYDRO = NTHYDRO + 1
  !
      If (NTRECNC .Gt. hydro_Nrec) Then
     !
         NTRECNC = 1 !reset time index for file to be 1
         IFNC = IFNC + 1 !A new file
     !
         NCFILE_NUMBER = ''
         Write (NCFILE_NUMBER(1:hydro_filenumwidth), TRIM(FORMAT_STR)) &
        & (IFNC+hydro_filenumstart-1)
         NCFILE = TRIM (hydro_dir) // TRIM (NCFILE_PREFIX) // TRIM &
        & (NCFILE_NUMBER) // TRIM (NCFILE_SUFFIX)
     !
         If (MSR) WRITE (*,*) TRIM (NCFILE)
     !
         Call NCD_READ_OPEN (NCFILE, UNC2, VNC2, WTSNC2, UARD_OBCNNC2, &
        & XFLUX_OBCNC2, DTFANC2, KHNC2, ELNC2, TNC2, SNC2, NTRECNC)
     !
      Else
     !
         Call NCD_READ (NCFILE, UNC2, VNC2, WTSNC2, UARD_OBCNNC2, &
        & XFLUX_OBCNC2, DTFANC2, KHNC2, ELNC2, TNC2, SNC2, NTRECNC)
     !
      End If
  !
  !***** Suspension feeders
      If (SFEEDER) Then
         Call SF_ALLOC
      End If
  !
  !
      If (SFEEDER) Then
     !
         Call SF_READ
     !
     ! initialize arrays that are passed in case susp. feeder model not used
         Call SF_INIT
     !
      End If
  !
  !
  !allocate Depoition Feeder arrays
      If (DFEEDER) Then
         Call DF_ALLOC
         Call DF_READ
         Call DF_INIT
      End If
  !
  !
  !
  !************************************************************************
  !**                   Parameters for SAV Submodel                      **
  !************************************************************************
  !
  !
      If (SAV_CALC) Then

     !
     !Use the maximum number of SAV species to allocate SAV arrays
         Call SAV_ALLOC
     !read SAV parameters
         Open (SVI, File=SVIFN(1), Status='OLD')
         Call SAV_READ
         Close (SVI)
		 Call SAV_INIT !AN
      End If
  !
  !: Open SVI for tvds.F to read if SAV_LOADS is ON
      If (SAV_LOADS) Then
         Call SAV_LOADS_ALLOC
         Open (SVI, File=SVIFN(SAVPTR), Status='OLD')
      End If
  !
  !******* Initial conditions
  !
  !******* Uniform constituent initial concentrations
  !
	  If (MSR) write (*,*) '  --> UNI_ICON_IN =', UNI_ICON_IN
   !
      If (UNI_ICON_IN) Then 
         Open (ICI, File=ICIFN, Status='OLD')
         Read (ICI,*)
     !
     !
     ! The above now in parallel, 
         Allocate (RTMP(0:MGL, KBM1, NCP))
         RTMP = 0.0 !, note index starts from zero
         Do I = 1, MGL
            Do K = 1, KBM1
               Read (ICI,*) (RTMP(I, K, JCON), JCON=1, NCP)
            End Do
         End Do
     !
         If (SERIAL) C1 = RTMP
     !
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do JCON = 1, NCP
               Do K = 1, KBM1
                  Do I = 1, MLOC
                     C1 (I, K, JCON) = RTMP (NGID(I), K, JCON)
                     C2 (I, K, JCON) = RTMP (NGID(I), K, JCON)
                  End Do
                  Do I = 1, NHN
                     C1 (I+MLOC, K, JCON) = RTMP (HN_LST(I), K, JCON)
                     C2 (I+MLOC, K, JCON) = RTMP (HN_LST(I), K, JCON)
                  End Do
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP)
         C2 = C1
         C1MIN = 1.E10
         C1MAX = 0.

         T = TNC1
         SALT = SNC1

         Call DENS2
     !
#if defined TA_OBC_INI_FROM_S || defined TA_INI_FROM_S
     !---------Make TALK (variable=34) a function of Salinity at initial conditions - LB Jan 8 2016-----!
         Write (*,*) '!  TALK at Initial time :    function of salinity&
        &'
!
         Do K = 1, KBM1
            Do I = 1, MTLOC
               TAasFsalt = (47.7_SP*SALT(I, K)+646.7_SP)!umol/kg-SW
           !If(MSR) Write (*,*) 'LBnote I2. I= ',I,' (47.7_SP * SALT(I, K) + 646.7_SP)= ',TAasFsalt
               TAasFsalt = TAasFsalt * RHO (I, K) / 1000.0_SP !mmol/m3
           !
               C1 (I, K, 34) = TAasFsalt
               C2 (I, K, 34) = TAasFsalt
!
           !If(MSR) Write (*,*) 'LBnote I2.      above * RHO(I, K) / 1000.0_SP= ',TAasFsalt
           !If(MSR) Write (*,*) 'LBnote I2.      SALT(I, K),RHO(I, K)= ',SALT(I, K),RHO(I, K)
           !If(MSR) Write (*,*) 'LBnote I2.      SALT(I, K),RHO(I, K)= ',SALT(I, K),RHO(I, K)
            End Do
         End Do
     
#endif
     !Go back to T and SALT equal to zero
         T = 0.0_SP
         SALT = 0.0_SP

         If (SEDIMENT_CALC) Then
            If (UNI_ICON_IN_SED_VAR) Then ! added: Sediment initial conditions are spatially variable, taken from previous run
           !
               Call SED_INIT ()!  Setup sediment diagenesis module (including SED_ALLOC)
           !
               Read (ICI, 1029) (CTEMP_GL(I), I=1, MGL)!Temperature (deg-C)
               Read (ICI, 1029) ((CPOP_GL(I, JG), I=1, MGL), JG=1, 3)!3G POP (mg-P/m^3)
               Read (ICI, 1029) ((CPON_GL(I, JG), I=1, MGL), JG=1, 3)!3G PON (mg-N/m^3)
               Read (ICI, 1029) ((CPOC_GL(I, JG), I=1, MGL), JG=1, 3)!3G POC (mg-C/m^3)
               Read (ICI, 1029) (CPOS_GL(I), I=1, MGL)!CPOS_GL  is global array of sediment diagenesis particulate biogenic silica
               Read (ICI, 1029) (PO4T2TM1S_GL(I), I=1, MGL)!PO4T2TM1S_GL is global array of PO4 in sediment diagenesis model (mg-P/m^3)
               Read (ICI, 1029) (NH4T2TM1S_GL(I), I=1, MGL)!NH4T2TM1S_GL is global array of NH4 in sediment diagenesis model (mg-N/m^3)
               Read (ICI, 1029) (NO3T2TM1S_GL(I), I=1, MGL)!NO3T2TM1S_GL is global array of NO3 in sediment diagenesis model (mg-N/m^3)
               Read (ICI, 1029) (HST2TM1S_GL(I), I=1, MGL)!HST2TM1S_GL  is global array of H2S in sediment diagenesis model (mg-OxygenEq/m^3)
               Read (ICI, 1029) (CH4T2TM1S_GL(I), I=1, MGL)!CH4T2TM1S_GL is global array of CH4 in 2nd (lower) layer of sediment
               Read (ICI, 1029) (CH41TM1S_GL(I), I=1, MGL)!CH41TM1S_GL  is global array of CH4 in 1st (upper) layer of sediment
               Read (ICI, 1029) (SO4T2TM1S_GL(I), I=1, MGL)!SO4T2TM1S_GL is global array of SO4 in sediment diagenesis model
               Read (ICI, 1029) (SIT2TM1S_GL(I), I=1, MGL)!SIT2TM1S_GL  is global array of sediment diagenesis dissolved silica
               Read (ICI, 1029) (BENSTRTM1S_GL(I), I=1, MGL)!BENSTR1S_GL  is global array of benthic stress (Unit ???)
               Read (ICI, 1029) (SODTM1S_GL(I), I=1, MGL)
               Read (ICI, 1029) (PO41TM1S_GL(I), I=1, MGL)
               Read (ICI, 1029) (NH41TM1S_GL(I), I=1, MGL)
               Read (ICI, 1029) (NO31TM1S_GL(I), I=1, MGL)
               Read (ICI, 1029) (HS1TM1S_GL(I), I=1, MGL)
               Read (ICI, 1029) (SI1TM1S_GL(I), I=1, MGL)
           !
           !
               If (BALGAE_CALC) Then
                  Read (ICI, 1029) (BBM_GL(I), I=1, MGL)
               End If

               If (DFEEDER) Then
                  Read (ICI, 1029) (DFEEDM1S_GL(I), I=1, MGL)!DFEEDSM1S_GL is global array of deposition feeder biomass (mgC/m^2)
               End If
               If (SFEEDER) Then
                  Read (ICI, 1029) ((SFEED_GL(I, K), I=1, MGL), K=1, &
                 & NSPECIES)!SFEED_GL is global array of suspension feeder biomass (mg/m^2)
               End If

!-----------------------------------------------------------------------------------------
        IF (SAV_CALC) THEN
                  Read (ICI, 1029) ((LEAF_GL(I, K), I=1, MGL), K=1, &
                  & NSAVM)! AN
                  Read (ICI, 1029) ((STEM_GL(I, K), I=1, MGL), K=1, &
                  & NSAVM)! AN
                  Read (ICI, 1029) ((ROOT_GL(I, K), I=1, MGL), K=1, &
                  & NSAVM)! AN
                  Read (ICI, 1029) ((TUBER_GL(I, K), I=1, MGL), K=1, &
                  & NSAVM)! AN
                  Read (ICI, 1029) ((EP_GL(I, K), I=1, MGL), K=1, &
                  & NSAVM)! AN
         END IF
!---------------------------------------------------------------------------------------
		    !
               If (SERIAL) Then
                  CTEMP = CTEMP_GL
                  CPOP = CPOP_GL
                  CPON = CPON_GL
                  CPOC = CPOC_GL
                  CPOS = CPOS_GL
              !
                  PO4T2TM1S = PO4T2TM1S_GL
                  NH4T2TM1S = NH4T2TM1S_GL
                  NO3T2TM1S = NO3T2TM1S_GL
                  HST2TM1S = HST2TM1S_GL
                  CH4T2TM1S = CH4T2TM1S_GL
                  CH41TM1S = CH41TM1S_GL
                  SO4T2TM1S = SO4T2TM1S_GL
                  SIT2TM1S = SIT2TM1S_GL
                  BENSTRTM1S = BENSTRTM1S_GL
                  SODTM1S = SODTM1S_GL
              !
              !: save other sediment variables needed for restarts
                  PO41TM1S = PO41TM1S_GL
                  NH41TM1S = NH41TM1S_GL
                  NO31TM1S = NO31TM1S_GL
                  HS1TM1S = HS1TM1S_GL
                  SI1TM1S = SI1TM1S_GL
              !
                  If (BALGAE_CALC) Then
                     BBM = BBM_GL
                  End If
                  If (DFEEDER) Then
                     DFEEDM1S = DFEEDM1S_GL
                  End If
                  If (SFEEDER) Then
                     SFEED = SFEED_GL
                  End If
				  If (SAV_CALC) Then
                     LEAF = LEAF_GL
					 STEM = STEM_GL
					 ROOT = ROOT_GL
					 TUBER = TUBER_GL
					 EP = EP_GL
                  End If
               End If
           !
# if defined (MULTIPROCESSOR)
               If (PAR) Then
                  Do I = 1, MLOC
                     CTEMP (I) = CTEMP_GL (NGID(I))
                     CPOP (I, :) = CPOP_GL (NGID(I), :)
                     CPON (I, :) = CPON_GL (NGID(I), :)
                     CPOC (I, :) = CPOC_GL (NGID(I), :)
                     CPOS (I) = CPOS_GL (NGID(I))
                     PO4T2TM1S (I) = PO4T2TM1S_GL (NGID(I))
                     NH4T2TM1S (I) = NH4T2TM1S_GL (NGID(I))
                     NO3T2TM1S (I) = NO3T2TM1S_GL (NGID(I))
                     HST2TM1S (I) = HST2TM1S_GL (NGID(I))
                     CH4T2TM1S (I) = CH4T2TM1S_GL (NGID(I))
                     CH41TM1S (I) = CH41TM1S_GL (NGID(I))
                     SO4T2TM1S (I) = SO4T2TM1S_GL (NGID(I))
                     SIT2TM1S (I) = SIT2TM1S_GL (NGID(I))
                     BENSTRTM1S (I) = BENSTRTM1S_GL (NGID(I))
                     SODTM1S (I) = SODTM1S_GL (NGID(I))
                 !
                 !: save other sediment variables needed for restarts
                     PO41TM1S (I) = PO41TM1S_GL (NGID(I))
                     NH41TM1S (I) = NH41TM1S_GL (NGID(I))
                     NO31TM1S (I) = NO31TM1S_GL (NGID(I))
                     HS1TM1S (I) = HS1TM1S_GL (NGID(I))
                     SI1TM1S (I) = SI1TM1S_GL (NGID(I))
                 !
                     If (BALGAE_CALC) Then
                        BBM (I) = BBM_GL (NGID(I))
                     End If
                     If (DFEEDER) Then
                        DFEEDM1S (I) = DFEEDM1S_GL (NGID(I))
                     End If
                     If (SFEEDER) Then
                        Do K = 1, NSPECIES
                           SFEED (I, K) = SFEED_GL (NGID(I), K)
                        End Do
                     End If
					 If (SAV_CALC) Then
                        Do K = 1, NSAVM
                           LEAF (I, K) = LEAF_GL (NGID(I), K)
						   STEM (I, K) = STEM_GL (NGID(I), K)
						   ROOT (I, K) = ROOT_GL (NGID(I), K)
						   TUBER (I, K) = TUBER_GL (NGID(I), K)
						   EP (I, K) = EP_GL (NGID(I), K)
                        End Do
                     End If
                  End Do
                  Do I = 1, NHN
                     CTEMP (I+MLOC) = CTEMP_GL (HN_LST(I))
                     CPOP (I+MLOC, :) = CPOP_GL (HN_LST(I), :)
                     CPON (I+MLOC, :) = CPON_GL (HN_LST(I), :)
                     CPOC (I+MLOC, :) = CPOC_GL (HN_LST(I), :)
                     CPOS (I+MLOC) = CPOS_GL (HN_LST(I))
                     PO4T2TM1S (I+MLOC) = PO4T2TM1S_GL (HN_LST(I))
                     NH4T2TM1S (I+MLOC) = NH4T2TM1S_GL (HN_LST(I))
                     NO3T2TM1S (I+MLOC) = NO3T2TM1S_GL (HN_LST(I))
                     HST2TM1S (I+MLOC) = HST2TM1S_GL (HN_LST(I))
                     CH4T2TM1S (I+MLOC) = CH4T2TM1S_GL (HN_LST(I))
                     CH41TM1S (I+MLOC) = CH41TM1S_GL (HN_LST(I))
                     SO4T2TM1S (I+MLOC) = SO4T2TM1S_GL (HN_LST(I))
                     SIT2TM1S (I+MLOC) = SIT2TM1S_GL (HN_LST(I))
                     BENSTRTM1S (I+MLOC) = BENSTRTM1S_GL (HN_LST(I))
                     SODTM1S (I+MLOC) = SODTM1S_GL (HN_LST(I))
                 !: save other sediment variables needed for restarts
                     PO41TM1S (I+MLOC) = PO41TM1S_GL (HN_LST(I))
                     NH41TM1S (I+MLOC) = NH41TM1S_GL (HN_LST(I))
                     NO31TM1S (I+MLOC) = NO31TM1S_GL (HN_LST(I))
                     HS1TM1S (I+MLOC) = HS1TM1S_GL (HN_LST(I))
                     SI1TM1S (I+MLOC) = SI1TM1S_GL (HN_LST(I))
                     If (BALGAE_CALC) Then
                        BBM (I+MLOC) = BBM_GL (HN_LST(I))
                     End If
                     If (DFEEDER) Then
                        DFEEDM1S (I+MLOC) = DFEEDM1S_GL (HN_LST(I))
                     End If
                     If (SFEEDER) Then
                        Do K = 1, NSPECIES
                           SFEED (I+MLOC, K) = SFEED_GL (HN_LST(I), K)
                        End Do
                     End If
					 If (SAV_CALC) Then
                        Do K = 1, NSAVM
                           LEAF (I+MLOC, K) = LEAF_GL (HN_LST(I), K)
						   STEM (I+MLOC, K) = STEM_GL (HN_LST(I), K)
						   ROOT (I+MLOC, K) = ROOT_GL (HN_LST(I), K)
						   TUBER (I+MLOC, K) = TUBER_GL (HN_LST(I), K)
						   EP (I+MLOC, K) = EP_GL (HN_LST(I), K)
                        End Do
                     End If
                  End Do
               End If
# endif
           !
           !
            Else !Sediment initial conditions are spatially uniform
               Read (ICI, 1030) CTEMPI
           !
               Read (ICI, 1030) (CPOPI(JG), JG=1, 3)
               Read (ICI, 1030) (CPONI(JG), JG=1, 3)
               Read (ICI, 1030) (CPOCI(JG), JG=1, 3)
               Read (ICI, 1030) CPOSI
           !
               Read (ICI, 1030) PO41TI, PO4T2I ! changed input file reading sequence here
               Read (ICI, 1030) NH41TI, NH4T2I
               Read (ICI, 1030) NO31TI, NO3T2I
               Read (ICI, 1030) CH41TI, CH4T2I
               Read (ICI, 1030) HS1TI, HST2I
               Read (ICI, 1030) SO41TI, SO4T2I
               Read (ICI, 1030) SI1TI, SIT2I
               Read (ICI, 1030) BENSTI, SODI !Input initial SOD
               If (BALGAE_CALC) Then
                  Read (ICI, 1030) BBMI !benthic algae initial condition
               End If

		       Read (ICI, 1032)
               If (SAV_CALC) Then
		          READ (ICI, *) (LEAFI(N), N=1, NSAVM)
		          READ (ICI, *) (STEMI(N), N=1, NSAVM)
                  READ (ICI, *) (ROOTI(N), N=1, NSAVM)
		          READ (ICI, *) (TUBERI(N), N=1, NSAVM)
		          READ (ICI, *) (EPI(N), N=1, NSAVM)
	           End If

               If (DFEEDER) READ (ICI, 1030) DFEEDI
               If (SFEEDER) READ (ICI, 1030) (SFEEDI(N), N=1, NSPECIES)
            End If
         End If
     !
     !***** Constituent concentrations

         Close (ICI)

     !******* RESTART constituent initial concentrations
     !
      Else If (RES_ICON_IN) Then
         Open (RSI, File=RSIFN, Status='OLD', Form='UNFORMATTED')
         Read (RSI) OLDTITLE
     !
         Allocate (RTMP(0:MGL, KBM1, NCP))
         RTMP = 0.0 !: note index from zero
         Read (RSI) (((RTMP(I, K, JCON), I=1, MGL), K=1, KBM1), JCON=1, &
        & NCP)
     !
         If (SERIAL) C1 = RTMP
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do JCON = 1, NCP
               Do K = 1, KBM1
                  Do I = 1, MLOC
                     C1 (I, K, JCON) = RTMP (NGID(I), K, JCON)
                  End Do
                  Do I = 1, NHN
                     C1 (I+MLOC, K, JCON) = RTMP (HN_LST(I), K, JCON)
                  End Do
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP)
     !
         Do JCON = 1, NAC
            Do K = 1, KBM1
               Do I = 1, MLOC
                  C1 (I, K, AC(JCON)) = Max (C1(I, K, AC(JCON)), 0.0)
                  C2 (I, K, AC(JCON)) = C1 (I, K, AC(JCON))
               End Do
            End Do
         End Do
     !
         If (SEDIMENT_CALC) Then
        !
        !
            Read (RSI) (CTEMP_GL(I), I=1, MGL), ((CPOP_GL(I, J), I=1, &
           & MGL), J=1, 3), ((CPON_GL(I, J), I=1, MGL), J=1, 3), &
           & ((CPOC_GL(I, J), I=1, MGL), J=1, 3), (CPOS_GL(I), I=1, &
           & MGL), (PO4T2TM1S_GL(I), I=1, MGL), (NH4T2TM1S_GL(I), I=1, &
           & MGL), (NO3T2TM1S_GL(I), I=1, MGL), (HST2TM1S_GL(I), I=1, &
           & MGL), (CH4T2TM1S_GL(I), I=1, MGL), (CH41TM1S_GL(I), I=1, &
           & MGL), (SO4T2TM1S_GL(I), I=1, MGL), (SIT2TM1S_GL(I), I=1, &
           & MGL), (BENSTRTM1S_GL(I), I=1, MGL), (SODTM1S_GL(I), I=1, &
           & MGL), (PO41TM1S_GL(I), I=1, MGL), (NH41TM1S_GL(I), I=1, &
           & MGL), (NO31TM1S_GL(I), I=1, MGL), (HS1TM1S_GL(I), I=1, &
           & MGL), (SI1TM1S_GL(I), I=1, MGL)
        !
            If (BALGAE_CALC) Then
               Read (RSI) (BBM_GL(I), I=1, MGL)
            End If
        !
            If (SERIAL) Then
           !
               CTEMP = CTEMP_GL
               CPOP = CPOP_GL
               CPON = CPON_GL
               CPOC = CPOC_GL
               CPOS = CPOS_GL
           !
               PO4T2TM1S = PO4T2TM1S_GL
               NH4T2TM1S = NH4T2TM1S_GL
               NO3T2TM1S = NO3T2TM1S_GL
               HST2TM1S = HST2TM1S_GL
               CH4T2TM1S = CH4T2TM1S_GL
               CH41TM1S = CH41TM1S_GL
               SO4T2TM1S = SO4T2TM1S_GL
               SIT2TM1S = SIT2TM1S_GL
               BENSTRTM1S = BENSTRTM1S_GL
               SODTM1S = SODTM1S_GL
           !
           !: save other sediment variables needed for restarts
               PO41TM1S = PO41TM1S_GL
               NH41TM1S = NH41TM1S_GL
               NO31TM1S = NO31TM1S_GL
               HS1TM1S = HS1TM1S_GL
               SI1TM1S = SI1TM1S_GL
           !
               If (BALGAE_CALC) Then
                  BBM = BBM_GL
               End If
            End If
        !
# if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do I = 1, MLOC
                  CTEMP (I) = CTEMP_GL (NGID(I))
                  CPOP (I, :) = CPOP_GL (NGID(I), :)
                  CPON (I, :) = CPON_GL (NGID(I), :)
                  CPOC (I, :) = CPOC_GL (NGID(I), :)
                  CPOS (I) = CPOS_GL (NGID(I))
                  PO4T2TM1S (I) = PO4T2TM1S_GL (NGID(I))
                  NH4T2TM1S (I) = NH4T2TM1S_GL (NGID(I))
                  NO3T2TM1S (I) = NO3T2TM1S_GL (NGID(I))
                  HST2TM1S (I) = HST2TM1S_GL (NGID(I))
                  CH4T2TM1S (I) = CH4T2TM1S_GL (NGID(I))
                  CH41TM1S (I) = CH41TM1S_GL (NGID(I))
                  SO4T2TM1S (I) = SO4T2TM1S_GL (NGID(I))
                  SIT2TM1S (I) = SIT2TM1S_GL (NGID(I))
                  BENSTRTM1S (I) = BENSTRTM1S_GL (NGID(I))
                  SODTM1S (I) = SODTM1S_GL (NGID(I))
              !
              !: save other sediment variables needed for restarts
                  PO41TM1S (I) = PO41TM1S_GL (NGID(I))
                  NH41TM1S (I) = NH41TM1S_GL (NGID(I))
                  NO31TM1S (I) = NO31TM1S_GL (NGID(I))
                  HS1TM1S (I) = HS1TM1S_GL (NGID(I))
                  SI1TM1S (I) = SI1TM1S_GL (NGID(I))
              !
                  If (BALGAE_CALC) Then
                     BBM (I) = BBM_GL (NGID(I))
                  End If
               End Do
               Do I = 1, NHN
                  CTEMP (I+MLOC) = CTEMP_GL (HN_LST(I))
                  CPOP (I+MLOC, :) = CPOP_GL (HN_LST(I), :)
                  CPON (I+MLOC, :) = CPON_GL (HN_LST(I), :)
                  CPOC (I+MLOC, :) = CPOC_GL (HN_LST(I), :)
                  CPOS (I+MLOC) = CPOS_GL (HN_LST(I))
                  PO4T2TM1S (I+MLOC) = PO4T2TM1S_GL (HN_LST(I))
                  NH4T2TM1S (I+MLOC) = NH4T2TM1S_GL (HN_LST(I))
                  NO3T2TM1S (I+MLOC) = NO3T2TM1S_GL (HN_LST(I))
                  HST2TM1S (I+MLOC) = HST2TM1S_GL (HN_LST(I))
                  CH4T2TM1S (I+MLOC) = CH4T2TM1S_GL (HN_LST(I))
                  CH41TM1S (I+MLOC) = CH41TM1S_GL (HN_LST(I))
                  SO4T2TM1S (I+MLOC) = SO4T2TM1S_GL (HN_LST(I))
                  SIT2TM1S (I+MLOC) = SIT2TM1S_GL (HN_LST(I))
                  BENSTRTM1S (I+MLOC) = BENSTRTM1S_GL (HN_LST(I))
                  SODTM1S (I+MLOC) = SODTM1S_GL (HN_LST(I))
              !: save other sediment variables needed for restarts
                  PO41TM1S (I+MLOC) = PO41TM1S_GL (HN_LST(I))
                  NH41TM1S (I+MLOC) = NH41TM1S_GL (HN_LST(I))
                  NO31TM1S (I+MLOC) = NO31TM1S_GL (HN_LST(I))
                  HS1TM1S (I+MLOC) = HS1TM1S_GL (HN_LST(I))
                  SI1TM1S (I+MLOC) = SI1TM1S_GL (HN_LST(I))
                  If (BALGAE_CALC) Then
                     BBM (I+MLOC) = BBM_GL (HN_LST(I))
                  End If
               End Do
            End If
# endif
        !
        !
         End If
     !
         If (DFEEDER) Then

            Read (RSI) (DFEEDM1S_GL(I), I=1, MGL)
            If (SERIAL) DFEEDM1S = DFEEDM1S_GL
#        if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do I = 1, MLOC
                  DFEEDM1S (I) = DFEEDM1S_GL (NGID(I))
               End Do
               Do I = 1, NHN
                  DFEEDM1S (I+MLOC) = DFEEDM1S_GL (HN_LST(I))
               End Do
            End If
#        endif
        !         DEALLOCATE(DFEEDM1S_GL)
         End If
     !
         If (SFEEDER) Then
        !                                                               !
            Read (RSI) ((SFEED_GL(I, K), I=1, MGL), K=1, NSPECIES)
            If (SERIAL) SFEED = SFEED_GL
#        if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do K = 1, NSPECIES
                  Do I = 1, MLOC
                     SFEED (I, K) = SFEED_GL (NGID(I), K)
                  End Do
                  Do I = 1, NHN
                     SFEED (I+MLOC, K) = SFEED_GL (HN_LST(I), K)
                  End Do
               End Do
            End If
#        endif

         End If
		 !
		 If (SAV_CALC) Then
             Read (RSI) ((LEAF_GL(I, K), I=1, MGL), K=1, NSAVM), &
		   &((STEM_GL(I, K), I=1, MGL), K=1, NSAVM), &
		   &((ROOT_GL(I, K), I=1, MGL), K=1, NSAVM), &
		   &((TUBER_GL(I, K), I=1, MGL), K=1, NSAVM), &
		   &((EP_GL(I, K), I=1, MGL), K=1, NSAVM)
            If (SERIAL) Then
			  LEAF = LEAF_GL
			  STEM = STEM_GL
			  ROOT = ROOT_GL
			  TUBER = TUBER_GL
			  EP = EP_GL
			 End If
#        if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do K = 1, NSAVM
                  Do I = 1, MLOC
                     LEAF (I, K) = LEAF_GL (NGID(I), K)
					 STEM (I, K) = STEM_GL (NGID(I), K)
					 ROOT (I, K) = ROOT_GL (NGID(I), K)
					 TUBER (I, K) = TUBER_GL (NGID(I), K)
					 EP (I, K) = EP_GL (NGID(I), K)
                  End Do
                  Do I = 1, NHN
                     LEAF (I+MLOC, K) = LEAF_GL (HN_LST(I), K)
					 STEM (I+MLOC, K) = STEM_GL (HN_LST(I), K)
					 ROOT (I+MLOC, K) = ROOT_GL (HN_LST(I), K)
					 TUBER (I+MLOC, K) = TUBER_GL (HN_LST(I), K)
					 EP (I+MLOC, K) = EP_GL (HN_LST(I), K)
                  End Do
               End Do
            End If
#        endif
         End If
         Close (RSI)

     !
      Else
         Write (*,*) 'initial conditions file specified incorrectly'
         Stop
      End If
  !
  !***** Mineralization rates
  !
      Open (MRL, File=MRLFN, Status='OLD')
      Read (MRL, 1032)
  !
  !***** Spatially-invariant kinetics parameters
  !
      Read (MRL, 1080) KHONT, KHNNT, KHOCOD, KHODOC, KHNDN
      Read (MRL, 1080) AOCR, AONT, KHCOAG
      Read (MRL, 1080) TRCOD, TRMNL, TRHDR, TRSUA
      Read (MRL, 1080) KTCOD, KTMNL, KTHDR, KTSUA
      Read (MRL, 1080) KTNT1, KTNT2, TMNT
      Read (MRL, 1080) KADPO4, KADSA
      Read (MRL, 1080) AANOX, ANDC 
      Read (MRL, 1080) AREAR, BREAR, CREAR
  !
  !***** Spatially-varying kinetics parameters
  !
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KLDC (1, 1)
         Do K = 1, KBM1

            Do I = 1, MTLOC
               KLDC (I, K) = KLDC (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KLDC = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KLDC (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KLDC (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KRDC (1, 1)
         Do K = 1, KBM1
            Do I = 1, MLOC
               KRDC (I, K) = KRDC (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KRDC = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KRDC (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KRDC (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KLPC (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KLPC (I, K) = KLPC (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KLPC = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KLPC (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KLPC (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KRPC (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KRPC (I, K) = KRPC (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KRPC = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KRPC (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KRPC (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KLDN (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KLDN (I, K) = KLDN (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KLDN = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KLDN (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KLDN (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KRDN (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KRDN (I, K) = KRDN (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KRDN = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KRDN (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KRDN (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KLPN (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KLPN (I, K) = KLPN (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KLPN = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KLPN (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KLPN (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KRPN (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KRPN (I, K) = KRPN (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KRPN = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KRPN (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KRPN (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KLDP (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KLDP (I, K) = KLDP (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KLDP = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KLDP (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KLDP (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KRDP (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KRDP (I, K) = KRDP (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KRDP = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KRDP (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KRDP (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KLPP (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KLPP (I, K) = KLPP (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KLPP = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KLPP (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KLPP (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KRPP (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KRPP (I, K) = KRPP (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KRPP = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KRPP (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KRPP (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KSUA (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KSUA (I, K) = KSUA (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KSUA = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KSUA (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KSUA (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KCOD (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KCOD (I, K) = KCOD (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KCOD = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KCOD (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KCOD (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KDCALG (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KDCALG (I, K) = KDCALG (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KDCALG = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KDCALG (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KDCALG (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KLCALG (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KLCALG (I, K) = KLCALG (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KLCALG = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KLCALG (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KLCALG (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM .Eq. 'CONSTANT') Then
         Read (MRL, 1033) KRCOAG (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KRCOAG (I, K) = KRCOAG (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KRCOAG = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KRCOAG (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KRCOAG (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KDNALG (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KDNALG (I, K) = KDNALG (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KDNALG = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KDNALG (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KDNALG (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KLNALG (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KLNALG (I, K) = KLNALG (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KLNALG = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KLNALG (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KLNALG (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KDPALG (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KDPALG (I, K) = KDPALG (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KDPALG = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KDPALG (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KDPALG (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) KLPALG (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               KLPALG (I, K) = KLPALG (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) KLPALG = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  KLPALG (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  KLPALG (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Read (MRL, 1060) SPVARM, PRINTM
      If (SPVARM == 'CONSTANT') Then
         Read (MRL, 1033) NTM (1, 1)
         Do K = 1, KBM1
            Do I = 1, MTLOC
               NTM (I, K) = NTM (1, 1)
            End Do
         End Do
      Else
         Allocate (RTMP2(MGL, KBM1))
         RTMP2 = 0.0
         Do K = 1, KBM1
            Read (MRL, 1033) (RTMP2(I, K), I=1, MGL)
         End Do
         If (SERIAL) NTM = RTMP2
#      if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  NTM (I, K) = RTMP2 (NGID(I), K)
               End Do
               Do I = 1, NHN
                  NTM (I+MLOC, K) = RTMP2 (HN_LST(I), K)
               End Do
            End Do
         End If
#      endif
         Deallocate (RTMP2)
      End If
      Close (MRL)
  !
  !***** Boundary concentrations
  !
      If (BOUNDARY_CONC) Then
         Read (CBC, 1100)
         Read (CBC, 1020) (NCB(JCON), JCON=1, NCP)
         Read (CBC, 1100)
      End If
  !
  !***** Sources One

      If (SOURCE_ONE) Then !
         Read (S1, 1100)
         Read (S1, 1020) (S1LN(JCON), JCON=1, NCP)
         Do JCON = 1, NCP
            Read (S1, 1020) (S1LB(JP, JCON), JP=1, S1LN(JCON))
         End Do
         Read (S1, 1100)
      End If
  !
  !***** Sources Two
  !
      If (SOURCE_TWO) Then !
         Read (S2, 1100)
         Read (S2, 1020) (S2LN(JCON), JCON=1, NCP)
         Do JCON = 1, NCP
            Read (S2, 1020) (S2LB(JNP, JCON), JNP=1, S2LN(JCON))
         End Do
         Read (S2, 1100)
      End If
  !
  !***** Sources Three
  !
      If (SOURCE_THR) Then
         Read (S3, 1100)
         Read (S3, 1020) (S3LN(JCON), JCON=1, NCP)
         Do JCON = 1, NCP
            Read (S3, 1020) (S3LB(JNP, JCON), JNP=1, S3LN(JCON))
         End Do
         Read (S3, 1100)
      End If
  !
  !***** Light extinction and algal growth
  !
      If (LIGHT_EXTINCTION) Then
         Open (KEI, File=KEIFN, Status='OLD')
         Call OWQ_ALLOC
         Call OWQ_READ
         Close (KEI)
      End If
  !
  !***** Algae
  !
      If (ALGAE_CALC) Then ! 
         Call ALG_READ
      End If
  !
  !***** Zooplankton
  !
      If (ZOO_CALC) Call ZOOP_READ
  !
  !***** Settling rates
  !
      If (SETTLING) Then
         Read (STL, 1032)
         Read (STL, 1060) SPVARM, PRINTM
         If (SPVARM == 'CONSTANT') Then
            Read (STL, 1033) WSS (1, 1), WSL (1, 1), WSR (1, 1), WS1 &
           & (1, 1), WS2 (1, 1), WS3 (1, 1), WSU (1, 1)
            Do K = 1, KBM1
               Do I = 1, MLOC
                  WSS (I, K) = WSS (1, 1)
                  WSL (I, K) = WSL (1, 1)
                  WSR (I, K) = WSR (1, 1)
                  WS1 (I, K) = WS1 (1, 1)
                  WS2 (I, K) = WS2 (1, 1)
                  WS3 (I, K) = WS3 (1, 1)
                  WSU (I, K) = WSU (1, 1)
               End Do
            End Do
         Else
            Allocate (RTMP2(0:MGL, KBM1))
            RTMP2 = 0.0 !: note index from zero because WSS index is from zero
            Allocate (RTMP21(0:MGL, KBM1))
            RTMP21 = 0.0 !: note index from zero because WSL index is from zero
            Allocate (RTMP22(0:MGL, KBM1))
            RTMP22 = 0.0 !: note index from zero because WSR index is from zero
            Allocate (RTMP23(0:MGL, KBM1))
            RTMP23 = 0.0 !: note index from zero because WS1 index is from zero
            Allocate (RTMP24(0:MGL, KBM1))
            RTMP24 = 0.0 !: note index from zero because WS2 index is from zero
            Allocate (RTMP25(0:MGL, KBM1))
            RTMP25 = 0.0 !: note index from zero becasue WS3 index is from zero
            Allocate (RTMP26(0:MGL, KBM1))
            RTMP26 = 0.0 !: note index from zero because WSU index is from zero
            Do K = 1, KBM1
               Do I = 1, MGL
!!! commented the above and fixed with the line below:
                  Read (STL, 1085) RTMP2 (I, K), RTMP21 (I, K), RTMP22 &
                 & (I, K), RTMP23 (I, K), RTMP24 (I, K), RTMP25 (I, K), &
                 & RTMP26 (I, K)
               End Do
            End Do
            If (SERIAL) Then
               WSS = RTMP2
               WSL = RTMP21
               WSR = RTMP22
               WS1 = RTMP23
               WS2 = RTMP24
               WS3 = RTMP25
               WSU = RTMP26
            End If
#        if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do K = 1, KBM1
                  Do I = 1, MLOC
                     WSS (I, K) = RTMP2 (NGID(I), K)
                     WSL (I, K) = RTMP21 (NGID(I), K)
                     WSR (I, K) = RTMP22 (NGID(I), K)
                     WS1 (I, K) = RTMP23 (NGID(I), K)
                     WS2 (I, K) = RTMP24 (NGID(I), K)
                     WS3 (I, K) = RTMP25 (NGID(I), K)
                     WSU (I, K) = RTMP26 (NGID(I), K)
                  End Do
                  Do I = 1, NHN
                     WSS (I+MLOC, K) = RTMP2 (HN_LST(I), K)
                     WSL (I+MLOC, K) = RTMP21 (HN_LST(I), K)
                     WSR (I+MLOC, K) = RTMP22 (HN_LST(I), K)
                     WS1 (I+MLOC, K) = RTMP23 (HN_LST(I), K)
                     WS2 (I+MLOC, K) = RTMP24 (HN_LST(I), K)
                     WS3 (I+MLOC, K) = RTMP25 (HN_LST(I), K)
                     WSU (I+MLOC, K) = RTMP26 (HN_LST(I), K)
                  End Do
               End Do
            End If
#        endif
            Deallocate (RTMP2, RTMP21, RTMP22, RTMP23, RTMP24, RTMP25, &
           & RTMP26)
         End If
         Close (STL)
      End If
  !
  !***** Atmospheric loads
  !
      If (ATMOS_LOADS) READ (ATM, 1000)
  !
  !***** Submerged aquatic vegetation
  !
  !***** Benthic fluxes
      If (BENTHIC_FLUXES) Then
     !       WRITE(*,*)'reading data from BFI file'
         Read (BFI, 11040)!skip 3 title lines
     !(the format skips two lines, the read itself skips one line)
         Read (BFI, 1080) KSDOC, KSNH4, KSNO3, KSPO4, KSO, KSSA !skip 2 lines and read
     !       WRITE(*,11051) KSDOC, KSNH4, KSNO3, KSPO4, KSO, KSSA
         Read (BFI, 1080) TRSDOC, TRSNH4, TRSNO3, TRSPO4, TRSO, TRSSA !skip 2 lines and read
     !       WRITE(*,11051) TRSDOC,TRSNH4,TRSNO3,TRSPO4,TRSO,TRSSA
         Read (BFI, 1080) MTCNO3, SEDNO3, KHSO !skip 2 lines and read
     !       WRITE(*,11051) MTCNO3, SEDNO3, KHSO
         Read (BFI, 1100)!skip 2 lines
     !note here 1100 format skips one line
     !but the read itself skips one line
      End If
  !***** Meteorologic data
  !
      Read (MET, 1000)
  !
  !***** atmospheric pCO2  !
      If (CARBONATE_CALC) Then
     !Get constants K0, FugFac needed in PCO2AIR_READ_INITIAL to turn pCO2atm to CO2star_sat
     ! and also CO2star_surf  -     Needs input file inputs/input_K1_K2_KSO4.dat
         TDIC = C1 (:, :, 33)
         TALK = C1 (:, :, 34)
         Call CALC_PH_PCO2
     !Read atmospheric pCO2 or fCO2 data (pCO2atm) from inputs/input_pCO2atm.dat
         Call PCO2AIR_READ_INITIAL
      End If
!
  !***** Input FORMAT statements
1000  Format (/ //)
1001  Format (I6, 10(I4, F10.1))
1002  Format (7 F10.3)
1003  Format (8 F10.3)
1005  Format (/ / F10.1)
1010  Format (A72)
1020  Format (/ / (8 X, 9 I8))
1025  Format (8 X, 9 I8)
1027  Format (/ / (11 X, 8 I8))
1029  Format (10(:F12.4, 1 X))
1030  Format (/ / (8 X, 9 F8.0))
1031  Format (/ / (8 X, 2 F8.0, I8, (5X,A3)))
1032  Format (/)
1033  Format (/ / : (8 X, 9 F8.0))
1035  Format (:///10(10 F8.0: /))
  !1036 FORMAT(//8X,2F8.1,I8,3F8.1)
  !1037 FORMAT(//8X,8F8.1)
  !1038 FORMAT(//8X,8F8.1)
1040  Format (/ / 8 X, I8, 8 F8.0)
1050  Format (/ / 13 X, A3, I8, 5 X, A3)
1052  Format (/ / 13 X, A3, 9 I8)
1055  Format (/ / 8 X, 4 (5 X, A3), I8)
1060  Format (/ / 8 X, 2 A8, 2 F8.0)
1065  Format (://(8 X, 9 A8))
1070  Format (/ / (8 X, 9(5 X, A3)))
1071  Format (/ / (8 X, 5(5 X, A3), A8, 5 X, A3, 5 X, A3, 5 X, A3,5 X, A3))
1072  Format (/ / (8 X, 6(5 X, A3)))
1080  Format (://8 X, 9 F8.0)
  !1082 FORMAT(//8X,2A8//)
1085  Format (8 X, 9 F8.0)
1087  Format (8 X, 4 F8.0)
1090  Format (/ / (8 X, A72))
1100  Format (/)
1110  Format (:////////(8 X, 5 I8))
1120  Format (/ // / (8 X, I8))
1130  Format (8 X, 9 I8)
1140  Format (://(2 I8))
1150  Format (5 X, 3 F15.0, F18.0, F12.0, I10)
1160  Format (://(13 X, F13.0))
1170  Format (://(2 I10))
1190  Format (I8, 9 F8.2)
1200  Format (/ (10 F12.0))
1220  Format (/ / (:8 X, 6 F8.0))
1230  Format (8 X, 1 F8.0)
11103 Format (10 F8.0)
11040 Format (/ /)
11051 Format (://8 X, 9 F8.4)!go down 2 lines and write, but do not
  
  !
  !***** Time variables
  !
      NIT = 0
  !     JDAY     = TMSTRT
      JDAYMBL = TMSTRT
      DLTDP = 1
      TFLDP = 1
      KFLDP = 1
      RSODP = 1
      SNPDP = 1
      PLTDP = 1
      APLDP = 1
      OPLDP = 1
      MBLDP = 1
      DIADP = 1
      NWQMR = 0
      DLT = DLTVAL (DLTDP)
	  DTI = DLT			!, use DLT for DTI all the time
      ELTMSPLT = JDAY * 86400.
      ELTMSTFL = JDAY * 86400.
      ELTMSKFL = JDAY * 86400.
  !
      NHMR = 0
      NXTVD = JDAY
      ELTMS = JDAY * 86400.

      MXDLT = DLTMAX (DLTDP)
      FNDLT = DLTFTN (DLTDP)
      NXSNP = SNPD (SNPDP)
      NXPLT = PLTD (PLTDP)
      NXTFL = TFLD (TFLDP)
      NXKFL = KFLD (KFLDP)
      NXOPL = OPLD (OPLDP)
      NXMBL = MBLD (MBLDP)
      NXDIA = DIAD (DIADP)
      NXAPL = APLTD (APLDP)
  !
  !***** Output control variables
  !
      SNPD (NSNP+1) = TMEND + 1.
      PLTD (NPLT+1) = TMEND + 1.
      OPLD (NOPL+1) = TMEND + 1.
      MBLD (NMBL+1) = TMEND + 1.
      TFLD (NTFL+1) = TMEND + 1.
      KFLD (NKFL+1) = TMEND + 1.
      RSOD (NRSO+1) = TMEND + 1.
      DLTD (NDLT+1) = TMEND + 1.
      DIAD (NDIA+1) = TMEND + 1.
      APLTD (NAPL+1) = TMEND + 1.
  !
      NXCBC = 0.
  ! 
      NXMET = 0.
      NXS1 = 0.
      NXS2 = 0.
      NXS3 = 0.
      NXBFI = 0.
      NXATM = 0.
      NXSAV = 0.
  !
      AC1 = 0.0
  !
      AKE = 0.0
      ACCHL1 = 0.0
      ACCHL2 = 0.0
      ACCHL3 = 0.0
      AFI1 = 0.0
      ANL1 = 0.0
      APL1 = 0.0
      ASL1 = 0.0
      AFI2 = 0.0
      ANL2 = 0.0
      APL2 = 0.0
      ASL2 = 0.0
      AFI3 = 0.0
      ANL3 = 0.0
      APL3 = 0.0
      ASL3 = 0.0
      ANPP = 0.0
      AGPP = 0.0
      ARESP = 0.0
  !
      AASRAT = 0.0
      ACFIX = 0.0
  !
  !Accumulation of stuff for diagonostics output
      ACLSZ = 0.0
      ACLLZ = 0.0
      ARSZ = 0.0
      ARLZ = 0.0
      ABMSZ = 0.0
      ABMLZ = 0.0
      AMSZ = 0.0
      AMLZ = 0.0
      APRSZLZ = 0.0
      AGSZ = 0.0
      AGLZ = 0.0
      ADOCSZ = 0.0
      APOCSZ = 0.0
      ADOCLZ = 0.0
      APOCLZ = 0.0
      ANH4SZ = 0.0
      ADONSZ = 0.0
      APONSZ = 0.0
      ANH4LZ = 0.0
      ADONLZ = 0.0
      APONLZ = 0.0
      APO4SZ = 0.0
      ADOPSZ = 0.0
      APOPSZ = 0.0
      APO4LZ = 0.0
      ADOPLZ = 0.0
      APOPLZ = 0.0
      APRSZ = 0.0
      APRLZ = 0.0
  !
      AB1SZ = 0.0
      AB2SZ = 0.0
      AB3SZ = 0.0
      AB1LZ = 0.0
      AB2LZ = 0.0
      AB3LZ = 0.0
      ADOSZ = 0.0
      ADOLZ = 0.0
      ASASZ = 0.0
      ASUSZ = 0.0
      ASALZ = 0.0
      ASULZ = 0.0
  !
      ACPOS = 0.0
      ACPIP = 0.0
      ASSFWS = 0.0
      APCFWS = 0.0
      APNFWS = 0.0
      APPFWS = 0.0
      APSFWS = 0.0
      ABENDO = 0.0
      ABENSA = 0.0
      ABENDOC = 0.0
      ABENNH4 = 0.0
      ABENNO3 = 0.0
      ABENPO4 = 0.0
      ABENCOD = 0.0
      ABENCH4G = 0.0
      ABENCH4A = 0.0
      AFIB = 0.0
      ANLB = 0.0
      APLB = 0.0
      ANPPB = 0.0
      ABBM = 0.0
      ABLITE = 0.0
      ADFEED = 0.0
  !
      ASFEED = 0.0
  !
      AJNSF = 0.0
      AJPSF = 0.0
      ASODSF = 0.0
      ASASF = 0.0
      ASUSF = 0.0
      ASFGCIN = 0.0
      ASFCFEC = 0.0
      ASFCPSF = 0.0
      AFLXCSF = 0.0
      AFLXNSF = 0.0
      AFLXPSF = 0.0
      ARPOCSF = 0.0
      ARPONSF = 0.0
      ARPOPSF = 0.0
      ASSISF = 0.0
      ASSISASF = 0.0
      ASSISUSF = 0.0
      ASSIPSF = 0.0
  !
      ACPOC = 0.0
      ACPON = 0.0
      ACPOP = 0.0
  !
  !
      A_T = 0.0
      AP1 = 0.0
      ABM1 = 0.0
      APR1 = 0.0
      AP2 = 0.0
      ABM2 = 0.0
      APR2 = 0.0
      AP3 = 0.0
      ABM3 = 0.0
      APR3 = 0.0
      AALGDOC = 0.0
      AALGPOC = 0.0
      ADENIT = 0.0
      AMNLDOC = 0.0
      AHDRPOC = 0.0
      AALGNH4 = 0.0
      AALGNO3 = 0.0
      AALGDON = 0.0
      AALGPON = 0.0
      ANT = 0.0
      ANFIX = 0.0
      ADENNO3 = 0.0
      AMNLDON = 0.0
      AHDRPON = 0.0
      AALGPO4 = 0.0
      AALGDOP = 0.0
      AALGPOP = 0.0
      AMNLDOP = 0.0
      AHDRPOP = 0.0
      APSD = 0.0
      ASAP = 0.0
      AALGUP = 0.0
      AALGRES = 0.0
      ADO = 0.0
      ADORALG = 0.0
      ADOPR = 0.0
      ADCOD = 0.0
      ADDOC = 0.0
      ANITRIF = 0.0
  !
  !
      ABMB = 0.0
      APB = 0.0
      APRB = 0.0
      ABADOC = 0.0
      ABAPOC = 0.0
      ABANH4 = 0.0
      ABANO3 = 0.0
      ABAPON = 0.0
      ABAPO4 = 0.0
      ABAPOP = 0.0
      ABADO = 0.0
  !
     
  !***** Initial water column mass
  !
      If (MASS_BALANCE) Then
         Do JCON = 1, NAC
            CMASS (AC(JCON)) = 0.0
            Do K = 1, KBM1
               Do I = 1, MTLOC
                  CMASS (AC(JCON)) = CMASS (AC(JCON)) + C1 (I, K, &
                 & AC(JCON)) * ART1 (I) * D (I) * DZ2D (I,K) / 1000.0
               End Do
            End Do
         End Do
         IWCMN = ANC1 * CMASS (4) + ANC2 * CMASS (5) + ANC3 * CMASS (6) &
        & + CMASS (12) + ANCSZ * CMASS (7) + ANCLZ * CMASS (8) + CMASS &
        & (13) + CMASS (14) + CMASS (15) + CMASS (16) + CMASS (17) + &
        & CMASS (18) + CMASS (19)
         IWCMP = CMASS (20) + CMASS (21) + CMASS (22) + CMASS (23) + &
        & CMASS (24) + CMASS (25) + APCSZ * CMASS (7) + APCLZ * CMASS &
        & (8) + CMASS (30) + CMASS (31) + CMASS (32)
         IWCMC = CMASS (4) + CMASS (5) + CMASS (6) + CMASS (7) + CMASS &
        & (8) + CMASS (9) + CMASS (10) + CMASS (11) + CMASS (12)
         IWCMS = ASC1 * CMASS (4) * ASC2 * CMASS (5) + ASC3 * CMASS (6) &
        & + CMASS (28) + CMASS (29)
      End If
  !
  !***** Change in concentrations
  !
      DTC = 0.0
  !
  !
  !    *** DO days
  !
      DOVDAYS = 0.0
  !
  !***** Maximum number of Sources One through Three
  !
      Do JCON = 1, NCP
         If (S1LN(JCON) > S1LNMAX) Then
            S1LNMAX = S1LN (JCON)
            JCS1MAX = JCON
         End If
         If (S2LN(JCON) > S2LNMAX) Then
            S2LNMAX = S2LN (JCON)
            JCS2MAX = JCON
         End If
         If (S3LN(JCON) > S3LNMAX) Then
            S3LNMAX = S3LN (JCON)
            JCS3MAX = JCON
         End If
      End Do
  !
  !***** Horizontal diffusion
  !
      VISCOFH = 0.0
      If (XY_DIFFUSION) Then
         Do K = 1, KBM1
            Do I = 1, MTLOC ! 9may2016 changed from MLOC
               VISCOFH (I, K) = XYDF
            End Do
         End Do
      End If
  !
  !
  !************************************************************************
  !**                   Parameters for Sediment Submodels                **
  !************************************************************************
  !
      If (SEDIMENT_CALC) Then
         If ( .Not. UNI_ICON_IN_SED_VAR) Then !: if UNI_ICON_IN_SED_VAR is true, SED_INIT called earlier
            Call SED_INIT ()!  Setup sediment diagenesis module (including SED_ALLOC and SED_INIT_ICI)
			If (SAV_CALC) Call SAV_INIT_ICI() ! AN
         End If
         Call SED_READ !
         Call SED_INIT2 ()!  initialize the values from input read by INPUTS() above
     !
      End If
  !
  !************************************************************************
  !**                              Outputs                               **
  !************************************************************************
  !
  !***** Open output files
  !
  !  - only output on MSR
      If (MSR) Then
# if defined KILLTREES
         If (SNAPSHOTS) OPEN (SNP, FILE=SNPFN)
         If (DIAGNOSTICS) OPEN (DIA, FILE=DIAFN)
         If (KINETIC_FLUXES) OPEN (KFL, FILE=KFLFN, FORM='UNFORMATTED')!MNOEL
         If (MASS_BALANCE) OPEN (MBL, FILE=MBLFN, FORM='UNFORMATTED')
         If (PLOTS) OPEN (PLT, FILE=PLTFN)
         If (AVERAGE_PLOTS) OPEN (APL, FILE=APLFN, FORM='UNFORMATTED')
         If (TRANSPORT_FLUXES) OPEN (TFL, FILE=TFLFN, FORM='UNFORMATTED&
        &')
         If (OXYGEN_PLOTS) OPEN (OPL, FILE=OPLFN, FORM='UNFORMATTED')
         If (BENTHIC_FLUXES .And. BENTHIC_OUTPUT) OPEN (BFO, &
        & FILE=BFOFN)
         If (SFEEDER) OPEN (SFO, FILE=SFOFN, FORM='UNFORMATTED')
         If (DFEEDER .And. DFEEDER_OUTPUT) OPEN (DFO, FILE=DFOFN)
# endif
      End If
  !
  !***** Snapshots
  !
# if defined KILLTREES
  !      If (SNAPSHOTS) Then
  !         If (MSR) WRITE (SNP, 2000) 'CE-QUAL-IC - Version 5.0 w LRDOM, &
  !        &Internal P'
  !         If (MSR) WRITE (SNP, 2010) TITLE
  !         If (MSR) WRITE (SNP, 2020) MAPFN, GEOFN, ICIFN, AGRFN, ZOOFN, &
  !        & SFIFN
  !         If (MSR) WRITE (SNP, 2027) STLFN, MRLFN, KEIFN
  !         If (MSR) WRITE (SNP, 2021) (HYDFN(JF), JF=1, NHYDF)
  !         If (MSR) WRITE (SNP, 2022) (METFN(JF), JF=1, NTVDF)
  !         If (MSR) WRITE (SNP, 2023) (S1FN(JF), JF=1, NTVDF)
  !         If (MSR) WRITE (SNP, 2024) (S2FN(JF), JF=1, NTVDF)
  !         If (MSR) WRITE (SNP, 2029) (S3FN(JF), JF=1, NTVDF)
  !         If (MSR) WRITE (SNP, 2031) (ATMFN(JF), JF=1, NTVDF)
  !         If (MSR) WRITE (SNP, 2028) (SVIFN(JF), JF=1, NTVDF)
  !         If (MSR) WRITE (SNP, 2025) (CBCFN(JF), JF=1, NTVDF)! JLM 2/23/01
  !         If (MSR) WRITE (SNP, 2026) (BFIFN(JF), JF=1, NTVDF)
  !         If (MSR) WRITE (SNP, 2030) ICOFN, SNPFN, RSOFN, PLTFN, APLFN, &
  !        & DIAFN, TFLFN, KFLFN, OPLFN, MBLFN, ALOFN, ZFOFN, BFOFN, &
  !        & SVOFN, SFOFN, BAOFN, DFOFN
  !         If (MSR) WRITE (SNP, 2040) NB, NSB, NQF, NL
  !         If (MSR) WRITE (SNP, 2050) TMSTRT, TMEND
  !         If (MSR) WRITE (SNP, 2060) (DLTD(J), J=1, NDLT)
  !         If (MSR) WRITE (SNP, 2070) (Int(DLTVAL(J)), J=1, NDLT)
  !         If (MSR) WRITE (SNP, 2071) AUTOC
  !         If (MSR) WRITE (SNP, 2072) (Int(DLTMAX(J)), J=1, NDLT)
  !         If (MSR) WRITE (SNP, 2073) (DLTFTN(J), J=1, NDLT)
  !         If (MSR) WRITE (SNP, 2080) RSIC, BCC, S1C, S2C, S3C, BFC, MDC, &
  !        & ICIC, ATMC, SAVLC, STLC
  !         If (MSR) WRITE (SNP, 2090) ICOC, SNPC, PLTC, QPLTC, SPLTC, &
  !        & APLTC, DIAC, TFLC, KFLC, OPLC, RSOC
  !         If (MSR) WRITE (SNP, 2100) (SNPD(K), K=1, NSNP)
  !         If (MSR) WRITE (SNP, 2110) (PLTD(J), J=1, NPLT)
  !         If (MSR) WRITE (SNP, 2115) (APLTD(J), J=1, NAPL)
  !         If (MSR) WRITE (SNP, 2116) (TFLD(J), J=1, NTFL)
  !         If (MSR) WRITE (SNP, 2118) (OPLD(J), J=1, NOPL)
  !         If (MSR) WRITE (SNP, 2120) (RSOD(J), J=1, NRSO)
  !         If (MSR) WRITE (SNP, 2130) (SNPF(K), K=1, NSNP)
  !         If (MSR) WRITE (SNP, 2135) (PLTF(J), J=1, NPLT)
  !         If (MSR) WRITE (SNP, 2136) (APLF(J), J=1, NAPL)
  !         If (MSR) WRITE (SNP, 2137) (TFLF(J), J=1, NTFL)
  !         If (MSR) WRITE (SNP, 2138) (OPLF(J), J=1, NOPL)
  !         If (MSR) WRITE (SNP, 2141) MBLC, VBC
  !         If (MSR) WRITE (SNP, 2150) HYDC, Int (AHMDLT), SLC, CONSC, TH, &
  !        & MINSTEP, FILGTH
  !         If (MSR) WRITE (SNP, 2152) SEDC
  !         If (MSR) WRITE (SNP, 2154) SAVMC
  !         If (MSR) WRITE (SNP, 2160) BNDTC
  !         If (MSR) WRITE (SNP, 2165) XYDF, ZDFMUL
  !         If (MSR) WRITE (SNP, 2170) FLC, XYDFC, ZDFC
  !         If (MSR) WRITE (SNP, 2180) (CNAME(JCON), ACC(JCON), CIC(JCON), &
  !        & JCON=1, NCP)
  !         If (MSR) WRITE (SNP, 2190) REDS1N, REDS2N, REDS3N, REDCBN, &
  !        & REDS1P, REDS2P, REDS3P, REDCBP, REDS1C, REDS2C, REDS3C, &
  !        & REDCBC
  !         If (MSR) WRITE (SNP, 2200) KHONT, KHNNT, KHODOC, KHOCOD, KHNDN
  !         If (MSR) WRITE (SNP, 2210) AOCR, AONT
  !         If (MSR) WRITE (SNP, 2260) KTNT1, KTNT2, TMNT
  !         If (MSR) WRITE (SNP, 2280) TRCOD, KTCOD, TRMNL, KTMNL, TRHDR, &
  !        & KTHDR, TRHDR, KTHDR, TRMNL, KTMNL, TRHDR, KTHDR, TRHDR, &
  !        & KTHDR, TRMNL, KTMNL, TRHDR, KTHDR, TRHDR, KTHDR, TRSUA, KTSUA
  !         If (MSR) WRITE (SNP, 2317) KADPO4, KADSA
  !         If (MSR) WRITE (SNP, 2330) AANOX, ANDC, KHCOAG
  !         If (MSR) WRITE (SNP, 2332) AREAR, BREAR, CREAR
  !         If (MSR) WRITE (SNP, 2310) KE (1, 1), KEISS (1, 1), KEVSS (1, &
  !        & 1), KEDOC (1, 1)
  !     End If
!
  !***** Benthic flux parameters
!
  !      If (BENTHIC_FLUXES .And. BENTHIC_OUTPUT) Then
  !         If (MSR) WRITE (BFO, 2501) KSDOC, KSNH4, KSNO3, KSPO4, KSO, &
  !        & KSSA
  !         If (MSR) WRITE (BFO, 2503) TRSDOC, TRSNH4, TRSNO3, TRSPO4, &
  !        & TRSO, TRSSA
  !         If (MSR) WRITE (BFO, 2502) MTCNO3, SEDNO3, KHSO
  !      End If
!
  !
  !***** Diagnostics
  !
  !      If (DIAGNOSTICS) Then !MNOEL 2-5-93
  !         If (MSR) WRITE (DIA, 2010) TITLE
  !         If (UNI_ICON_IN) Then
  !            If (MSR) WRITE (DIA, 2012)
  !            If (MSR) WRITE (DIA, 2011) ICIC
  !         End If
  !         If (BIN_ICON_IN) Then
  !            If (MSR) WRITE (DIA, 2012)
  !            If (MSR) WRITE (DIA, 2011) ICIC
  !            If (MSR) WRITE (DIA, 2013)
  !            If (MSR) WRITE (DIA, 2010) OLDTITLE
  !         End If
  !      End If
!
!
  !***** Plots
  !
  !      If (MSR .And. PLOTS) WRITE (PLT,*) TITLE, NAC, AC, NB, NSPECIES, &
  !     & QUALITY_DIAG, SEDIMENT_DIAG, SAV_PLOTS, ANC1, ANC2, ANC3, ASC1, &
  !     & ASC2, ASC3, ANCSZ, ANCLZ, APCSZ, APCLZ, KADPO4, KADSA, M2, &
  !     & PIE2HS, PIE2PO4, PIE2SI, PIENH4, ADWCEPI !CFC
  !      If (MSR) WRITE (PLT,*) NSB,(ART1(I), I=1, MTLOC)!Replaced SFA by ART1
!
  !***** Average plots
!
  !      If (AVERAGE_PLOTS) Then
  !         If (MSR) WRITE (APL) TITLE, NAC, AC, NB, NSPECIES, &
  !        & QUALITY_DIAG, SEDIMENT_DIAG, SAV_PLOTS, ANC1, ANC2, ANC3, &
  !        & ASC1, ASC2, ASC3, ANCSZ, ANCLZ, APCSZ, APCLZ, KADPO4, KADSA, &
  !        & ADWCEPI
  !         If (MSR) WRITE (APL) NSB, (ART1(I), I=1, MTLOC)!Replaced SFA by ART1
  !      End If
!
  !***** Constituent fluxes
!
  !JQI     IF (TRANSPORT_FLUXES) WRITE (TFL) TITLE,NQF
!
  !***** Kinetics fluxes
!
  !      If (MSR .And. KINETIC_FLUXES) WRITE (KFL) TITLE, MTLOC, KBM1, &
  !     & (SBN(I), I=1, MTLOC), (BBN(I), I=1, MTLOC), &
  !	 & (ART1(I), I=1, MTLOC), SAV_CALC, &
  !     & BALGAE_CALC !: replaced SFA by ART1
!
  !***** Oxygen plots
!
  !      If (OXYGEN_PLOTS) Then
  !         If (MSR) WRITE (OPL) TITLE, MTLOC, KBM1, NOINT, (OINT(I), I=1, &
  !        & NOINT)
  !      End If
!
  !***** Mass balances
!
  !      If (MASS_BALANCE) Then
  !         If (MSR) WRITE (MBL) NSB, NB, (SBN(B), B=1, MTLOC), (BBN(B), &
  !        & B=1, MTLOC)
  !      End If
!
#endif
  !
  !***** Output FORMAT statements
  !
2000  Format (1 X, A24/'+', 10('_')//)
2010  Format (1 X, A72)
2012  Format (/ 1 X, 79 ('*') /)
2011  Format (3 X, 'Type of initial conditions file ', T35, ' is ', A8)
2013  Format (/ / ' Title lines read from binary initial conditions fil&
     &e' /)
2020  Format (/ / 1 X, 'Input filenames' / '+', 15 ('_') // 3 X, 'Map', T28, '= ', A72 / 3 X, 'Geometry', T28, '= ', A72 / 3 X, 'In&
     &itial conditions in', T28, '= ', A72 / 3 X, 'Algal parameters', T28, '= ', A72 / 3 X, 'Zooplankton', T28, '= ', A72 / 3 X, 'S&
     &uspension Feeders', T28, '= ', A72)
2021  Format (3 X, 'Hydrodynamic', T28, '= ', A72:/(T30, A72))
2022  Format (3 X, 'Meteorologic', T28, '= ', A72:/(T30, A72))
2023  Format (3 X, 'Source One', T28, '= ', A72:/(T30, A72))
2024  Format (3 X, 'Source Two', T28, '= ', A72:/(T30, A72))
2025  Format (3 X, 'Boundary concentrations', T28, '= ', A72:/(T30, &
     & A72))
2026  Format (3 X, 'Benthic fluxes', T28, '= ', A72:/(T30, A72))
2027  Format (3 X, 'Settling', T28, '= ', A72/3 X, 'Mineralization', &
     & T28, '= ', A72/3 X, 'Light Extinction', T28, '= ', A72)
2028  Format (3 X, 'Aquatic vegetation', T28, '= ', A72:/(T30, A72))
2029  Format (3 X, 'Source Three', T28, '= ', A72:/(T30, A72))
2031  Format (3 X, 'Atmospheric Loads', T28, '= ', A72:/(T30, A72))
2030  Format (/ / '1', 'Output filenames' / '+', 16 ('_') // 3 X, 'Init&
     &ial conditions out', T28, '= ', A72 / 3 X, 'Snapshot', T28, '= ', &
     & A72 / 3 X, 'Restart output', T28, '= ', A72 / 3 X, 'Plot', T28, &
     & '= ', A72 / 3 X, 'Average plot', T28, '= ', A72 / 3 X, 'Diagnost&
     &ics', T28, '= ', A72 / 3 X, 'Transport fluxes', T28, '= ', A72 / &
     & 3 X, 'Kinetic fluxes', T28, '= ', A72 / 3 X, 'Oxygen plot', T28, &
     & '= ', A72 / 3 X, 'Mass balance', T28, '= ', A72 / 3 X, 'Algal pa&
     &rameters', T28, '= ', A72 / 3 X, 'Zooplankton', T28, '= ', A72 / &
     & 3 X, 'Sediment model', T28, '= ', A72 / 3 X, 'Aquatic vegetation&
     &', T28, '= ', A72 / 3 X, 'Suspension feeders', T28, '= ', A72)
2040  Format (/ / 1 X, 'Geometry' / '+', 8 ('_') // 3 X, 'Total boxes', &
     & T21, '= ', I5 / 3 X, 'Surface boxes', T21, '= ', I5 / 3 X, 'Tota&
     &l faces', T21, '= ', I5 / 3 X, 'Horizontal faces', T21, '= ', I5 &
     & / 3 X, 'Layers', T21, '= ', I5)
2050  Format (/ / 1 X, 'Time control' / '+', 12 ('_') // 3 X, 'Starting&
     & time = ', F8.2, ' Julian day' / 3 X, 'Ending time   = ', F8.2, '&
     & Julian day')
2060  Format (/ / 1 X, 'Timestep control' / '+', 16 ('_') // 3 X, 'Fixe&
     &d' / (5 X, 'Julian Day ', T39, '= ', 8 F8.2) : /)
2070  Format (5 X, 'Timestep (sec)', T39, '= ', 8 I8/(T49, 8 I8))
2071  Format (3 X, 'Variable'/5 X, 'Autostepping', T39, '= ', 5 X, A3)
2072  Format (5 X, 'Maximum allowable timestep (sec)', T39, '= ', 8 &
     & I8/(T41, 8 I8))
2073  Format (5 X, 'Fraction of calculated timestep', T39, '= ', 8 &
     & F8.2/(T41, 8 F8.2))
2080  Format (/ / 1 X, 'Input controls' / '+', 14 ('_') // 3 X, 'Restar&
     &t ', T33, '= ', A8 / 3 X, 'Boundary constituents', T33, '= ', A8 &
     & / 3 X, 'Source One loadings', T33, '= ', A8 / 3 X, 'Source Two l&
     &oadings', T33, '= ', A8 / 3 X, 'Source Three loadings', T33, '= '&
     & , A8 / 3 X, 'Benthic fluxes', T33, '= ', A8 / 3 X, 'Constituent &
     &modifications', T33, '= ', A8 / 3 X, 'Initial conditions in', &
     & T33, '= ', A8 / 3 X, 'Atmospheric loadings', T33, '= ', A8 / 3 &
     & X, 'SAV loads', T33, '= ', A8 / 3 X, 'Settling', T33, '= ', A8)
2090  Format (/ / '1', 'Output controls' / '+', 15 ('_') // 3 X, 'Initi&
     &al conditions out', T32, '= ', A3 / 3 X, 'Snapshot', T32, '= ', &
     & A3 / 3 X, 'Plot', T32, '= ', A3 / 5 X, 'Water quality diagnostic&
     &s', T32, '= ', A3 / 5 X, 'Sediment diagnostics', T32, '= ', A3 / &
     & 3 X, 'Average plot', T32, '= ', A3 / 3 X, 'Diagnostics', T32, '=&
     & ', A3 / 3 X, 'Transport fluxes', T32, '= ', A3 / 3 X, 'Kinetic f&
     &luxes', T32, '= ', A3 / 3 X, 'Oxygen plots', T32, '= ', A3 / 3 X, &
     & 'Restart', T32, '= ', A3)
2100  Format (/ / 3 X, 'Output dates (Julian day)' / '+', 2 X, 12 ('_') &
     & // 5 X, 'Snapshots', T24, '= ', 9 F8.2, (:/T26, 9 F8.2))
2110  Format (5 X, 'Plots', T24, '= ', 9 F8.2, (:/T26, 9 F8.2))
2115  Format (5 X, 'Average plots', T24, '= ', 9 F8.2, (:/T26, 9 F8.2))
2116  Format (5 X, 'Transport fluxes', T24, '= ', 9 F8.2, (:/T26, 9 &
     & F8.2))
2118  Format (5 X, 'Oxygen plots', T24, '= ', 9 F8.2, (:/T26, 9 F8.2))
2120  Format (5 X, 'Restarts', T24, '= ', 9 F8.2, (:/T26, 9 F8.2))
2130  Format (/ 3 X, 'Output frequencies (days)' / '+', 2 X, 18 ('_') &
     & // 5 X, 'Snapshots', T24, '= ', 9 F8.2, (:/T26, 9 F8.2))
2135  Format (5 X, 'Plots', T24, '= ', 9 F8.2, (:/T26, 9 F8.2))
2136  Format (5 X, 'Average plots', T24, '= ', 9 F8.2, (:/T26, 9 F8.2))
2137  Format (5 X, 'Transport fluxes', T24, '= ', 9 F8.2, (:/T26, 9 &
     & F8.2))
2138  Format (5 X, 'Oxygen plots', T24, '= ', 9 F8.2, (:/T26, 9 F8.2))
2141  Format (/ / 1 X, 'Balance calculations' / '+', 20 ('_') // 3 X, '&
     &Mass balance', T21, '= ', A3 / 3 X, 'Volume balance', T21, '= ', &
     & A3)
2150  Format (/ / 1 X, 'Hydrodynamics' / '+', 13 ('_') // 3 X, 'Origin',&
     &  T40, '= ', A8 / 3 X, 'Interval between updates', T40, '= ', I8, &
     & ' sec' / 3 X, 'Horizontal transport', T40, '= ', A8 / 3 X, 'Cons&
     &ervation type', T40, '= ', A8 / 3 X, 'Theta for vertical solution&
     &', T40, '= ', F8.2 / 3 X, 'Minimum time step', T40, '= ', F8.2, '&
     & sec' / 3 X, 'Days in Time-variable input files', T40, '= ', F8.2 &
     & //)
2152  Format ('1', 'Sediment calculations'/'+', 21('_')//3 X, 'Model = &
     &', A3)
2154  Format ('1', 'SAV calculations'/'+', 21('_')//3 X, 'Model = ', &
     & A3)
2160  Format (/ / 1 X, 'Boundary interpolation' / '+', 22 ('_') // 3 X, &
     & 'Temporal', T13, '=', A8)
2165  Format (/ / 1 X, 'Diffusion' / '+', 9 ('_') // 3 X, 'Horizontal', &
     & T24, '=', F6.1 / 3 X, 'Vertical multiplier', T24, '=', F6.1)
2170  Format (/ / 1 X, 'Dead sea case' / '+', 13 ('_') // 3 X, 'Flow', &
     & T25, '= ', A3 / 3 X, 'Horizontal diffusion', T25, '= ', A3 / 3 &
     & X, 'Vertical diffusion', T25, '= ', A3)
2180  Format (/ / 1 X, 'Constituent controls' / '+', 20 ('_') // 5 X, '&
     &Constituent', T28, 'Computation', T43, 'Initial concen', 'tration&
     &' / '+', 4 X, 11 ('_'), T28, 11 ('_'), T43, 7 ('_'), 1 X, 13 ('_'&
     & ) / T49, '(g/m**3)' // (3 X, A24, T32, A3, T47, F8.3))
  !           (3X,A24,T32,A3,T47,F8.3)//)
2190  Format (/ / '1', 2 X, 'Nutrient reductions' / '+', 2 X, 19 ('_') &
     & // 5 X, 'Nitrogen' / 7 X, 'Fall line', T24, '=', F5.2 / 7 X, 'Be&
     &low fall line', T24, '=', F5.2 / 7 X, 'Atmospheric', T24, '=', &
     & F5.2 / 7 X, 'Boundaries', T24, '=', F5.2 / 5 X, 'Phosphorus' / 7 &
     & X, 'Fall line', T24, '=', F5.2 / 7 X, 'Below fall line', T24, '=&
     &', F5.2 / 7 X, 'Atmospheric', T24, '=', F5.2 / 7 X, 'Boundaries', &
     & T24, '=', F5.2 / 5 X, 'Carbon' / 7 X, 'Fall line', T24, '=', &
     & F5.2 / 7 X, 'Below fall line', T24, '=', F5.2 / 7 X, 'Atmospheri&
     &c', T24, '=', F5.2 / 7 X, 'Boundaries', T24, '=', F5.2)
2200  Format (/ / 3 X, 'Half-saturation coefficients' / '+', 2 X, 28 ('&
     &_') // T25, 'Oxygen', T38, 'Nitrogen', T51, 'Phosphorus', T67, 'C&
     &arbon', T80, 'Silica', T93, 'Respiration' / '+', T25, 6 ('_'), &
     & T38, 8 ('_'), T51, 10 ('_'), T67, 6 ('_'), T80, 6 ('_'), T93, 11 &
     & ('_') / T23, '(g DO/m**3)', T37, '(g N/m**3)', T51, '(g P/m**3)',&
     &  T65, '(g C/m**3)', T78, '(g Si/m**3)', T93, '(g DO/m**3)' // 5 &
     & X, 'Nitrification', T24, F6.3, T38, F6.3 / 5 X, 'COD oxidation', &
     & T24, F6.3 / 5 X, 'DOC oxidation', T24, F6.3 / 5 X, 'Denitrificat&
     &ion', T24, F6.3 //)
2210  Format (3 X, 'Nutrient to carbon ratios'/'+', 2 X, 25('_')//T26, &
     & 'Oxygen', T37, 'Nitrogen', T50, 'Phosphorus', T66, 'Silica'/'+', &
     & T26, 6('_'), T37, 8('_'), T50, 10('_'), T66, 6('_')//5 X, 'Respi&
     &ration', T25, F6.3/5 X, 'Nitrification', T25, F6.3//)
2260  Format (3 X, 'Temperature/algal effects'/'+', 2 X, 25('_')//T28, &
     & 'Suboptimal', T50, 'Superoptimal', T77, 'Maximum'/'+', T28, 10('&
     &_'), T50, 12('_'), T77, 7('_')/T27, '(/degrees C)', T50, '(/degre&
     &es C)', T75, '(degrees C)'//5 X, 'Nitrification', T29, F6.3, T52, &
     & F6.3, T77, F6.1 /)
2270  Format (T26, 'Reference rate', T46, 'Reference temperature', T72, &
     & 'Temperature effect', T98, 'Algal effect'/'+', T26, 14('_'), &
     & T46, 21('_'), T72, 18('_'), T98, 12('_')/T30, '(/day)', T51, '(d&
     &egrees C)', T75, '(/degrees C)', T97, '(m**3/day/g C)'//5 X, 'Bas&
     &e metabolism'/7 X, 'Cyanobacteria', T52, F6.1, T77, F6.3/7 X, 'Di&
     &atoms', T52, F6.1, T77, F6.3/7 X, 'Greens', T52, F6.1, T77, F6.3)
2280  Format (5 X, 'COD oxidation', T52, F6.1, T77, F6.3/5 X, 'DOC oxid&
     &ation', T52, F6.1, T77, F6.3/5 X, 'LPOC hydrolysis', T52, F6.1, &
     & T77, F6.3/5 X, 'RPOC hydrolysis', T52, F6.1, T77, F6.3/5 X, 'DON&
     & oxidation', T52, F6.1, T77, F6.3/5 X, 'LPON hydrolysis', T52, &
     & F6.1, T77, F6.3/5 X, 'RPON hydrolysis', T52, F6.1, T77, F6.3/5 &
     & X, 'DOP oxidation', T52, F6.1, T77, F6.3/5 X, 'LPOP hydrolysis', &
     & T52, F6.1, T77, F6.3/5 X, 'RPOP hydrolysis', T52, F6.1, T77, &
     & F6.3/5 X, 'Silica dissolution', T52, F6.1, T77, F6.3//)
2317  Format (3 X, 'Sorption coefficients'/'+', 2 X, 21('_')//5 X, 'Pho&
     &sphorus', T17, '=', F6.2, ' m**3/gm'/5 X, 'Silica', T17, '=', &
     & F6.2, ' m**3/gm'//)
2330  Format (3 X, 'Miscellaneous coefficients'/'+', 2 X, 26('_')//5 X, &
     & 'Ratio of anoxic to oxic metabolism', T41, '=', F6.2/5 X, 'N/C r&
     &atio for anoxic metabolism', T41, '=', F6.3//, 5 X, 'Half-saturat&
     &ion salt for coagulation', T41, '=', F6.3)
2332  Format (3 X, 'Reaeration (m/d) = ', F8.3, ' * Rnu * (', F6.3, ' *&
     & WMS) **', F6.3 /)
2310  Format (3 X, 'Extinction coefficients'/'+', 2 X, 26('_')//5 X, 'B&
     &ackground', T41, '=', F6.2, ' 1/m'/5 X, 'Inorganic Solids', T41, &
     & '=', F6.2, ' m**2/gm'/5 X, 'Volatile Solids', T41, '=', F6.2, ' &
     &m**2/gm'/5 X, 'Dissolved Organic C', T41, '=', F6.2, ' m**2/gm'//&
     & )
2501  Format (' Temperature effects'//' DOC = ', F8.4, ' per Degree C'/&
     & ' NH4 = ', F8.4, ' per Degree C'/' NO3 = ', F8.4, ' per Degree C&
     &'/' PO4 = ', F8.4, ' per Degree C'/' SOD = ', F8.4, ' per Degree &
     &C'/' SIAT= ', F8.4, ' per Degree C')
2502  Format (/ / ' Other parameters' // ' MTCNO3   = ', F8.4, ' m/day' &
     & / ' SEDNO3   = ', F8.4, ' gm/m3' / ' KHSO     = ', F8.4, ' gm/m3&
     &')
2503  Format (/ / ' Reference Temperatures' // ' DOC = ', F8.4, ' Degre&
     &es C' / ' NH4 = ', F8.4, ' Degrees C' / ' NO3 = ', F8.4, ' Degree&
     &s C' / ' PO4 = ', F8.4, ' Degrees C' / ' SOD = ', F8.4, ' Degrees&
     & C' / ' SIAT= ', F8.4, ' Degrees C')
  !
  !
      NXAPL = APLTD (APLDP)
  !
      INFLOW = 0
      Do F = 1, NHQF !
         If (RIGHT_FLOWB(F) .Or. LEFT_FLOWB(F)) Then
            INFLOW = INFLOW + 1
            IFLOWP (INFLOW) = F
         End If
      End Do
  !
  !
      Return
End Subroutine INPUTS
!
Subroutine ERROR_CHECK (HEADER, FLAG)

      Implicit None
      Character (Len=120) HEADER
      Integer FLAG
      If (FLAG .Eq.-1) Then
         Write (*, "(a120)") HEADER !print error message in report.txt
         Stop
      Else
         Write (*, "('***** WARNING *****')")!print warning message on the screen
         Write (*, "(a120)") HEADER
      End If
      Return
End Subroutine ERROR_CHECK
