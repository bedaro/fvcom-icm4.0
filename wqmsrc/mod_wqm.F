!mod_wqm.F
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
Module MOD_WQM
  !
      Use MOD_SIZES, Only: NCP, NQFP, NS1P, NS2P, NS3P, NBCP, NFLP, &
     & NSSFP, MGL, NGL, NOBTY, NSAVM
  !
      Use MOD_PREC, Only: SP, CDF_PREC
  !
      Use MOD_LIMS, Only: MTLOC, KBM1, KB, MLOC
  !
  !
  !
#if defined (MULTIPROCESSOR)
      Use MOD_BUFFERS, Only: SNDBUF, RCVBUF
#endif
      Implicit None
      Save
  !********************************************************************************
  !**                                                                            **
  !**                        COMMON file for CE-QUAL-IC                          **
  !**                                                                            **
  !**                                Version 1.1                                 **
  !**                               July 28, 1992                                **
  !**                                                                            **
  !**                        Water Quality Modeling Group                        **
  !**                        U.S. Army Corps of Engineers                        **
  !**                        Waterways Experiment Station                        **
  !**                        Vicksburg, Mississippi 39180                        **
  !**                                                                            **
  !********************************************************************************
  !
  !********************************************************************************
  !**                                                                            **
  !**                           PARAMETER definitions                            **
  !**                                                                            **
  !**    NBP   - Number of boxes                                                 **
  !**    NQFP  - Number of horizontal and vertical flow faces                    **
  !**    NHQP  - Number of horizontal flow faces                                 **
  !**    NSBP  - Number of boxes in the surface layer                            **
  !**    NLP   - Maximum number of layers                                        **
  !**    NS1P  - Number of source one inputs                                     **
  !**    NS2P  - Number of source two inputs                                     **
  !**    NS3P  - Number of source three inputs                                   **
  !**    NBCP  - Number of boundary concentration inputs                         **
  !**    NMP   - Number of modifications to initial box concentrations           **
  !**    NDP   - Maximum number of days for any output                           **
  !**    NFLP  - Number of files for each type of time-varying data              **
  !**    NCP   - Number of constituent state variables                           **
  !**    NSAVM - Maximum Number of SAV Species                                   ** !:
  !**    NSSFP - Number of suspension-feeding species                            **
  !**                                                                            **
  !********************************************************************************
  !
  !
  !********************************************************************************
  !**                        Water Quality Model Setup                           **
  !********************************************************************************
  !
  !
  !
  ! index variables (typically used for do loops)
  !
  !
      Integer :: B !index for overall  boxes (i.e. nodes in fvcom)
      Integer :: BB !index for surface or bottom boxes (now equivalent to B)
  !
      Integer :: SAVB !array storing boxes indices where there is SAV in it in SAV model
  !this should be globalized and then broadcast to individual processors
  !******!!Never used!!*****
  !
  ! water column constituents
  !
      Real (SP), Allocatable, Dimension (:, :) :: T, SALT, SSI, B1, B2, &
     & B3, LDOC, RDOC, LPOC, RPOC, NH4, NO3, UREA, LDON, RDON, LPON, &
     & RPON, PO4, LDOP, RDOP, LPOP, RPOP, PIP, COD, DOXG, SIUPB, SIAT, &
     & PIB1, PIB2, PIB3, TDIC, TALK, pH, pCO2
  !
  !T,    & !1 water column temperature (degC)
  !SALT, & !2 salinity in water column (psu)
  !SSI,  & !3 inorganic solids (g/m^3)
  !B1,   & !4 algal group 1 biomass concentration (gC/m^3)
  !B2,   & !5 algal group 2 biomass concentration (gC/m^3)
  !B3,   & !6 algal group 3 biomass concentration (gC/m^3)
  !                                    SZ,   & !7 micro(small) zooplankton concentration (gC/m^3)  !  moved to mod_zoop
  !                                    LZ,   & !8 meso(large) zooplankton concentration (gC/m^3)     !  moved to mod_zoop
  !LDOC, & !9 labile dissolved organic carbon concentration (gC/m^3)
  !RDOC, & !10 refractory dissolved organic carbon concentration (gC/m^3)
  !LPOC, & !11 labile particulate organic carbon concentration (gC/m^3)
  !RPOC, & !12 refractory particulate organic carbon concentration (gC/m^3)
  !NH4,  & !13 ammonia concentration in water column (gN/m^3)
  !NO3,  & !14 nitrate concentration in water column (gN/m^3)
  !UREA, & !15 urea concentration (gN/m^3)    !: this term is never used !!!
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
  !SIAT, & !29 total silica (available) concentration in water column (gSi/m^3) (Changed SA to DSI and to SIAT)
  !PIB1, & !Phosphorus in algae 1 (gP/m^3)
  !PIB2, & !Phosphorus in algae 2 (gP/m^3)
  !PIB3    !Phosphorus in algae 3 (gP/m^3)
  !TDIC,  & !33 Dissolved inorg carbon (mmolC/m^3)    					   !
  !TALK, & !34 Total alkalinity (mmol/m^3) 		   		  			   !
  !pH,   & !-- pH (seawater scale)    									   !
  !pCO2    !-- partial pressure of CO2 in seawater (ppmv or microatm)     !
!
  !Water column settling rates
  !3D Spatially varying (horizontally and vertically) settling velocities (read from STLFN)
      Real (SP), Allocatable, Dimension (:, :) :: WSS, WSL, WSR, WS1, &
     & WS2, WS3, WSU, WSSHI !fixed solids settling rate (m/day) when solids concentration is high
  ! Set to WSS when SSI (fixed solids concentration ) <100 g/m^3
  ! Set to 5 (m/day) when fixed solids concentration > 100 g/m^3
  !
      Real (SP) :: KHR1, KHR2, KHR3, KHN1, KHN2, KHN3, KHNAVG, KHP1, &
     & KHP2, KHP3, KHPAVG, KHONT, KHOCOD, KHODOC, KHNNT, KHNDN, KHNFIX, &
     & FNFIX, KHS1, KHS2, KHS3, KHST1, KHST2, KHST3 !Half saturation constant of salinity at which mortality is half the maximum value for algal group 3 (ppt)
  !
  !
  !parameters for specifying benthic flux for the water column model
  !original water column model was not coupled with sediment flux model
  !and benthic fluxes are specified using these parameters
  !
      Real (SP) :: KSO, KSNH4, KSNO3, KSPO4, KSSA, KSDOC, TRSO, TRSNH4, &
     & TRSNO3, TRSPO4, TRSSA, TRSDOC, MTCNO3, SEDNO3, KHSO !half saturation constant of SOD for oxygen (gO2/m^3)
  !
      Real (SP), Allocatable :: BENDOCB (:), BENNH4B (:), BENNO3B (:), &
     & BENPO4B (:), BENCODB (:), BENDOB (:), BENSAB (:)!benthic dissolved silica flux (gSi/m^2/day) at reference temperature
  !
  !Temperature control on nitriciation
      Real (SP) :: TMNT, KTNT1, KTNT2 !effect of temperature above TMNT on nitrifcation (1/degC^2)
  !
      Real (SP) :: KTCOD, KTMNL, KTHDR, KTSUA !effect of temperature on particulate biogenic silica dissolution rates (1/degC)
  !
      Real (SP) :: KHCOAG !salinity at which coagulation of RDOC and RDON is half the maxium value (ppt)
  !
      Integer :: HYDPTR, METPTR, CBCPTR, S1PTR, S2PTR, S3PTR, BFIPTR, &
     & ATMPTR, SAVPTR, BAIPTR, WCLPTR, NTVWCLF !number of sediment flux model water column input files
  !
      Real (SP) :: KHTIS ! Never used
  !
      Real (SP) :: NXCBC, JDAY, JDAYMBL
  !
      Real (SP) :: KT !Coef of heat exchange (watts/m^2/C)
  !
      Real (SP) :: NXMET, NXS1, NXS2, NXS3, NXBFI, NXATM, NXSAV
  !
      Real (SP) :: I0, IT, ITNX, ISMIN !NOT used in the code  !IT is in E/m2/day unit for light
  !
      Real (SP) :: KADPO4, KADSA !Partition coefficient for sorption of silica on inorganic solids (m^3/gSolids)
  !SA-- available silica
  !SUA -- unavailable particulate biogenic silica
  !
  !
  !Moved all these to the mod_zoop.F
  !!Zooplankton parameters
  !
  !     REAL(SP) :: CTSZ,    &  !carbon threshold for grazing (gC/m^3) for microzooplankton
  !             CTLZ,    &  !carbon threshold for grazing (gC/m^3) for mesozooplankton
       !             KHCSZ,   &  !prey density at which grazing is halved (gC/m^3) for microzooplankton
  !             KHCLZ,   &  !prey density at which grazing is halved (gC/m^3) for mesozooplankton
       !             MZEROSZ, &  !mortality at zero dissolved oxygen (1/day) for microzooplankton
  !             MZEROLZ     !mortality at zero dissolved oxygen (1/day) for mesozooplankton
  !
  !     REAL(SP) :: UB1SZ, &   !Utilization of algal group 1 by microzooplankton, range [0,1]
  !             UB2SZ, &   !Utilization of algal group 2 by microzooplankton, range [0,1]
       !             UB3SZ, &   !Utilization of algal group 3 by microzooplankton, range [0,1]
  !             UB1LZ, &   !Utilization of algal group 1 by mesozooplankton, range [0,1]
       !             UB2LZ, &   !Utilization of algal group 2 by mesozooplankton, range [0,1]
  !             UB3LZ, &   !Utilization of algal group 3 by mesozooplankton, range [0,1]
       !             UDSZ,  &   !Utilization of dissolved organic carbon by microzooplankton, range [0,1]
  !             ULDSZ, &   !Utilization of labile dissolved organic carbon by microzooplankton, range [0,1]
       !             URDSZ, &   !Utilization of refractory dissolved organic carbon by microzooplankton, range [0,1]
  !             ULPSZ, &   !microzooplankton utilization of labile particulate organic carbon, range [0,1]
       !             URPSZ, &   !microzooplankton utilization of refractory particulate organic carbon, range [0,1]
  !             ULLZ,  &   !mesozooplankton utilization of labile particulate organic carbon, range [0,1]
       !             URLZ,  &   !mesozooplankton utilization of refractory particulate organic carbon, range [0,1]
  !             USZLZ, &   !mesozooplankton utilization of microzooplankton, range [0,1], (~1.0)
       !             TRSZ,  &   !reference temperature for microzooplankton metabolism (degC)
  !             TRLZ,  &   !reference temperature for mesozooplankton metabolism (degC)
       !             DOCRITSZ, &!dissolved oxygen below which mortality occurs for microzooplankton (gO2/m^3)
  !             DOCRITLZ, &!dissolved oxygen below which mortality occurs for mesozooplankton (gO2/m^3)
       !             ANCSZ,  &  !microzooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
  !             ANCLZ,  &  !mesozooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
       !             APCSZ,  &  !microzooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
  !             APCLZ,  &  !mesozooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
       !             AOCRSZ, &  !ratio of oxygen consumed to microzooplankton carbon metabolized (gO2/gC) (~2.67)
  !             AOCRLZ, &  !ratio of oxygen consumed to mesozooplankton carbon metabolized (gO2/gC) (~2.67)
       !             FRSASZ, &  !fraction of microzooplankton silica recycled to dissolved silica pool, range [0,1]
  !             FRSALZ     !fraction of mesozooplankton silica recycled to dissolved silica pool, range [0,1]
  !
  !     REAL(SP) :: FLDOCSZ,&    !fraction of microzooplankton carbon released to LDOC, range [0,1]
  !             FRDOCSZ,&    !fraction of microzooplankton carbon released to RDOC, range [0,1]
       !             FLPOCSZ,&    !fraction of microzooplankton carbon released to LPOC, range [0,1]
  !             FRPOCSZ,&    !fraction of microzooplankton carbon released to RPOC, range [0,1]
       !             FLDONSZ,&    !fraction of microzooplankton nitrogen released to LDON, range [0,1]
  !             FRDONSZ,&    !fraction of microzooplankton nitrogen released to RDON, range [0,1]
       !             FLPONSZ,&    !fraction of microzooplankton nitrogen released to LPON, range [0,1]
  !             FRPONSZ,&    !fraction of microzooplankton nitrogen released to RPON, range [0,1]
       !             FLDOPSZ,&    !fraction of microzooplankton phosphorus released to LDOP range [0,1]
  !             FRDOPSZ,&    !fraction of microzooplankton phosphorus released to RDOP range [0,1]
       !             FLPOPSZ,&    !fraction of microzooplankton phosphorus released to LPOP range [0,1]
  !             FRPOPSZ,&    !fraction of microzooplankton phosphorus released to RPOP range [0,1]
       !             FNH4SZ, &    !fraction of microzooplankton nitrogen recycled to DIN as NH4, range [0,1]
  !             FPO4SZ, &    !fraction of microzooplankton phosphorus recycled to DIP as PO4, range [0,1]
       !             FUREASZ,&    !Never used !!!
  !             FLDOCLZ,&    !fraction of mesozooplankton carbon released to LDOC, range [0,1]
       !             FRDOCLZ,&    !fraction of mesozooplankton carbon released to RDOC, range [0,1]
  !             FLPOCLZ,&    !fraction of mesozooplankton carbon released to LPOC, range [0,1]
       !             FRPOCLZ,&    !fraction of mesozooplankton carbon released to RPOC, range [0,1]
  !             FLDONLZ,&    !fraction of mesozooplankton nitrogen released to LDON, range [0,1]
       !             FRDONLZ,&    !fraction of mesozooplankton nitrogen released to RDON, range [0,1]
  !             FLPONLZ,&    !fraction of mesozooplankton nitrogen released to LPON, range [0,1]
       !             FRPONLZ,&    !fraction of mesozooplankton nitrogen released to RPON, range [0,1]
  !             FLDOPLZ,&    !fraction of mesozooplankton phosphorus released to LDOP range [0,1]
       !             FRDOPLZ,&    !fraction of mesozooplankton phosphorus released to RDOP range [0,1]
  !             FLPOPLZ,&    !fraction of mesozooplankton phosphorus released to LPOP range [0,1]
       !             FRPOPLZ,&    !fraction of mesozooplankton phosphorus released to RPOP range [0,1]
  !             FNH4LZ, &    !fraction of mesozooplankton nitrogen recycled to DIN as NH4, range [0,1]
       !             FPO4LZ, &    !fraction of mesozooplankton phosphorus recycled to DIP as PO4, range [0,1]
  !             FUREALZ      !Never used!!!
  !
  !     REAL(SP),ALLOCATABLE::  B1ASZ(:,:),   &
  !                         B2ASZ(:,:),   &
       !                         B3ASZ(:,:),   &
  !                         LPOCASZ(:,:), &
       !                         RPOCASZ(:,:), &
  !                         PRASZ(:,:),   &
       !                         B1ALZ(:,:),   &
  !                         B2ALZ(:,:),   &
       !                         B3ALZ(:,:),   &
  !                         SZALZ(:,:),   &
       !                         LPOCALZ(:,:), &
  !                         RPOCALZ(:,:), &
       !                         PRALZ(:,:),   &
  !                         CLSZ(:,:),    &
       !                         CLLZ(:,:),    &
  !                         RSZ(:,:),     &
       !                         RLZ(:,:),     &
  !                         RMAXSZ(:,:),  &
       !                         RMAXLZ(:,:),  &
  !                         BMSZ(:,:),    &
       !                         BMLZ(:,:),    &
  !                         BMRSZ(:,:),   &
       !                         BMRLZ(:,:),   &
  !                         MSZ(:,:),     &
       !                         MLZ(:,:),     &
  !                         PRSZLZ(:,:),  &
       !                         GSZ(:,:),     &
  !                         GLZ(:,:),     &
       !                         ESZ(:,:),     &
  !                         ELZ(:,:),     &
       !                         RFSZ(:,:),    &
  !                         RFLZ(:,:),    &
       !                         PRSZ(:,:),    &
  !                         PRLZ(:,:),    &
       !                         LDOCASZ(:,:), &
  !                         BPRSZ(:,:),   &
       !                         BPRLZ(:,:),   &
  !                         RDOCASZ(:,:), &
       !                         LDOCSZ(:,:),  &
  !                         LPOCSZ(:,:),  &
       !                         RPOCSZ(:,:),  &
  !                         LDOCLZ(:,:),  &
       !                         LPOCLZ(:,:),  &
  !                         RPOCLZ(:,:),  &
       !                         NH4SZ(:,:),   &
  !                         LDONSZ(:,:),  &
       !                         LPONSZ(:,:),  &
  !                         RPONSZ(:,:),  &
       !                         NH4LZ(:,:),   &
  !                         LDONLZ(:,:),  &
       !                         LPONLZ(:,:),  &
  !                         RPONLZ(:,:),  &
       !                         PO4SZ(:,:),   &
  !                         LDOPSZ(:,:),  &
       !                         LPOPSZ(:,:),  &
  !                         RPOPSZ(:,:),  &
       !                         PO4LZ(:,:),   &
  !                         LDOPLZ(:,:),  &
       !                         LPOPLZ(:,:),  &
  !                         RPOPLZ(:,:),  &
       !                         RDOCSZ(:,:),  &
  !                         RDONSZ(:,:),  &
       !                         RDOPSZ(:,:),  &
  !                         RDOCLZ(:,:),  &
       !                         RDONLZ(:,:),  &
  !                         RDOPLZ(:,:),  &
       !                         PIB1SZ(:,:),  &
  !                         PIB2SZ(:,:),  &
       !                         PIB3SZ(:,:),  &
  !                         PIB1LZ(:,:),  &
       !                         PIB2LZ(:,:),  &
  !                         PIB3LZ(:,:),  &
       !                         B1SZ(:,:),    &
  !                         B2SZ(:,:),    &
       !                         B3SZ(:,:),    &
  !                         B1LZ(:,:),    &
       !                         B2LZ(:,:),    &
  !                         B3LZ(:,:),    &
       !                         DOSZ(:,:),    &
  !                         DOLZ(:,:),    &
       !                         SASZ(:,:),    &
  !                         SUSZ(:,:),    &
       !                         SALZ(:,:),    &
  !                         SULZ(:,:)
  !     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: ACLSZ,   ACLLZ,    ARSZ,    ARLZ,     &
  !               ABMSZ,   ABMLZ,    AMSZ,    AMLZ,     &
       !               APRSZLZ, AGSZ,     AGLZ,    ADOCSZ,   &
  !               APOCSZ,  ADOCLZ,   APOCLZ,                 &
       !               ANH4SZ,  ADONSZ,   APONSZ,                 &
  !               ANH4LZ,  ADONLZ,   APONLZ,                 &
       !               APO4SZ,  ADOPSZ,   APOPSZ,                 &
  !               APO4LZ,  ADOPLZ,   APOPLZ,                 &
       !               APRSZ,   APRLZ,    APISZ,  APILZ
  !
  !     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: AB1SZ,   AB2SZ,    AB3SZ,  AB1LZ,     &
  !               AB2LZ,   AB3LZ,    ADOSZ,  ADOLZ,     &
       !               ASASZ,   ASUSZ,    ASALZ,  ASULZ
  !
  !     REAL(SP),DIMENSION(-50:400) :: FTLZ,  FTSZ,  FTBMSZ,          &
  !               FTBMLZ,FTPRSZ,FTPRLZ,
  !
  !
  !!Tempreature control on algae, should move to mod_algal.F   !!: commented because we now use get_ft1,get_ftbm1,etc
  !       REAL(SP),DIMENSION(-50:400) ::  FT1,            &
  !                                   FT2,            &
       !                                   FT3,         &
  !                                FTBM1,            &
       !                                FTBM2,            &
  !                                FTBM3,          &
       !                                FTPR
  !
  !
  !  variables needed by autostep and main
      Integer, Save :: NXHYD = 0
      Real (SP), Save :: COURMX, COURQS, COURVS
      Real (SP), Save :: DIFFMX, DIFFAS, DIFFDS, FNDLT
  !
!
  !
      Real (SP), Save :: DLT8, MAXDLTDP, MXDLTDP, NXHYDDP, ELTMS1
  !
      Real (SP) :: ELTMS
  !
      Logical :: BOUNDARY_CONC, SOURCE_ONE, SOURCE_TWO, SOURCE_THR, &
     & SOLIDS_CALC, PIP_CALC, BENTHIC_FLUXES, SEDIMENT_CALC, &
     & LIGHT_EXTINCTION, ATMOS_LOADS, TRANSPORT_FLUXES, AVERAGE_PLOTS, &
     & QUALITY_DIAG, KINETIC_FLUXES, SEDIMENT_DIAG, DIAGNOSTICS, &
     & CONSERVE_MASS, SETTLING, STEP_BOUNDARY, MASS_BALANCE,SEDTEMP_CALC
  !
      Logical :: SAV_LOADS, SAV_CALC, SAV_PLOTS
  !
      Logical :: FLOW, XY_DIFFUSION, Z_DIFFUSION, BINARY_HYDRO, &
     & ASCII_HYDRO, DEPTH_AVG_HYDRO, SIGMA_HYDRO
      Character (Len=8) :: SPVARKE, PRINTKE
      Character (Len=72) :: KEIFN, SVOFN, BFOFN, ZOOFN, ZFOFN, ALOFN, &
     & BAOFN, DFOFN ! , Deposition feeder output
  !
      Character (Len=8) :: MINPROD
  !
  !***** Dimension declarations
  !
      Real (SP), Allocatable :: DTM (:, :, :)
  !
      Real (SP), Allocatable :: ADVFLUX (:, :, :)!: 
      Real (SP), Allocatable :: VDIFFLUX (:, :, :)!: 
      Real (SP), Allocatable :: ADVFLUX_GL (:, :, :)!: 
      Real (SP), Allocatable :: VDIFFLUX_GL (:, :, :)!:
!
  !
  !
      Real (SP), Allocatable :: DEPTHR (:, :)
      Real (SP), Allocatable, Dimension (:, :, :) :: C1, C2, C2F, &
     & CSTAR, DTC, AC1
  ! Added C2_GL
      Real (SP), Allocatable, Dimension (:, :, :) :: C2_GL
      Real (SP), Allocatable, Dimension (:, :) :: CCHL1_GL, T_GL, S_GL !
      Real (SP), Allocatable, Dimension (:) :: total_netPP_GL
      Real (SP), Allocatable, Dimension (:) :: D_GL, H_GL, EL_GL !, D_GL : total depth (global)
  !          H_GL : still depth (global)
  !         EL_GL : surface elevation (global)
      Real (SP), Allocatable, Dimension (:, :) :: pH_GL, pCO2_GL !
  !: to save GPP and NPP in history output
      Real (SP), Allocatable, Dimension (:, :) :: GPP_GL, NPP_GL
  !
  !Moved following BFLUX terms from tvds.F to here
      Real (SP), Allocatable, Dimension (:, :) :: BFLUX_GL
      Real (SP), Allocatable, Dimension (:, :) :: BFLUXB_GL
      Real (SP), Allocatable, Dimension (:, :) :: BFLUXNX, BFLUXNX_GL, &
     & BFLUX
  !
  !
  !Moved the following variables from ncdio.F to here
  !
      Real (CDF_PREC), Allocatable, Dimension (:, :) :: UL_GL, VL_GL !u and v
      Real (CDF_PREC), Allocatable, Dimension (:, :) :: WTSL_GL, KHL_GL !wts, kh
      Real (CDF_PREC), Allocatable, Dimension (:, :) :: SL_GL, TL_GL !salinity,temp
      Real (CDF_PREC), Allocatable, Dimension (:) :: ELL_GL, DTFAL_GL !zeta,dtfa
      Real (CDF_PREC), Allocatable, Dimension (:) :: UARD_OBCN_GL !uard_obcn
      Real (CDF_PREC), Allocatable, Dimension (:, :) :: XFLUX_OBC_GL !xflux_obc


      Real (SP), Allocatable, Dimension (:, :) :: TSS_GL, KESS_GL, SSI_GL, TURB_GL, PARAD_GL

	  Real (SP), Allocatable, Dimension (:, :) :: LEAF_GL, ROOT_GL, STEM_GL, &
	 & TUBER_GL, EP_GL
  !
      Real (SP), Dimension (NCP) :: CMASS !Total mass of all constiutents (kg)
  !
      Real (SP) :: S1L (NS1P, NCP), S2L (NS2P, NCP), S3L (NS3P, NCP)
      Integer :: S1LB (NS1P, NCP), S2LB (NS2P, NCP), S3LB (NS3P, NCP)
  !
      Real (SP) :: CB (NBCP, NCP), CBNX (NBCP, NCP), CBOLD (NBCP, NCP)
  !
      Real (SP), Allocatable, Dimension (:) :: CTEMP !sediment temperature (degC)
  !
  !!The following needs to be moved to SED_INIT

      Real (SP), Allocatable, Dimension (:) :: BSVOL !
  !
  !water column nitrogen balance check  (positive into water column)
      Real (SP), Allocatable, Dimension (:) :: ATMFLXNB, ATMFLXPB, &
     & ATMFLXCB !atmospheric flux of carbon
!
      Real (SP), Allocatable, Dimension (:) :: BENFLXPNB, BENFLXDNB, &
     & BENFLXPPB, BENFLXDPB, BENFLXPCB !benthic POC flux (kgC/day)

  !sediment nitrogen check (positive into sediments)
      Real (SP), Allocatable, Dimension (:) :: DLSEDKNB, DLSEDKCB, &
     & BURIALFLXNB, BURIALFLXPB, BURIALFLXCB
      Real (SP) :: FLUXT (0:NQFP, NCP), AFLUX (NQFP, 13)!fluxes through all faces of 13 chosen species
  !
      Real (SP), Dimension (NQFP) :: FLXTTEM, FLXTSAL, FLXTSSI, FLXT1, &
     & FLXT2, FLXT3, FLXTLDOC, FLXTLPOC, FLXTRPOC, FLXTNH4, FLXTNO3, &
     & FLXTLDON, FLXTLPON, FLXTRPON, FLXTPO4, FLXTLDOP, FLXTLPOP, &
     & FLXTRPOP, FLXTCOD, FLXTDO, FLXTSU, FLXTSA, FLXTSZ, FLXTLZ, &
     & FLXTRDOC, FLXTRDON, FLXTRDOP, FLXTUREA, FLXTPIP, FLXTPIB1, &
     & FLXTPIB2, FLXTPIB3 !32
  !
  !settling flux of all particulate constituents in water column through
  !bottom face of each layer for each TCE
      Real (SP), Allocatable :: FLUXS (:, :, :)!size MTLOC x KBM1 x 13

      Real (SP), Allocatable, Dimension (:, :) :: FLXSPIB1, FLXSPIB2, &
     & FLXSPIB3 !13 of FLUXS above
  !
      Real (SP), Allocatable, Dimension (:) :: ABENDOC, ABENNH4, &
     & ABENNO3, ABENPO4, ABENCOD, ABENDO, ABENSA, ABENCH4G, ABENCH4A, &
     & ASSFWS, APCFWS, APNFWS, APPFWS, APSFWS, ACPIP, ACPOS, ADFEED
  !
      Real (SP), Allocatable, Dimension (:, :) :: ACPOC, ACPON, ACPOP, &
     & ASFEED
  !
      Real (SP), Allocatable, Dimension (:) :: AJNSF, AJPSF, ASODSF, &
     & ASASF, ASUSF, ASFGCIN, ASFCFEC, ASFCPSF, AFLXCSF, AFLXNSF, &
     & AFLXPSF, ARPOCSF, ARPONSF, ARPOPSF, ASSISF, ASSISASF, ASSISUSF, &
     & ASSIPSF
  !
      Real (SP), Allocatable, Dimension (:, :) :: ACCHL1, ACCHL2, &
     & ACCHL3
  !
      Integer :: S1LN (NCP), S2LN (NCP), S3LN (NCP), AC (NCP), NCB &
     & (NCP)!J    ,IBC(NCP,NBCP)
  !
      Real (SP), Dimension (0:NQFP) :: Q, A, DIFF
  !
      Real (SP), Allocatable, Dimension (:, :) :: V2 !,V1?? LB - lost in cleanup or not needed?
  !
      Real (SP), Allocatable, Dimension (:, :) :: PN1, PN2, PN3, BM1, &
     & BM2, BM3, FTCOD, NPP, P1, P2, P3, NT, PR1, PR2, PR3, RATOX, GPP, &
     & P1NNF
  !
      Real (SP), Allocatable, Dimension (:) :: ASRAT, total_netPP
  !
      Real (SP), Allocatable, Dimension (:) :: CFIX, SNFIX
  !
      Real (SP), Allocatable, Dimension (:, :) :: AFI1, ANL1, APL1, &
     & AFI2, ANL2, APL2, ASL2, AFI3, ANL3, APL3, ANPP, ARESP, AKE, &
     & ASL1, ASL3, AGPP
  !
      Real (SP), Allocatable, Dimension (:) :: AASRAT, ACFIX
  !
      Real (SP), Allocatable, Dimension (:) :: AFIB, ANLB, APLB, ANPPB, &
     & ABBM, ABLITE, ABMB, APB, APRB, ABADOC, ABAPOC, ABANH4, ABANO3, &
     & ABAPON, ABAPO4, ABAPOP, ABADO
  !
  !
     Real (SP), Allocatable, Dimension (:, :) :: PM1
	 Real (SP), Allocatable, Dimension (:, :) :: PM2
	 Real (SP), Allocatable, Dimension (:, :) :: PM3
	 Real (SP), Allocatable, Dimension (:, :) :: BMR1

     Real (SP), Allocatable, Dimension (:, :) :: BMR2
	 Real (SP), Allocatable, Dimension (:, :) :: BMR3
	 Real (SP), Allocatable, Dimension (:, :) :: BPR1	!Holling Type III predation, has unit [1/day]/[gC/m^3]
	 Real (SP), Allocatable, Dimension (:, :) :: BPR2   !Holling Type III predation, has unit [1/day]/[gC/m^3]
	 Real (SP), Allocatable, Dimension (:, :) :: BPR3   !Holling Type III predation, has unit [1/day]/[gC/m^3]
	 Real (SP), Allocatable, Dimension (:, :) :: CCHL1
	 Real (SP), Allocatable, Dimension (:, :) :: CCHL2
	 Real (SP), Allocatable, Dimension (:, :) :: CCHL3
	 Real (SP), Allocatable, Dimension (:, :) :: Q1
	 Real (SP), Allocatable, Dimension (:, :) :: Q2
	 Real (SP), Allocatable, Dimension (:, :) :: Q3 !Phosphorus to Carbon ratio (gP/gC) for algal group 3

  !
  !     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: DTT, DTSSI, DTB1, DTB2, DTB3,      &
  !      REAL(SP),ALLOCATABLE,DIMENSION(:,:) ::   &
       !               DTLDOC, DTLPOC, DTRPOC, DTNH4, DTNO3, DTLDON, DTLPON,       &
  !               DTRPON, DTPO4, DTLDOP, DTLPOP, DTRPOP, DTCOD, DTDO, DTSIUPB,   &
       !               DTSIAT, DTSZ, DTLZ, DTRDOC, DTRDON, DTRDOP, DTUREA, DTPIP,    &
  !                DTUREA
  !               DTPIB1, DTPIB2, DTPIB3
  !
      Real (SP), Allocatable, Dimension (:, :) :: FI1, FI2, FI3, NL1, &
     & NL2, NL3, PL1, PL2, PL3, SL2, RESP, SL1, SL3, KE, KEISS, KEVSS, &
     & KEDOC !  
  !
      Real (SP), Allocatable, Dimension (:, :) :: FTMNL, FTHDR
  !
      Real (SP), Allocatable, Dimension (:, :) :: RESPC, DLALGC
  !
      Real (SP), Allocatable, Dimension (:, :) :: IWCMNB, IWCMPB, &
     & IWCMCB, IWCMSB, WCMNB, WCMPB, WCMCB, WCMSB, DLWCMNB, DLWCMPB, &
     & DLWCKMNB, DLWCKMCB, S1FLXNB, S1FLXPB, S1FLXCB, S2FLXNB, S2FLXPB, &
     & S2FLXCB, S3FLXNB, S3FLXPB, S3FLXCB
  !
      Real (SP), Dimension (NCP) :: DLWCMCB
  !
      Real (SP), Allocatable, Dimension (:) :: BENDOC, BENNH4, BENNO3, &
     & BENDON, BENDOP, BENPO4, BENCOD, BENDO, BENSA, BENDEN, BENCH4G, &
     & BENCH4A
#if defined(CSOD_CALC)
      Real (SP), Allocatable, Dimension (:) :: CSODflux !: CSODflux=CSOD gO2/m^2/day
#endif
  !
      Integer, Allocatable, Dimension (:) :: SBN, BBN, HMSBV
  !SFA, replaced by ART1
  !
      Real (SP), Allocatable, Dimension (:) :: WSSNET, WSLNET, WSRNET, &
     & WS1NET, WS2NET, WS3NET, WSUNET !settling rate of particulate biogenic (unavaiable) silica (m/d)
!
      Real (SP), Allocatable, Dimension (:) :: PPFWS, PNFWS, PCFWS, &
     & PSFWS, SSFWS !suspended solids flux (gSolids/m^2/dat)
  !
      Character (Len=72), Dimension (NFLP) :: METFN, S1FN, S2FN, S3FN, &
     & HYDFN, CBCFN, ATMFN, SVIFN, BFIFN !Benthic Flux data file or Sediment Diagenesis module input file
  !
      Character (Len=72) BAIFN !benthic algae input file name
  !
      Character (Len=72) DFIFN !deposition feeder input file name
!
      Real (SP), Allocatable, Dimension (:, :) :: KLDC, KRDC, KLDN, &
     & KRDN, KLDP, KRDP, KLPC, KRPC, KLPN, KRPN, KLPP, KRPP, KDCALG, &
     & KDNALG, KDPALG, KLCALG, KLNALG, KLPALG, KSUA, KCOD, KRCOAG, NTM
  ! KLDC,    &  !labile dissolved organic carbin mineralization rate (1/day)
  !KRDC,    &  !refractory dissolved organic carbin mineralization rate (1/day)
  !KLDN,    &  !labile dissolved organic nitrogen mineralization rate (1/day)
  !KRDN,    &  !refractory dissolved organic nitrogen mineralization rate (1/day)
  !KLDP,    &  !labile dissolved organic phosphorus mineralization rate (1/day)
  !KRDP,    &  !refractory dissolved organic phosphorus mineralization rate (1/day)
  !KLPC,    &  !labile particulate organic carbon mineralization rate (1/day)
  !KRPC,    &  !refractory particulate organic carbon mineralization rate (1/day)
  !KLPN,    &  !labile particulate organic nitrogen remin rate (1/day)
  !KRPN,    &  !refractory particulate organic nitrogen remin rate (1/day)
  !KLPP,    &  !labile particulate organic phosphorus remin rate (1/day)
  !KRPP,    &  !refractory particulate organic phosphorus remin rate (1/day)
  !KDCALG,  &  !constant relating DOC respiration rate to total algae biomass C (1/day)/(gC/m^3)
  !KDNALG,  &  !constant relating DON respiration rate to total algae biomass N (1/day)/(gN/m^3)
  !KDPALG,  &  !constant relating DOP respiration rate to total algae biomass P (1/day)/(gP/m^3)
  !KLCALG,  &  !constant relating LPOC hydrolysis rate to total algae biomass C (1/day)/(gC/m^3)
  !KLNALG,  &  !constant relating LPON hydrolysis rate to total algae biomass N (1/day)/(gN/m^3)
  !KLPALG,  &  !constant relating LPOP hydrolysis rate to total algae biomass P (1/day)/(gP/m^3)
  !KSUA,    &  !particulate silica dissolution rate (1/day)
  !KCOD,    &  !oxidation rate of chemicao oxygen demand (1/day)
  !KRCOAG,  &  !dissolved organic carbon coagulation rate (1/day)
  !			 Coagulation rate increases as salinity increases
  !NTM         !maximum nitrication rate at optimal temperature (gN/m^3/day from mannual) !Laura: should be 1/d
  !			 nitrification is limited by NH4 and DOXG availability
  !
      Real (SP), Allocatable, Dimension (:, :) :: ALGDOC, ALGPOC, &
     & ALGDON, ALGPON, ALGDOP, ALGPOP, ALGNH4, ALGNO3, ALGPO4, MNLLDOC, &
     & MNLRDOC, MNLLDON, MNLRDON, MNLLDOP, MNLRDOP, HDRLPOC, HDRRPOC, &
     & HDRLPON, HDRRPON, HDRLPOP, HDRRPOP, DENNO3, DENIT, NFIX !Nitrogen fixation rate (gN/m^3/day)

  !ALGDOC,  &      !Dissolvoed organic carbon produced due to algae (all groups)
  !!mortality (respiration, metabolism and predation)
  !!labile and refractory DOCs together (gC/m^3/day)
  !ALGPOC,  &      !Particulate organic carbon produced due to algae (all groups)
  !mortaility (respiration, metabolism and predation)
  !(gC/m^3/day)
  !ALGDON,  &      !Dissolvoed organic nitrogen produced due to algae (all groups)
  !mortality (respiration, metabolism and predation)
  !labile and refractory DONs together
  !(gN/m^3/day)
  !ALGPON,  &      !Particulate organic nitrogen produced due to algae (all groups)
  !mortaility (respiration, metabolism and predation)
  !(gN/m^3/day)
  !ALGDOP,  &      !Dissolvoed organic phosphorus produced due to algae (all groups)
  !mortality (respiration, metabolism and predation)
  !labile and refractory DONs together
  !(gP/m^3/day)
  !ALGPOP,  &      !Particulate organic phosphorus produced due to algae (all groups)
  !mortality (respiration, metabolism and predation)
  !(gP/m^3/day)
  !ALGNH4,  &      !Soure+sink term for water column NH4 due to algae uptake (gN/m^3/day)
  !positive as source of NH4
  !ALGNO3,  &      !Source+sink term for water column NO3 due to algae uptake (gN/m^3/day)
  ! positive as source of NO3
  ! ALGPO4,  &      !Source+sink term for water column PO4 due to algae uptake (gP/m^3/day)
  ! positive as source of PO4
  ! MNLLDOC, &      !Mineralization rate of labile DOC (positive losing LDOC)
  ! gC/m^3/day
  ! MNLRDOC, &      !Mineralization rate of refractory DOC (positive losing RDOC)
  ! gC/m^3/day
  ! MNLLDON, &      !Mineralization rate of labile DON (positive losing LDON)
  ! gN/m^3/day
  ! MNLRDON, &      !Mineralization rate of refractory DON (positive losing RDON)
  ! gN/m^3/day
  ! MNLLDOP, &      !Mineralization rate of labile DOP (positive losing LDOP)
  ! gP/m^3/day
  ! MNLRDOP, &      !Mineralization rate of refractory DOC (positive losing RDOP)
  ! gP/m^3/day
  ! HDRLPOC, &      !Hydrolysis loss rate of labile POC (positive losing LPOC)
  ! gaining LDOC (gC/m^3/day)
  ! HDRRPOC, &      !Hydrolysis loss rate of refractory POC (positive losing RPOC)
  ! gaining RDOC (gC/m^3/day)
  ! HDRLPON, &      !Hydrolysis loss rate of labile PON (positive losing LPON)
  ! gaining LDON (gN/m^3/day)
  ! HDRRPON, &      !Hydrolysis loss rate of refractory PON (positive losing RPON)
  ! gaining RDOC (gN/m^3/day)
  ! HDRLPOP, &      !Hydrolysis loss rate of labile POP (positive losing LPOP)
  ! gaining LDOP (gP/m^3/day)
  ! HDRRPOP, &      !Hydrolysis loss rate of refractory POP (positive losing RPOP)
  ! gaining LDOC (gP/m^3/day)
  ! DENNO3,  &      !Denitrification loss rate of NO3 (gN/m^3/day)
  ! Positive losing NO3
  ! DENIT,   &      !Denitrication loss rate ot DOC (gC/m^3/day)
  ! positive losing DOC (nitrate limited reaction)
  ! NFIX            !Nitrogen fixation rate (gN/m^3/day)
!
      Real (SP), Allocatable, Dimension (:, :) :: PSD, SAP, ALGUP, &
     & ALGRES
  !
      Real (SP), Allocatable, Dimension (:, :) :: DORALG, DOPR, DCOD, &
     & DDOC, NITRIF
  !
      Real (SP), Allocatable, Dimension (:, :) :: A_T, AP1, ABM1, APR1, &
     & AP2, ABM2, APR2, AP3, ABM3, APR3, AALGDOC, AALGPOC, ADENIT, &
     & AMNLDOC, AHDRPOC, AALGNH4, AALGNO3, AALGDON, AALGPON, ANT, &
     & ADENNO3, AMNLDON, AHDRPON, ANFIX, AALGPO4, AALGDOP, AALGPOP, &
     & AMNLDOP, AHDRPOP, APSD, ASAP, AALGUP, AALGRES, ADO, ADORALG, &
     & ADOPR, ADCOD, ADDOC, ANITRIF
!
  !
      Real (SP), Dimension (366) :: TVPR, B2GR
  !
      Real (SP), Allocatable :: QLIT (:, :), CONLIT (:, :, :)
  !
      Integer, Save :: LENGTH, IJDAY, LENGTH2, IWRIT
  !
      Integer, Save :: NAPL, NDLT, NSNP, NPLT, NTFL
  !
      Integer, Save :: NKFL, NOPL, NOINT, NMBL
  !
      Integer, Save :: KFLDP, MBLDP, INFLOW
  !
      Integer, Save :: NDIA, NRSO, NHYDF, NTVDF, KEIPTR !, NUMBOX
  !
      Integer, Save :: JCB, IT1, IT2, IT3, JT1, IVTMP, JVTMP
  !
      Integer, Save :: JCS1MAX, JCS2MAX, JCS3MAX, JDAYVB
  !
      Integer :: BAOPTR
  !
      Real (SP), Save :: TH, TMSTRT, XYDF
  !
      Real (SP), Save :: BBMI, DFEEDI, DUMMY
  !
      Real (SP), Save :: ELTMSPLT, ELTMSTFL, ELTMSKFL
  !
      Real (SP), Save :: RTI, SFATOT, SHI
      Real (SP), Save :: RATIO, RELTMS, TODS, DSSR, SAVENT, SEDENT, &
     & SFENT
  !
      Real (SP), Save :: SFATMP, ALFAS, CONC1, CONC2, CONC3, CADV, OMTH
  !
      Real (SP), Save :: BNDC, BNDN, BNDP
  !
  !
      Real (SP), Save :: DOS, DLTAV, ELTMJD, HOURVB, HMTV2
  !
      Real (SP), Save :: WCMN, WCMP, WCMC, WCMS, WQMDLV, WQMTV2
  !
      Real (SP), Save :: DLWCMC, DLWCMN, DLWCMP, DLWCSSMC, DLWCSSMN, &
     & DLWCSSMP
  !
      Real (SP), Save :: DTMP1
  !
      Real (SP), Save :: AVGINT, DTMP2, DFTMP
  !
      Real (SP), Save :: DLSEDMC, DLSEDMN, DLSEDMP, DLSEDSSMC, &
     & DLSEDSSMN, DLSEDSSMP
  !
      Real (SP), Save :: ERRWCC, ERRWCN, ERRWCP, ERRSEDC, ERRSEDN, &
     & ERRSEDP
  !
      Real (SP), Save :: ELTMMBL
  !
  !
  !***** Common declarations
  !
      Real (SP) :: TMEND, OLDNXCBC, TE, WMS
  !
      Real (SP) :: AOCR, AONT, FDOP, AREAR, BREAR, CREAR
  !
      Integer :: NAC, NB, NQF, NL, NSB, NBB, NHQF, NSQF
  !
      Real (SP) :: TRCOD
  !
      Real (SP) :: ASC1, ASC2, ASC3, FSAP, TRSUA
  !
      Real (SP) :: FR, SCTOX
      Real (SP) :: PRSP1, PRSP2, PRSP3, STF1, STF2, STF3, VMAX1, VMAX2, &
     & VMAX3 !maximum P:C ratio for algae 3? (gP/gC)
  !
      Real (SP) :: ANC1, ANC2, ANC3, FNI1, FNI2, FNI3, FNLD1, FNLD2, &
     & FNLD3, FNRD1, FNRD2, FNRD3, FNLP1, FNLP2, FNLP3, FNRP1, FNRP2, &
     & FNRP3, FNIP, FNLDP, FNRDP, FNLPP, FNRPP, FNUP
  !
      Real (SP) :: ALPHMIN1, ALPHMIN2, ALPHMIN3, ALPHRAT1, ALPHRAT2, &
     & ALPHRAT3, CHLCMN1, CHLCMN2, CHLCMN3, ACHLC1, ACHLC2, ACHLC3, &
     & BCHLC1, BCHLC2, BCHLC3, CCHLC1, CCHLC2, CCHLC3
  !
      Real (SP) :: ANDC
      Real (SP) :: Q01, Q02, Q03 !Previously called APC3, Algal group 3 Phosphorus-to-carbon ratio (gP/gC)
  !
      Real (SP) :: FPI1, FPI2, FPI3, FPLD1, FPLD2, FPLD3, FPRD1, FPRD2, &
     & FPRD3, FPLP1, FPLP2, FPLP3, FPRP1, FPRP2, FPRP3, FPIP, FPLDP, &
     & FPRDP, FPLPP, FPRPP
  !
      Real (SP) :: FCLD1, FCLD2, FCLD3, FCRD1, FCRD2, FCRD3, FCLP1, &
     & FCLP2, FCLP3, FCRP1, FCRP2, FCRP3, FCLDP, FCRDP, FCLPP, FCRPP, &
     & TRMNL, TRHDR, AANOX
  !
      Real (SP) :: FD, TTSS, SREFLECT, PRECIP
  !
      Real (SP) :: ATMNH4, ATMNO3, ATMLDON, ATMRDON, ATMPO4, ATMLDOP, &
     & ATMRDOP
  !
      Real (SP) :: REDS1C, REDS1N, REDS1P, REDS2C, REDS2N, REDS2P, &
     & REDCBC, REDCBN, REDCBP, REDS3C, REDS3N, REDS3P
  !
      Real (SP) :: PS1, PS2, PS3
  !
      Real (SP) :: Z1TP, Z2TP
  !
      Real (SP) :: DLT, AHMDLT, FILGTH, ZDFMUL, ZDFBCK
  !
  !  moved to mod_zoop.F
  !INTEGER  ::     KTBSZ,        &
  !                KTBLZ
      Integer :: NHMDLT, NWQMR, NHMR, NIT !count of integration time steps
  !
  !
  !Added the following for debugging algae production and mortality
  !
      Real (SP), Allocatable :: P1_GL (:, :), P2_GL (:, :), BM1_GL (:, &
     & :), BM2_GL (:, :), NL1_GL (:, :), NL2_GL (:, :), PL1_GL (:, :), &
     & PL2_GL (:, :), FI1_GL (:, :), FI2_GL (:, :), B1SZ_GL (:, :), &
     & B1LZ_GL (:, :), B2SZ_GL (:, :), B2LZ_GL (:, :), PR1_GL (:, :), &
     & PR2_GL (:, :), IAVG_GL (:, :), IK1_GL (:, :), IK2_GL (:, :), IK1 &
     & (:, :), IK2 (:, :)! local IK for alg 2
  !
  !Added the following for debugging DOXG source and sinks in water column
  !
      Real (SP), Allocatable, Dimension (:, :) :: DDOC_GL, DCOD_GL, &
     & NITRIF_GL, DOSZ_GL, DOLZ_GL, DOPR_GL, DORALG_GL, RESP_GL, &
     & FTCOD_GL, COD_GL, REAERDO_GL, REAERDO !local variable of DOXG areation
!
  ! DDOC_GL,       &   !oxygen loss due to remin of DOC (LDOC and RDOC)          (mgO2/L/day)  (sink)
  ! DCOD_GL,       &   !oxygen loss due to oxidation of COD in water column      (mgO2/L/day)  (sink)
  ! NITRIF_GL,     &   !oxygen loss due to nitrification                         (mgO2/L/day)  (sink)
  ! DOSZ_GL,       &   !oxygen loss due to microzooplankton respiration          (mgO2/L/day)  (sink)
  ! DOLZ_GL,       &   !oxygen loss due to macrozooplankton respiration          (mgO2/L/day)  (sink)
  ! DOPR_GL,       &   !sum of DOP1 and DOP2 (total loss due to predation)       (mgO2/L/day)  (sink)
  ! DORALG_GL,     &   !sum of DOR1 and DOR2 (total production)                  (mgO2/L/day)  (source)
  ! RESP_GL,       &   !total oxygen loss due to respiration (ALGDO, DOPR, DDOC, DCOD, NITRIF, DOSZ, DOLZ)  (mgO2/L/day)
  ! FTCOD_GL,      &   !                  where ALGDO is oxygen loss due to algal resipration and basal metabolism  (sink)
  ! COD_GL,        &   !COD decay rate after applying temperature control (m/day)
  ! REAERDO_GL,     &      !COD (mgO2-equiv/L)!areation of DOXG (mg-O2/L/day)  (global)
  ! REAERDO            !local variable of DOXG areation
  !
  !********************************************************************************
  !**                           Sediment Model Setup                             **
  !********************************************************************************
!
      Real (SP) :: S1FLXN, S1FLXP, S1FLXC, S3FLXN, S3FLXP, S3FLXC, &
     & S2FLXN, S2FLXP, S2FLXC !source two load C (gC/day)
  !
      Real (SP) :: ATMFLXN, ATMFLXP !Atmospheric phosphorus flux  (kgP) in DLT (sec) time
  !
      Real (SP) :: BENFLXPN, BENFLXPP, BENFLXPC, BENFLXDN, BENFLXDP, &
     & DLWCKMN, DLWCKMC, BNDFLXN, BNDFLXP, BNDFLXC, DLSEDKN, DLSEDKC, &
     & BURIALFLXN, BURIALFLXP, BURIALFLXC
  !  BENFLXPN,   &  !benthic particulate organic nitrogen flux  (kgN) positive into water column, leaving sediments
  ! BENFLXPP,   &  !benthic particulate organic phosphorus flux (kgP)
  ! BENFLXPC,   &  !benthic particulate organic Carbon flx (kgC)
  ! BENFLXDN,   &  !benthic dissolved nitrogen (NH4+NO3) flux (kgN)
  ! BENFLXDP,   &  !benthic dissolved phosphorus (PO4) flux (kgP)
  ! DLWCKMN,    &  !nitogen change (kgN) due to denitrification kinetics in water column
  ! DLWCKMC,    &  !carbon change (kgC) due to kinetics in water column
  ! BNDFLXN,    &  !latetal boundary nitrogen flux due to flow (advection) and diffusion
  ! in water column  (kgN)
  ! BNDFLXP,    &  !lateral boundary phosphorus flux due to flow and diffusion
  ! in water colum (kgP)
  ! BNDFLXC,    &  !lateral boundary carbon flux due to flow (advction)
  ! and diffusion  in water column (kgC)
  ! DLSEDKN,    &  !nitrogen change (kgN) due to denitrification in
  ! sediments (layer 1 and layer 2)
  ! DLSEDKC,    &  !carbon change (kgC) due to sediment carbon
  ! diagenesis
  ! BURIALFLXN, &  !nitrogen change (kgN) due to sediment burial flux
  ! positive into sediments
  ! BURIALFLXP, &  !phosphorus change (kgP) due to sediment burial flux
  ! positive into sediments
  ! BURIALFLXC     !carbon change (kgC) due to sediment burial flux positive into sediments
  !
      Real (SP) :: IWCMS, IWCMN, IWCMP, IWCMC !total initial carbon in water column
  !
Contains
  !
  !Subroutine WQM_ALLOC()
  !Subroutine WQM_DEALLOC()
  !
      Subroutine WQM_ALLOC
    !
         Implicit None
    !
         Allocate (T(0:MTLOC, KBM1))
         T = 0.0
         Allocate (SSI(0:MTLOC, KBM1))
         SSI = 0.0
         Allocate (B1(0:MTLOC, KBM1))
         B1 = 0.0
         Allocate (B2(0:MTLOC, KBM1))
         B2 = 0.0
         Allocate (B3(0:MTLOC, KBM1))
         B3 = 0.0
         Allocate (LDOC(0:MTLOC, KBM1))
         LDOC = 0.0
         Allocate (LPOC(0:MTLOC, KBM1))
         LPOC = 0.0
         Allocate (RPOC(0:MTLOC, KBM1))
         RPOC = 0.0
         Allocate (NH4(0:MTLOC, KBM1))
         NH4 = 0.0
         Allocate (NO3(0:MTLOC, KBM1))
         NO3 = 0.0
         Allocate (LDON(0:MTLOC, KBM1))
         LDON = 0.0
         Allocate (LPON(0:MTLOC, KBM1))
         LPON = 0.0
         Allocate (RPON(0:MTLOC, KBM1))
         RPON = 0.0
         Allocate (PO4(0:MTLOC, KBM1))
         PO4 = 0.0
         Allocate (LDOP(0:MTLOC, KBM1))
         LDOP = 0.0
         Allocate (LPOP(0:MTLOC, KBM1))
         LPOP = 0.0
         Allocate (RPOP(0:MTLOC, KBM1))
         RPOP = 0.0
         Allocate (COD(0:MTLOC, KBM1))
         COD = 0.0
         Allocate (DOXG(0:MTLOC, KBM1))
         DOXG = 0.0
         Allocate (SIUPB(0:MTLOC, KBM1))
         SIUPB = 0.0
         Allocate (SIAT(0:MTLOC, KBM1))
         SIAT = 0.0
         Allocate (SALT(0:MTLOC, KBM1))
         SALT = 0.0
    !
         Allocate (RDOC(0:MTLOC, KBM1))
         RDOC = 0.0
         Allocate (RDON(0:MTLOC, KBM1))
         RDON = 0.0
         Allocate (RDOP(0:MTLOC, KBM1))
         RDOP = 0.0
         Allocate (UREA(0:MTLOC, KBM1))
         UREA = 0.0
         Allocate (PIP(0:MTLOC, KBM1))
         PIP = 0.0
         Allocate (PIB1(0:MTLOC, KBM1))
         PIB1 = 0.0
         Allocate (PIB2(0:MTLOC, KBM1))
         PIB2 = 0.0
         Allocate (PIB3(0:MTLOC, KBM1))
         PIB3 = 0.0
!
         Allocate (TDIC(0:MTLOC, KBM1))
         TDIC = 0.0 !
         Allocate (TALK(0:MTLOC, KBM1))
         TALK = 0.0 !
         Allocate (pH(0:MTLOC, KBM1))
         pH = 8.0 ! hopefully this will work to initialize the first pH guess as 8
         Allocate (pCO2(0:MTLOC, KBM1))
         pCO2 = 0.0 !
!
         Allocate (pH_GL(MGL, KBM1))
         pH_GL = 8.0 !
         Allocate (pCO2_GL(MGL, KBM1))
         pCO2_GL = 0.0 !
!
    !
         Allocate (PN1(MTLOC, KBM1))
         PN1 = 0.0
         Allocate (PN2(MTLOC, KBM1))
         PN2 = 0.0
         Allocate (PN3(MTLOC, KBM1))
         PN3 = 0.0
         Allocate (BM1(MTLOC, KBM1))
         BM1 = 0.0
         Allocate (BM2(MTLOC, KBM1))
         BM2 = 0.0
         Allocate (BM3(MTLOC, KBM1))
         BM3 = 0.0
         Allocate (MNLLDOC(MTLOC, KBM1))
         MNLLDOC = 0.0
         Allocate (FTCOD(MTLOC, KBM1))
         FTCOD = 0.0
         Allocate (NPP(MTLOC, KBM1))
         NPP = 0.0
         Allocate (P1(MTLOC, KBM1))
         P1 = 0.0
         Allocate (P2(MTLOC, KBM1))
         P2 = 0.0
         Allocate (P3(MTLOC, KBM1))
         P3 = 0.0
         Allocate (NT(MTLOC, KBM1))
         NT = 0.0
         Allocate (PR1(MTLOC, KBM1))
         PR1 = 0.0
         Allocate (PR2(MTLOC, KBM1))
         PR2 = 0.0
         Allocate (PR3(MTLOC, KBM1))
         PR3 = 0.0
         Allocate (DENIT(MTLOC, KBM1))
         DENIT = 0.0
         Allocate (RATOX(MTLOC, KBM1))
         RATOX = 0.0
         Allocate (GPP(MTLOC, KBM1))
         GPP = 0.0
         Allocate (MNLRDOC(MTLOC, KBM1))
         MNLRDOC = 0.0
         Allocate (P1NNF(MTLOC, KBM1))
         P1NNF = 0.0
    !
         Allocate (ASRAT(MTLOC))
         ASRAT = 0.0
         Allocate (total_netPP(0:MTLOC))
         total_netPP = 0.0 !: 
    !
         Allocate (KE(MTLOC, KBM1))
         KE = 0.0
         Allocate (FI1(MTLOC, KBM1))
         FI1 = 0.0
         Allocate (FI2(MTLOC, KBM1))
         FI2 = 0.0
         Allocate (FI3(MTLOC, KBM1))
         FI3 = 0.0
         Allocate (NL1(MTLOC, KBM1))
         NL1 = 0.0
         Allocate (NL2(MTLOC, KBM1))
         NL2 = 0.0
         Allocate (NL3(MTLOC, KBM1))
         NL3 = 0.0
         Allocate (PL1(MTLOC, KBM1))
         PL1 = 0.0
         Allocate (PL2(MTLOC, KBM1))
         PL2 = 0.0
         Allocate (PL3(MTLOC, KBM1))
         PL3 = 0.0
    !
         Allocate (RESP(MTLOC, KBM1))
         RESP = 0.0
    !
         Allocate (SL1(MTLOC, KBM1))
         SL1 = 0.0
         Allocate (SL2(MTLOC, KBM1))
         SL2 = 0.0
         Allocate (SL3(MTLOC, KBM1))
         SL3 = 0.0
    !
         Allocate (KEISS(MTLOC, KBM1))
         KEISS = 0.0
         Allocate (KEVSS(MTLOC, KBM1))
         KEVSS = 0.0
         Allocate (KEDOC(MTLOC, KBM1))
         KEDOC = 0.0
    !
         Allocate (PM1(MTLOC, KBM1))
         PM1 = 0.0
         Allocate (PM2(MTLOC, KBM1))
         PM2 = 0.0
         Allocate (PM3(MTLOC, KBM1))
         PM3 = 0.0
    !
         Allocate (BMR1(MTLOC, KBM1))
         BMR1 = 0.0
         Allocate (BMR2(MTLOC, KBM1))
         BMR2 = 0.0
         Allocate (BMR3(MTLOC, KBM1))
         BMR3 = 0.0
    !
         Allocate (BPR1(MTLOC, KBM1))
         BPR1 = 0.0
         Allocate (BPR2(MTLOC, KBM1))
         BPR2 = 0.0
         Allocate (BPR3(MTLOC, KBM1))
         BPR3 = 0.0
    !
         Allocate (DEPTHR(MTLOC, KBM1))
         DEPTHR = 0.0 !Never used
         Allocate (CCHL1(MTLOC, KBM1))
         CCHL1 = 0.0
         Allocate (CCHL2(MTLOC, KBM1))
         CCHL2 = 0.0
         Allocate (CCHL3(MTLOC, KBM1))
         CCHL3 = 0.0
         Allocate (Q1(MTLOC, KBM1))
         Q1 = 0.0
         Allocate (Q2(MTLOC, KBM1))
         Q2 = 0.0
         Allocate (Q3(MTLOC, KBM1))
         Q3 = 0.0
    !
    !           ALLOCATE(DTUREA(MTLOC,KBM1));       DTUREA = 0.0 !Never used
    !
    !
         Allocate (CFIX(MTLOC))
         CFIX = 0.0
         Allocate (SNFIX(MTLOC))
         SNFIX = 0.0
    !
         Allocate (WSS(0:MTLOC, KBM1))
         WSS = 0.0
         Allocate (WSL(0:MTLOC, KBM1))
         WSL = 0.0
         Allocate (WSR(0:MTLOC, KBM1))
         WSR = 0.0
         Allocate (WS1(0:MTLOC, KBM1))
         WS1 = 0.0
         Allocate (WS2(0:MTLOC, KBM1))
         WS2 = 0.0
         Allocate (WS3(0:MTLOC, KBM1))
         WS3 = 0.0
         Allocate (WSU(0:MTLOC, KBM1))
         WSU = 0.0
         Allocate (WSSHI(0:MTLOC, KBM1))
         WSSHI = 0.0
    !
         Allocate (FLXSPIB1(MTLOC, KBM1))
         FLXSPIB1 = 0.0
         Allocate (FLXSPIB2(MTLOC, KBM1))
         FLXSPIB2 = 0.0
         Allocate (FLXSPIB3(MTLOC, KBM1))
         FLXSPIB3 = 0.0
!
         Allocate (V2(0:MTLOC, KBM1))
         V2 = 0.0
!
         Allocate (WSSNET(MTLOC))
         WSSNET = 0.0
         Allocate (WSLNET(MTLOC))
         WSLNET = 0.0
         Allocate (WSRNET(MTLOC))
         WSRNET = 0.0
         Allocate (WS1NET(MTLOC))
         WS1NET = 0.0
         Allocate (WS2NET(MTLOC))
         WS2NET = 0.0
         Allocate (WS3NET(MTLOC))
         WS3NET = 0.0
    !
    !
         Allocate (WSUNET(MTLOC))
         WSUNET = 0.0
    !
    !
         Allocate (DTM(0:MTLOC, KBM1, NCP))
         DTM = 0.0
    !
         Allocate (ADVFLUX(0:MTLOC, KBM1, NCP))
         ADVFLUX = 0.0 !
         Allocate (VDIFFLUX(0:MTLOC, KBM1, NCP))!
         VDIFFLUX = 0.0 !
         Allocate (ADVFLUX_GL(MGL, KBM1, NCP))!
         ADVFLUX_GL = 0.0 !
         Allocate (VDIFFLUX_GL(MGL, KBM1, NCP))!
         VDIFFLUX_GL = 0.0 !
    !
         Allocate (ALGDOC(MTLOC, KBM1))
         ALGDOC = 0.0
         Allocate (ALGPOC(MTLOC, KBM1))
         ALGPOC = 0.0
         Allocate (HDRLPOC(MTLOC, KBM1))
         HDRLPOC = 0.0
         Allocate (HDRRPOC(MTLOC, KBM1))
         HDRRPOC = 0.0
    !
         Allocate (KLDC(MTLOC, KBM1))
         KLDC = 0.0
         Allocate (KRDC(MTLOC, KBM1))
         KRDC = 0.0
         Allocate (KDCALG(MTLOC, KBM1))
         KDCALG = 0.0
         Allocate (KLPC(MTLOC, KBM1))
         KLPC = 0.0
         Allocate (KLCALG(MTLOC, KBM1))
         KLCALG = 0.0
         Allocate (KRPC(MTLOC, KBM1))
         KRPC = 0.0
         Allocate (KSUA(MTLOC, KBM1))
         KSUA = 0.0
         Allocate (KCOD(MTLOC, KBM1))
         KCOD = 0.0
         Allocate (KLDN(MTLOC, KBM1))
         KLDN = 0.0
         Allocate (KRDN(MTLOC, KBM1))
         KRDN = 0.0
         Allocate (KDNALG(MTLOC, KBM1))
         KDNALG = 0.0
         Allocate (KLDP(MTLOC, KBM1))
         KLDP = 0.0
         Allocate (KRDP(MTLOC, KBM1))
         KRDP = 0.0
         Allocate (KDPALG(MTLOC, KBM1))
         KDPALG = 0.0
         Allocate (KLPN(MTLOC, KBM1))
         KLPN = 0.0
         Allocate (KLNALG(MTLOC, KBM1))
         KLNALG = 0.0
         Allocate (KLPP(MTLOC, KBM1))
         KLPP = 0.0
         Allocate (KLPALG(MTLOC, KBM1))
         KLPALG = 0.0
         Allocate (KRPN(MTLOC, KBM1))
         KRPN = 0.0
         Allocate (KRPP(MTLOC, KBM1))
         KRPP = 0.0
         Allocate (NTM(MTLOC, KBM1))
         NTM = 0.0
         Allocate (KRCOAG(MTLOC, KBM1))
         KRCOAG = 0.0
    !
         Allocate (FTMNL(MTLOC, KBM1))
         FTMNL = 0.0
         Allocate (FTHDR(MTLOC, KBM1))
         FTHDR = 0.0
         Allocate (RESPC(MTLOC, KBM1))
         RESPC = 0.0
         Allocate (DLALGC(MTLOC, KBM1))
         DLALGC = 0.0
    !
         Allocate (ALGNH4(MTLOC, KBM1))
         ALGNH4 = 0.0
         Allocate (ALGNO3(MTLOC, KBM1))
         ALGNO3 = 0.0
         Allocate (ALGDON(MTLOC, KBM1))
         ALGDON = 0.0
         Allocate (ALGPON(MTLOC, KBM1))
         ALGPON = 0.0
         Allocate (DENNO3(MTLOC, KBM1))
         DENNO3 = 0.0
         Allocate (MNLLDON(MTLOC, KBM1))
         MNLLDON = 0.0
         Allocate (HDRLPON(MTLOC, KBM1))
         HDRLPON = 0.0
         Allocate (HDRRPON(MTLOC, KBM1))
         HDRRPON = 0.0
         Allocate (MNLRDON(MTLOC, KBM1))
         MNLRDON = 0.0
         Allocate (NFIX(MTLOC, KBM1))
         NFIX = 0.0
    !
         Allocate (ALGPO4(MTLOC, KBM1))
         ALGPO4 = 0.0
         Allocate (ALGDOP(MTLOC, KBM1))
         ALGDOP = 0.0
         Allocate (ALGPOP(MTLOC, KBM1))
         ALGPOP = 0.0
         Allocate (MNLLDOP(MTLOC, KBM1))
         MNLLDOP = 0.0
         Allocate (HDRLPOP(MTLOC, KBM1))
         HDRLPOP = 0.0
         Allocate (HDRRPOP(MTLOC, KBM1))
         HDRRPOP = 0.0
         Allocate (MNLRDOP(MTLOC, KBM1))
         MNLRDOP = 0.0
    !
         Allocate (DORALG(MTLOC, KBM1))
         DORALG = 0.0
         Allocate (DOPR(MTLOC, KBM1))
         DOPR = 0.0
         Allocate (DCOD(MTLOC, KBM1))
         DCOD = 0.0
         Allocate (DDOC(MTLOC, KBM1))
         DDOC = 0.0
         Allocate (NITRIF(MTLOC, KBM1))
         NITRIF = 0.0
    !
         Allocate (PSD(MTLOC, KBM1))
         PSD = 0.0
         Allocate (SAP(MTLOC, KBM1))
         SAP = 0.0
         Allocate (ALGUP(MTLOC, KBM1))
         ALGUP = 0.0
         Allocate (ALGRES(MTLOC, KBM1))
         ALGRES = 0.0
    !
         Allocate (BENDOC(MTLOC))
         BENDOC = 0.0
         Allocate (BENNH4(MTLOC))
         BENNH4 = 0.0
         Allocate (BENNO3(MTLOC))
         BENNO3 = 0.0
         Allocate (BENDON(MTLOC))
         BENDON = 0.0
         Allocate (BENDOP(MTLOC))
         BENDOP = 0.0
         Allocate (BENPO4(MTLOC))
         BENPO4 = 0.0
         Allocate (BENCOD(MTLOC))
         BENCOD = 0.0
         Allocate (BENDO(MTLOC))
         BENDO = 0.0
         Allocate (BENSA(MTLOC))
         BENSA = 0.0
         Allocate (BENDEN(MTLOC))
         BENDEN = 0.0
         Allocate (BENCH4G(MTLOC))
         BENCH4G = 0.0
         Allocate (BENCH4A(MTLOC))
         BENCH4A = 0.0
    !
#if defined(CSOD_CALC)
         Allocate (CSODflux(MTLOC))
         CSODflux = 0.0 !
#endif
    !
         Allocate (C1(0:MTLOC, KBM1, NCP))
         C1 = 0.0
         Allocate (C2(0:MTLOC, KBM1, NCP))
         C2 = 0.0
         Allocate (C2F(0:MTLOC, KBM1, NCP))
         C2F = 0.0
         Allocate (CSTAR(0:MTLOC, KBM1, NCP))
         CSTAR = 0.0
         Allocate (DTC(0:MTLOC, KBM1, NCP))
         DTC = 0.0
         Allocate (AC1(0:MTLOC, KBM1, NCP))
         AC1 = 0.0
    !
         Allocate (QLIT(MTLOC, KBM1))
         QLIT = 0.0
         Allocate (CONLIT(MTLOC, KBM1, NCP))
         CONLIT = 0.0
    !
         Allocate (AFI1(MTLOC, KBM1))
         AFI1 = 0.0
         Allocate (ANL1(MTLOC, KBM1))
         ANL1 = 0.0
         Allocate (APL1(MTLOC, KBM1))
         APL1 = 0.0
         Allocate (AFI2(MTLOC, KBM1))
         AFI2 = 0.0
         Allocate (ANL2(MTLOC, KBM1))
         ANL2 = 0.0
         Allocate (APL2(MTLOC, KBM1))
         APL2 = 0.0
         Allocate (ASL2(MTLOC, KBM1))
         ASL2 = 0.0
         Allocate (AFI3(MTLOC, KBM1))
         AFI3 = 0.0
         Allocate (ANL3(MTLOC, KBM1))
         ANL3 = 0.0
         Allocate (APL3(MTLOC, KBM1))
         APL3 = 0.0
         Allocate (ANPP(MTLOC, KBM1))
         ANPP = 0.0
         Allocate (ARESP(MTLOC, KBM1))
         ARESP = 0.0
         Allocate (AKE(MTLOC, KBM1))
         AKE = 0.0
         Allocate (ASL1(MTLOC, KBM1))
         ASL1 = 0.0
         Allocate (ASL3(MTLOC, KBM1))
         ASL3 = 0.0
         Allocate (AGPP(MTLOC, KBM1))
         AGPP = 0.0
    !
         Allocate (ACCHL1(MTLOC, KBM1))
         ACCHL1 = 0.0
         Allocate (ACCHL2(MTLOC, KBM1))
         ACCHL2 = 0.0
         Allocate (ACCHL3(MTLOC, KBM1))
         ACCHL3 = 0.0
    !
    !
         Allocate (ABENDOC(MTLOC))
         ABENDOC = 0.0
         Allocate (ABENNH4(MTLOC))
         ABENNH4 = 0.0
         Allocate (ABENNO3(MTLOC))
         ABENNO3 = 0.0
         Allocate (ABENPO4(MTLOC))
         ABENPO4 = 0.0
         Allocate (ABENCOD(MTLOC))
         ABENCOD = 0.0
         Allocate (ABENDO(MTLOC))
         ABENDO = 0.0
         Allocate (ABENSA(MTLOC))
         ABENSA = 0.0
         Allocate (ASSFWS(MTLOC))
         ASSFWS = 0.0
         Allocate (APCFWS(MTLOC))
         APCFWS = 0.0
         Allocate (APNFWS(MTLOC))
         APNFWS = 0.0
         Allocate (APPFWS(MTLOC))
         APPFWS = 0.0
         Allocate (APSFWS(MTLOC))
         APSFWS = 0.0
         Allocate (ACPIP(MTLOC))
         ACPIP = 0.0
         Allocate (ACPOS(MTLOC))
         ACPOS = 0.0
         Allocate (ADFEED(MTLOC))
         ADFEED = 0.0
         Allocate (ABENCH4G(MTLOC))
         ABENCH4G = 0.0
         Allocate (ABENCH4A(MTLOC))
         ABENCH4A = 0.0
    !
         Allocate (AFIB(MTLOC))
         AFIB = 0.0
         Allocate (ANLB(MTLOC))
         ANLB = 0.0
         Allocate (APLB(MTLOC))
         APLB = 0.0
         Allocate (ANPPB(MTLOC))
         ANPPB = 0.0
         Allocate (ABBM(MTLOC))
         ABBM = 0.0
         Allocate (ABLITE(MTLOC))
         ABLITE = 0.0
         Allocate (ABMB(MTLOC))
         ABMB = 0.0
         Allocate (APB(MTLOC))
         APB = 0.0
         Allocate (APRB(MTLOC))
         APRB = 0.0
         Allocate (ABADOC(MTLOC))
         ABADOC = 0.0
         Allocate (ABAPOC(MTLOC))
         ABAPOC = 0.0
         Allocate (ABANH4(MTLOC))
         ABANH4 = 0.0
         Allocate (ABANO3(MTLOC))
         ABANO3 = 0.0
         Allocate (ABAPON(MTLOC))
         ABAPON = 0.0
         Allocate (ABAPO4(MTLOC))
         ABAPO4 = 0.0
         Allocate (ABAPOP(MTLOC))
         ABAPOP = 0.0
         Allocate (ABADO(MTLOC))
         ABADO = 0.0
    !
         Allocate (ACPOC(MTLOC, 3))
         ACPOC = 0.0
         Allocate (ACPON(MTLOC, 3))
         ACPON = 0.0
         Allocate (ACPOP(MTLOC, 3))
         ACPOP = 0.0
         Allocate (ASFEED(MTLOC, NSSFP))
         ASFEED = 0.0
    !
         Allocate (AJNSF(MTLOC))
         AJNSF = 0.0
         Allocate (AJPSF(MTLOC))
         AJPSF = 0.0
         Allocate (ASODSF(MTLOC))
         ASODSF = 0.0
         Allocate (ASASF(MTLOC))
         ASASF = 0.0
         Allocate (ASUSF(MTLOC))
         ASUSF = 0.0
         Allocate (ASFGCIN(MTLOC))
         ASFGCIN = 0.0
         Allocate (ASFCFEC(MTLOC))
         ASFCFEC = 0.0
         Allocate (ASFCPSF(MTLOC))
         ASFCPSF = 0.0
         Allocate (AFLXCSF(MTLOC))
         AFLXCSF = 0.0
         Allocate (AFLXNSF(MTLOC))
         AFLXNSF = 0.0
         Allocate (AFLXPSF(MTLOC))
         AFLXPSF = 0.0
         Allocate (ARPOCSF(MTLOC))
         ARPOCSF = 0.0
         Allocate (ARPONSF(MTLOC))
         ARPONSF = 0.0
         Allocate (ARPOPSF(MTLOC))
         ARPOPSF = 0.0
         Allocate (ASSISF(MTLOC))
         ASSISF = 0.0
         Allocate (ASSISASF(MTLOC))
         ASSISASF = 0.0
         Allocate (ASSISUSF(MTLOC))
         ASSISUSF = 0.0
         Allocate (ASSIPSF(MTLOC))
         ASSIPSF = 0.0
    !
         Allocate (A_T(MTLOC, KBM1))
         A_T = 0.0
         Allocate (AP1(MTLOC, KBM1))
         AP1 = 0.0
         Allocate (ABM1(MTLOC, KBM1))
         ABM1 = 0.0
         Allocate (APR1(MTLOC, KBM1))
         APR1 = 0.0
         Allocate (AP2(MTLOC, KBM1))
         AP2 = 0.0
         Allocate (ABM2(MTLOC, KBM1))
         ABM2 = 0.0
         Allocate (APR2(MTLOC, KBM1))
         APR2 = 0.0
         Allocate (AP3(MTLOC, KBM1))
         AP3 = 0.0
         Allocate (ABM3(MTLOC, KBM1))
         ABM3 = 0.0
         Allocate (APR3(MTLOC, KBM1))
         APR3 = 0.0
         Allocate (AALGDOC(MTLOC, KBM1))
         AALGDOC = 0.0
         Allocate (AALGPOC(MTLOC, KBM1))
         AALGPOC = 0.0
         Allocate (ADENIT(MTLOC, KBM1))
         ADENIT = 0.0
         Allocate (AMNLDOC(MTLOC, KBM1))
         AMNLDOC = 0.0
         Allocate (AHDRPOC(MTLOC, KBM1))
         AHDRPOC = 0.0
         Allocate (AALGNH4(MTLOC, KBM1))
         AALGNH4 = 0.0
         Allocate (AALGNO3(MTLOC, KBM1))
         AALGNO3 = 0.0
         Allocate (AALGDON(MTLOC, KBM1))
         AALGDON = 0.0
         Allocate (AALGPON(MTLOC, KBM1))
         AALGPON = 0.0
         Allocate (ANT(MTLOC, KBM1))
         ANT = 0.0
         Allocate (ADENNO3(MTLOC, KBM1))
         ADENNO3 = 0.0
         Allocate (AMNLDON(MTLOC, KBM1))
         AMNLDON = 0.0
         Allocate (AHDRPON(MTLOC, KBM1))
         AHDRPON = 0.0
    !
    !
         Allocate (AALGPO4(MTLOC, KBM1))
         AALGPO4 = 0.0
         Allocate (AALGDOP(MTLOC, KBM1))
         AALGDOP = 0.0
         Allocate (AALGPOP(MTLOC, KBM1))
         AALGPOP = 0.0
         Allocate (AMNLDOP(MTLOC, KBM1))
         AMNLDOP = 0.0
         Allocate (AHDRPOP(MTLOC, KBM1))
         AHDRPOP = 0.0
         Allocate (APSD(MTLOC, KBM1))
         APSD = 0.0
         Allocate (ASAP(MTLOC, KBM1))
         ASAP = 0.0
         Allocate (AALGUP(MTLOC, KBM1))
         AALGUP = 0.0
         Allocate (AALGRES(MTLOC, KBM1))
         AALGRES = 0.0
         Allocate (ADO(MTLOC, KBM1))
         ADO = 0.0
         Allocate (ADORALG(MTLOC, KBM1))
         ADORALG = 0.0
         Allocate (ADOPR(MTLOC, KBM1))
         ADOPR = 0.0
         Allocate (ADCOD(MTLOC, KBM1))
         ADCOD = 0.0
         Allocate (ADDOC(MTLOC, KBM1))
         ADDOC = 0.0
         Allocate (ANITRIF(MTLOC, KBM1))
         ANITRIF = 0.0
    !
!
         Allocate (SBN(MTLOC))
         SBN = 0
         Allocate (BBN(MTLOC))
         BBN = 0
         Allocate (HMSBV(MTLOC))
         HMSBV = 0
    !    ALLOCATE(SFA(MTLOC));          SFA = 0   !SFA, replaced by ART1
    !
         Allocate (FLUXS(0:MTLOC, KBM1, 13))
         FLUXS = 0.0
    !
         Allocate (IWCMNB(MTLOC, KBM1))
         IWCMNB = 0.0
         Allocate (IWCMPB(MTLOC, KBM1))
         IWCMPB = 0.0
         Allocate (IWCMCB(MTLOC, KBM1))
         IWCMCB = 0.0
         Allocate (IWCMSB(MTLOC, KBM1))
         IWCMSB = 0.0
         Allocate (WCMNB(MTLOC, KBM1))
         WCMNB = 0.0
         Allocate (WCMPB(MTLOC, KBM1))
         WCMPB = 0.0
         Allocate (WCMCB(MTLOC, KBM1))
         WCMCB = 0.0
         Allocate (WCMSB(MTLOC, KBM1))
         WCMSB = 0.0
         Allocate (DLWCMNB(MTLOC, KBM1))
         DLWCMNB = 0.0
         Allocate (DLWCMPB(MTLOC, KBM1))
         DLWCMPB = 0.0
         Allocate (DLWCKMNB(MTLOC, KBM1))
         DLWCKMNB = 0.0
         Allocate (DLWCKMCB(MTLOC, KBM1))
         DLWCKMCB = 0.0
         Allocate (S1FLXNB(MTLOC, KBM1))
         S1FLXNB = 0.0
         Allocate (S1FLXPB(MTLOC, KBM1))
         S1FLXPB = 0.0
         Allocate (S1FLXCB(MTLOC, KBM1))
         S1FLXCB = 0.0
         Allocate (S2FLXNB(MTLOC, KBM1))
         S2FLXNB = 0.0
         Allocate (S2FLXPB(MTLOC, KBM1))
         S2FLXPB = 0.0
         Allocate (S2FLXCB(MTLOC, KBM1))
         S2FLXCB = 0.0
         Allocate (S3FLXNB(MTLOC, KBM1))
         S3FLXNB = 0.0
         Allocate (S3FLXPB(MTLOC, KBM1))
         S3FLXPB = 0.0
         Allocate (S3FLXCB(MTLOC, KBM1))
         S3FLXCB = 0.0
    !
         Allocate (CTEMP(MTLOC))
         CTEMP = 0.0
    !
         Allocate (BSVOL(MTLOC))
         BSVOL = 0.0
    !
         Allocate (PCFWS(MTLOC))
         PCFWS = 0.0
         Allocate (PNFWS(MTLOC))
         PNFWS = 0.0
         Allocate (PPFWS(MTLOC))
         PPFWS = 0.0
         Allocate (PSFWS(MTLOC))
         PSFWS = 0.0
         Allocate (SSFWS(MTLOC))
         SSFWS = 0.0
    !
         Allocate (AASRAT(MTLOC))
         AASRAT = 0.0
    !
         Allocate (ACFIX(MTLOC))
         ACFIX = 0.0
         Allocate (ANFIX(MTLOC, KBM1))
         ANFIX = 0.0

         Allocate (BENDOCB(MTLOC))
         BENDOCB = 0.0
         Allocate (BENNH4B(MTLOC))
         BENNH4B = 0.0
         Allocate (BENNO3B(MTLOC))
         BENNO3B = 0.0
         Allocate (BENPO4B(MTLOC))
         BENPO4B = 0.0
         Allocate (BENCODB(MTLOC))
         BENCODB = 0.0
         Allocate (BENDOB(MTLOC))
         BENDOB = 0.0
         Allocate (BENSAB(MTLOC))
         BENSAB = 0.0
    !
         Allocate (ATMFLXNB(MTLOC))
         ATMFLXNB = 0.0
         Allocate (ATMFLXPB(MTLOC))
         ATMFLXPB = 0.0
         Allocate (ATMFLXCB(MTLOC))
         ATMFLXCB = 0.0
         Allocate (BENFLXPNB(MTLOC))
         BENFLXPNB = 0.0
         Allocate (BENFLXDNB(MTLOC))
         BENFLXDNB = 0.0
         Allocate (BENFLXPPB(MTLOC))
         BENFLXPPB = 0.0
         Allocate (BENFLXDPB(MTLOC))
         BENFLXDPB = 0.0
         Allocate (BENFLXPCB(MTLOC))
         BENFLXPCB = 0.0
         Allocate (DLSEDKNB(MTLOC))
         DLSEDKNB = 0.0
         Allocate (DLSEDKCB(MTLOC))
         DLSEDKCB = 0.0
         Allocate (BURIALFLXNB(MTLOC))
         BURIALFLXNB = 0.0
         Allocate (BURIALFLXPB(MTLOC))
         BURIALFLXPB = 0.0
         Allocate (BURIALFLXCB(MTLOC))
         BURIALFLXCB = 0.0
    !
    !
    ! Add C2_GL
         Allocate (C2_GL(MGL, KBM1, NCP))
         C2_GL = 0.0
    ! Added below for output of CCHL1_GL
         Allocate (CCHL1_GL(MGL, KBM1))
         CCHL1_GL = 0.0
         Allocate (T_GL(MGL, KBM1))
         T_GL = 0.0
         Allocate (S_GL(MGL, KBM1))
         S_GL = 0.0
         Allocate (total_netPP_GL(MGL))
         total_netPP_GL = 0.0
    !
    ! : to save GPP and NPP in history output
         Allocate (GPP_GL(MGL, KBM1))
         GPP_GL = 0.0
         Allocate (NPP_GL(MGL, KBM1))
         NPP_GL = 0.0
    !
    !Added D_GL, EL_GL, H_GL (total depth, surface elevation and bathymetric depth)
         Allocate (D_GL(MGL))
         D_GL = 0.0
         Allocate (H_GL(MGL))
         H_GL = 0.0
         Allocate (EL_GL(MGL))
         EL_GL = 0.0
    !
    !--debugging benthic fluxes ---
         Allocate (BFLUX_GL(MGL, 9))
         BFLUX_GL = 0.0
         Allocate (BFLUXB_GL(MGL, 9))
         BFLUXB_GL = 0.0
         Allocate (BFLUXNX(MTLOC, 9))
         BFLUXNX = 0.0
         Allocate (BFLUXNX_GL(MGL, 9))
         BFLUXNX_GL = 0.0
         Allocate (BFLUX(MTLOC, 9))
         BFLUX = 0.0
    !
    !--Allocating global hydro variables here ---
    !
         Allocate (UL_GL(NGL, KBM1))
         UL_GL = 0.0

         Allocate (VL_GL(NGL, KBM1))
         VL_GL = 0.0
         Allocate (WTSL_GL(MGL, KB))
         WTSL_GL = 0.0
         Allocate (KHL_GL(MGL, KB))
         KHL_GL = 0.0
         Allocate (SL_GL(MGL, KBM1))
         SL_GL = 0.0
         Allocate (TL_GL(MGL, KBM1))
         TL_GL = 0.0
         Allocate (ELL_GL(1:MGL))
         ELL_GL = 0.0
         Allocate (DTFAL_GL(MGL))
         DTFAL_GL = 0.0
         Allocate (UARD_OBCN_GL(NOBTY+1))
         UARD_OBCN_GL = 0.0
         Allocate (XFLUX_OBC_GL(NOBTY, KBM1))
         XFLUX_OBC_GL = 0.0
    !
    !
    !--Added the following algae debugging global variables here ---
    !
         Allocate (P1_GL(MGL, KBM1))
         P1_GL = 0.0
         Allocate (P2_GL(MGL, KBM1))
         P2_GL = 0.0
         Allocate (BM1_GL(MGL, KBM1))
         BM1_GL = 0.0
         Allocate (BM2_GL(MGL, KBM1))
         BM2_GL = 0.0
         Allocate (NL1_GL(MGL, KBM1))
         NL1_GL = 0.0
         Allocate (NL2_GL(MGL, KBM1))
         NL2_GL = 0.0
         Allocate (PL1_GL(MGL, KBM1))
         PL1_GL = 0.0
         Allocate (PL2_GL(MGL, KBM1))
         PL2_GL = 0.0
         Allocate (B1SZ_GL(MGL, KBM1))
         B1SZ_GL = 0.0
         Allocate (B2SZ_GL(MGL, KBM1))
         B2SZ_GL = 0.0
         Allocate (B1LZ_GL(MGL, KBM1))
         B1LZ_GL = 0.0
         Allocate (B2LZ_GL(MGL, KBM1))
         B2LZ_GL = 0.0
         Allocate (PR1_GL(MGL, KBM1))
         PR1_GL = 0.0
         Allocate (PR2_GL(MGL, KBM1))
         PR2_GL = 0.0
    !
         If (LIGHT_EXTINCTION) Then
            Allocate (IAVG_GL(MGL, KBM1))
            IAVG_GL = 0.0
            Allocate (FI1_GL(MGL, KBM1))
            FI1_GL = 0.0
            Allocate (FI2_GL(MGL, KBM1))
            FI2_GL = 0.0
         End If

		 !If (SAV_CALC) Then !Adi 4/23/18
		    Allocate (LEAF_GL(MGL,NSAVM))
            LEAF_GL = 0.0
            Allocate (STEM_GL(MGL,NSAVM))
            STEM_GL = 0.0
            Allocate (ROOT_GL(MGL,NSAVM))
            ROOT_GL = 0.0
            Allocate (TUBER_GL(MGL,NSAVM))
            TUBER_GL = 0.0
			Allocate (EP_GL(MGL,NSAVM))
            EP_GL = 0.0
		 !end if
         !
		 !If (SOLIDS_CALC) Then !Adi 2/13/18
		    Allocate (TSS_GL(MGL, KBM1))
            TSS_GL = 0.0
            Allocate (SSI_GL(MGL, KBM1))
            SSI_GL = 0.0
            Allocate (KESS_GL(MGL, KBM1))
            KESS_GL = 0.0
            Allocate (TURB_GL(MGL, KBM1))
            TURB_GL = 0.0
			Allocate (PARAD_GL(MGL, KBM1))
            PARAD_GL = 0.0
		 !end if


         Allocate (IK1_GL(MGL, KBM1))
         IK1_GL = 0.0
         Allocate (IK2_GL(MGL, KBM1))
         IK2_GL = 0.0
         Allocate (IK1(MTLOC, KBM1))
         IK1 = 0.0
         Allocate (IK2(MTLOC, KBM1))
         IK2 = 0.0
    !
    !---Added the following for debugging DOXG
    !
         Allocate (DDOC_GL(MGL, KBM1))
         DDOC_GL = 0.0
         Allocate (DCOD_GL(MGL, KBM1))
         DCOD_GL = 0.0
         Allocate (NITRIF_GL(MGL, KBM1))
         NITRIF_GL = 0.0
         Allocate (DOSZ_GL(MGL, KBM1))
         DOSZ_GL = 0.0
         Allocate (DOLZ_GL(MGL, KBM1))
         DOLZ_GL = 0.0
         Allocate (DOPR_GL(MGL, KBM1))
         DOPR_GL = 0.0
         Allocate (DORALG_GL(MGL, KBM1))
         DORALG_GL = 0.0
         Allocate (RESP_GL(MGL, KBM1))
         RESP_GL = 0.0
         Allocate (FTCOD_GL(MGL, KBM1))
         FTCOD_GL = 0.0
         Allocate (COD_GL(MGL, KBM1))
         COD_GL = 0.0
         Allocate (REAERDO_GL(MGL, KBM1))
         REAERDO_GL = 0.0
    !
         Allocate (REAERDO(MTLOC, KBM1))  !
										  !
         REAERDO = 0.0 ! local variabl
    !
    !---Declared a large buffer (so that we do not have to dynamically allocate and deallocate buffers
    !   all the time--
#if defined (MULTIPROCESSOR)
         If ( .Not. ALLOCATED(RCVBUF)) ALLOCATE (RCVBUF(Max(MGL, &
        & NGL)*3*KB))
         RCVBUF = 0.0
         If ( .Not. ALLOCATED(SNDBUF)) ALLOCATE (SNDBUF(Max(MGL, &
        & NGL)*3*KB))
         SNDBUF = 0.0
#endif
    !
         Return
    !
      End Subroutine WQM_ALLOC
  !
      Subroutine WQM_DEALLOC
    !
         Implicit None
    !
         If (ALLOCATED(T)) DEALLOCATE (T)
         If (ALLOCATED(SALT)) DEALLOCATE (SALT)
         If (ALLOCATED(SSI)) DEALLOCATE (SSI)
         If (ALLOCATED(B1)) DEALLOCATE (B1)
         If (ALLOCATED(B2)) DEALLOCATE (B2)
         If (ALLOCATED(B3)) DEALLOCATE (B3)
         If (ALLOCATED(LDOC)) DEALLOCATE (LDOC)
         If (ALLOCATED(LPOC)) DEALLOCATE (LPOC)
         If (ALLOCATED(RPOC)) DEALLOCATE (RPOC)
         If (ALLOCATED(NH4)) DEALLOCATE (NH4)
         If (ALLOCATED(NO3)) DEALLOCATE (NO3)
         If (ALLOCATED(LDON)) DEALLOCATE (LDON)
         If (ALLOCATED(LPON)) DEALLOCATE (LPON)
         If (ALLOCATED(RPON)) DEALLOCATE (RPON)
         If (ALLOCATED(PO4)) DEALLOCATE (PO4)
         If (ALLOCATED(LDOP)) DEALLOCATE (LDOP)
         If (ALLOCATED(LPOP)) DEALLOCATE (LPOP)
         If (ALLOCATED(RPOP)) DEALLOCATE (RPOP)
         If (ALLOCATED(COD)) DEALLOCATE (COD)
         If (ALLOCATED(DOXG)) DEALLOCATE (DOXG)
         If (ALLOCATED(SIUPB)) DEALLOCATE (SIUPB)
         If (ALLOCATED(SIAT)) DEALLOCATE (SIAT)
    !
    !
         If (ALLOCATED(RDOC)) DEALLOCATE (RDOC)
         If (ALLOCATED(RDON)) DEALLOCATE (RDON)
         If (ALLOCATED(RDOP)) DEALLOCATE (RDOP)
    !
         If (ALLOCATED(UREA)) DEALLOCATE (UREA)
    !
         If (ALLOCATED(PIP)) DEALLOCATE (PIP)
         If (ALLOCATED(PIB1)) DEALLOCATE (PIB1)
         If (ALLOCATED(PIB2)) DEALLOCATE (PIB2)
         If (ALLOCATED(PIB3)) DEALLOCATE (PIB3)
!
         If (ALLOCATED(TDIC)) DEALLOCATE (TDIC)!
         If (ALLOCATED(TALK)) DEALLOCATE (TALK)!
         If (ALLOCATED(pH)) DEALLOCATE (pH)!
         If (ALLOCATED(pCO2)) DEALLOCATE (pCO2)!
!
         If (ALLOCATED(pH_GL)) DEALLOCATE (pH_GL)!
         If (ALLOCATED(pCO2_GL)) DEALLOCATE (pCO2_GL)!
!
         If (ALLOCATED(PN1)) DEALLOCATE (PN1)
         If (ALLOCATED(PN2)) DEALLOCATE (PN2)
         If (ALLOCATED(PN3)) DEALLOCATE (PN3)
    !
         If (ALLOCATED(BM1)) DEALLOCATE (BM1)
         If (ALLOCATED(BM2)) DEALLOCATE (BM2)
         If (ALLOCATED(BM3)) DEALLOCATE (BM3)
    !
         If (ALLOCATED(MNLLDOC)) DEALLOCATE (MNLLDOC)
         If (ALLOCATED(MNLRDOC)) DEALLOCATE (MNLRDOC)
    !
         If (ALLOCATED(FTCOD)) DEALLOCATE (FTCOD)
    !
         If (ALLOCATED(NPP)) DEALLOCATE (NPP)
         If (ALLOCATED(P1)) DEALLOCATE (P1)
         If (ALLOCATED(P2)) DEALLOCATE (P2)
         If (ALLOCATED(P3)) DEALLOCATE (P3)
    !
         If (ALLOCATED(NT)) DEALLOCATE (NT)
         If (ALLOCATED(PR1)) DEALLOCATE (PR1)
         If (ALLOCATED(PR2)) DEALLOCATE (PR2)
         If (ALLOCATED(PR3)) DEALLOCATE (PR3)
    !
         If (ALLOCATED(DENIT)) DEALLOCATE (DENIT)
    !
         If (ALLOCATED(RATOX)) DEALLOCATE (RATOX)
         If (ALLOCATED(GPP)) DEALLOCATE (GPP)
    !
         If (ALLOCATED(P1NNF)) DEALLOCATE (P1NNF)
    !
         If (ALLOCATED(ASRAT)) DEALLOCATE (ASRAT)
         If (ALLOCATED(total_netPP)) DEALLOCATE (total_netPP)
    !
         If (ALLOCATED(KE)) DEALLOCATE (KE)
         If (ALLOCATED(KEISS)) DEALLOCATE (KEISS)
         If (ALLOCATED(KEVSS)) DEALLOCATE (KEVSS)
         If (ALLOCATED(KEDOC)) DEALLOCATE (KEDOC)
    !
         If (ALLOCATED(FI1)) DEALLOCATE (FI1)
         If (ALLOCATED(FI2)) DEALLOCATE (FI2)
         If (ALLOCATED(FI3)) DEALLOCATE (FI3)
         If (ALLOCATED(NL1)) DEALLOCATE (NL1)
         If (ALLOCATED(NL2)) DEALLOCATE (NL2)
         If (ALLOCATED(NL3)) DEALLOCATE (NL3)
         If (ALLOCATED(PL1)) DEALLOCATE (PL1)
         If (ALLOCATED(PL2)) DEALLOCATE (PL2)
         If (ALLOCATED(PL3)) DEALLOCATE (PL3)
    !
         If (ALLOCATED(RESP)) DEALLOCATE (RESP)
    !
         If (ALLOCATED(SL1)) DEALLOCATE (SL1)
         If (ALLOCATED(SL2)) DEALLOCATE (SL2)
         If (ALLOCATED(SL3)) DEALLOCATE (SL3)
    !
    !
         If (ALLOCATED(PM1)) DEALLOCATE (PM1)
         If (ALLOCATED(PM2)) DEALLOCATE (PM2)
         If (ALLOCATED(PM3)) DEALLOCATE (PM3)
    !
         If (ALLOCATED(BMR1)) DEALLOCATE (BMR1)
         If (ALLOCATED(BMR2)) DEALLOCATE (BMR2)
         If (ALLOCATED(BMR3)) DEALLOCATE (BMR3)
    !
         If (ALLOCATED(BPR1)) DEALLOCATE (BPR1)
         If (ALLOCATED(BPR2)) DEALLOCATE (BPR2)
         If (ALLOCATED(BPR3)) DEALLOCATE (BPR3)
    !
         If (ALLOCATED(DEPTHR)) DEALLOCATE (DEPTHR)
    !
         If (ALLOCATED(CCHL1)) DEALLOCATE (CCHL1)
         If (ALLOCATED(CCHL2)) DEALLOCATE (CCHL2)
         If (ALLOCATED(CCHL3)) DEALLOCATE (CCHL3)
    !
         If (ALLOCATED(Q1)) DEALLOCATE (Q1)
         If (ALLOCATED(Q2)) DEALLOCATE (Q2)
         If (ALLOCATED(Q3)) DEALLOCATE (Q3)
    !
    !          IF(ALLOCATED(DTUREA))DEALLOCATE()
    !
         If (ALLOCATED(CFIX)) DEALLOCATE (CFIX)
         If (ALLOCATED(SNFIX)) DEALLOCATE (SNFIX)
    !
         If (ALLOCATED(WSS)) DEALLOCATE (WSS)
         If (ALLOCATED(WSL)) DEALLOCATE (WSL)
         If (ALLOCATED(WSR)) DEALLOCATE (WSR)
    !
         If (ALLOCATED(WS1)) DEALLOCATE (WS1)
         If (ALLOCATED(WS2)) DEALLOCATE (WS2)
         If (ALLOCATED(WS3)) DEALLOCATE (WS3)
         If (ALLOCATED(WSU)) DEALLOCATE (WSU)
         If (ALLOCATED(WSSHI)) DEALLOCATE (WSSHI)
    !
         If (ALLOCATED(FLXSPIB1)) DEALLOCATE (FLXSPIB1)
         If (ALLOCATED(FLXSPIB2)) DEALLOCATE (FLXSPIB2)
         If (ALLOCATED(FLXSPIB3)) DEALLOCATE (FLXSPIB3)
    !
!
         If (ALLOCATED(V2)) DEALLOCATE (V2)
    !
         If (ALLOCATED(WSSNET)) DEALLOCATE (WSSNET)
         If (ALLOCATED(WSLNET)) DEALLOCATE (WSLNET)
         If (ALLOCATED(WSRNET)) DEALLOCATE (WSRNET)
         If (ALLOCATED(WS1NET)) DEALLOCATE (WS1NET)
         If (ALLOCATED(WS2NET)) DEALLOCATE (WS2NET)
         If (ALLOCATED(WS3NET)) DEALLOCATE (WS3NET)
         If (ALLOCATED(WSUNET)) DEALLOCATE (WSUNET)
    !
    !
         If (ALLOCATED(DTM)) DEALLOCATE (DTM)
    !
         If (ALLOCATED(ADVFLUX)) DEALLOCATE (ADVFLUX)!
         If (ALLOCATED(VDIFFLUX)) DEALLOCATE (VDIFFLUX)!
         If (ALLOCATED(ADVFLUX_GL)) DEALLOCATE (ADVFLUX_GL)!
         If (ALLOCATED(VDIFFLUX_GL)) DEALLOCATE (VDIFFLUX_GL)!
    !
         If (ALLOCATED(ALGDOC)) DEALLOCATE (ALGDOC)
         If (ALLOCATED(ALGPOC)) DEALLOCATE (ALGPOC)
    !
         If (ALLOCATED(HDRLPOC)) DEALLOCATE (HDRLPOC)
         If (ALLOCATED(HDRRPOC)) DEALLOCATE (HDRRPOC)
    !
         If (ALLOCATED(KLDC)) DEALLOCATE (KLDC)
         If (ALLOCATED(KRDC)) DEALLOCATE (KRDC)
    !
         If (ALLOCATED(KDCALG)) DEALLOCATE (KDCALG)
         If (ALLOCATED(KLCALG)) DEALLOCATE (KLCALG)
    !
         If (ALLOCATED(KLPC)) DEALLOCATE (KLPC)
         If (ALLOCATED(KRPC)) DEALLOCATE (KRPC)
    !
    !
    !
         If (ALLOCATED(KSUA)) DEALLOCATE (KSUA)
         If (ALLOCATED(KCOD)) DEALLOCATE (KCOD)
    !
         If (ALLOCATED(KLDN)) DEALLOCATE (KLDN)
         If (ALLOCATED(KRDN)) DEALLOCATE (KRDN)
    !
         If (ALLOCATED(KLDP)) DEALLOCATE (KLDP)
         If (ALLOCATED(KRDP)) DEALLOCATE (KRDP)
    !
         If (ALLOCATED(KDNALG)) DEALLOCATE (KDNALG)
         If (ALLOCATED(KLNALG)) DEALLOCATE (KLNALG)
    !
         If (ALLOCATED(KLPN)) DEALLOCATE (KLPN)
         If (ALLOCATED(KRPN)) DEALLOCATE (KRPN)
    !
    !
         If (ALLOCATED(KLPP)) DEALLOCATE (KLPP)
         If (ALLOCATED(KRPP)) DEALLOCATE (KRPP)
    !
         If (ALLOCATED(KLPALG)) DEALLOCATE (KLPALG)
         If (ALLOCATED(KDPALG)) DEALLOCATE (KDPALG)
    !
    !
         If (ALLOCATED(NTM)) DEALLOCATE (NTM)
    !
         If (ALLOCATED(KRCOAG)) DEALLOCATE (KRCOAG)!coagulation
    !
         If (ALLOCATED(FTMNL)) DEALLOCATE (FTMNL)
         If (ALLOCATED(FTHDR)) DEALLOCATE (FTHDR)
    !
         If (ALLOCATED(RESPC)) DEALLOCATE (RESPC)
    !
         If (ALLOCATED(DLALGC)) DEALLOCATE (DLALGC)
    !
         If (ALLOCATED(ALGNH4)) DEALLOCATE (ALGNH4)
         If (ALLOCATED(ALGNO3)) DEALLOCATE (ALGNO3)
         If (ALLOCATED(ALGDON)) DEALLOCATE (ALGDON)
         If (ALLOCATED(ALGPON)) DEALLOCATE (ALGPON)
         If (ALLOCATED(DENNO3)) DEALLOCATE (DENNO3)
         If (ALLOCATED(MNLLDON)) DEALLOCATE (MNLLDON)
         If (ALLOCATED(MNLRDON)) DEALLOCATE (MNLRDON)
         If (ALLOCATED(HDRLPON)) DEALLOCATE (HDRLPON)
         If (ALLOCATED(HDRRPON)) DEALLOCATE (HDRRPON)
    !
    !
         If (ALLOCATED(NFIX)) DEALLOCATE (NFIX)
    !
         If (ALLOCATED(ALGPO4)) DEALLOCATE (ALGPO4)
         If (ALLOCATED(ALGDOP)) DEALLOCATE (ALGDOP)
         If (ALLOCATED(ALGPOP)) DEALLOCATE (ALGPOP)
         If (ALLOCATED(MNLLDOP)) DEALLOCATE (MNLLDOP)
         If (ALLOCATED(MNLRDOP)) DEALLOCATE (MNLRDOP)
         If (ALLOCATED(HDRLPOP)) DEALLOCATE (HDRLPOP)
         If (ALLOCATED(HDRRPOP)) DEALLOCATE (HDRRPOP)
    !
         If (ALLOCATED(DORALG)) DEALLOCATE (DORALG)
         If (ALLOCATED(DOPR)) DEALLOCATE (DOPR)
         If (ALLOCATED(DCOD)) DEALLOCATE (DCOD)
         If (ALLOCATED(DDOC)) DEALLOCATE (DDOC)
    !
         If (ALLOCATED(NITRIF)) DEALLOCATE (NITRIF)
    !
         If (ALLOCATED(PSD)) DEALLOCATE (PSD)
         If (ALLOCATED(SAP)) DEALLOCATE (SAP)
    !
         If (ALLOCATED(ALGUP)) DEALLOCATE (ALGUP)
         If (ALLOCATED(ALGRES)) DEALLOCATE (ALGRES)
    !
         If (ALLOCATED(BENDOC)) DEALLOCATE (BENDOC)
         If (ALLOCATED(BENNH4)) DEALLOCATE (BENNH4)
         If (ALLOCATED(BENNO3)) DEALLOCATE (BENNO3)
         If (ALLOCATED(BENDON)) DEALLOCATE (BENDON)
         If (ALLOCATED(BENDOP)) DEALLOCATE (BENDOP)
         If (ALLOCATED(BENPO4)) DEALLOCATE (BENPO4)
         If (ALLOCATED(BENCOD)) DEALLOCATE (BENCOD)
         If (ALLOCATED(BENDO)) DEALLOCATE (BENDO)
         If (ALLOCATED(BENSA)) DEALLOCATE (BENSA)
         If (ALLOCATED(BENDEN)) DEALLOCATE (BENDEN)
         If (ALLOCATED(BENCH4G)) DEALLOCATE (BENCH4G)
         If (ALLOCATED(BENCH4A)) DEALLOCATE (BENCH4A)
#if defined(CSOD_CALC)
         If (ALLOCATED(CSODflux)) DEALLOCATE (CSODflux)!
#endif
!
         If (ALLOCATED(C1)) DEALLOCATE (C1)
         If (ALLOCATED(C2)) DEALLOCATE (C2)
         If (ALLOCATED(C2F)) DEALLOCATE (C2F)
         If (ALLOCATED(CSTAR)) DEALLOCATE (CSTAR)
         If (ALLOCATED(DTC)) DEALLOCATE (DTC)
         If (ALLOCATED(AC1)) DEALLOCATE (AC1)
    !
         If (ALLOCATED(QLIT)) DEALLOCATE (QLIT)
         If (ALLOCATED(CONLIT)) DEALLOCATE (CONLIT)
    !
         If (ALLOCATED(AFI1)) DEALLOCATE (AFI1)
         If (ALLOCATED(AFI2)) DEALLOCATE (AFI2)
         If (ALLOCATED(AFI3)) DEALLOCATE (AFI3)
    !
         If (ALLOCATED(ANL1)) DEALLOCATE (ANL1)
         If (ALLOCATED(ANL2)) DEALLOCATE (ANL2)
         If (ALLOCATED(ANL3)) DEALLOCATE (ANL3)
    !
         If (ALLOCATED(APL1)) DEALLOCATE (APL1)
         If (ALLOCATED(APL2)) DEALLOCATE (APL2)
         If (ALLOCATED(APL3)) DEALLOCATE (APL3)
    !
         If (ALLOCATED(ASL1)) DEALLOCATE (ASL1)
         If (ALLOCATED(ASL2)) DEALLOCATE (ASL2)
         If (ALLOCATED(ASL3)) DEALLOCATE (ASL3)
    !
         If (ALLOCATED(ACCHL1)) DEALLOCATE (ACCHL1)
         If (ALLOCATED(ACCHL2)) DEALLOCATE (ACCHL2)
         If (ALLOCATED(ACCHL3)) DEALLOCATE (ACCHL3)
    !
         If (ALLOCATED(ANPP)) DEALLOCATE (ANPP)
    !
         If (ALLOCATED(ARESP)) DEALLOCATE (ARESP)
         If (ALLOCATED(AKE)) DEALLOCATE (AKE)
         If (ALLOCATED(AGPP)) DEALLOCATE (AGPP)
    !
    !
    !
         If (ALLOCATED(ABENDOC)) DEALLOCATE (ABENDOC)
         If (ALLOCATED(ABENNH4)) DEALLOCATE (ABENNH4)
         If (ALLOCATED(ABENNO3)) DEALLOCATE (ABENNO3)
         If (ALLOCATED(ABENPO4)) DEALLOCATE (ABENPO4)
         If (ALLOCATED(ABENCOD)) DEALLOCATE (ABENCOD)
         If (ALLOCATED(ABENDO)) DEALLOCATE (ABENDO)
         If (ALLOCATED(ABENSA)) DEALLOCATE (ABENSA)
    !
         If (ALLOCATED(ASSFWS)) DEALLOCATE (ASSFWS)
         If (ALLOCATED(APCFWS)) DEALLOCATE (APCFWS)
         If (ALLOCATED(APNFWS)) DEALLOCATE (APNFWS)
         If (ALLOCATED(APPFWS)) DEALLOCATE (APPFWS)
         If (ALLOCATED(APSFWS)) DEALLOCATE (APSFWS)
         If (ALLOCATED(ACPIP)) DEALLOCATE (ACPIP)
         If (ALLOCATED(ACPOS)) DEALLOCATE (ACPOS)
    !
         If (ALLOCATED(ADFEED)) DEALLOCATE (ADFEED)
    !
         If (ALLOCATED(ABENCH4G)) DEALLOCATE (ABENCH4G)
         If (ALLOCATED(ABENCH4A)) DEALLOCATE (ABENCH4A)
    !
    !benthic algae
         If (ALLOCATED(AFIB)) DEALLOCATE (AFIB)
         If (ALLOCATED(ANLB)) DEALLOCATE (ANLB)
         If (ALLOCATED(APLB)) DEALLOCATE (APLB)
    !
         If (ALLOCATED(ANPPB)) DEALLOCATE (ANPPB)
    !
         If (ALLOCATED(ABBM)) DEALLOCATE (ABBM)!benthic algae
         If (ALLOCATED(ABLITE)) DEALLOCATE (ABLITE)
         If (ALLOCATED(ABMB)) DEALLOCATE (ABMB)
         If (ALLOCATED(APB)) DEALLOCATE (APB)
         If (ALLOCATED(APRB)) DEALLOCATE (APRB)
    !
    !absolute values
         If (ALLOCATED(ABADOC)) DEALLOCATE (ABADOC)
         If (ALLOCATED(ABAPOC)) DEALLOCATE (ABAPOC)
         If (ALLOCATED(ABANH4)) DEALLOCATE (ABANH4)
         If (ALLOCATED(ABANO3)) DEALLOCATE (ABANO3)
         If (ALLOCATED(ABAPON)) DEALLOCATE (ABAPON)
         If (ALLOCATED(ABAPO4)) DEALLOCATE (ABAPO4)
         If (ALLOCATED(ABAPOP)) DEALLOCATE (ABAPOP)
         If (ALLOCATED(ABADO)) DEALLOCATE (ABADO)
    !
         If (ALLOCATED(ACPOC)) DEALLOCATE (ACPOC)
         If (ALLOCATED(ACPON)) DEALLOCATE (ACPON)
         If (ALLOCATED(ACPOP)) DEALLOCATE (ACPOP)
    !
    !Suspension feeder
         If (ALLOCATED(ASFEED)) DEALLOCATE (ASFEED)
         If (ALLOCATED(AJNSF)) DEALLOCATE (AJNSF)
         If (ALLOCATED(AJPSF)) DEALLOCATE (AJPSF)
         If (ALLOCATED(ASODSF)) DEALLOCATE (ASODSF)
         If (ALLOCATED(ASASF)) DEALLOCATE (ASASF)
         If (ALLOCATED(ASUSF)) DEALLOCATE (ASUSF)
    !
         If (ALLOCATED(ASFGCIN)) DEALLOCATE (ASFGCIN)
         If (ALLOCATED(ASFCFEC)) DEALLOCATE (ASFCFEC)
         If (ALLOCATED(ASFCPSF)) DEALLOCATE (ASFCPSF)
    !
         If (ALLOCATED(AFLXCSF)) DEALLOCATE (AFLXCSF)
         If (ALLOCATED(AFLXNSF)) DEALLOCATE (AFLXNSF)
         If (ALLOCATED(AFLXPSF)) DEALLOCATE (AFLXPSF)
    !
         If (ALLOCATED(ARPOCSF)) DEALLOCATE (ARPOCSF)
         If (ALLOCATED(ARPONSF)) DEALLOCATE (ARPONSF)
         If (ALLOCATED(ARPOPSF)) DEALLOCATE (ARPOPSF)
    !
         If (ALLOCATED(ASSISF)) DEALLOCATE (ASSISF)
         If (ALLOCATED(ASSISASF)) DEALLOCATE (ASSISASF)
         If (ALLOCATED(ASSISUSF)) DEALLOCATE (ASSISUSF)
         If (ALLOCATED(ASSIPSF)) DEALLOCATE (ASSIPSF)
    !
    !
         If (ALLOCATED(A_T)) DEALLOCATE (A_T)
    !
         If (ALLOCATED(AP1)) DEALLOCATE (AP1)
         If (ALLOCATED(ABM1)) DEALLOCATE (ABM1)
         If (ALLOCATED(APR1)) DEALLOCATE (APR1)
    !
         If (ALLOCATED(AP2)) DEALLOCATE (AP2)
         If (ALLOCATED(ABM2)) DEALLOCATE (ABM2)
         If (ALLOCATED(APR2)) DEALLOCATE (APR2)
    !
         If (ALLOCATED(AP3)) DEALLOCATE (AP3)
         If (ALLOCATED(ABM3)) DEALLOCATE (ABM3)
         If (ALLOCATED(APR3)) DEALLOCATE (APR3)
    !
         If (ALLOCATED(AALGDOC)) DEALLOCATE (AALGDOC)
         If (ALLOCATED(AALGPOC)) DEALLOCATE (AALGPOC)
         If (ALLOCATED(ADENIT)) DEALLOCATE (ADENIT)
         If (ALLOCATED(AMNLDOC)) DEALLOCATE (AMNLDOC)
         If (ALLOCATED(AHDRPOC)) DEALLOCATE (AHDRPOC)
         If (ALLOCATED(AALGNH4)) DEALLOCATE (AALGNH4)
         If (ALLOCATED(AALGNO3)) DEALLOCATE (AALGNO3)
         If (ALLOCATED(AALGDON)) DEALLOCATE (AALGDON)
         If (ALLOCATED(AALGPON)) DEALLOCATE (AALGPON)
         If (ALLOCATED(ANT)) DEALLOCATE (ANT)
         If (ALLOCATED(ADENNO3)) DEALLOCATE (ADENNO3)
         If (ALLOCATED(AMNLDON)) DEALLOCATE (AMNLDON)
         If (ALLOCATED(AHDRPON)) DEALLOCATE (AHDRPON)
    !
         If (ALLOCATED(ACFIX)) DEALLOCATE (ACFIX)
         If (ALLOCATED(ANFIX)) DEALLOCATE (ANFIX)
    !
         If (ALLOCATED(AALGPO4)) DEALLOCATE (AALGPO4)
         If (ALLOCATED(AALGDOP)) DEALLOCATE (AALGDOP)
         If (ALLOCATED(AALGPOP)) DEALLOCATE (AALGPOP)
         If (ALLOCATED(AMNLDOP)) DEALLOCATE (AMNLDOP)
         If (ALLOCATED(AHDRPOP)) DEALLOCATE (AHDRPOP)
    !
         If (ALLOCATED(APSD)) DEALLOCATE (APSD)
         If (ALLOCATED(ASAP)) DEALLOCATE (ASAP)
    !
         If (ALLOCATED(AALGUP)) DEALLOCATE (AALGUP)
         If (ALLOCATED(AALGRES)) DEALLOCATE (AALGRES)
    !
         If (ALLOCATED(ADO)) DEALLOCATE (ADO)
         If (ALLOCATED(ADORALG)) DEALLOCATE (ADORALG)
         If (ALLOCATED(ADOPR)) DEALLOCATE (ADOPR)
    !
         If (ALLOCATED(ADCOD)) DEALLOCATE (ADCOD)
         If (ALLOCATED(ADDOC)) DEALLOCATE (ADDOC)
    !
         If (ALLOCATED(ANITRIF)) DEALLOCATE (ANITRIF)
    !
    !
         If (ALLOCATED(SBN)) DEALLOCATE (SBN)
         If (ALLOCATED(BBN)) DEALLOCATE (BBN)
         If (ALLOCATED(HMSBV)) DEALLOCATE (HMSBV)
    !
    !    IF(ALLOCATED(SFA))DEALLOCATE(SFA)
    !
         If (ALLOCATED(IWCMCB)) DEALLOCATE (IWCMCB)
         If (ALLOCATED(IWCMNB)) DEALLOCATE (IWCMNB)
         If (ALLOCATED(IWCMPB)) DEALLOCATE (IWCMPB)
         If (ALLOCATED(IWCMSB)) DEALLOCATE (IWCMSB)
    !
         If (ALLOCATED(WCMNB)) DEALLOCATE (WCMNB)
         If (ALLOCATED(WCMPB)) DEALLOCATE (WCMPB)
         If (ALLOCATED(WCMCB)) DEALLOCATE (WCMCB)
         If (ALLOCATED(WCMSB)) DEALLOCATE (WCMSB)
    !
         If (ALLOCATED(DLWCMNB)) DEALLOCATE (DLWCMNB)
         If (ALLOCATED(DLWCMPB)) DEALLOCATE (DLWCMPB)
         If (ALLOCATED(DLWCKMNB)) DEALLOCATE (DLWCKMNB)
         If (ALLOCATED(DLWCKMCB)) DEALLOCATE (DLWCKMCB)
    !
         If (ALLOCATED(FLUXS)) DEALLOCATE (FLUXS)
    !
         If (ALLOCATED(S1FLXCB)) DEALLOCATE (S1FLXCB)
         If (ALLOCATED(S1FLXNB)) DEALLOCATE (S1FLXNB)
         If (ALLOCATED(S1FLXPB)) DEALLOCATE (S1FLXPB)
    !
         If (ALLOCATED(S2FLXCB)) DEALLOCATE (S2FLXCB)
         If (ALLOCATED(S2FLXNB)) DEALLOCATE (S2FLXNB)
         If (ALLOCATED(S2FLXPB)) DEALLOCATE (S2FLXPB)
    !
         If (ALLOCATED(S3FLXCB)) DEALLOCATE (S3FLXCB)
         If (ALLOCATED(S3FLXNB)) DEALLOCATE (S3FLXNB)
         If (ALLOCATED(S3FLXPB)) DEALLOCATE (S3FLXPB)
    !
    !

    !
         If (ALLOCATED(BSVOL)) DEALLOCATE (BSVOL)
    !           IF(ALLOCATED(HSED))DEALLOCATE(HSED)
    !
         If (ALLOCATED(CTEMP)) DEALLOCATE (CTEMP)
    !
         If (ALLOCATED(PCFWS)) DEALLOCATE (PCFWS)
         If (ALLOCATED(PNFWS)) DEALLOCATE (PNFWS)
         If (ALLOCATED(PPFWS)) DEALLOCATE (PPFWS)
         If (ALLOCATED(PSFWS)) DEALLOCATE (PSFWS)
         If (ALLOCATED(SSFWS)) DEALLOCATE (SSFWS)
    !
         If (ALLOCATED(AASRAT)) DEALLOCATE (AASRAT)
    !
    !
    !----
    !
         If (ALLOCATED(BENDOCB)) DEALLOCATE (BENDOCB)
         If (ALLOCATED(BENNH4B)) DEALLOCATE (BENNH4B)
         If (ALLOCATED(BENNO3B)) DEALLOCATE (BENNO3B)
         If (ALLOCATED(BENPO4B)) DEALLOCATE (BENPO4B)
         If (ALLOCATED(BENCODB)) DEALLOCATE (BENCODB)
         If (ALLOCATED(BENDOB)) DEALLOCATE (BENDOB)
         If (ALLOCATED(BENSAB)) DEALLOCATE (BENSAB)
    !
         If (ALLOCATED(ATMFLXCB)) DEALLOCATE (ATMFLXCB)
         If (ALLOCATED(ATMFLXNB)) DEALLOCATE (ATMFLXNB)
         If (ALLOCATED(ATMFLXPB)) DEALLOCATE (ATMFLXPB)
    !
         If (ALLOCATED(BENFLXPCB)) DEALLOCATE (BENFLXPCB)
    !
         If (ALLOCATED(BENFLXDNB)) DEALLOCATE (BENFLXDNB)
         If (ALLOCATED(BENFLXPNB)) DEALLOCATE (BENFLXPNB)
    !
    !
         If (ALLOCATED(BENFLXDPB)) DEALLOCATE (BENFLXDPB)
         If (ALLOCATED(BENFLXPPB)) DEALLOCATE (BENFLXPPB)
    !
         If (ALLOCATED(DLSEDKNB)) DEALLOCATE (DLSEDKNB)
         If (ALLOCATED(DLSEDKCB)) DEALLOCATE (DLSEDKCB)
    !
         If (ALLOCATED(BURIALFLXCB)) DEALLOCATE (BURIALFLXCB)
         If (ALLOCATED(BURIALFLXNB)) DEALLOCATE (BURIALFLXNB)
         If (ALLOCATED(BURIALFLXPB)) DEALLOCATE (BURIALFLXPB)
    !
    !
    ! Added C2_GL
         If (ALLOCATED(C2_GL)) DEALLOCATE (C2_GL)
    ! Added below for output of CCHL1_GL
         If (ALLOCATED(CCHL1_GL)) DEALLOCATE (CCHL1_GL)
         If (ALLOCATED(T_GL)) DEALLOCATE (T_GL)
         If (ALLOCATED(S_GL)) DEALLOCATE (S_GL)
         If (ALLOCATED(total_netPP_GL)) DEALLOCATE (total_netPP_GL)
    !
    !  GPP and NPP to be saved in history output
         If (ALLOCATED(GPP_GL)) DEALLOCATE (GPP_GL)
         If (ALLOCATED(NPP_GL)) DEALLOCATE (NPP_GL)
    !
    !Added D_GL, EL_GL, H_GL (total depth, surface elevation and bathymetric depth)
         If (ALLOCATED(D_GL)) DEALLOCATE (D_GL)
         If (ALLOCATED(H_GL)) DEALLOCATE (H_GL)
         If (ALLOCATED(EL_GL)) DEALLOCATE (EL_GL)
    !
    !--Debugging benthic fluxes ---
         If (ALLOCATED(BFLUX)) DEALLOCATE (BFLUX)
         If (ALLOCATED(BFLUX_GL)) DEALLOCATE (BFLUX_GL)
    !
         If (ALLOCATED(BFLUXNX)) DEALLOCATE (BFLUXNX)
         If (ALLOCATED(BFLUXNX_GL)) DEALLOCATE (BFLUXNX_GL)
    !
         If (ALLOCATED(BFLUXB_GL)) DEALLOCATE (BFLUXB_GL)
    !
    !
    !--Deallocating global hydro variables here ---
    !
         If (ALLOCATED(UL_GL)) DEALLOCATE (UL_GL)
         If (ALLOCATED(VL_GL)) DEALLOCATE (VL_GL)
         If (ALLOCATED(WTSL_GL)) DEALLOCATE (WTSL_GL)
         If (ALLOCATED(KHL_GL)) DEALLOCATE (KHL_GL)
         If (ALLOCATED(SL_GL)) DEALLOCATE (SL_GL)
         If (ALLOCATED(TL_GL)) DEALLOCATE (TL_GL)
         If (ALLOCATED(ELL_GL)) DEALLOCATE (ELL_GL)
         If (ALLOCATED(DTFAL_GL)) DEALLOCATE (DTFAL_GL)
         If (ALLOCATED(UARD_OBCN_GL)) DEALLOCATE (UARD_OBCN_GL)
         If (ALLOCATED(XFLUX_OBC_GL)) DEALLOCATE (XFLUX_OBC_GL)
    !
    !
    !--Added the following algae debugging global variables here ---
    !
         If (ALLOCATED(P1_GL)) DEALLOCATE (P1_GL)
         If (ALLOCATED(P2_GL)) DEALLOCATE (P2_GL)
         If (ALLOCATED(BM1_GL)) DEALLOCATE (BM1_GL)
         If (ALLOCATED(BM2_GL)) DEALLOCATE (BM2_GL)
         If (ALLOCATED(NL1_GL)) DEALLOCATE (NL1_GL)
         If (ALLOCATED(NL2_GL)) DEALLOCATE (NL2_GL)
         If (ALLOCATED(PL1_GL)) DEALLOCATE (PL1_GL)
         If (ALLOCATED(PL2_GL)) DEALLOCATE (PL2_GL)
         If (ALLOCATED(B1SZ_GL)) DEALLOCATE (B1SZ_GL)
         If (ALLOCATED(B2SZ_GL)) DEALLOCATE (B2SZ_GL)
         If (ALLOCATED(B1LZ_GL)) DEALLOCATE (B1LZ_GL)
         If (ALLOCATED(B2LZ_GL)) DEALLOCATE (B2LZ_GL)
         If (ALLOCATED(PR1_GL)) DEALLOCATE (PR1_GL)
         If (ALLOCATED(PR2_GL)) DEALLOCATE (PR2_GL)
    !
         If (LIGHT_EXTINCTION) Then
            If (ALLOCATED(IAVG_GL)) DEALLOCATE (IAVG_GL)
            If (ALLOCATED(FI1_GL)) DEALLOCATE (FI1_GL)
            If (ALLOCATED(FI2_GL)) DEALLOCATE (FI2_GL)
         End If
    !
         If (ALLOCATED(IK1_GL)) DEALLOCATE (IK1_GL)
         If (ALLOCATED(IK2_GL)) DEALLOCATE (IK2_GL)
         If (ALLOCATED(IK1)) DEALLOCATE (IK1)
         If (ALLOCATED(IK2)) DEALLOCATE (IK2)
    !
    !---Added the following for debugging DOXG
    !
         If (ALLOCATED(DDOC_GL)) DEALLOCATE (DDOC_GL)
         If (ALLOCATED(DCOD_GL)) DEALLOCATE (DCOD_GL)
    !
         If (ALLOCATED(NITRIF_GL)) DEALLOCATE (NITRIF_GL)
         If (ALLOCATED(DOSZ_GL)) DEALLOCATE (DOSZ_GL)
         If (ALLOCATED(DOLZ_GL)) DEALLOCATE (DOLZ_GL)
         If (ALLOCATED(DOPR_GL)) DEALLOCATE (DOPR_GL)
         If (ALLOCATED(DORALG_GL)) DEALLOCATE (DORALG_GL)
         If (ALLOCATED(RESP_GL)) DEALLOCATE (RESP_GL)
         If (ALLOCATED(FTCOD_GL)) DEALLOCATE (FTCOD_GL)
         If (ALLOCATED(COD_GL)) DEALLOCATE (COD_GL)
         If (ALLOCATED(REAERDO_GL)) DEALLOCATE (REAERDO_GL)
         If (ALLOCATED(REAERDO)) DEALLOCATE (REAERDO)

		 !If (SOLIDS_CALC) Then !Adi 2/13/18
            If (ALLOCATED(TSS_GL)) DEALLOCATE (TSS_GL)
            If (ALLOCATED(SSI_GL)) DEALLOCATE (SSI_GL)
            If (ALLOCATED(TURB_GL)) DEALLOCATE (TURB_GL)
			If (ALLOCATED(KESS_GL)) DEALLOCATE (KESS_GL)
			If (ALLOCATED(PARAD_GL)) DEALLOCATE (PARAD_GL)
         !End If


		 !If (SAV_CALC) Then !Adi 4/23/18
            If (ALLOCATED(LEAF_GL)) DEALLOCATE (LEAF_GL)
            If (ALLOCATED(STEM_GL)) DEALLOCATE (STEM_GL)
            If (ALLOCATED(ROOT_GL)) DEALLOCATE (ROOT_GL)
			If (ALLOCATED(TUBER_GL)) DEALLOCATE (TUBER_GL)
			If (ALLOCATED(EP_GL)) DEALLOCATE (EP_GL)
         !End If
    !
    !---Declared a large buffer (so that we do not have to dynamically allocate and deallocate buffers
    !   all the time--
#if defined (MULTIPROCESSOR)
         If (ALLOCATED(RCVBUF)) DEALLOCATE (RCVBUF)
         If (ALLOCATED(SNDBUF)) DEALLOCATE (SNDBUF)
#endif
    !
    !
         Return
    !
      End Subroutine WQM_DEALLOC
  !
End Module MOD_WQM
!
!
