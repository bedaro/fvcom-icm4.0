!mod_sav.F
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
Module MOD_SAV
  !
  !************************************************************************
  !**                                                                    **
  !**       Submerged Aquatic Vegetation Subroutines for FVCOM-ICM       **
  !**               				   **
  !**                                                                    **
  !************************************************************************
  !
      Use MOD_PREC, Only: SP
  !
	  Use MOD_SIZES, Only: NSAVM ! AN
  !
	  Use MOD_LIMS, Only: MTLOC !AN
  !
      Real (SP) :: LEAFI (10), STEMI (10), ROOTI (10), &
	 & TUBERI (10), EPI (10) !AN
  !
      Character (Len=8) SPNAM, EPINAM !name of the epiphyte
  !
  !   Integer :: NSAVM !global maximum number of SAV species in all of tracer control elements (TCEs) - polygon
  !surrounding a given node ! Move to MOD_SIZES.F - AN
  !
      Real (SP) :: KHSTOX, KESAV, KEEPI, KHEP, NLIMEPI !Never used ???
  !
      Real (SP) :: KHNLEAF, KHNROOT, KHPLEAF, KHPROOT, KHNEPI, KHPEPI, &
     & KHNPSAV, KHNPEPI !Threshold for NH4 in calculating NH4 preference over NO3 in N uptake(gN/m^3) for epiphytes
  !
      Real (SP) :: PMSAV, PRSPSAVB, BMSAV, BMTBRREF, FDOSR, SALMAX, &
     & ALPHSAV, ALAC, SLSAV, ACDWSAV, ANDWSAV, APDWSAV, ANCSAV, APCSAV, &
     & HCAN, ACAN, BCAN !coefficient for calculating leaf (gC/m^2) and stem (gC/m^2)
  !PMSAV,       &   !maximum production rate of SAV, gC/gDW/day  where DW is dry weight
  !PRSPSAVB,    &   !fraction of leaf production rate that is counted as photochemical respiration at salt << saltmax
  !					at salt >> saltmax, PRSPSAV approaches 1, and there will be no production of SAV leaf hence
  !no photosynthesis for SAV growth
  !BMSAV,       &   !Basal metabolism rate of SAV (1/day)
  !BMTBRREF,    &   !Reference basal metabolism rate of tube [1/day] BMTUBE=BMTBRREF*f(T), i.e. function of temperature
  !FDOSR,       &   !Fraction of DO production discounted due to route of Leaf production to root
  !SALMAX,      &   !Salinity at which PRSPSAV is approaching 1 (maximum) from PRSPSAVB, [psu]
  !ALPHSAV,     &   !Initial slope of P-I curve for SAV [ (gC/gDW)/(E/m^2) ]
  !ALAC,        &   !m^2-leaf-area/gC-leaf (leaf area to leaf carbon ratio)
  !Note in 2002 version, it is the inverse ACLA: gC-leaf/m^2-leaf-area
  !SLSAV,       &   ! 1/day
  !ACDWSAV,     &   !Carbon to Dry Weight ratio of SAV (gC/gDW)
  !ANDWSAV,     &   !Nitrogen to Dry Weight ratio of SAV (gN/gDW)
  !APDWSAV,     &   !Phosphorus to Dry Weight ratio of Sav (gP/gDW)
  !ANCSAV,      &   !Nitrogen to Carbon ratio of SAV (gN/gC)
  !APCSAV,      &   !P to C ratio of SAV (gP/gC)
  !HCAN,        &   !Canopy height (m), HCAN = ACAN + BCAN *( leaf+stem )
  !ACAN,        &   !basic canopy height (m)
  !BCAN             !coefficient for calculating leaf (gC/m^2) and stem (gC/m^2)
  !contribution to canopy height,  [m^3/gC]
!
  !Distribution of C, N,P on SAV mortality
      Real (SP) :: FNISAV, FNLDSAV, FNRDSAV, FNLPSAV, FNRPSAV, FPISAV, &
     & FPLDSAV, FPRDSAV, FPLPSAV, FPRPSAV, FCLDSAV, FCRDSAV, FCLPSAV, &
     & FCRPSAV, FDOSAV !Fraction of Carbon routed to CO2 (oxygen consumption) in water column due to SAV mortality
  !Note FDOSAV + FCLDSAV+ FCRDSAV + FCLPSAV + FCRPSAV = 1
  !
  !
      Real (SP) :: WSSSAV, WSLSAV, WSRSAV, WS1SAV, WS2SAV, WS3SAV, &
     & WSUSAV !Increase of particulate biogenic silica settling rate by SAV leaf and stem, m^3/gC/day
  !
      Integer :: NDDI, NSAVCELL, NSAVSPCM !maximum number of species in a cell
  !
      Real (SP) :: PMEPI, BMEPI, PREPI, PRSPEPI !Photochemical respiration rate of epiphyte --> BMEP
  !
      Real (SP) :: ALPHEPI, CCHLEPI, ANCEPI, APCEPI, ADWCEPI !ratio of total epiphyte dry weight to viable carbon (gDW/gC)
  !
      Real (SP) :: FCLDEPI, FCRDEPI, FNIEPI, FNLDEPI, FNRDEPI, FNLPEPI, &
     & FNRPEPI, FPIEPI, FPLDEPI, FPRDEPI, FPLPEPI, FPRPEPI
  !FCLDEPI         !Fraction of C routed to LDOC due to basal metabolism of epiphyte
  !FCRDEPI	   !Fraction of C routed to RDOC due to basal metabolism of epiphyte
  !				Missing FCLPEPI and FCRPEPI to routt to LPOC and RPOC 
  !				FCLPEPI and FCRPEPI should be added 
  !				Missing FDOEPI compared to FDOSAV for SAV
  !				Code assumes that no carbon will be routed to LPOP and RPOP due to metabolism
  !				and also FDOEPI=1-FCLDEPI-FCRDEPI, and FDOEPI is only implicitly used in calculating DOEPIW
  !FNIEPI         !Fraction of N routed to NH4  due to epiphyte basal metabolism
  !FNLDEPI        !Fraction of N routed to LDON due to epiphyte basal metabolism
  !FNRDEPI        !Fraction of N routed to RDON due to epiphyte basal metabolism
  !FNLPEPI        !Fraction of N routed to LPON due to epiphyte basal metabolism
  !FNRPEPI        !Fraction of N routed to RPON due to epiphyte basal metabolism
  !FPIEPI         !Fraction of P routed to  PO4 due to basal metabolism of epiphyte
  !FPLDEPI        !Fraction of P routed to LDOP due to basal metabolism of epiphyte
  !FPRDEPI        !Fraction of P routed to RPOP due to basal metabolism of epiphyte
  !FPLPEPI        !Fraction of P routed to LPOP due to basal metabolism of epiphyte
  !FPRPEPI        !Fraction of P routed to RPOP due to basal metabolism of epiphyte
!
  !Distribution of C, N, P upon epiphyte mortality (predation)
      Real (SP) :: FCLDPEP, FCRDPEP, FCLPPEP, FCRPPEP, FNIPEP, FNLDPEP, &
     & FNRDPEP, FNLPPEP, FNRPPEP, FPIPEP, FPLDPEP, FPRDPEP, FPLPPEP, &
     & FPRPPEP !Fraction of P routed to RPOP due to predation on epiphytes and slough
  !
      Real (SP) :: FPSR 
  !
  !: these fractions should change due to aging (senescence)
      Real (SP), Dimension (0:1000) :: FPLEAF, FPROOT, FPSTEM, FPTUBER, &
     & TRTBRLF !transfer rate from tuber to leaf (gC-Leaf/day/gC-Tuber)
  !FPLEAF, & !fraction of production assigned to leaf growth
  !FPROOT, & !fraction of production assigned to root growth
  !FPSTEM, & !fraction of production assigned to stem growth
  !FPTUBER,& !fraction of production assigned to tuber growth
  !			:  FPLEAF +FPROOT+FPSTEM+FPTUBER = 1

  !TRTBRLF   !transfer rate from tuber to leaf (gC-Leaf/day/gC-Tuber)
  !
  !Look up table of temperature control in a year on growth and mortality rates of SAV and epiphytes
  !    Real (SP), Dimension (-50:400) :: FTPSAV, FTRSAV, FTPEP, FTREP, &
  !   & FTPREP !temperature control on predation on epiphyte 
  !
  !fluxes to sediments by SAV mortality
      Real (SP), Dimension (3) :: FRPOCSAV, FRPONSAV, FRPOPSAV !3G partitioning of POC, PON, POP sources to sediments from SAV (non-dimensional)
  !
  !State variables of SAV model
      Real (SP), Allocatable, Dimension (:, :) :: LEAF, ROOT, STEM, &
     & TUBER, EP !Epiphytes biomass (gC-epiphytes/gm^2-leaf)
  
      REAL (SP), ALLOCATABLE, DIMENSION(:) :: BM_LEAF
  !
      Real (SP), Allocatable, Dimension (:) :: SAVAREA !area of SAV coverage for each surface cell (m^2) --> needs to be changed for FVCOM grid
      Integer, Allocatable, Dimension (:) :: SAVCELL !grid cell number for all cells that have SAV in it
      INTEGER, DIMENSION(:), ALLOCATABLE :: SAVNDS_GL   ! array for global sav nodes
  !size 1xNSAVCELL --> need to change to element number or node number!!
      Integer, Allocatable :: NSAVSPC (:)!number of SAV species in each grid cell
  !size 1xNSAVCELL
  !
      Real (SP), Allocatable, Dimension (:, :) :: SAVDPH, SAVFRAC !coverage for each SAV cell and each vertical deph increment
  !REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: SAVDPH,  &  !Depth of SAV (mean depth of the plant) for all species in each grid cell
  !size, NSAVCELL x NSAVSPC(B) where B is cell number, B=1,...,NSAVCELL
  !SAVFRAC     !coverage for each SAV cell and each vertical deph increment
  !i.e. coverage is function of grid cell and SAV species
  !
  !Growth and mortality rates of SAV
      Real (SP), Allocatable, Dimension (:, :) :: PLEAF !Leaf growth rate (1/day), P (SAV shoots growth) in C-M 2001 papers SH equation
      Real (SP), Allocatable, Dimension (:) :: BMLEAF, BMSTEM, BMROOT, &
     & BMTUBER !respiration rate of tuber of SAV (1/day) with temperature control
  !
  !Growth and mortality rates of SAV
      Real (SP), Allocatable, Dimension (:, :) :: PEP, BMEP, PREP !predation rate on epiphytes (1/day) = PREPI*EP
  !where PREPI has unit gC-shoot/gC-epiphyte/day
  !
      Real (SP), Allocatable :: SLSH (:)!epiphyte sloughing rate (1/day)
  !
  !Nutrient limitations on SAV and epiphytes growth
      Real (SP), Allocatable, Dimension (:) :: NLSAV, PLSAV, FNSEDSAV, &
     & FPSEDSAV, NLEPI, PLEPI, FHS !fresh water SAV species growth limitation due to sulfide HS^-
  !NLSAV,    &  !Nitrogen limiation on SAV growth (non-dimensional)
  !PLSAV,    &  !Phosphorus limitation on SAV growth (non-dimensional)
  !FNSEDSAV, &  !fraction of nitrogen uptake by SAV due to nitrogen in sediment pore water
  !				: further split into NH4 and NO3 in future versions
  !          uptake by SAV from sediments
  !FPSEDSAV, &  !fraction of phosphorus uptake by SAV due to P in sediment pore water
  !NLEPI,    &  !Nitrogen limitation on epiphytes growth (non-dimensional)
  !PLEPI,    &  !Phosphorus limitation on epiphytes growth (non-dimensional)
  !FHS          !fresh water SAV species growth limitation due to sulfide HS^-
  !
  !(in sediment anerobic layer) toxicity (dimensionless)
  !
      Real (SP), Allocatable, Dimension (:, :) :: EPATN, SAVATN, &
     & WATATN, FISH, FIEP, NPPSAV, NPPEPI !net epiphyte production per day (growth - basal metabolism) (gC-epi/gC-leaf/day)
  !REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: EPATN,  & !Light attenuation by epiphyte
  !							e^(-Kep*Acla*Adwcep*EP) in Cerco-Moore 2001 paper for Ish eqn
  !	e^(-Kep*Adwcep*EP) in new version as EP changed to gC-Epi/m^2-leaf unit
  !SAVATN, & !Light attenuation by SAV (shoots)
  ![1-e^(-Ksh*SH+KESS*HCAN)]/(Ksh*SH+KESS*HCAN) in C-M 2001 Iwc eqn
  !WATATN, & !water column attenuation of light(dimensionless) at canopy top
  != exp(-KESS*ZTOCNPY) where ZTOCNPY is distance from surface
  !                       to canopy top
  !FISH,   & !light limitation function for LEAF (shoots) f(I)=Ish/sqrt(Ish^2+Ik^2)
  !FIEP,   & !light limitation function for epiphytes, f(I)=Iep/sqrt(Iep^2+Ik^2)
  !NPPSAV, & !net leaf production per day (growth-basal metabolism) (gC/m^2/day)
  !NPPEPI    !net epiphyte production per day (growth - basal metabolism) (gC-epi/gC-leaf/day)
  !
  !fluxes to water column
      Real (SP), Allocatable, Dimension (:) :: DOSAVW, LDOCSAVW, &
     & RDOCSAVW, LPOCSAVW, RPOCSAVW, NH4SAVW, NO3SAVW, LDONSAVW, &
     & RDONSAVW, LPONSAVW, RPONSAVW, PO4SAVW, LDOPSAVW, RDOPSAVW, &
     & LPOPSAVW, RPOPSAVW, DOEPIW, LDOCEPIW, RDOCEPIW, LPOCEPIW, &
     & RPOCEPIW, NH4EPIW, NO3EPIW, LDONEPIW, RDONEPIW, LPONEPIW, &
     & RPONEPIW, PO4EPIW, LDOPEPIW, RDOPEPIW, LPOPEPIW, RPOPEPIW
  !DOSAVW, &  !DO source (flux) due to SAV (source-sink) (gO2/m^2/day )
  !DOSAVW --> DTDO (gO2/m^3/sec)--> DTC (gO2/m^3/sec)
  !-->DTM(gO2/m^3/sec)-->DTM*DLT-->C1(:,:,27)(gO2/m^3) !27th constituents is DO
  !LDOCSAVW, &  !LDOC flux due to SAV (gC/m^2/day)
  !WL                                              !***ALSO used by TVDS.F ****
  !LDOCSAVW --> DTDOC (gC/m^3/sec) -> DTC --> DTM --> DTM*DLT --> C1(:,:,9)
  !RDOCSAVW, &  !RDOC flux due to SAV (gC/m^2/day)
  !LPOCSAVW, &  !LPOC flux due to SAV (gC/m^2/day)
  !WL                                              !***ALSO used by TVDS.F ****
  !RPOCSAVW, &  !RPOC flux due to SAV (gC/m^2/day)
  !WL                                              !***ALSO used by TVDS.F ****
  !NH4SAVW,  &  !NH4 flux to water column (respiration source - uptake sink by SAV growth) (gN/m^2/day)
  !NH4SAVW (gN/m^2/day) --> DTNH4  (gN/m^3/sec) --> DTC --> DTM --> C1(:,:,13)
  !NO3SAVW,  &  !NO3 flux to water column (due to uptake by SAV, positive increasing NO3 in water) (gN/m^2/day)
  !LDONSAVW, &  !LDON flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf, basal metabolism of
  !RDONSAVW, &  !RDON flux
  !LPONSAVW, &  !LPON flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem (gN/m^2/day)
  !RPONSAVW, &  !RPON flux
  !PO4SAVW,  &  !PO4 flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem and uptake by
  !LDOPSAVW, &  !LDOP flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem) (gP/m^2/day
  !RDOPSAVW, &  !RDOP flux (gP/m^2/day)
  !LPOPSAVW, &  !LPOP flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem) (gP/m^2/day
  !RPOPSAVW, &  !RPOP flux (gP/m^2/day)
  !DOEPIW, &  !DO   flux (positive into water column) due to epiphytes (gO2/m^2/day)
  !LDOCEPIW, &  !LDOC flux (positive into water column) due to epiphytes (gC/m^2/day)
  !RDOCEPIW, &  !RDOC flux (positive into water column) due to epiphytes (gC/m^2/day)
  !LPOCEPIW, &  !LPOC flux (positive into water column) due to epiphytes (gC/m^2/day)
  !RPOCEPIW, &  !RPOC flux (positive into water column) due to epiphytes (gC/m^2/day)
  !NH4EPIW, &  !NH4  flux (positive into water column) due to epiphytes (gN/m^2/day) (basal meabolism + predation + sloughing - upta
  !NO3EPIW, &  !NO3  flux
  !LDONEPIW, &  !LDON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/
  !RDONEPIW, &  !RDON flux (positive into water column) (gN/m^2/day)
  !LPONEPIW, &  !LPON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/
  !RPONEPIW, &  !RPON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/
  !PO4EPIW,  &  !PO4  flux (positive into water column) due to epiphytes mortality (basal metabolism, salt toxicity, sloughing) and
  !LDOPEPIW, &  !LDOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/
  !RDOPEPIW, &  !RPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/
  !LPOPEPIW, &  !LPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/
  !RPOPEPIW     !RPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/
  !fluxes to sediments
!
      Real (SP), Allocatable, Dimension (:) :: SEDPOCSAV, SEDPONSAV, &
     & SEDPOPSAV, SEDDOSAV, SEDNH4SAV, SEDNO3SAV, SEDPO4SAV
  !SEDPOCSAV, &  !POC flux to sediments (positive into sediments) due to SAV (gC/m^2/day)
  !SEDPOCSAV (gC/m^2/day) --> JPOC (mgC/m^2/day) -> POC1, POC2, POC3 (source of particulate organic carbon in sediments)
  !SEDPONSAV, &  !PON flux to sediments (positive into sediments) due to SAV (gN/m^2/day)
  !SEDPONSAV (gN/m^2/day) --> JPON (mgN/m^2/day) -> PON1, PON2, PON3 (source of particulate organic nitrogen in sediments)
  !SEDPOPSAV, &  !POP flux to sedimetns (positive into sediments) due to SAV (gP/m^2/day)
  !SEDPOPSAV (gP/m^2/day) --> JPON (mgP/m^2/day) -> POP1, POP2, POP3 (source of particulate organic phosphorus in sediments)
  !SEDDOSAV,  &  !DO flux to water column (positive into water) due to SAV growth in sediments (source) and mortality of SAV in sedi
  !(gO2/m^2/day) 
  !SEDNH4SAV, &  !NH4 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
  !(gN/m^2/day), SEDNH4SAV --> NH4T2TM1 (mg/N/m^3) decerase in sediments
  !SEDNO3SAV, &  !NO3 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
  !(gN/m^2/day), SEDNO3SAV --> NO3T2TM1 (mg/N/m^3) decerase in sediments
  !SEDPO4SAV     !PO4 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
  !
!
      Real (SP), Allocatable, Dimension (:) :: DGRDAYS !Degree day calculated from temperature of water
!
  !
  !time average values
  !
      Real (SP), Allocatable, Dimension (:, :) :: ALEAF, AROOT, ASTEM, &
     & ATUBER, AEP !Average epiphyte biomass (gC-epi/gC-leaf)
  !
      Real (SP), Allocatable, Dimension (:, :) :: APLEAF !Averaged leaf growth rate (1/day)
      Real (SP), Allocatable, Dimension (:) :: ABMLEAF, ABMTUBER !Average basal metabolism rate of tuber (1/day)
  !
      Real (SP), Allocatable, Dimension (:, :) :: APEP, ABMEP, APREP !Average predation rate on epiphytes (1/day) = APREP*EP
  !
      Real (SP), Allocatable, Dimension (:) :: ASLSH !Average epiphute sloughing rate (1/day)
  !
      Real (SP), Allocatable, Dimension (:) :: ANLSAV, APLSAV, ANLEPI, &
     & APLEPI, AFNSED, AFPSED, AFHS !Averaged SAV growth limiation due to sulfide HS^-
  !
  !attenuation and light limitation
      Real (SP), Allocatable, Dimension (:, :) :: AEPATN, AWATATN, &
     & AFISH, AFIEP !Averaged light limitation for Epiphytes
  !net production
      Real (SP), Allocatable, Dimension (:, :) :: ANPPSAV, ANPPEPI !Averaged net epiphyte production per day (gC-epi/gC-leaf/day)
  !
  !time averaged flux to water column
      Real (SP), Allocatable, Dimension (:) :: ADOSAVW, ADOCSAVW, &
     & APOCSAVW, ANH4SAVW, ANO3SAVW, ADONSAVW, APONSAVW, APO4SAVW, &
     & ADOPSAVW, APOPSAVW, ADOEPIW, ANH4EPIW, ANO3EPIW, APO4EPIW, &
     & ADOCEPIW, APOCEPIW, ADONEPIW, APONEPIW, ADOPEPIW, APOPEPIW
  !ADOCSAVW,  &  !average DOC flux due to SAV
  !APOCSAVW,  &  !average POC flux due to SAV
  !ANH4SAVW,  &  !average NH4 flux due to SAV
  !ANO3SAVW,  &  !average NO3 flux due to SAV
  !ADONSAVW,  &  !average DON flux due to SAV
  !APONSAVW,  &  !average PON flux due to SAV
  !APO4SAVW,  &  !average PO4 flux due to SAV
  !ADOPSAVW,  &  !average DOP flux due to SAV
  !APOPSAVW,  &  !average POL flux due to SAV
  !ADOEPIW,  &   !DO  flux (positive into water column ) due to epiphytes (gO2/m2/day)
  !ANH4EPIW,  &  !NH4 flux (positive into water column) due to epiphytes (gN/m2/day)
  !ANO3EPIW,  &  !NO3 flux (positive into water column) due to epiphytes (gN/m2/day)
  !APO4EPIW,  &  !PO4 flux (positive into water column) due to epiphytes (gP/m2/day)

  !ADOCEPIW,  &  ! DOC flux to water column due to epiphytes (gC/m^2/day)
  !APOCEPIW,  &  ! POC flux to water column due to epiphytes (gC/m^2/day)
  !ADONEPIW,  &  ! DON flux to water column due to epiphytes (gN/m^2/day)
  !APONEPIW,  &  ! PON flux to water column due to epiphytes (gN/m^2/day)
  !ADOPEPIW,  &  ! DOP flux to water column due to epiphytes (gP/m^2/day)
  !APOPEPIW      ! POP flux to water column due to epiphytes (gP/m^2/day)
!
  ! time averaged flux to sediments
      Real (SP), Allocatable, Dimension (:) :: ASEDDOSAV, ASEDPOCSAV, &
     & ASEDPONSAV, ASEDPOPSAV, ASEDNH4SAV, ASEDPO4SAV !PO4 flux due to SAV, positive leaving sediments
!
	  Real (SP) :: TRPMSAV, TRBMSAV, TRPMEPI, TRBMEPI, TRPREPI  ! AN: moved for calculating f(T)
	  Real (SP) :: KTPS1, KTPS2, KTBMSAV, KTPE1, KTPE2, KTBME, KTPRE !temperature control coefficient for epiphyte predation loss rate [1/degC]
!

	  Real (SP) :: AveT, AveSalt, AveNO3, AvePO4, AveNH4, AveSSI !AN
    INTEGER :: NSAVCELL_GL


Contains
  !Subroutines:
  !	subroutine SAV_INIT()
  !	subroutine SAV_ALLOC()
  !	subroutine SAV_DEALLOC()
  !	subroutine SAV_LOADS_ALLOC()
  !	subroutine SAV_LOADS_DEALLOC()
  !	subroutine SAV_READ()
  !	subroutine SAV_COMP()
  ! subroutine SAV_INIT_ICI() !AN
!
  !************************************************************************
  !**             S U B R O U T I N E   SAV_INIT                        **
  !************************************************************************

      Subroutine SAV_INIT
         Implicit None
    !WLong: moved these here wqm_inputs.F
         ALEAF = 0.0
         AROOT = 0.0
         ASTEM = 0.0
         ATUBER = 0.0
         AEP = 0.0
     !
         APLEAF = 0.0
         ABMLEAF = 0.0
         ABMTUBER = 0.0
     !
         APEP = 0.0
         ABMEP = 0.0
         APREP = 0.0
         ASLSH = 0.0
     !
         ANLSAV = 0.0
         APLSAV = 0.0
         ANLEPI = 0.0
         APLEPI = 0.0
         AFNSED = 0.0
         AFPSED = 0.0
         AFHS = 0.0
     !
         AEPATN = 0.0
         AWATATN = 0.0
         AFISH = 0.0
         AFIEP = 0.0
         ANPPSAV = 0.0
         ANPPEPI = 0.0
     !
         ADOCSAVW = 0.0
         APOCSAVW = 0.0
         ADOCEPIW = 0.0
         APOCEPIW = 0.0
     !
         ANH4SAVW = 0.0
         ANO3SAVW = 0.0
         ADONSAVW = 0.0
         APONSAVW = 0.0
         ANH4EPIW = 0.0
         ANO3EPIW = 0.0
         ADONEPIW = 0.0
         APONEPIW = 0.0
     !
         APO4SAVW = 0.0
         ADOPSAVW = 0.0
         APOPSAVW = 0.0
         APO4EPIW = 0.0
         ADOPEPIW = 0.0
         APOPEPIW = 0.0
     !
         ADOSAVW = 0.0
         ADOEPIW = 0.0
     !
         ASEDDOSAV = 0.0
         ASEDPOCSAV = 0.0
         ASEDPONSAV = 0.0
         ASEDNH4SAV = 0.0
         ASEDPOPSAV = 0.0
         ASEDPO4SAV = 0.0
    !
      End Subroutine SAV_INIT





  !************************************************************************
  !**             S U B R O U T I N E   SAV_ALLOC                        **
  !************************************************************************
!
      Subroutine SAV_ALLOC
    !Allocate all variables related to SAV
    !much of this section is moved from wqm_modules.F
!         Use MOD_LIMS, Only: MTLOC 
!
         Allocate (LEAF(MTLOC, NSAVM))
         LEAF = 0.0 !1
         Allocate (ROOT(MTLOC, NSAVM))
         ROOT = 0.0 !2
         Allocate (STEM(MTLOC, NSAVM))
         STEM = 0.0 !3
         Allocate (TUBER(MTLOC, NSAVM))
         TUBER = 0.0 !4
         Allocate (EP(MTLOC, NSAVM))
         EP = 0.0 !5
    !
         Allocate (SAVAREA(MTLOC))
         SAVAREA = 0.0 !6
         Allocate (SAVCELL(MTLOC))
         SAVCELL = 0 !7
    !
         Allocate (NSAVSPC(MTLOC))
         NSAVSPC = 0.0 !8
    !  ALLOCATE(SAVWDH(MTLOC,NSAVM));       SAVWDH  = 0.0         !9
         Allocate (SAVDPH(MTLOC, NSAVM))
         SAVDPH = 0.0 !10
         Allocate (SAVFRAC(MTLOC, NSAVM))
         SAVFRAC = 0.0 !11
    !
    !Production of SAV and Basal Metabolism of SAV
         Allocate (PLEAF(MTLOC, NSAVM))
         PLEAF = 0.0 !12
         Allocate (BMLEAF(MTLOC))
         BMLEAF = 0.0 !13
         Allocate (BMSTEM(MTLOC))
         BMSTEM = 0.0 !14
         Allocate (BMROOT(MTLOC))
         BMROOT = 0.0 !15
         Allocate (BMTUBER(MTLOC))
         BMTUBER = 0.0 !16
         ALLOCATE (BM_LEAF(MTLOC))    ! 
         BM_LEAF = 0.0
    !
    !Production, Basal Metabolism and Predation of Epiphytes
         Allocate (PEP(MTLOC, NSAVM))
         PEP = 0.0 !17
         Allocate (BMEP(MTLOC, NSAVM))
         BMEP = 0.0 !18
         Allocate (PREP(MTLOC, NSAVM))
         PREP = 0.0 !19
    !
    !Slough rate of Epiphytes
         Allocate (SLSH(MTLOC))
         SLSH = 0.0 !20
    !
    !Nitrogen and Phosphorus limitation to SAV growth, Fraction Nitrogen and Phosphorus uptake fron sediments by SAV
    !Nitrogen and Phosphorus limitation on Epiphyte growth
    !Sulfide HS^- in sediments toxicity to limited SAV growth
         Allocate (NLSAV(MTLOC))
         NLSAV = 0.0 !21
         Allocate (PLSAV(MTLOC))
         PLSAV = 0.0 !22
         Allocate (FNSEDSAV(MTLOC))
         FNSEDSAV = 0.0 !23
         Allocate (FPSEDSAV(MTLOC))
         FPSEDSAV = 0.0 !24
         Allocate (NLEPI(MTLOC))
         NLEPI = 0.0 !25
         Allocate (PLEPI(MTLOC))
         PLEPI = 0.0 !26
         Allocate (FHS(MTLOC))
         FHS = 0.0 !27
    !
    !Light attenuation by SAV, Epiphytes, Water
    !Light limitation to SAV (shoots) growth, Epiphyte growth
    !Net leaf, epiphytes production per day
         Allocate (SAVATN(MTLOC, NSAVM))
         SAVATN = 0.0 !28
         Allocate (EPATN(MTLOC, NSAVM))
         EPATN = 0.0 !29
         Allocate (WATATN(MTLOC, NSAVM))
         WATATN = 0.0 !30
         Allocate (FISH(MTLOC, NSAVM))
         FISH = 0.0 !31
         Allocate (FIEP(MTLOC, NSAVM))
         FIEP = 0.0 !32
         Allocate (NPPSAV(MTLOC, NSAVM))
         NPPSAV = 0.0 !33
         Allocate (NPPEPI(MTLOC, NSAVM))
         NPPEPI = 0.0 !34
    !
    !DO, LDOC, RDOC, LPOC, RPOC flux to water column due to SAV
         Allocate (DOSAVW(MTLOC))
         DOSAVW = 0.0 !35
         Allocate (LDOCSAVW(MTLOC))
         LDOCSAVW = 0.0 !36
         Allocate (RDOCSAVW(MTLOC))
         RDOCSAVW = 0.0 !37
         Allocate (LPOCSAVW(MTLOC))
         LPOCSAVW = 0.0 !38
         Allocate (RPOCSAVW(MTLOC))
         RPOCSAVW = 0.0 !39
    !NH4, NO4, LDON, RDON, LPON, RPON flux to water column due to SAV
         Allocate (NH4SAVW(MTLOC))
         NH4SAVW = 0.0 !40
         Allocate (NO3SAVW(MTLOC))
         NO3SAVW = 0.0 !41
         Allocate (LDONSAVW(MTLOC))
         LDONSAVW = 0.0 !42
         Allocate (RDONSAVW(MTLOC))
         RDONSAVW = 0.0 !43
         Allocate (LPONSAVW(MTLOC))
         LPONSAVW = 0.0 !44
         Allocate (RPONSAVW(MTLOC))
         RPONSAVW = 0.0 !45
    !PO4, LDOP, RDOP, LPOP, RPOP flux to water column due to SAV
         Allocate (PO4SAVW(MTLOC))
         PO4SAVW = 0.0 !46
         Allocate (LDOPSAVW(MTLOC))
         LDOPSAVW = 0.0 !47
         Allocate (RDOPSAVW(MTLOC))
         RDOPSAVW = 0.0 !48
         Allocate (LPOPSAVW(MTLOC))
         LPOPSAVW = 0.0 !49
         Allocate (RPOPSAVW(MTLOC))
         RPOPSAVW = 0.0 !50
    !
    !DO, LDOC, RDOC, RPOC flux to water colum due to epiphytes
         Allocate (DOEPIW(MTLOC))
         DOEPIW = 0.0 !51
         Allocate (LDOCEPIW(MTLOC))
         LDOCEPIW = 0.0 !52
         Allocate (RDOCEPIW(MTLOC))
         RDOCEPIW = 0.0 !53
         Allocate (LPOCEPIW(MTLOC))
         LPOCEPIW = 0.0 !54
         Allocate (RPOCEPIW(MTLOC))
         RPOCEPIW = 0.0 !55
    !NH4, NO3, LDON, RDON, LPON, RPON flux to water column due to epiphytes
         Allocate (NH4EPIW(MTLOC))
         NH4EPIW = 0.0 !56
         Allocate (NO3EPIW(MTLOC))
         NO3EPIW = 0.0 !57
         Allocate (LDONEPIW(MTLOC))
         LDONEPIW = 0.0 !58
         Allocate (RDONEPIW(MTLOC))
         RDONEPIW = 0.0 !59
         Allocate (LPONEPIW(MTLOC))
         LPONEPIW = 0.0 !60
         Allocate (RPONEPIW(MTLOC))
         RPONEPIW = 0.0 !61
    !PO4, LDOP, RDOP, LPOP, RPOP flux to water column due to epiphytes
         Allocate (PO4EPIW(MTLOC))
         PO4EPIW = 0.0 !62
         Allocate (LDOPEPIW(MTLOC))
         LDOPEPIW = 0.0 !63
         Allocate (RDOPEPIW(MTLOC))
         RDOPEPIW = 0.0 !64
         Allocate (LPOPEPIW(MTLOC))
         LPOPEPIW = 0.0 !65
         Allocate (RPOPEPIW(MTLOC))
         RPOPEPIW = 0.0 !66
    !
    !loading flux : POC, PON, POP (positive into sediments)      to sediments
    !return  flux : DO, NH4, PO4  (positive leaving sediments) from sediments
    !applied to sediments due to SAV
    !
         Allocate (SEDPOCSAV(MTLOC))
         SEDPOCSAV = 0.0 !67
         Allocate (SEDPONSAV(MTLOC))
         SEDPONSAV = 0.0 !68
         Allocate (SEDPOPSAV(MTLOC))
         SEDPOPSAV = 0.0 !69
         Allocate (SEDDOSAV(MTLOC))
         SEDDOSAV = 0.0 !70
         Allocate (SEDNH4SAV(MTLOC))
         SEDNH4SAV = 0.0 !71
         Allocate (SEDNO3SAV(MTLOC))
         SEDNO3SAV = 0.0 !71-1
         Allocate (SEDPO4SAV(MTLOC))
         SEDPO4SAV = 0.0 !72
    !
    !misc
         Allocate (DGRDAYS(MTLOC))
         DGRDAYS = 0.0 !73
    !
    !
    !SAV model time averaged quantities
    !
    !state varaiables
         Allocate (ALEAF(MTLOC, NSAVM))
         ALEAF = 0.0
         Allocate (AROOT(MTLOC, NSAVM))
         AROOT = 0.0
         Allocate (ASTEM(MTLOC, NSAVM))
         ASTEM = 0.0
         Allocate (ATUBER(MTLOC, NSAVM))
         ATUBER = 0.0
         Allocate (AEP(MTLOC, NSAVM))
         AEP = 0.0
    !
    !growth and mortality of SAV
         Allocate (APLEAF(MTLOC, NSAVM))
         APLEAF = 0.0
         Allocate (ABMLEAF(MTLOC))
         ABMLEAF = 0.0
         Allocate (ABMTUBER(MTLOC))
         ABMTUBER = 0.0
    !
    !growth and mortality of Epiphytes
         Allocate (APEP(MTLOC, NSAVM))
         APEP = 0.0
         Allocate (ABMEP(MTLOC, NSAVM))
         ABMEP = 0.0
         Allocate (APREP(MTLOC, NSAVM))
         APREP = 0.0
         Allocate (ASLSH(MTLOC))
         ASLSH = 0.0
    !
    !nutrient, light, HS^- limitation on SAV and Epiphytes growth
         Allocate (ANLSAV(MTLOC))
         ANLSAV = 0.0
         Allocate (APLSAV(MTLOC))
         APLSAV = 0.0
         Allocate (ANLEPI(MTLOC))
         ANLEPI = 0.0
         Allocate (APLEPI(MTLOC))
         APLEPI = 0.0
         Allocate (AFNSED(MTLOC))
         AFNSED = 0.0 !fraction uptake of N from sediment
         Allocate (AFPSED(MTLOC))
         AFPSED = 0.0 !fraction uptake of P from sediment
         Allocate (AFHS(MTLOC))
         AFHS = 0.0
    !
    !Attenuation, Light limitation and net production
         Allocate (AEPATN(MTLOC, NSAVM))
         AEPATN = 0.0
         Allocate (AWATATN(MTLOC, NSAVM))
         AWATATN = 0.0
         Allocate (AFISH(MTLOC, NSAVM))
         AFISH = 0.0
         Allocate (AFIEP(MTLOC, NSAVM))
         AFIEP = 0.0
         Allocate (ANPPSAV(MTLOC, NSAVM))
         ANPPSAV = 0.0 !net production of SAV
         Allocate (ANPPEPI(MTLOC, NSAVM))
         ANPPEPI = 0.0 !net production of Epiphytes
    !
    !time averaged flux to water column due to SAV and epiphytes
         Allocate (ADOCSAVW(MTLOC))
         ADOCSAVW = 0.0 !average of LDOCSAVW and RDOCSAVW accumulated over time
         Allocate (APOCSAVW(MTLOC))
         APOCSAVW = 0.0
         Allocate (ADOCEPIW(MTLOC))
         ADOCEPIW = 0.0
         Allocate (APOCEPIW(MTLOC))
         APOCEPIW = 0.0
    !
         Allocate (ANH4SAVW(MTLOC))
         ANH4SAVW = 0.0
         Allocate (ANO3SAVW(MTLOC))
         ANO3SAVW = 0.0
         Allocate (ADONSAVW(MTLOC))
         ADONSAVW = 0.0
         Allocate (APONSAVW(MTLOC))
         APONSAVW = 0.0
         Allocate (ANH4EPIW(MTLOC))
         ANH4EPIW = 0.0
         Allocate (ANO3EPIW(MTLOC))
         ANO3EPIW = 0.0
         Allocate (ADONEPIW(MTLOC))
         ADONEPIW = 0.0
         Allocate (APONEPIW(MTLOC))
         APONEPIW = 0.0
    !
         Allocate (APO4SAVW(MTLOC))
         APO4SAVW = 0.0
         Allocate (ADOPSAVW(MTLOC))
         ADOPSAVW = 0.0
         Allocate (APOPSAVW(MTLOC))
         APOPSAVW = 0.0
         Allocate (APO4EPIW(MTLOC))
         APO4EPIW = 0.0
         Allocate (ADOPEPIW(MTLOC))
         ADOPEPIW = 0.0
         Allocate (APOPEPIW(MTLOC))
         APOPEPIW = 0.0
    !
         Allocate (ADOSAVW(MTLOC))
         ADOSAVW = 0.0
         Allocate (ADOEPIW(MTLOC))
         ADOEPIW = 0.0
    !
    !time averaged flux to sediment due to SAV and epiphytes
         Allocate (ASEDDOSAV(MTLOC))
         ASEDDOSAV = 0.0
         Allocate (ASEDPOCSAV(MTLOC))
         ASEDPOCSAV = 0.0
         Allocate (ASEDPONSAV(MTLOC))
         ASEDPONSAV = 0.0
         Allocate (ASEDNH4SAV(MTLOC))
         ASEDNH4SAV = 0.0
         Allocate (ASEDPOPSAV(MTLOC))
         ASEDPOPSAV = 0.0
         Allocate (ASEDPO4SAV(MTLOC))
         ASEDPO4SAV = 0.0
    !
      End Subroutine SAV_ALLOC
  !
  !
  !************************************************************************
  !**             S U B R O U T I N E   SAV_LOADS_ALLOC                  **
  !************************************************************************
  !
      Subroutine SAV_LOADS_ALLOC
    !Allocate all variables related to SAV
    !much of this section is moved from wqm_modules.F
     !    Use MOD_LIMS, Only: MTLOC !AN
    !
         Allocate (DOSAVW(MTLOC))
         DOSAVW = 0.0 !35
         Allocate (LDOCSAVW(MTLOC))
         LDOCSAVW = 0.0 !36
         Allocate (RDOCSAVW(MTLOC))
         RDOCSAVW = 0.0 !37
         Allocate (LPOCSAVW(MTLOC))
         LPOCSAVW = 0.0 !38
         Allocate (RPOCSAVW(MTLOC))
         RPOCSAVW = 0.0 !39
      End Subroutine SAV_LOADS_ALLOC
  !
  !************************************************************************
  !**             S U B R O U T I N E   SAV_LOADS_DEALLOC                **
  !************************************************************************
  !
      Subroutine SAV_LOADS_DEALLOC
    !Allocate all variables related to SAV
    !much of this section is moved from wqm_modules.F
         If (ALLOCATED(DOSAVW)) DEALLOCATE (DOSAVW)!35
         If (ALLOCATED(LDOCSAVW)) DEALLOCATE (LDOCSAVW)!36
         If (ALLOCATED(RDOCSAVW)) DEALLOCATE (RDOCSAVW)!37
         If (ALLOCATED(LPOCSAVW)) DEALLOCATE (LPOCSAVW)!38
         If (ALLOCATED(RPOCSAVW)) DEALLOCATE (RPOCSAVW)!39
      End Subroutine SAV_LOADS_DEALLOC
  !
  !
  !************************************************************************
  !**             S U B R O U T I N E   SAV_DEALLOC                      **
  !************************************************************************
  !
      Subroutine SAV_DEALLOC
    !
    !Deallocation of all variables related to SAV
    !much of this section is moved from wqm_main.F
    !
         If (ALLOCATED(LEAF)) DEALLOCATE (LEAF)!1
         If (ALLOCATED(ROOT)) DEALLOCATE (ROOT)!2
         If (ALLOCATED(STEM)) DEALLOCATE (STEM)!3
         If (ALLOCATED(TUBER)) DEALLOCATE (TUBER)!4
         If (ALLOCATED(EP)) DEALLOCATE (EP)!5
    !
         If (ALLOCATED(SAVAREA)) DEALLOCATE (SAVAREA)!6
         If (ALLOCATED(SAVCELL)) DEALLOCATE (SAVCELL)!7
    !
         If (ALLOCATED(NSAVSPC)) DEALLOCATE (NSAVSPC)!8
    !IF(ALLOCATED(SAVWDH)) DEALLOCATE (SAVWDH)        !9 !this is now completely deprecated and replaced by SAVFRAC
         If (ALLOCATED(SAVDPH)) DEALLOCATE (SAVDPH)!10
         If (ALLOCATED(SAVFRAC)) DEALLOCATE (SAVFRAC)!11
    !
         If (ALLOCATED(PLEAF)) DEALLOCATE (PLEAF)!12
         If (ALLOCATED(BMLEAF)) DEALLOCATE (BMLEAF)!13
         If (ALLOCATED(BMSTEM)) DEALLOCATE (BMSTEM)!14
         If (ALLOCATED(BMROOT)) DEALLOCATE (BMROOT)!15
         If (ALLOCATED(BMTUBER)) DEALLOCATE (BMTUBER)!16
    !
         If (ALLOCATED(PEP)) DEALLOCATE (PEP)!17
         If (ALLOCATED(BMEP)) DEALLOCATE (BMEP)!18
         If (ALLOCATED(PREP)) DEALLOCATE (PREP)!19
         If (ALLOCATED(SLSH)) DEALLOCATE (SLSH)!20
    !
         If (ALLOCATED(NLSAV)) DEALLOCATE (NLSAV)!21
         If (ALLOCATED(PLSAV)) DEALLOCATE (PLSAV)!22
         If (ALLOCATED(FNSEDSAV)) DEALLOCATE (FNSEDSAV)!23
         If (ALLOCATED(FPSEDSAV)) DEALLOCATE (FPSEDSAV)!24
         If (ALLOCATED(NLEPI)) DEALLOCATE (NLEPI)!25
         If (ALLOCATED(PLEPI)) DEALLOCATE (PLEPI)!26
         If (ALLOCATED(FHS)) DEALLOCATE (FHS)!27
    !
         If (ALLOCATED(EPATN)) DEALLOCATE (EPATN)!28
         If (ALLOCATED(SAVATN)) DEALLOCATE (SAVATN)!29
         If (ALLOCATED(WATATN)) DEALLOCATE (WATATN)!30
         If (ALLOCATED(FISH)) DEALLOCATE (FISH)!31
         If (ALLOCATED(FIEP)) DEALLOCATE (FIEP)!32
         If (ALLOCATED(NPPSAV)) DEALLOCATE (NPPSAV)!33
         If (ALLOCATED(NPPEPI)) DEALLOCATE (NPPEPI)!34
    !
         If (ALLOCATED(DOSAVW)) DEALLOCATE (DOSAVW)!35
         If (ALLOCATED(LDOCSAVW)) DEALLOCATE (LDOCSAVW)!36
         If (ALLOCATED(RDOCSAVW)) DEALLOCATE (RDOCSAVW)!37
         If (ALLOCATED(LPOCSAVW)) DEALLOCATE (LPOCSAVW)!38
         If (ALLOCATED(RPOCSAVW)) DEALLOCATE (RPOCSAVW)!39
         If (ALLOCATED(NH4SAVW)) DEALLOCATE (NH4SAVW)!40
         If (ALLOCATED(NO3SAVW)) DEALLOCATE (NO3SAVW)!41
         If (ALLOCATED(LDONSAVW)) DEALLOCATE (LDONSAVW)!42
         If (ALLOCATED(RDONSAVW)) DEALLOCATE (RDONSAVW)!43
         If (ALLOCATED(LPONSAVW)) DEALLOCATE (LPONSAVW)!44
         If (ALLOCATED(RPONSAVW)) DEALLOCATE (RPONSAVW)!45
         If (ALLOCATED(PO4SAVW)) DEALLOCATE (PO4SAVW)!46
         If (ALLOCATED(LDOPSAVW)) DEALLOCATE (LDOPSAVW)!47
         If (ALLOCATED(RDOPSAVW)) DEALLOCATE (RDOPSAVW)!48
         If (ALLOCATED(LPOPSAVW)) DEALLOCATE (LPOPSAVW)!49
         If (ALLOCATED(RPOPSAVW)) DEALLOCATE (RPOPSAVW)!50
    !
         If (ALLOCATED(DOEPIW)) DEALLOCATE (DOEPIW)!51
         If (ALLOCATED(LDOCEPIW)) DEALLOCATE (LDOCEPIW)!52
         If (ALLOCATED(RDOCEPIW)) DEALLOCATE (RDOCEPIW)!53
         If (ALLOCATED(LPOCEPIW)) DEALLOCATE (LPOCEPIW)!54
         If (ALLOCATED(RPOCEPIW)) DEALLOCATE (RPOCEPIW)!55
         If (ALLOCATED(NH4EPIW)) DEALLOCATE (NH4EPIW)!56
         If (ALLOCATED(NO3EPIW)) DEALLOCATE (NO3EPIW)!57
         If (ALLOCATED(LDONEPIW)) DEALLOCATE (LDONEPIW)!58
         If (ALLOCATED(RDONEPIW)) DEALLOCATE (RDONEPIW)!59
         If (ALLOCATED(LPONEPIW)) DEALLOCATE (LPONEPIW)!60
         If (ALLOCATED(RPONEPIW)) DEALLOCATE (RPONEPIW)!61
         If (ALLOCATED(PO4EPIW)) DEALLOCATE (PO4EPIW)!62
         If (ALLOCATED(LDOPEPIW)) DEALLOCATE (LDOPEPIW)!63
         If (ALLOCATED(RDOPEPIW)) DEALLOCATE (RDOPEPIW)!64
         If (ALLOCATED(LPOPEPIW)) DEALLOCATE (LPOPEPIW)!65
         If (ALLOCATED(RPOPEPIW)) DEALLOCATE (RPOPEPIW)!66
    !
    !
         If (ALLOCATED(SEDPOCSAV)) DEALLOCATE (SEDPOCSAV)!67
         If (ALLOCATED(SEDDOSAV)) DEALLOCATE (SEDDOSAV)!68
         If (ALLOCATED(SEDNH4SAV)) DEALLOCATE (SEDNH4SAV)!69
         If (ALLOCATED(SEDNO3SAV)) DEALLOCATE (SEDNO3SAV)!69-1
         If (ALLOCATED(SEDPO4SAV)) DEALLOCATE (SEDPO4SAV)!70
         If (ALLOCATED(SEDPONSAV)) DEALLOCATE (SEDPONSAV)!71
         If (ALLOCATED(SEDPOPSAV)) DEALLOCATE (SEDPOPSAV)!72
    !
         If (ALLOCATED(DGRDAYS)) DEALLOCATE (DGRDAYS)!73
    !
    ! SAV time averaged quantities
    !
         If (ALLOCATED(AEPATN)) DEALLOCATE (AEPATN)
         If (ALLOCATED(AWATATN)) DEALLOCATE (AWATATN)
         If (ALLOCATED(AFISH)) DEALLOCATE (AFISH)
         If (ALLOCATED(AFIEP)) DEALLOCATE (AFIEP)
         If (ALLOCATED(ANPPSAV)) DEALLOCATE (ANPPSAV)
         If (ALLOCATED(ANPPEPI)) DEALLOCATE (ANPPEPI)
    !
         If (ALLOCATED(ADOSAVW)) DEALLOCATE (ADOSAVW)
         If (ALLOCATED(ADOCSAVW)) DEALLOCATE (ADOCSAVW)
         If (ALLOCATED(APOCSAVW)) DEALLOCATE (APOCSAVW)
         If (ALLOCATED(ANH4SAVW)) DEALLOCATE (ANH4SAVW)
         If (ALLOCATED(ANO3SAVW)) DEALLOCATE (ANO3SAVW)
         If (ALLOCATED(ADONSAVW)) DEALLOCATE (ADONSAVW)
         If (ALLOCATED(APONSAVW)) DEALLOCATE (APONSAVW)
         If (ALLOCATED(APO4SAVW)) DEALLOCATE (APO4SAVW)
         If (ALLOCATED(ADOPSAVW)) DEALLOCATE (ADOPSAVW)
         If (ALLOCATED(APOPSAVW)) DEALLOCATE (APOPSAVW)
    !
         If (ALLOCATED(ADOEPIW)) DEALLOCATE (ADOEPIW)
         If (ALLOCATED(ADOCEPIW)) DEALLOCATE (ADOCEPIW)
         If (ALLOCATED(APOCEPIW)) DEALLOCATE (APOCEPIW)
         If (ALLOCATED(ANH4EPIW)) DEALLOCATE (ANH4EPIW)
         If (ALLOCATED(ANO3EPIW)) DEALLOCATE (ANO3EPIW)
         If (ALLOCATED(ADONEPIW)) DEALLOCATE (ADONEPIW)
         If (ALLOCATED(APONEPIW)) DEALLOCATE (APONEPIW)
         If (ALLOCATED(APO4EPIW)) DEALLOCATE (APO4EPIW)
         If (ALLOCATED(ADOPEPIW)) DEALLOCATE (ADOPEPIW)
         If (ALLOCATED(APOPEPIW)) DEALLOCATE (APOPEPIW)
    !
         If (ALLOCATED(ASEDDOSAV)) DEALLOCATE (ASEDDOSAV)
         If (ALLOCATED(ASEDPOCSAV)) DEALLOCATE (ASEDPOCSAV)
         If (ALLOCATED(ASEDPONSAV)) DEALLOCATE (ASEDPONSAV)
         If (ALLOCATED(ASEDPOPSAV)) DEALLOCATE (ASEDPOPSAV)
         If (ALLOCATED(ASEDNH4SAV)) DEALLOCATE (ASEDNH4SAV)
         If (ALLOCATED(ASEDPO4SAV)) DEALLOCATE (ASEDPO4SAV)
    !
    !
      End Subroutine SAV_DEALLOC
  !
  !
  !************************************************************************
  !**             S U B R O U T I N E   S A V _ R E A D                  **
  !************************************************************************
  !
      Subroutine SAV_READ
    !
         Use MOD_SIZES, Only: MGL
         Use MOD_WQM, Only: B, SVOFN
    !
         Use MOD_FILEINFO, Only: SVI, SVO
    !
         Use MOD_LIMS, Only: MLOC
    !
    !Took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
         Use MOD_CONTROL, Only: SERIAL, MSR, PAR
    !
#  if defined (MULTIPROCESSOR)
         Use MOD_PAR, Only: NLID
#  endif
         Implicit None
         Character (Len=72) :: TITLE (6)
         Character (Len=8) :: TIMVAR, PRINTIT
         ! Real (SP) :: KTPS1, KTPS2, KTBMSAV, KTPE1, KTPE2, KTBME, KTPRE !temperature control coefficient for epiphyte predation loss rate [1/degC]

         Integer :: I, J, N, K ! : N is never used here
    !
         Integer :: NSAV !Number of SAV species
    !
    !     Real (SP) :: TRPMSAV, TRBMSAV, TRPMEPI, TRBMEPI, TRPREPI, &
    !    & TLOOK !Tempoary variable for temperature
    !
    !Added following temporary variables
         Integer, Allocatable :: SAVCELL_GL (:), NSAVSPC_GL (:)!number of species in each global TCE
    !
         Real (SP), Allocatable :: SAVFRAC_GL (:, :)!global SAV coverage in global TCE
         Real (SP), Allocatable :: SAVDPH_GL (:, :)!global SAV depth in global TCE
    !
         !Integer :: NSAVCELL_GL !global total number of SAV TCE's

		 Character (Len=120) :: ERROR_MSG

         Integer NSAVM_TMP
    !
    ! TITLE CARDS
    !
         Read (SVI, 1000)
    !1000 FORMAT(///)          !go to the end of 4th line
    !
         Read (SVI, 1010) (TITLE(J), J=1, 6)!read one line
    !1010 FORMAT(A72)
    ! READ SAV PROPERTIES
    !
         Read (SVI, '(//8X,I8)') NSAVM_TMP !go down two lines and read the third line
    !
    !1040 FORMAT(//8X,2A8)    !go down two lines and read the third line
         Read (SVI, 1040) SPNAM !read name of SAV species
    !growth rate, respiration ratel, basal metabolism, DO consumption fraction, sulfide toxicity coeffient on SAV
    !cut-off salinity for salinity toxicity on freshwater species of SAV. (ramp up respiration (mortality) with high salt)
         Read (SVI, 1030) PMSAV, PRSPSAVB, BMSAV, BMTBRREF, FDOSR, &
        & KHSTOX, SALMAX
    !C, N, P to SAV dryweight ratios (gC/gDW), (gN/gDW), (gP/gDW)
         Read (SVI, 1030) ACDWSAV, ANDWSAV, APDWSAV
    !1030 FORMAT(//(8X,9F8.0))  !go down two lines and read the third line
    !
    ! COMPUTE ADDITIONAL STOICHIOMETRY OF SAV
    !
         ANCSAV = ANDWSAV / ACDWSAV !N to C ratio  gN/gC in SAV (WL: should it be different for leaf, root and epiphytes ?)
         APCSAV = APDWSAV / ACDWSAV !P to C ratio  gP/gC in SAV (WL: should it be different for leaf, root and epiphytes ?)
    !
    !Half saturation constant for N uptake by leaf (shoots) (gN/m^3)
    !Threshold for NH4 in calculating NH4 preference over NO3 in N uptake(gN/m^3)
    !Half saturation constant for N uptake by root (gN/m^3)
    !Half saturation constant for P uptake by leaf (KHPLEAF) and root (KHPROOT) (gP/m^3)
    !
         Read (SVI, 1030) KHNLEAF, KHNPSAV, KHNROOT, KHPLEAF, KHPROOT !go down two lines and read the third line
    !
    !PI curve slope alpha for SAV, ALPHSAV, [(gC/gDW)/(Einstein/m^2)]
    !Light attenuation due to SAV leaf and stem, KESAV (gC/m^3)
    !ALAC = carbon weight per unit leaf area (m^2-leaf-area/gC-leaf)
    !ACAN = basic canopy height (m)
    !BCAN = effect of SAV leaf and stem biomass on canopy height (m^3/gC)
    !
         Read (SVI, 1030) ALPHSAV, KESAV, ALAC, ACAN, BCAN !go down two lines and read the third line
    !
    !Optimal temperature for SAV growth (degC) and reference temperature for SAV basal metabolism
         Read (SVI, 1030) TRPMSAV, TRBMSAV !go down two lines and read the third line
    !
    !effect of temperature below optimal temperature and above optimal temperature for SAV growth
    !effect of temperature on SAV basal metabolism
         Read (SVI, 1030) KTPS1, KTPS2, KTBMSAV !go down two lines and read the third line
    !
    ! READ PIECEWISE FUNCTIONS DESCRIBING PARTITIONING OF PRODUCTION,
    ! TRANSFER FROM TUBERS TO LEAVES
    ! VALUES ARE ENTERED AT INTERVALS OF 10 DEGREE DAYS FROM ZERO TO
    ! SENESCENCE
    !
    !number of data entries for FPLEAF, FPSTEM, FPROOT, FPTUBER, TRTBRLF
    !if these change over time of the year (temperature). It is specified every 1/10 of DGRDAYS
    !where DGRDAYS is caluclated as TDIFF*DLTDY, and TDIFF = max(T-3.0,0)  !
    !                                                    !conversion of temperature to days
    !                                                    !spring and autumn have same temperature
    !                                                    !but may not have same entries for FPLEAF etc here
    !
         Read (SVI, 1042) NDDI, PRINTIT
    !1042 FORMAT(//8X,I8,A8)                   !go down two lines and read the third
         Read (SVI, 1032) (FPLEAF(J), FPSTEM(J), FPROOT(J), FPTUBER(J), &
        & TRTBRLF(J), J=0, NDDI-1)
    !1032 FORMAT(//(8X,5F8.0))  !go down two lines and read the third on ward until finish
    ! READ DISTRIBUTION OF SAV BIOMASS AFTER MORTALITY
    !
    !Nitrogen distribution to NH4, LDON, RDON, LPON, RPON
         Read (SVI, 1030) FNISAV, FNLDSAV, FNRDSAV, FNLPSAV, FNRPSAV ! go down two lines and read the third line
    !
    !3G partitioning of PON flux into sediments
         Read (SVI, 1030) FRPONSAV (1), FRPONSAV (2), FRPONSAV (3)! go down two lines and read the third line
    !
    !Phosphorus distribution to PO4, LDOP, RDOP, LPOL, RPOP
         Read (SVI, 1030) FPISAV, FPLDSAV, FPRDSAV, FPLPSAV, FPRPSAV ! go down two lines and read the third line
    !
    !3G partitioning of POP flux into sediments
         Read (SVI, 1030) FRPOPSAV (1), FRPOPSAV (2), FRPOPSAV (3)! go down two lines and read the third line
    !
    !carbon distribution into DO consumption (CO2), LDOC, RDOC, LPOC, RPOC
         Read (SVI, 1030) FDOSAV, FCLDSAV, FCRDSAV, FCLPSAV, FCRPSAV ! go down two lines and read the third line
    !
    !3G partitioning of POC flux into sediments
         Read (SVI, 1030) FRPOCSAV (1), FRPOCSAV (2), FRPOCSAV (3)! go down two lines and read the third line
    !
    ! READ EFFECT OF SAV ON NET SETTLING TO SEDIMENTS
    !
    !fixed solids, LPOM, RPOM, algae1, algae 2, algae 3, and particulate biogenic silica
         Read (SVI, 1030) WSSSAV, WSLSAV, WSRSAV, WS1SAV, WS2SAV, &
        & WS3SAV, WSUSAV ! go down two lines and read the third line
    !
    ! READ EPIPHYTE NAME AND PROPERTIES
    !Name of epipytes
    !1040 FORMAT(//8X,2A8)
         Read (SVI, 1040) EPINAM !go down two lines and read the third
    !
    !growth rate, photoresprication rate, basal metabolism and predation rate
         Read (SVI, 1030) PMEPI, PRSPEPI, BMEPI, PREPI ! go down two lines and read the third
    !
    !C to CHLA ratio, DW to C ratio, N to C ratio, P to C ratio of Epiphytes
         Read (SVI, 1030) CCHLEPI, ADWCEPI, ANCEPI, APCEPI !go down two lines and read the third
    !
    !Nitrogen uptake half constant, NH4 preference constant, Phosphors uptake half constant
         Read (SVI, 1030) KHNEPI, KHNPEPI, KHPEPI !go down two lines and read the third
    !
    !PI curve slope for Epiphyte, Light attenuation of epiphyte
         Read (SVI, 1030) ALPHEPI, KEEPI !go down two lines and read the third
    !
    !primary production optimal temperature, basal metabolism reference termperature, predation for epiphytes
         Read (SVI, 1030) TRPMEPI, TRBMEPI, TRPREPI !go down two lines and read the third
    !
    !temperature effect of growth (production) below and above optimal temperature TRPMEPI,
    !temperature effect on basal metabolism and predation for epiphytes
         Read (SVI, 1030) KTPE1, KTPE2, KTBME, KTPRE !go down two lines and read the third
    !
    !half saturation density of epiphyte for epiphyte growth
         Read (SVI, 1030) KHEP !go down two lines and read the third
    !
    ! READ DISTRIBUTION OF EPIPHYTE BIOMASS AFTER MORTALITY
    !
    !nitrogen distribution to inoragnic N (NH4), LDON, RDON, LPON, RPON upon basal metabolism loss
         Read (SVI, 1030) FNIEPI, FNLDEPI, FNRDEPI, FNLPEPI, FNRPEPI !go down two lines and read the third
    !
    !nitrogen distribution to inorganic N (NH4), LDON, RDON, LPON, RPON upon predation loss
         Read (SVI, 1030) FNIPEP, FNLDPEP, FNRDPEP, FNLPPEP, FNRPPEP !go down two lines and read the third
    !
    !phosphorus distribution to inorganic P (PO4), LDOP, RDOP, LPOP, RPOP upon basal metabolism loss
         Read (SVI, 1030) FPIEPI, FPLDEPI, FPRDEPI, FPLPEPI, FPRPEPI !go down two lines and read the third
    !
    !phosphorus distribution to inorganic P (PO4), LDOP, RDOP, LPOP, RPOP upon predation loss
         Read (SVI, 1030) FPIPEP, FPLDPEP, FPRDPEP, FPLPPEP, FPRPPEP !go down two lines and read the third
    !
    !Carbon distribution to LDOC, RDOC up on basal metabolism
         Read (SVI, 1030) FCLDEPI, FCRDEPI !go down two lines and read the third
    !Potential future improvements
    !Missing FCLPEPI and FCRPEPI to routt to LPOC and RPOC (not coded)
    !FCLPEPI and FCRPEPI should be added (),as these are coded for nitrogen
    !and there is no reason not to do the same for carbon. Also below LPOP and RPOP
    !have their share of carbon from predation loss. No reason that metabolism should be
    !different
    !
    !also missing FDOEPI compared to FDOSAV for SAV
    !code assumes that no carbon will be routed to LPOP and RPOP due to metabolism
    !and also FDOEPI=1-FCLDEPI-FCRDEPI, and FDOEPI is only implicitly used in calculating
    !DOEPIW
    !
    !Carbon distribution to LDOP, RDOP, LPOP, RPOP up on predation loss for epiphytes
         Read (SVI, 1030) FCLDPEP, FCRDPEP, FCLPPEP, FCRPPEP !go down two lines and read the third
    !
    ! READ NUMBER OF CELLS WITH SAV, MAXIMUM NUMBER OF DEPTH INCREMENTS
    !
         Read (SVI, 1020) NSAVCELL_GL, NSAVSPCM !go down two lines and read the third
    !: Note that NSAVSPCM is never used
    !1020 FORMAT(//(8X,9I8))
    !
    ! READ NUMBER OF DEPTH INCREMENTS, CHARACTERISTIC WIDTHS, AND DEPTHS
    !
         Do I = 1, MLOC
            SAVCELL (I) = 0
            NSAVSPC (I) = 0
         End Do
         Do J = 1, NSAVM
            Do I = 1, MLOC
          !       SAVWDH(I,J) = 0.0
               SAVFRAC (I, J) = 0.0 !: replace SAVWDH with SAVFRAC
               SAVDPH (I, J) = 0.0
            End Do
         End Do
    !
         Read (SVI, 1006)
    !1006 FORMAT(/)
    !
    !
    !
         Allocate (SAVCELL_GL(NSAVCELL_GL))
         SAVCELL_GL = 0
         ALLOCATE(SAVNDS_GL(NSAVCELL_GL))
         SAVNDS_GL = 0
    !global TCE ID
         Allocate (NSAVSPC_GL(MGL))
         NSAVSPC_GL = 0
    !number of species in each global TCE
         Allocate (SAVFRAC_GL(MGL, NSAVM))
         SAVFRAC_GL = 0.0
    !global SAV coverage in global TCE
         Allocate (SAVDPH_GL(MGL, NSAVM))
         SAVDPH_GL = 0.0
    !global SAV depth in global TCE
    !
         Do I = 1, NSAVCELL_GL !Loop through all SAV cells
            Read (SVI,*) NSAV, B
			If (NSAV > NSAVM) Then
               ERROR_MSG = 'Oops, NSAV has to be less than&
              & NSAVM!!'
               Call ERROR_CHECK (ERROR_MSG,-1)
            End If
            SAVCELL_GL (I) = B !The overall box number of I'th SAV cell
            NSAVSPC_GL (B) = NSAV !Number of SAV species in box B 
            SAVNDS_GL (I) = B  ! Storing the global SAV-nodes for initialization
       !
            Read (SVI,*) (SAVFRAC_GL(B, J), J=1, NSAVSPC_GL(B))!Fraction of SAV coverage in the projected surface area of the
       !bottom grid cell
            Read (SVI,*) (SAVDPH_GL(B, J), J=1, NSAVSPC_GL(B))!Depth of SAV (mean depth of the plant) for all
       !species in each grid cell (should replace with
       !still water depth of the SAV cell as we no longer
       !use fixed z grid that was in the old ICM)
         End Do
    !
         If (SERIAL) Then
            NSAVCELL = NSAVCELL_GL
            NSAVSPC = NSAVSPC_GL
            SAVFRAC = SAVFRAC_GL
            SAVDPH = SAVDPH_GL
         End If
    !
# if defined (MULTIPROCESSOR)
         If (PAR) Then
            NSAVCELL = 0
            Do I = 1, NSAVCELL_GL
          !find global node
               B = SAVCELL_GL (I)!global TCE (node) ID
               J = NLID (B)!local TCE (node) ID
               If (J .Ne. 0) Then !find if B is in local domain
             !accumulate number of local NSAVCELL
                  NSAVCELL = NSAVCELL + 1
             !local TCE ID
                  SAVCELL (NSAVCELL) = J
             !set the local SAVFRAC and SAVDPH from the global value
                  NSAVSPC (J) = NSAVSPC_GL (B)!number of species
                  SAVFRAC (J, 1:NSAVSPC(J)) = SAVFRAC_GL (B, &
                 & 1:NSAVSPC_GL(B))!fraction of coverage in TCE J
                  SAVDPH (J, 1:NSAVSPC(J)) = SAVDPH_GL (B, &
                 & 1:NSAVSPC_GL(B))!depth of SAV in TCE J
               End If
            End Do
       !
       !
         End If
# endif
    !
    !***** Input FORMAT statements
    !
1000     Format (/ //)
1006     Format (/)
1010     Format (A72)
1020     Format (/ / (8 X, 9 I8))
1030     Format (/ / (8 X, 9 F8.0))
1032     Format (/ / (8 X, 5 F8.0))
1040     Format (/ / 8 X, 2 A8)
1042     Format (/ / 8 X, I8, A8)
1050     Format (/ / (8 X, A8, F8.0))
    !
    ! OUTPUT WHAT WAS INPUT
    !
         If (MSR) Then
       !
            Open (SVO, File=SVOFN)
            Write (SVO, 2010) (TITLE(J), J=1, 6)
       !
       ! WRITE SAV SPECIES NAME AND PROPERTIES
       !
            Write (SVO, 2040) SPNAM
            Write (SVO, 2050) PMSAV, TRPMSAV, KTPS1, KTPS2
            Write (SVO, 2052) PMSAV / ACDWSAV, TRPMSAV
            Write (SVO, 2056) ACDWSAV, ANDWSAV, APDWSAV
            Write (SVO, 4010) KHNLEAF, KHNPSAV, KHNROOT, KHPLEAF, &
           & KHPROOT
            Write (SVO, 2060) BMSAV, TRBMSAV, KTBMSAV
            Write (SVO, 2058) PRSPSAVB, SALMAX
            Write (SVO, 3070) FDOSR
            Write (SVO, 2090) KHSTOX
            Write (SVO, 3000) ALPHSAV
            Write (SVO, 3010) KESAV, ALAC, ACAN, BCAN
			write (SVO, 2030) PRINTIT
            If (PRINTIT == '     ALL') Then
               Write (SVO, 3080)
               Write (SVO, 3090) (J, FPLEAF(J), FPSTEM(J), FPROOT(J), &
              & FPTUBER(J), TRTBRLF(J), J=1, NDDI)
            End If
       !
       ! WRITE DISTRIBUTION OF BIOMASS UPON MORTALITY
       !
            Write (SVO, 4020)
            Write (SVO, 4030) FNISAV, FNLDSAV, FNRDSAV, FNLPSAV, &
           & FNRPSAV
            Write (SVO, 4040) FPISAV, FPLDSAV, FPRDSAV, FPLPSAV, &
           & FPRPSAV
            Write (SVO, 4050) FDOSAV, FCLDSAV, FCRDSAV, FCLPSAV, &
           & FCRPSAV
            Write (SVO, 4080)
            Write (SVO, 4090) FRPONSAV (1), FRPONSAV (2), FRPONSAV (3)
            Write (SVO, 5000) FRPOPSAV (1), FRPOPSAV (2), FRPOPSAV (3)
            Write (SVO, 5010) FRPOCSAV (1), FRPOCSAV (2), FRPOCSAV (3)
       !
       ! WRITE EFFECT OF SAV ON NET SETTLING TO SEDIMENTS
            Write (SVO, 5060) WSSSAV, WSLSAV, WSRSAV, WS1SAV, WS2SAV, &
           & WS3SAV, WSUSAV
       !
       ! WRITE EPIPHYTE PROPERTIES
       !
            Write (SVO, 2040) EPINAM
            Write (SVO, 2051) PMEPI, TRPMEPI, KTPE1, KTPE2
            Write (SVO, 2052) PMEPI / CCHLEPI, TRPMEPI
            Write (SVO, 2062) CCHLEPI, ADWCEPI, ANCEPI, APCEPI
            Write (SVO, 5020) KHNEPI, KHNPEPI, KHPEPI
            Write (SVO, 2065) BMEPI, TRBMEPI, KTBME
            Write (SVO, 3020) PREPI, TRPREPI, KTPRE
            Write (SVO, 3002) ALPHEPI
            Write (SVO, 3030) KEEPI
            Write (SVO, 4000) KHEP
       !
       ! WRITE DISTRIBUTION OF BIOMASS UPON MORTALITY
       !
            Write (SVO, 5030)
            Write (SVO, 4030) FNIEPI, FNLDEPI, FNRDEPI, FNLPEPI, &
           & FNRPEPI
            Write (SVO, 4040) FPIEPI, FPLDEPI, FPRDEPI, FPLPEPI, &
           & FPRPEPI
            Write (SVO, 5050) FCLDEPI, FCRDEPI
            Write (SVO, 5040)
            Write (SVO, 4030) FNIPEP, FNLDPEP, FNRDPEP, FNLPPEP, &
           & FNRPPEP
            Write (SVO, 4040) FPIPEP, FPLDPEP, FPRDPEP, FPLPPEP, &
           & FPRPPEP
            Write (SVO, 4070) FCLDPEP, FCRDPEP, FCLPPEP, FCRPPEP
       !
       ! WRITE NUMBER OF SAV CELLS, SAV WIDTH AND DEPTH INCREMENTS
       !
            Write (SVO, 5070) NSAVCELL_GL, NSAVSPCM !: note that  NSAVSPCM is never used
            If (PRINTIT == '     ALL') Then
               Do I = 1, NSAVCELL_GL
                  B = SAVCELL_GL (I)
                  Write (SVO, 5080) I, B, NSAVSPC_GL (B)
                  Write (SVO, 5082) (SAVFRAC_GL(B, J), J=1, &
                 & NSAVSPC_GL(B))!
                  Write (SVO, 5090) (SAVDPH_GL(B, J), J=1, &
                 & NSAVSPC_GL(B))
               End Do
            End If
       !
         End If !END of MSR block
    !
         Deallocate (SAVCELL_GL, NSAVSPC_GL, SAVFRAC_GL, SAVDPH_GL)
    !
    ! COMPUTE FRACTIONAL AREA OF EACH SAV DEPTH INCREMENT
    !
    ! CREATE LOOKUP TABLE OF TEMPERATURE EFFECTS
    
    ! INITIALIZE DEGREE DAYS TO ZERO
    !
         Do I = 1, MLOC
            DGRDAYS (I) = 0
         End Do
    !
    ! ZERO OUT ALL COUNTERS FOR EFECTS ON WATER COLUMN AND SEDIMENTS
    !
         Do I = 1, MLOC
       !
            DOSAVW (I) = 0.
            LDOCSAVW (I) = 0.
            RDOCSAVW (I) = 0.
            LPOCSAVW (I) = 0.
            RPOCSAVW (I) = 0.
       !
            DOEPIW (I) = 0.
            LDOCEPIW (I) = 0.
            RDOCEPIW (I) = 0.
            LPOCEPIW (I) = 0.
            RPOCEPIW (I) = 0.
       !
            NH4SAVW (I) = 0.
            NO3SAVW (I) = 0.
            LDONSAVW (I) = 0.
            RDONSAVW (I) = 0.
            LPONSAVW (I) = 0.
            RPONSAVW (I) = 0.
       !
            NH4EPIW (I) = 0.
            NO3EPIW (I) = 0.
            RDONEPIW (I) = 0.
            LDONEPIW (I) = 0.
            LPONEPIW (I) = 0.
            RPONEPIW (I) = 0.
       !
            PO4SAVW (I) = 0.
            LDOPSAVW (I) = 0.
            RDOPSAVW (I) = 0.
            LPOPSAVW (I) = 0.
            RPOPSAVW (I) = 0.
       !
            PO4EPIW (I) = 0.
            LDOPEPIW (I) = 0.
            RDOPEPIW (I) = 0.
            LPOPEPIW (I) = 0.
            RPOPEPIW (I) = 0.
       !
         End Do
    !
         Do I = 1, MLOC
       !
            SEDDOSAV (I) = 0. !Dissolved oxygen flux leaving sediments due to SAV and epiphytes
            SEDNH4SAV (I) = 0. !NH4 flux leaving sediments due to SAV and epiphytes
            SEDNO3SAV (I) = 0. !NO3 flux (Added NO3 here) leaving sediments due to SAV and epiphytes
            SEDPO4SAV (I) = 0. !PO4 flux leaving sediments due to SAV
       !
            SEDPOCSAV (I) = 0. !POC flux into sediments due to SAV and epiphytes
            SEDPONSAV (I) = 0. !PON flux into sediments due to SAV and epiphytes
            SEDPOPSAV (I) = 0. !POP flux into sediments due to SAV and epiphytes
       !
         End Do
    !
    !***** Output FORMAT statements
    !
2010     Format (1 X, A72)
2030     Format (/ ' PRINTIT---', A8,'----')
2040     Format (/ 1 X, A8)
2050     Format (' PRODUCTION = ', F8.3, ' GM C/GM DW/DAY AT ', F8.2, '&
        & C.'/' KT1 = ', F8.3, ' KT2 = ', F8.3, ' PER DEGREE**2')
2051     Format (' PRODUCTION = ', F8.3, ' GM C/GM CHL/DAY AT ', F8.2, &
        & ' C.'/' KT1 = ', F8.3, ' KT2 = ', F8.3, ' PER DEGREE**2')
2052     Format (' CARBON SPECIFIC GROWTH RATE = ', F8.3, ' PER DAY AT &
        &', F8.2, ' C.')
2056     Format (' CARBON TO DRY WEIGHT RATIO = ', F8.3/' NITROGEN TO D&
        &RY WEIGHT RATIO = ', F8.3/' PHOSPHORUS TO DRY WEIGHT RATIO = '&
        & , F8.3)
2058     Format (' PHOTORESPIRATION = ', F8.3, ' * PRODUCTION @ ZERO SA&
        &LT'/' PHOTORESPIRATION = 1.0 WHEN SALINITY > ', F8.3)
2060     Format (/ ' PLANT RESPIRATION = ', F8.3, ' PER DAY AT ', F8.2, &
        & ' C.  KT = ', F8.3, ' PER DEGREE')
2062     Format (' CARBON TO CHLOROPHYLL RATIO = ', F8.3/' DRY WEIGHT T&
        &O CARBON RATIO = ', F8.3/' NITROGEN TO CARBON RATIO = ', &
        & F8.3/' PHOSPHORUS TO CARBON RATIO = ', F8.3)
2065     Format (/ ' RESPIRATION = ', F8.3, ' PER DAY AT ', F8.2, ' C. &
        & KT = ', F8.3, ' PER DEGREE')
2090     Format (' SULFIDE TOXICITY EFFECT = ', F8.3, ' PER MOL HS')
3000     Format (/ ' ALPHA = ', F8.3, ' (GM C/GM DW) / (E/M**2)')
3002     Format (/ ' ALPHA = ', F8.3, ' (GM C/GM CHL) / (E/M**2)')
3010     Format (' LIGHT ATTENUATION = ', F8.3, ' M**2/GM SHOOT C',/' M&
        &**2 LEAF AREA PER GM LEAF CARBON  = ', F8.3,/' CANOPY HEIGHT =&
        & ', F8.3, ' + ', F8.3, ' * SHOOTS (M)')
3020     Format (' PREDATION = ', F8.3, ' PER DAY AT ', F8.2, ' C.  KT &
        &= ', F8.3, ' PER DEGREE')
3030     Format (' LIGHT ATTENUATION = ', F8.3, ' M**2 LEAF AREA PER GM&
        & EPIPHYTE C' /)
3040     Format (' SPECIES SPECIFIED INCORRECTLY IN CELL ', I5, ' AS ', &
        & A8)
3050     Format (/ 3 (' CELL SPECIES CODE BLSAV'))
3060     Format (3(I5, 1 X, A8, I4, F8.2))
3070     Format (' FRACTION DO TRANSFERRED FROM SHOOTS TO ROOTS =', &
        & F8.3)
3080     Format (/ '  DGR DAY  FPLEAF   FPSTEM   FPROOT  FPTUBER TRTBRL&
        &F' / '           0<  <1   0<  <1   0<  <1  0<  <1  1/DAY')
3090     Format (I7, 3 X, 5 F10.3)
4000     Format (/ ' EPIPHYTE DENSITY AT WHICH GROWTH IS HALVED = ', &
        & F8.3, ' GM C/GM C')
4010     Format (/ ' NITROGEN HALF-SATURATION CONC (SHOOTS) = ', F8.3, &
        & ' GM/M**3' / ' NH4 HALF-SATURATION CONC (PREFERENCE) = ', &
        & F8.3, ' GM/M**3' / ' NITROGEN HALF-SATURATION CONC (ROOTS) = &
        &', F8.3, ' GM/M**3' / ' PHOSPHORUS HALF-SATURATION CONC (SHOOT&
        &S) = ', F8.3, ' GM/M**3' / ' PHOSPHORUS HALF-SATURATION CONC (&
        &ROOTS) = ', F8.3, ' GM/M**3')
4020     Format (/ ' DISTRIBUTION OF SAV UPON MORTALITY' / ' SHOOT RESP&
        &IRATION  DIS INORG  LAB DISS  REF DISS  ', 'LAB PART  REF PART&
        &')
4030     Format (' NITROGEN          ', 5 F10.3)
4040     Format (' PHOSPHORUS        ', 5 F10.3)
4050     Format (' CARBON            ', 5 F10.3)
4070     Format (' CARBON            ', 10 X, 4 F10.3)
4080     Format (/ ' ROOT MORTALITY   LAB PART  REF PART     INERT')
4090     Format (' NITROGEN       ', 3 F10.3)
5000     Format (' PHOSPHORUS     ', 3 F10.3)
5010     Format (' CARBON         ', 3 F10.3)
5020     Format (/ ' NITROGEN HALF-SATURATION CONC = ', F8.3, ' GM/M**3&
        &' / ' NH4 HALF-SATURATION CONC (PREFERENCE) = ', F8.3, ' GM/M*&
        &*3' / ' PHOSPHORUS HALF-SATURATION CONC = ', F8.3, ' GM/M**3')
5030     Format (/ ' DISTRIBUTION OF EPIPHYTES UPON MORTALITY' / ' RESP&
        &IRATION        DIS INORG  LAB DISS  REF DISS  ', 'LAB PART  RE&
        &F PART')
5040     Format (' PREDATION')
5050     Format (' CARBON            ', 10 X, 2 F10.3)
5060     Format (/ 'EFFECT OF SAV ON NET SETTLING TO SEDIMENTS ', '(M**&
        &3/GM C/DAY)' / '  SOLIDS  LABILE REFRACT      B1      B2', '  &
        &    B3     PBS' / 7 F8.2)
5070     Format (/ ' NUMBER OF SAV CELLS = ', I8, / ' MAXIMUM NUMBER OF&
        & DEPTHS = ', I8)
    !5080 FORMAT(/' SAV CELL ',I6,', WQM CELL ',I6, ', # DEPTH INCR ',I6)
5080     Format (/ ' SAV CELL ', I6, ', WQM CELL ', I6, ', # Species ', &
        & I6)
5082     Format (' Coverage Fraction ', 10 F8.1)!Changed from Width (m) to Coverage (dimensionless)
5090     Format (' DEPTH (M) ', 10 F8.2)
    !
         Return
      End Subroutine SAV_READ
  !
  !************************************************************************
  !**                    S U B R O U T I N E   S A V                     **
  !************************************************************************
  !
      Subroutine SAV_COMP
         Use MOD_LIMS, Only: MLOC, KBM1
         Use MOD_WQM, Only: AOCR, B, DLT, I0, JDAY, KADPO4, NH4, NO3, &
        & PO4, SALT, SREFLECT, SSI, T !,             &!
    !
         Use MOD_OWQ, Only: KESS
         Use MOD_SED_SAV_EXCHANGE_VARS, Only: NH4T2TM1S_SHARE, &
        & PO4T2TM1S_SHARE, HST2TM1S_SHARE, M2_SHARE, PIE2HS_SHARE
         Use MOD_CONTROL, Only: PAR, MSR
         Use MOD_HYDROVARS, Only: EL 
         Implicit None
         Integer :: IDGD (MLOC)!Degree day (basically infering day of year from temperature)
    !
         Real (SP) :: DENLIM (MLOC, NSAVM), PRSPSAV (MLOC)!Photochemical respiration and salinity toxicity to salinity
    !
         Real (SP) :: IK, IATCNPY, IAVEP, IAVSH, NSED, NWAT, PSED, &
        & PWAT, LFOLD !SAV leaf biomass (gC/m^2) of previous time step
    !
         Real (SP) :: NPRSAV (MLOC), NPREPI (MLOC)!NH4 preference on nitrogen uptake by Epiphytes (dimensionless)
!
         !Integer :: ITEMP, I, J, N !AN
		 Integer :: I, J, N, K
         Real (SP) :: ATN, DF, HSMOL, STOXI, GLIMSAV, GLIMEPI, DELLF, TFTBM
!
    !  ATN,     &
    !  DF,      & !fraction of dissolved PO4 concentration relative to total PO4 cocentration in water column
    !  HSMOL,   & !HS (sulfide) concentration in sediment pore water (mole-S2/m^3)
    !             STEST,   & !salinity used in calculating salt toxicity (Changed this to STOXI below) (psu)
    !  STOXI,   & !salinity used in calculating salt toxicity (psu)
    !  GLIMSAV, & !overall (global) limitation to SAV growth (N, P, light and sulfide HS toxicity)
    !  GLIMEPI, & !overall (global) limitation to SAV growth (N, P, light)
    !  DELLF      !delta leaf, i.e. leaf biomass change (gC/m^2/day)
    !
         Real (SP) :: DLTDY, TDIFF, ZMEAN, ZTOCNPY
    !  DLTDY,  &  !time step in days (day) for updating SAV and epiphyte biomass
    !  TDIFF,  &  !temerature - 3 (degC), used to calculate DGRDAYS and then IDGD for discretizing
    !  !leaf,stem,root,tuber share of the SAV growth etc as a function of temperature to
    !  !account for change of functions due to aging of SAV (senescence)
    !  ZMEAN,  &  !SAVDPH (vertical coordinate of SAV distribution discretized vertically) (distance from surface) (m)
    !  ZTOCNPY    !Ztc in Cerco and Moore 2001 paper, depth (m) at the top of the canopy (distance from surface to
    !              canopy top
!
         Integer :: KWC !water column layer index
		 Integer :: KWL !water column layer number
    !
    ! CONVERT MODEL TIME STEP TO DAYS, DETERMINE JULIAN DAY
    !
         DLTDY = DLT / 86400.
         J = AINT(1.0 + AMOD (JDAY, 365.25))
    !
    ! COMPUTE DEGREE DAYS FOR EACH CELL, DISCRETE ARRAY ELEMENT FOR
    ! TRANSFERS ETC
    !
	     KWL=6

         Do I = 1, NSAVCELL
            B = SAVCELL (I)

			AveT=0.0_SP
			DO K=1, KWL
			   AveT=AveT+T(B, K)
			END DO
			AveT=AveT/KWL

            ! TDIFF = Max (T(B, KBM1)-3., 0.)!Changed to KBM1
			  TDIFF = Max (AveT-3., 0.)!Wen Long changed to KBM1 !AN change to surface
       !when T< 3, TDIFF is bounded from below by 0
       !so that when T changes over time, DGRDAYS keeps increasing
       !
            DGRDAYS (B) = DGRDAYS (B) + TDIFF * DLTDY
            IDGD (B) = DGRDAYS (B) / 10. !Index in FPLEAF, FPSTEM, FPROOT, FPTUBER, TRTBRLF
       !which are functions of DGRDAYS/10
       !Here DGRDAYS mean integration of tempreature in time
       ! i.e. the accumulative history of temperature that the
       ! plant has experienced and endured/enjoyed, which units [degC day]

       !
       !The look up table FPLEAF, FPSTEM, FPROOT, FPTUBER, TRTBRLF
       !has to be large enough to contain the whole life history
       !in [degC*day] units of the plant
       !when DGRDAYS --> infinity as time increases, the look up table
       !should asymptotically approach a set of limiting values corresponding
       !to the dying phase of the plant.

       !Laugier T. et. al., Seasonal dynamics in mixed eelgrass beds,
       !Zostera marina L. and Z. noltii Hornem.,in a Mediterranean coastal lagoon
       !(Thau lagoon, France) Aquatic Botany, 63, pp 51-69, 1999 gives log distribution
       !of mortality as function of shoot and root density
       !
       !Shoot recruitment should also be a function of age and a parameter
       !to the model.  see Laugier paper Table 5 above.
       !Shoot turn over (1/shoot-life [1/year]) should also be modeled
       !separately from stem and root
       !
       !Tuber is an interesting component, it stores nutrients for next
       !year or when nutrients are depleted. So its function should
       !aslo be modeled differently.  See reference : Harwell and Orth, Influence
       ! of a tube-dwelling polychaete on the dispersal of fragmented
       ! reproductive shoots of eelgrass, Aquatic Botany, 70, pp. 1-7, 2001
       !
       !When rhizomes become over crowded, the plant stops flowering too.
       !
       !Alternatively, we model mortality difference through size of patches
       !
       !
         End Do
    !
    ! ZERO OUT AFFECTS OF SAV CELLS ON WATER COLUMN AND SEDIMENTS
    !
         Do I = 1, NSAVCELL
            B = SAVCELL (I)
            DOSAVW (B) = 0.
            LDOCSAVW (B) = 0.
            RDOCSAVW (B) = 0.
            LPOCSAVW (B) = 0.
            RPOCSAVW (B) = 0.
       !
            DOEPIW (B) = 0.
            LDOCEPIW (B) = 0.
            RDOCEPIW (B) = 0.
            LPOCEPIW (B) = 0.
            RPOCEPIW (B) = 0.
       !
            NH4SAVW (B) = 0.
            NO3SAVW (B) = 0.
            LDONSAVW (B) = 0.
            RDONSAVW (B) = 0.
            LPONSAVW (B) = 0.
            RPONSAVW (B) = 0.
       !
            NH4EPIW (B) = 0.
            NO3EPIW (B) = 0.
            RDONEPIW (B) = 0.
            LDONEPIW (B) = 0.
            LPONEPIW (B) = 0.
            RPONEPIW (B) = 0.
       !
            PO4SAVW (B) = 0.
            LDOPSAVW (B) = 0.
            RDOPSAVW (B) = 0.
            LPOPSAVW (B) = 0.
            RPOPSAVW (B) = 0.
       !
            PO4EPIW (B) = 0.
            LDOPEPIW (B) = 0.
            RDOPEPIW (B) = 0.
            LPOPEPIW (B) = 0.
            RPOPEPIW (B) = 0.
       !
         End Do
    !
         Do I = 1, NSAVCELL
            B = SAVCELL (I)
            SEDDOSAV (B) = 0.
            SEDPOCSAV (B) = 0.
            SEDNH4SAV (B) = 0.
            SEDNO3SAV (B) = 0.
            SEDPO4SAV (B) = 0.
            SEDPONSAV (B) = 0.
            SEDPOPSAV (B) = 0.
         End Do
    !
    ! LIGHT EFFECTS
    !
    !     KWC = KBM1 !water column layer (bottom layer)
		 KWL = 6

    !
         Do I = 1, NSAVCELL
            B = SAVCELL (I)

			!AveT=0.0_SP

			!DO K=1, KWL
			   !AveT=AveT+T(B, K)
			!END DO
		!	AveT=AveT/KWL

      AveT = SUM(T(B,1:6))/6.0

       !
            !ITEMP = 10. * T (B, KBM1) + 0.05 !Here T should be T(B,KBM1), i.e. bottom temperature !AN
       !
            Do N = 1, NSAVSPC (B)
          !
          ! COMPUTE IRRADIANCE AT CANOPY HEIGHT
          !
          !    Titus and Adams 1979, "Coexistence and comparative light relations of
          !          the submersed macrophytes Myriophyllum spicatum L. and Vallisneria americana Michx",
          !          Oecologia (Berl.) 40, 273-286)
          !    Titus, J.E., Goldstein, R.A., Adams, M.S., Mankin, J.B., O'Neill, R.V., Weiler, P.R. Jr., Shugart,
          !           H.H., Booth, R.S. : A production model for Myriophyllum spicatum L. Ecology 56, 1129-1138
          !
          !     ZMEAN = SAVDPH (B, N)!depth of SAV veritcal increment (mean depth level of canopy (m from surface)?)
               ZMEAN = SAVDPH (B, N)+EL (B)-1.5_SP
          !     HCAN = ACAN + BCAN * (LEAF(B, N)+STEM(B, N))!canopy height increases (m) with leaf and stem growth
          !
               HCAN = (ACAN + BCAN * (LEAF(B, N)+STEM(B, N)))*0.01 ! cm to m - AN
               HCAN = Min (ZMEAN, HCAN)!make sure ZMEAN > HCAN
          !so that ZTOCNPY will be positive below
          !     ZTOCNPY = ZMEAN - HCAN !depth level of canopy top (m from surface), always positive
               ZTOCNPY = Max(ZMEAN - HCAN, 0.0_SP) !depth level of canopy top (m from surface), always positive
          !     WATATN(B,N) = EXP(-KESS(B,1)*ZTOCNPY)       !attenuation due to water column

          !     WATATN (B, N) = Exp (-KESS(B, KBM1)*ZTOCNPY)!attenuation due to water column
               WATATN (B, N) = Exp (-KESS(B, 1)*ZTOCNPY)!attenuation due to water column
          !
               IATCNPY = I0 * (1.-SREFLECT/100.) * WATATN (B, N)!Ic in Cerco and Moore 2001, light at canopy top
          !
          ! COMPUTE ATTENUATION BY EPIPHYTES AND SELF SHADING
          !
          ! EPATN(B,N) = EXP(-KEEPI*ACLA*ADWCEPI*EP(B,N))  !Old version with ACLA= gC-leaf/m^2-leaf, EP=gC-epi/gC-leaf
          !Cerco & Moore,
          !Estuaries Vol. 24, No. 4, p. 522-534 August 2001
          !equation for Ish
          !KEEPI unit:           m^2 leaf/gDW-epi
          !ADWCEPI unit:         gDW-epi/gC-epi
          !EP :                   gC-epi/m^2-leaf
          !---------------------------------------------------
          !KEEPI*ADWCEP*EP :      1    (nondimensional)  !check !
          !
               EPATN (B, N) = Exp (-KEEPI*ADWCEPI*EP(B, N))!New version now EP is in unit gC-epi/m^2-leaf
          !
          !     ATN = KESS (B, KBM1) * HCAN + KESAV * (LEAF(B, &
          !    & N)+STEM(B, N))
               ATN = KESS (B, 1) * HCAN + KESAV * (LEAF(B, &
              & N)+STEM(B, N)) !AN
               SAVATN (B, N) = (1.-Exp(-ATN)) / (ATN+1.0E-6)
          !
          !mean light within canopy
               IAVEP = IATCNPY * SAVATN (B, N)!Iwc in Cerco and Moore 2001 (mean light within canopy)
          !
          !mean light on shoots
               IAVSH = IAVEP * EPATN (B, N)!Ish in Cerco and Moore 2001
          !
          ! COMPUTE LIGHT LIMITATIONS TO GROWTH
          !     ALPHSAV = 0.002 * SAVDPH (B, N) + 0.0016 
          !     ALPHSAV = Max (0.0,ALPHSAV)
		  !
          !     IK = PMSAV * GET_FT (T(B, KWC), TRPMSAV, KTPS1, KTPS2) / ALPHSAV
          !     FISH (B, N) = IAVSH / (Sqrt(IK*IK+IAVSH*IAVSH)+1.0E-30)
          !     IK = PMEPI * GET_FT (T(B, KWC), TRPMSAV, KTPE1, KTPE2) / ALPHEPI
          !     FIEP (B, N) = IAVEP / (Sqrt(IK*IK+IAVEP*IAVEP)+1.0E-30)
			   IK = PMSAV * GET_FT (AveT, TRPMSAV, KTPS1, KTPS2) / ALPHSAV
               FISH (B, N) = IAVSH / (Sqrt(IK*IK+IAVSH*IAVSH)+1.0E-30)
               IK = PMEPI * GET_FT (AveT, TRPMSAV, KTPE1, KTPE2) / ALPHEPI
               FIEP (B, N) = IAVEP / (Sqrt(IK*IK+IAVEP*IAVEP)+1.0E-30)
          !
            End Do
       !
         End Do
    !
    ! COMPUTE NUTRIENT LIMITATIONS TO GROWTH, SAV FIRST
    !
    !     KWC = KBM1 !water column layer (bottom layer)
	!	 KWC = 1 !water column layer (surface layer) 
		 KWL = 6
    !
         Do I = 1, NSAVCELL
            B = SAVCELL (I)
       !
		    AveNH4=0.0_SP
			AveNO3=0.0_SP
			AvePO4=0.0_SP
			AveSSI=0.0_SP
			DO K=1, KWL
			   AveNH4=AveNH4+NH4(B, K)
			   AveNO3=AveNO3+NO3(B, K)
			   AvePO4=AvePO4+PO4(B, K)
			   AveSSI=AveSSI+SSI(B, K)
			END DO
			AveNH4=AveNH4/KWL
			AveNO3=AveNO3/KWL
			AvePO4=AvePO4/KWL
			AveSSI=AveSSI/KWL




            NSED = Max (0., NH4T2TM1S_SHARE(B)/1000.)!NH4 in sediment anaerobic (2nd) layer  mgN/m^3 --> gN/m^3  ==  mgN/L

       !     NWAT = Max (0., NH4(B, KWC)+NO3(B, KWC))!NH4+NO3 in bottom layer of water column (gN/m^3)  (WAT means water)
			NWAT = Max (0., AveNH4+AveNO3)!NH4+NO3 in bottom layer of water column (gN/m^3)  (WAT means water)
            PSED = Max (0., PO4T2TM1S_SHARE(B)/1000.)!PO4 in pore water of sediment anaerobic layer (gP/m^3)

       !     DF = 1. / (1.+KADPO4*SSI(B, KWC))!fraction of dissolved PO4 in water column
			DF = 1. / (1.+KADPO4*AveSSI)!fraction of dissolved PO4 in water column
       ! where KADPO4 is PO4 partitioning coefficient m^3/gSolids, as SSI has unit gSolids/m^3
       !     PWAT = Max (0., DF*PO4(B, KWC))!Dissolved PO4 (gP/m^3), where PO4 is total concentration (gP/m^3)
	        PWAT = Max (0., DF*AvePO4)!Dissolved PO4 (gP/m^3), where PO4 is total concentration (gP/m^3)
       !
       !Nitrogen limitation on SAV growth
            NLSAV (B) = (KHNROOT*NWAT+KHNLEAF*NSED) / &
           & (KHNLEAF*KHNROOT+(KHNROOT*NWAT+KHNLEAF*NSED))
       !
       !fraction of nitrogen uptake by SAV due to inorganic nitrogen in sediment pore water
            FNSEDSAV (B) = (KHNLEAF*NSED) / &
           & (KHNROOT*NWAT+KHNLEAF*NSED+1.0E-6)
       !
       !preference of NH4 over NO3 on N uptake by SAV, controled by threshould KHNPSAV (gN/m^3)
       !
       !     NPRSAV (B) = NH4 (B, KWC) / (KHNPSAV+NO3(B, KWC)+1.E-6) * &
       !    & (NO3(B, KWC)/(KHNPSAV+NH4(B, &
       !    & KWC)+1.E-6)+KHNPSAV/(1.E-6+NH4(B, KWC)+NO3(B, KWC)))
		    NPRSAV (B) = AveNH4 / (KHNPSAV+AveNO3+1.E-6) * &
           & (AveNO3/(KHNPSAV+AveNH4 &
           & +1.E-6)+KHNPSAV/(1.E-6+AveNH4+AveNO3))
       !
       !Phosphorus limitation on SAV growth
            PLSAV (B) = (KHPROOT*PWAT+KHPLEAF*PSED) / &
           & (KHPLEAF*KHPROOT+KHPROOT*PWAT+KHPLEAF*PSED)
       !
       !fraction of PO4 uptake by SAV due to PO4 in sediment pore water
            FPSEDSAV (B) = (KHPLEAF*PSED) / &
           & (KHPROOT*PWAT+KHPLEAF*PSED+1.0E-6)
       !
       ! EPIPHYTES
       !
            NLEPI (B) = NWAT / (KHNEPI+NWAT+1.0E-6)
       !     NPREPI (B) = NH4 (B, KWC) / (KHNPEPI+NO3(B, KWC)+1.E-6) * &
       !    & (NO3(B, KWC)/(KHNPEPI+NH4(B, &
       !    & KWC)+1.E-6)+KHNPEPI/(1.E-6+NH4(B, KWC)+NO3(B, KWC)))!Same as eqn (6-17) of 2002 Chesapeake Bay Model report
            NPREPI (B) = AveNH4 / (KHNPEPI+AveNO3+1.E-6) * &
           & (AveNO3 /(KHNPEPI+AveNH4 &
           & +1.E-6)+KHNPEPI/(1.E-6+AveNH4+AveNO3))!Same as eqn (6-17) of 2002 Chesapeake Bay Model report
       !"The 2002 Chesapeake Bay Eutrophication Model",
       !Cerco C.F. and Noel M. R., July 2004, EPA 903-R-04-004
            PLEPI (B) = PWAT / (KHPEPI+PWAT+1.0E-6)
       !
       ! COMPUTE SULFIDE TOXICITY EFECT ON SAV, BASICALLY LIMITING GROWTH ON SAV
       !
            HSMOL = HST2TM1S_SHARE (B) / (1.+M2_SHARE*PIE2HS_SHARE) / &
           & 64000. !HS concetration in sediment pore water
       !where PIE2HS is Sulfide partitiong coefficient  (L/Kg)
       !M2 is solids concentration (kg/L)
       !HST2TM1S is total concentation of HS in sediments (mg/m^3)
       !==> HSMOL is in units of mole-S2/m^3
       !
            FHS (B) = Exp (-KHSTOX*HSMOL)!KHSTOX is of unit (m^3/mole-S2)
       !
       ! COMPUTE EPIPHYTE DENSITY LIMITATION
       !
            Do N = 1, NSAVSPC (B)
               DENLIM (B, N) = KHEP / (KHEP+EP(B, N)+1.0E-6)!density limitation of species N
            End Do
         End Do
    !
    !     KWC = KBM1
	!	 KWC = 1 !water column layer (surface layer) 
		 KWL = 6
    !
         Do I = 1, NSAVCELL
            B = SAVCELL (I)

		    AveSalt=0.0_SP
			  AveT=0.0_SP

			! DO K=1, KWL
			 !  AveSalt=AveSalt+SALT(B, K)
			 !  AveT=AveT+T(B, K)
		!	END DO
		!	AveSalt=AveSalt/KWL
		!	AveT=AveT/KWL

      AveT = SUM(T(B,1:6))/6.0
      AveSalt = SUM(SALT(B,1:6))/6.0

       ! COMPUTE RESPIRATION
       !
       !Photosynshetic respiration with toxicity to SAV enhanced by salinity
       !    STOXI = Max (SALT(B, KWC), 0.0)
		   STOXI = Max (AveSalt, 0.0)
            PRSPSAV (B) = PRSPSAVB + (1.-PRSPSAVB) * 0.5 * &
           & (1.+Tanh(STOXI-SALMAX))!Transition from PRSPSAVB (value < 1) to 1
       !using a tanh function with transition
       !centered at salinity = SALMAX
	   !     TFTBM = GET_FTBM (T(B, KWC), TRBMSAV, KTBMSAV)
		    TFTBM = GET_FTBM (AveT, TRBMSAV, KTBMSAV)
            BMLEAF (B) = BMSAV * TFTBM
            BMSTEM (B) = BMSAV * TFTBM
            BMROOT (B) = BMSAV * TFTBM
            BMTUBER (B) = BMTBRREF * TFTBM
            SLSH (B) = BMSAV * TFTBM * (FCLPSAV+FCRPSAV)
       !
       ! DETERMINE GROWTH FOR INDIVIDUAL DEPTH INCREMENTS
       !
            Do N = 1, NSAVSPC (B)!Loop over all species
          !
               GLIMSAV = Min (NLSAV(B), PLSAV(B), FISH(B, N), FHS(B))!global limitation due to N, P, light and HS toxicity
               GLIMEPI = Min (NLEPI(B), PLEPI(B), FIEP(B, N))
          !
          !leaf growth rate (1/day)
          !     PLEAF (B, N) = PMSAV * GET_FT (T(B, KWC), TRPMSAV, &
		  !   & KTPS1, KTPS2) * GLIMSAV / ACDWSAV !carbon to dry weight ratio gC/gDW
			    PLEAF (B, N) = PMSAV * GET_FT (AveT, TRPMSAV, &
			  & KTPS1, KTPS2) * GLIMSAV / ACDWSAV !carbon to dry weight ratio gC/gDW

          !BM_LEAF(B) = PMSAV * GET_FT (AveT, TRPMSAV, KTPS1, KTPS2)   ! added for DIC  interactions

          !==> PLEAF has unit 1/day
          !net production of SAV (gC/m^2/day)
               NPPSAV (B, N) = (PLEAF(B, N)-BMLEAF(B)) * LEAF (B, N)!!net production of SAV (gC/m^2/day)
          !
          ! EPIPHYTE GROWTH AND RESPIRATION
          ! PHOTORESPIRATION IS NOT CONSISTENT WITH SAV
          ! FOR NOW LEAVE IT ALONE
          !
          !epiphytes growth rate (1/day)
          !     PEP (B, N) = PMEPI * GET_FT (T(B, KWC), TRPMSAV, &
		  !	  & KTPE1, KTPE2) * GLIMEPI * DENLIM (B, N) / CCHLEPI !gC/gCHLA
          !!==> PEP has unit 1/day

          !     BMEP (B, N) = BMEPI * GET_FTBM (T(B, KWC), TRBMEPI, &
		  !	  & KTBME) + PRSPEPI * PEP (B, N)!basal metabolism and photochemical respiration (1/day)

          !     PREP (B, N) = PREPI * GET_FTPR (T(B, KWC), TRPREPI, &
		  !	  & KTPRE) * EP (B, N)!predation rate m^2-leaf/gC-Epi/day * (gC-Epi/m^2-leaf) ==> (1/day)
            !net production gC-Epi/gC-leaf/day (growth - respiration)
			    PEP (B, N) = PMEPI * GET_FT (AveT, TRPMSAV, &
			  & KTPE1, KTPE2) * GLIMEPI * DENLIM (B, N) / CCHLEPI !gC/gCHLA
          !==> PEP has unit 1/day
             !respiration rate 1/day
               BMEP (B, N) = BMEPI * GET_FTBM (AveT, TRBMEPI, &
			  & KTBME) + PRSPEPI * PEP (B, N)!basal metabolism and photochemical respiration (1/day)
             !predation 1/day
               PREP (B, N) = PREPI * GET_FTPR (AveT, TRPREPI, &
			  & KTPRE) * EP (B, N)!predation rate m^2-leaf/gC-Epi/day * (gC-Epi/m^2-leaf) ==> (1/day)
             !net production gC-Epi/gC-leaf/day (growth - respiration)

               NPPEPI (B, N) = (PEP(B, N)-BMEP(B, N)) * EP (B, N)!net epiphytes growth [gC-epiphyte/m^2-leaf)/day]
          !
            End Do
       !
         End Do
    !
    ! INTERACTIONS WITH WATER COLUMN
    !
    !
         Do I = 1, NSAVCELL
            B = SAVCELL (I)
            Do N = 1, NSAVSPC (B)!Loop over all species
          !
          !DISSOLVED OXYGEN AND CARBON
          !DO flux due to SAV (gO2/m^2/day), positive generating DO in the bottom layer (SAV cell)
          !
               DOSAVW (B) = DOSAVW (B) + SAVFRAC (B, N) * &
              & (AOCR*(1.-FDOSR)*PLEAF(B, N)*LEAF(B, &
              & N)-AOCR*((BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))*LEAF(B, &
              & N)+BMSTEM(B)*STEM(B, N))*FDOSAV)!Discount (reduce) by (1-FDOSAV) because
          !oxygen consumption of  STEM and LEAF are
          !routed to ROOT and Tuber ??
          !
               LDOCSAVW (B) = LDOCSAVW (B) + SAVFRAC (B, N) * ((LEAF(B, &
              & N)*(BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))+STEM(B, &
              & N)*BMSTEM(B))*FCLDSAV)!WLong: PLEAF has unit 1/day, PRSPSAV has no units
          !
               RDOCSAVW (B) = RDOCSAVW (B) + SAVFRAC (B, N) * ((LEAF(B, &
              & N)*(BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))+STEM(B, &
              & N)*BMSTEM(B))*FCRDSAV)
               LPOCSAVW (B) = LPOCSAVW (B) + SAVFRAC (B, N) * ((LEAF(B, &
              & N)*(BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))+STEM(B, &
              & N)*BMSTEM(B))*FCLPSAV)
               RPOCSAVW (B) = RPOCSAVW (B) + SAVFRAC (B, N) * ((LEAF(B, &
              & N)*(BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))+STEM(B, &
              & N)*BMSTEM(B))*FCRPSAV)
          !
          !flux of oxygen to water column due to epiphytes growth (photosynthesis) - respiration consumption of O2
          !gO2/m^2/day (positive into water column)
               DOEPIW (B) = DOEPIW (B) + SAVFRAC (B, N) * EP (B, N) * &
              & LEAF (B, N) * ALAC * AOCR * (PEP(B, N)-BMEP(B, &
              & N)*(1.-FCLDEPI-FCRDEPI))!             1
          !---------------------------
          !          gO2/m^2/day
          !
          !flux to LDOC by epi mortality (gC/m^2/day), positive to water column
               LDOCEPIW (B) = LDOCEPIW (B) + SAVFRAC (B, N) * (EP(B, &
              & N)*LEAF(B, N)*ALAC*BMEP(B, N)*FCLDEPI+EP(B, N)*LEAF(B, &
              & N)*ALAC*(PREP(B, N)+SLSH(B))*FCLDPEP)
          !
          !flux to RDOC by epi mortality (gC/m^2/day), positive to water column
               RDOCEPIW (B) = RDOCEPIW (B) + SAVFRAC (B, N) * (EP(B, &
              & N)*LEAF(B, N)*ALAC*BMEP(B, N)*FCRDEPI+EP(B, N)*LEAF(B, &
              & N)*ALAC*(PREP(B, N)+SLSH(B))*FCRDPEP)
          !flux to LPOC by epi mortality (gC/m^2/day), positive to water column
               LPOCEPIW (B) = LPOCEPIW (B) + SAVFRAC (B, N) * (EP(B, &
              & N)*LEAF(B, N)*ALAC*(PREP(B, N)+SLSH(B))*FCLPPEP)
          !
          !flux to RPOC by epi mortality (gC/m^2/day), positive to water column
               RPOCEPIW (B) = RPOCEPIW (B) + SAVFRAC (B, N) * (EP(B, &
              & N)*LEAF(B, N)*ALAC*(PREP(B, &
              & N)+SLSH(B))*(1.-FCLDPEP-FCRDPEP-FCLPPEP))!Note that 1 - FCLDPEP - FCRDPEP - FCLPPEP <==> FCRPPEP
          !
          ! NITROGEN AND PHOSPHORUS

               NH4SAVW (B) = NH4SAVW (B) + SAVFRAC (B, N) * &
              & (ANCSAV*FNISAV*((BMLEAF(B)+PLEAF(B, &
              & N)*PRSPSAV(B))*LEAF(B, N)+BMSTEM(B)*STEM(B, &
              & N))-(1.-FNSEDSAV(B))*NPRSAV(B)*PLEAF(B, &
              & N)*ANCSAV*LEAF(B, N))
          !
               NO3SAVW (B) = NO3SAVW (B) - SAVFRAC (B, N) * &
              & (1.-FNSEDSAV(B)) * (1.-NPRSAV(B)) * PLEAF (B, N) * &
              & ANCSAV * LEAF (B, N)
          !
          !
          !          I think it should also play a role in SEDNH4SAV and SEDNO3SAV
          !          It seems that we are missing the split to SEDNO3SAV
          !
          !
          !LDON flux (gN/m^2/day)
               LDONSAVW (B) = LDONSAVW (B) + SAVFRAC (B, N) * FNLDSAV * &
              & ANCSAV * ((BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))*LEAF(B, &
              & N)+BMSTEM(B)*STEM(B, N))
          !RDON flux (gN/m^2/day)
               RDONSAVW (B) = RDONSAVW (B) + SAVFRAC (B, N) * FNRDSAV * &
              & ANCSAV * ((BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))*LEAF(B, &
              & N)+BMSTEM(B)*STEM(B, N))
          !
          !LPON flux (gN/m^2/day)
               LPONSAVW (B) = LPONSAVW (B) + SAVFRAC (B, N) * FNLPSAV * &
              & ANCSAV * ((BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))*LEAF(B, &
              & N)+BMSTEM(B)*STEM(B, N))
          !
          !RPON flux (gN/m^2/day)
               RPONSAVW (B) = RPONSAVW (B) + SAVFRAC (B, N) * FNRPSAV * &
              & ANCSAV * ((BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))*LEAF(B, &
              & N)+BMSTEM(B)*STEM(B, N))
          !
          !NH4 flux (gN/m^2/day)
               NH4EPIW (B) = NH4EPIW (B) + SAVFRAC (B, N) * &
              & (FNIEPI*BMEP(B, N)+FNIPEP*(PREP(B, &
              & N)+SLSH(B))-NPREPI(B)*PEP(B, N)) * ANCEPI * EP (B, N) * &
              & LEAF (B, N) * ALAC
          !
          !NO3 flux (gN/m^2/day)
               NO3EPIW (B) = NO3EPIW (B) - SAVFRAC (B, N) * &
              & (1.-NPREPI(B)) * PEP (B, N) * ANCEPI * EP (B, N) * LEAF &
              & (B, N) * ALAC
          !
          !LDON flux due to Epiphyte (gN/m^2/day)
               LDONEPIW (B) = LDONEPIW (B) + SAVFRAC (B, N) * &
              & (FNLDEPI*BMEP(B, N)+FNLDPEP*(PREP(B, N)+SLSH(B))) * &
              & ANCEPI * EP (B, N) * LEAF (B, N) * ALAC
          !
          !RDON flux due to Epiphyte (gN/m^2/day)
               RDONEPIW (B) = RDONEPIW (B) + SAVFRAC (B, N) * &
              & (FNRDEPI*BMEP(B, N)+FNRDPEP*(PREP(B, N)+SLSH(B))) * &
              & ANCEPI * EP (B, N) * LEAF (B, N) * ALAC
          !
          !LPON flux (gN/m^2/day)
               LPONEPIW (B) = LPONEPIW (B) + SAVFRAC (B, N) * &
              & (FNLPEPI*BMEP(B, N)+FNLPPEP*(PREP(B, N)+SLSH(B))) * &
              & ANCEPI * EP (B, N) * LEAF (B, N) * ALAC
          !
          !RPON flux (gN/m^2/day)
               RPONEPIW (B) = RPONEPIW (B) + SAVFRAC (B, N) * &
              & (FNRPEPI*BMEP(B, N)+FNRPPEP*(PREP(B, N)+SLSH(B))) * &
              & ANCEPI * EP (B, N) * LEAF (B, N) * ALAC
          !
          !P fluxes to water column (gP/m^2/day)
          !
          !PO4 flux (gP/m^2/day)
               PO4SAVW (B) = PO4SAVW (B) + SAVFRAC (B, N) * &
              & (APCSAV*FPISAV*((BMLEAF(B)+PLEAF(B, &
              & N)*PRSPSAV(B))*LEAF(B, N)+BMSTEM(B)*STEM(B, &
              & N))-(1.-FPSEDSAV(B))*PLEAF(B, N)*APCSAV*LEAF(B, N))
          !
          !LDOP flux (gP/m^2/day)
               LDOPSAVW (B) = LDOPSAVW (B) + SAVFRAC (B, N) * FPLDSAV * &
              & APCSAV * ((BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))*LEAF(B, &
              & N)+BMSTEM(B)*STEM(B, N))
          !
          !RDOP flux (gP/m^2/day)
               RDOPSAVW (B) = RDOPSAVW (B) + SAVFRAC (B, N) * FPRDSAV * &
              & APCSAV * ((BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))*LEAF(B, &
              & N)+BMSTEM(B)*STEM(B, N))
          !
          !LPOP flux (gP/m^2/day)
               LPOPSAVW (B) = LPOPSAVW (B) + SAVFRAC (B, N) * FPLPSAV * &
              & APCSAV * ((BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))*LEAF(B, &
              & N)+BMSTEM(B)*STEM(B, N))
          !
          !RPOP flux (gP/m^2/day)
               RPOPSAVW (B) = RPOPSAVW (B) + SAVFRAC (B, N) * FPRPSAV * &
              & APCSAV * ((BMLEAF(B)+PLEAF(B, N)*PRSPSAV(B))*LEAF(B, &
              & N)+BMSTEM(B)*STEM(B, N))
          !
          !PO4 flux due to Epiphyte
               PO4EPIW (B) = PO4EPIW (B) + SAVFRAC (B, N) * &
              & (FPIEPI*BMEP(B, N)+FPIPEP*(PREP(B, N)+SLSH(B))-PEP(B, &
              & N)) * APCEPI * EP (B, N) * LEAF (B, N) * ALAC
          !
          !
               LDOPEPIW (B) = LDOPEPIW (B) + SAVFRAC (B, N) * &
              & (FPLDEPI*BMEP(B, N)+FPLDPEP*(PREP(B, N)+SLSH(B))) * &
              & APCEPI * EP (B, N) * LEAF (B, N) * ALAC
               RDOPEPIW (B) = RDOPEPIW (B) + SAVFRAC (B, N) * &
              & (FPRDEPI*BMEP(B, N)+FPRDPEP*(PREP(B, N)+SLSH(B))) * &
              & APCEPI * EP (B, N) * LEAF (B, N) * ALAC
               LPOPEPIW (B) = LPOPEPIW (B) + SAVFRAC (B, N) * &
              & (FPLPEPI*BMEP(B, N)+FPLPPEP*(PREP(B, N)+SLSH(B))) * &
              & APCEPI * EP (B, N) * LEAF (B, N) * ALAC
               RPOPEPIW (B) = RPOPEPIW (B) + SAVFRAC (B, N) * &
              & (FPRPEPI*BMEP(B, N)+FPRPPEP*(PREP(B, N)+SLSH(B))) * &
              & APCEPI * EP (B, N) * LEAF (B, N) * ALAC
          !
            End Do
         End Do
    !
    ! INTERACTIONS WITH SEDIMENTS
    !
         Do I = 1, NSAVCELL
            B = SAVCELL (I)
            Do N = 1, NSAVSPC (B)!loop over all species
          !
          ! DISSOLVED OXYGEN AND CARBON
          !
               SEDDOSAV (B) = SEDDOSAV (B) + SAVFRAC (B, N) * (LEAF(B, &
              & N)*AOCR*FDOSR*PLEAF(B, N)-(ROOT(B, &
              & N)*BMROOT(B)+TUBER(B, N)*BMTUBER(B))*AOCR*FDOSAV)
!
          !   SEDDOSAV(B) = SEDDOSAV(B)+SAVFRAC(B,N)*                    &    !Wen Long: why positive here? Should be negative
          !        (LEAF(B,N)*AOCR*FDOSR*PLEAF(B,N)              &    !Leaf mortality related DO consumption
          !        !(oxidation of C to CO2)
          !                !
          !   -( ROOT(B,N)*BMROOT(B)                        &    !Root respiration and tuber respiration
          !        +TUBER(B,N)*BMTUBER(B)                      &    !cost of oxygen (mortality loss of
          !        )*AOCR*FDOSAV                                &    !Root and Tuber carbon converted to oxygen
          !        )                                                  !consumption and oxydation of C to CO2 )
          !FDOSAV is the portion of carbon loss
          !assocatied with oxydation of Organic Carbon
          !to CO2 (TDIC)
          !
          !SEDDOSAV should be added to SOD directly and let it affect sediment
          !structure (such as first layer thickness)
          !
          !
               SEDPOCSAV (B) = SEDPOCSAV (B) + SAVFRAC (B, N) * &
              & (ROOT(B, N)*BMROOT(B)+TUBER(B, N)*BMTUBER(B)) * &
              & (1.-FDOSAV)!Root and Tuber mortality
          !
          ! NITROGEN AND PHOSPHORUS
          !
               SEDNH4SAV (B) = SEDNH4SAV (B) + SAVFRAC (B, N) * &
              & FNSEDSAV (B) * NPRSAV (B) * PLEAF (B, N) * ANCSAV * &
              & LEAF (B, N)
          !
               SEDNO3SAV (B) = SEDNO3SAV (B) + SAVFRAC (B, N) * &
              & FNSEDSAV (B) * (1.-NPRSAV(B)) * PLEAF (B, N) * ANCSAV * &
              & LEAF (B, N)
          !
               SEDPO4SAV (B) = SEDPO4SAV (B) + SAVFRAC (B, N) * &
              & FPSEDSAV (B) * PLEAF (B, N) * APCSAV * LEAF (B, N)
          !
               SEDPONSAV (B) = SEDPONSAV (B) + SAVFRAC (B, N) * &
              & (ROOT(B, N)*BMROOT(B)+TUBER(B, N)*BMTUBER(B)) * ANCSAV
          !
               SEDPOPSAV (B) = SEDPOPSAV (B) + SAVFRAC (B, N) * &
              & (ROOT(B, N)*BMROOT(B)+TUBER(B, N)*BMTUBER(B)) * APCSAV
          !
            End Do
         End Do
    !
    ! COMPUTE NEW SHOOT (leaf, tube, stem), ROOT, AND EPIPHYTE BIOMASS
    !
         Do I = 1, NSAVCELL
            B = SAVCELL (I)
       !
       ! DETERMINE DEGREE DAY INTERVAL.  SET MAXIMUM TO OCCURANCE OF SENESCENCE
       ! TO PREVENT ARRAY OVERFLOW
            J = Min (IDGD(B), NDDI-1)!entry to the loop up table in FPROOT, FPLEAF, etc for senescence of SAV (aging)
       !
            Do N = 1, NSAVSPC (B)!loop over all species
          !
               LFOLD = LEAF (B, N)!previous leaf biomass gC/m^2
          !
          !leaf biomass change rate(gC/m^2/day)
               DELLF = (PLEAF(B, &
              & N)*(1.-PRSPSAV(B))*FPLEAF(J)-BMLEAF(B)) * LFOLD + &
              & TRTBRLF (J) * TUBER (B, N)
!
          !   DELLF = (  PLEAF(B,N)                      & !leaf growth 1/day
          !        *(1.-PRSPSAV(B))                    & !discount mortality due to salt toxicity
          !        !when salt> SALTMAX, PRSPSAV --> 1 and there will be no growth
          !   *FPLEAF(J)                          & !
          !        -BMLEAF(B)                          & !discount basal metabolism
          !        )*LFOLD                              &
          !        + TRTBRLF(J)*TUBER(B,N)                 !transfer from tuber to leaf
          !
          !update root biomass (gC/m^2)
               ROOT (B, N) = ROOT (B, N) + (FPROOT(J)*PLEAF(B, &
              & N)*(1.-PRSPSAV(B))*LFOLD-BMROOT(B)*ROOT(B, N)) * DLTDY
          !
          !update stem biomass (gC/m^2)
               STEM (B, N) = STEM (B, N) + (FPSTEM(J)*PLEAF(B, &
              & N)*(1.-PRSPSAV(B))*LFOLD-BMSTEM(B)*STEM(B, N)) * DLTDY
          !
               TUBER (B, N) = TUBER (B, N) + (FPTUBER(J)*PLEAF(B, &
              & N)*(1.-PRSPSAV(B))*LFOLD-(BMTUBER(B)+&
              & TRTBRLF(J))*TUBER(B, N)) * DLTDY
          !
          !update leaf biomass (gC/m^2)
               LEAF (B, N) = LEAF (B, N) + DLTDY * DELLF !leaf mass increase (gC/m^2) in DLTDY (days) time step
          !Note LEAF is updated before EP!
          !
          !update epiphyte biomass (gC-Epiphyte)/(gm^2-Leaf)
               EP (B, N) = EP (B, N) * LFOLD * (1.0+DLTDY*(PEP(B, &
              & N)-BMEP(B, N)-PREP(B, N)-SLSH(B))) / (LEAF(B, &
              & N)+1.0E-20)!normalize by new LEAF to get EP
          !Note LFOLD is LEAF mass in previous time step
            End Do
         End Do
    !
    ! TEMPORARILY SET MINIMUM BIOMASS
    !
         Do I = 1, NSAVCELL
            B = SAVCELL (I)
            Do N = 1, NSAVSPC (B)!loop over all species
               ROOT (B, N) = Max (ROOT(B, N), 0.01)
               LEAF (B, N) = Max (LEAF(B, N), 0.01)
               STEM (B, N) = Max (STEM(B, N), 0.01)
               TUBER (B, N) = Max (TUBER(B, N), 0.01)
               EP (B, N) = Max (EP(B, N), 0.01)
            End Do
         End Do
    !
         Return
      End Subroutine SAV_COMP
  !
  !

!set the initial condition to sav if it is read from ICI file ! AN
      Subroutine SAV_INIT_ICI
    !
    !     Use MOD_WQM, Only: C2
	     Use MOD_LIMS, Only: MLOC
    !
         Implicit None
         Save
         Integer :: I, N
    !
    !if(MSR)write(*,*)'Init temp sedim = ',CTEMPI,'*******************************************'
         Do I = 1, MLOC

               DO N=1,NSAVM
                   LEAF(I,N)      = LEAFI(N)
				   STEM(I,N)      = STEMI(N)
				   ROOT(I,N)      = ROOTI(N)
				   TUBER(I,N)     = TUBERI(N)
				   EP(I,N)        = EPI(N)
               ENDDO

       !PIE1=PIENH4
       !PIE2=PIENH4
       !FD1=1/(1.0+m1*PIE1)
       !FD2=1/(1.0+m2*PIE2)
       !NH42I=NH4T2I*FD2
       !
       !NO32I=  NO3T2I
       !
       !IF(SALT(I,KBM1)<SALTND)THEN
       !        DPIE1PO4=DPIE1PO4F
       !ELSE
       !        DPIE1PO4=DPIE1PO4S
       !ENDIF
       !IF(O20>O2CRITPO4)THEN
       !        PIE1=PIE2PO4*DPIE1PO4
       !ELSE
       !        PIE1=PIE2PO4*DPIE1PO4**(O20/O2CRITPO4)
       !ENDIF
       !PIE2=PIE2PO4
       !FD1=1/(1.0+m1*PIE1)
       !FD2=1/(1.0+m2*PIE2)
       !PO42I=PO4T2I*FD2
       !
       !IF(O20>O2CRITSI)THEN
       !        PIE1=PIE2SI*DPIE1SI
       !ELSE
       !        PIE1=PIE2SI*DPIE1SI**(O20/O2CRITSI)
       !ENDIF
       !PIE2=PIE2SI
       !FD1=1/(1.0+m1*PIE1)
       !FD2=1/(1.0+m2*PIE2)
       !SI2I= SIT2I*FD2
       !
       !PIE1=PIE1HS
       !PIE2=PIE2HS
       !FD1=1/(1.0+m1*PIE1)
       !FD2=1/(1.0+m2*PIE2)
       !HS2I= HST2I *FD2
       !
         End Do
    !
    !
      End Subroutine SAV_INIT_ICI
	!
	!
      Function GET_FT (TEMPVAL, TMP, KTG1, KTG2)!m/d
       Real (SP) :: TEMPVAL, TMP, KTG1, KTG2 !temperature (degC)
       Real (SP) :: GET_FT
       If (TEMPVAL < TMP) Then !shape is different for < TMP1 and > TMP1
          GET_FT = Exp (-KTG1*(TEMPVAL-TMP)**2)!TMP1 is optimal temperature
       Else !KTG11 is the shape parameter (bell shape) for T < TMP1
          GET_FT = Exp (-KTG2*(TMP-TEMPVAL)**2)!KTG12 is for T > TMP1
       End If
       Return
      End Function GET_FT
!
      Function GET_FTBM (TEMPVAL, TR, KTB)!m/d
       Real (SP) :: TEMPVAL, TR, KTB !temperature (degC)
	   Real (SP) :: GET_FTBM
       GET_FTBM = Exp (KTB*(TEMPVAL-TR))
       Return
      End Function GET_FTBM
!
      Function GET_FTPR (TEMPVAL, TRPR, KTPR)!m/d
       Real (SP) :: TEMPVAL, TRPR, KTPR !temperature (degC)
       Real (SP) :: GET_FTPR
       GET_FTPR = Exp (KTPR*(TEMPVAL-TRPR))
       Return
      End Function GET_FTPR
!
End Module MOD_SAV
