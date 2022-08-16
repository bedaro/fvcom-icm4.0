!mod_sed.F
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
!************************************************************************
!**                                                                    **
!**                                                                    **
!************************************************************************
!**                                                                    **
!** Inputs                                                             **
!**                                                                    **
!**            Required inputs for sediment sub-model                  **
!**                                                                    **
!** A.  Passed to sediment subroutine from water quality subroutine    **
!**                                                                    **
!**   3.  Overlying water column depth                                 **
!**   4.  Overlying water column segment temperature and salinity      **
!**   5.  Overlying water column ammonia, nitrate, phosphate, silica   **
!**       and dissolved oxygen concentrations                          **
!**                                                                    **
!** B. Inputs supplied via direct input to the sediment subroutine     **
!**                                                                    **
!**  Variable names        Description                         Units   **
!**                                                                    **
!**                                                                    **
!**     HSEDALL      Depth of sediment layer  (h2)               cm    **
!**      DIFFT       Water column-sediment layer diffusion             **
!**                  coefficient for temperature                       **
!**                  cm^2/sec on input, model converts it to m^/sec    **
!**     SALTSW       Salinity concentration for determining      ppt   **
!**                  whether methane or sulfide SOD formulation        **
!**                  is to be used.  Also determines PO4 sorption.     **
!**     SALTND       Determines whether fresh or saltwater             **
!**                  nitrification/denitrification rates are used      **
!**                                                                    **
!**   Diagenesis stoichiometry                                         **
!**                                                                    **
!**                 Fractions of G1, G2, and G3 contained in ...       **
!**                                                                    **
!**   FRPALG1(3)     Algal group no.1 phosphorus                       **
!**   FRPALG2(3)     Algal group no.2 phosphorus                       **
!**   FRPALG3(3)     Algal group no.3 phosphorus                       **
!**   FRPOP(MTLOC,3) Non-algal particulate organic phosphorus          **
!**   FRNALG1(3)     Algal group no.1 nitrogen                         **
!**   FRNALG2(3)     Algal group no.2 nitrogen                         **
!**   FRNALG3(3)     Algal group no.3 nitrogen                         **
!**   FRPON(MTLOC,3) Non-algal particulate organic nitrogen            **
!**   FRCALG1(3)     Algal group no.1 carbon                           **
!**   FRCALG2(3)     Algal group no.2 carbon                           **
!**   FRCALG3(3)     Algal group no.3 carbon                           **
!**   FRPOC(MTLOC,3) Non-algal particulate organic carbon              **
!**                                                                    **
!**   Diagenesis kinetics                                              **
!**                                                                    **
!**   KPDIAG(3)     Reaction rates for POP G1, G2, and G3        /day  **
!**   DPTHTA(3)     Temperature thetas for POP G1, G2, and G3          **
!**   KNDIAG(3)     Reaction rates for PON G1, G2, and G3        /day  **
!**   DNTHTA(3)     Temperature thetas for PON G1, G2, and G3          **
!**   KCDIAG(3)     Reaction rates for POC G1, G2, and G3        /day  **
!**   DCTHTA(3)     Temperature thetas for POC G1, G2, and G3          **
!**   KSI           Reaction rate for Part. Biogenic Si (PSISED) /day  **
!**   THTASI        Temperature theta for PSISED                       **
!**   THTASISAT     Temperature theta for Si saturation concentration  **
!**                                                                    **
!**   Solids and transport                                             **
!**                                                                    **
!**   VPMIX(MTLOC)   Particulate diffusion coefficient (Dpp)   m**2/day**
!**   THTADP        Temperature theta for Dpp                          **
!**   VDMIX(MTLOC)   Porewater diffusion coefficient (Ddp)     m**2/day**
!**   THTADD        Temperature theta for Ddp                          **
!**      M1         Concentration of solids in layer 1       kg/l      **
!**      M2         Concentration of solids in layer 2       kg/l      **
!**                                                                    **
!**   Reaction kinetics                                                **
!**                                                                    **
!**    KAPPNH4F     Nitrification reaction velocity                    **
!**                 for freshwater in layer 1                m/day     **
!**    KAPPNH4S     Nitrification reaction velocity                    **
!**                 for saltwater  in layer 1                m/day     **
!**    PIENH4       Ammonia partition coefficient            L/kg      **
!**    THTANH4      Theta for nitrification reaction velicities        **
!**    KMNH4        Nitrification half saturation constant             **
!**                 for ammonia                              mg N/m3   **
!**    KMNH4O2      Nitrification half saturation constant             **
!**                 for oxygen                               mg O2/L   **
!**    KAPPNO3F     Denitrification reaction velocity                  **
!**                 for freshwater in layer 1                m/day     **
!**    KAPPNO3S     Denitrification reaction velocity                  **
!**                 for saltwater  in layer 1                m/day     **
!**    K2NO3        Denitrification reaction velocity                  **
!**                 in layer 2                               m/day     **
!**    THTANO3      Theta for denitrification                          **
!**    KAPP1HSD       Reaction velocity for dissolved sulfide          **
!**                 oxidation in layer 1                     m/day     **
!**    KAPP1HSP       Reaction velocity for particulate sulfide        **
!**                 oxidation in layer 1                     m/day     **
!**    PIE1HS        Partition coefficient for sulfide                 **
!**                 in layer 1                               L/kg      **
!**    PIE2HS        Partition coefficient for sulfide                 **
!**                 in layer 2                               L/kg      **
!**    THTAH2S      Theta for both dissolved and particulate           **
!**                 sulfide oxidation                                  **
!**    KMHSO2       Sulfide oxidation normalization constant           **
!**                 for oxygen                               mg O2/L   **
!**    CSISAT       Saturation concentration for pore water            **
!**                 silica                                     mg Si/m3**
!**    DPIE1SI      Incremental partition coefficient for              **
!**                 silica in layer 1                        L/kg      **
!**    PIE2SI       Partition coefficient for silica in                **
!**                 layer 2                                  L/kg      **
!**    KMPSI        Particulate biogenic silica half saturation        **
!**                 constant for dissolution                 mg Si/m3  **
!**    O2CRITSI     Critical dissolved oxygen concentration            **
!**                 for layer 1 incremental silica sorption  mg O2/L   **
!**    JSIDETR      Detrital biogenic silica source to                 **
!**                 sediment                                 mg Si/m2-d**
!**    DPIE1PO4F    Incremental partition coefficient for              **
!**                 phosphate in layer 1 (freshwater)        unitless  **
!**    DPIE1PO4S    Incremental partition coefficient for              **
!**                 phosphate in layer 1 (saltwater)         L/kg      **
!**    PIE2PO4      Partition coefficient for phosphate                **
!**                 in layer 2                               L/kg      **
!**    O2CRITPO4       Critical dissolved oxygen concentration for     **
!**                 layer 1 incremental phosphate sorption   mg O2/L   **
!**    KMO2DP       Particle mixing half saturation constant           **
!**                 for oxygen                               mg O2/L   **
!**    TEMPBEN      Temperature at which benthic stress                **
!**                 accumulation is reset to zero            deg C     **
!**    KBENSTR      Decay constant for benthic stress        /day      **
!**    KLBNTH       Ratio of bio-irrigation to bioturbation            **
!**    DPMIN        Minimum particle diffusion coefficient   m2/day    **
!**    KAPPCH4      methane oxidation reaction velocity      m/day     **
!**    THTACH4      theta for methane oxidation                        **
!**                                                                    **
!** Output                                                             **
!**                                                                    **
!**    The subroutine returns fluxes for                               **
!**                                                                    **
!**     JSOD, JAQSOD, JCH4AQ and JCH4G  [g o2*/m2-day]                 **
!**     JNH4, JPO4, JNO3 and JSI  [mg/m2-day]                          **
!**                                                                    **
!************************************************************************
!
Module MOD_SED
  !
      Use MOD_PREC, Only: SP
  !
      Use MOD_SIZES, Only: MGL, NFLP
  !
  Use MOD_WQM, Only: BFOFN, WCLPTR, NTVWCLF, DLT, TMSTRT, SETTLING, &
       & B, B1, B2, B3, LPOP, Q1, Q2, Q3, SAV_CALC, WS1NET, WS2NET, &
       & WS3NET, WSLNET, WSRNET, WSSNET, WSUNET, ANC1, ANC2, ANC3, ASC1, &
       & ASC2, ASC3, CTEMP, LPOC, LPON, PCFWS, PNFWS, PPFWS, PSFWS, &
       & SIUPB, RPOC, RPON, RPOP, SSFWS, SSI, BENCH4A, BENCH4G, BENCOD, &
       & BENDEN, BENDO, BENDOC, BENNH4, BENNO3, BENPO4, BENSA, COD, DOXG, &
#if defined(CSOD_CALC)
  & CSODflux,       &!: keep CSOD for later use in TOTINC (mod_DIC_TALK.F)
#endif
  & SIAT, JDAY, KADPO4, KADSA, NH4, NO3, PIP, PO4, SALT, T, SEDTEMP_CALC
  !
      Use MOD_FILEINFO, Only: BFI, WCL, BFO
  !
      Use MOD_LIMS, Only: MLOC, MTLOC, KBM1
  !
      Use MOD_HYDROVARS, Only: ART1, DZ, D
  !! took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: SERIAL, MSR, PAR
      Use MOD_SF, Only: SFEEDER
      Use MOD_SED_SF_EXCHANGE_VARS, Only: JLPOC_SED_SF, JLPOC_SED_SF, &
     & JLPON_SED_SF, JLPOP_SED_SF, JRPOC_SED_SF, JRPON_SED_SF, &
     & JRPOP_SED_SF, JNH4_SED_SF, JPO4_SED_SF, SOD_SED_SF, JSA_SED_SF, &
     & JSU_SED_SF, SSI_SED_SF, SU_SED_SF, SA_SED_SF, PIP_SED_SF
  !
  !
      Use MOD_BA, Only: BALGAE_CALC, BA_CALC, BANH4, BANO3, BAPO4, &
     & BADOC, BADO, DTPOC1_BA, DTPOC2_BA, DTPOC3_BA, DTPON1_BA, &
     & DTPON2_BA, DTPON3_BA, DTPOP1_BA, DTPOP2_BA, DTPOP3_BA, &
     & DTAG3CFL_BA, DTAG3NFL_BA, DTAG3PFL_BA
  !
      Use MOD_DF, Only: DFEEDER, DF_GROW_POC1, DF_GROW_POC2, DF_MORT, &
     & DF_RESP, DF_PRED, DF_SOD, DFEEDM1S, AMCN, AMCP, DF_CALC
      Implicit None
      Save
  !
  !
  !variables that are used for as flags of model simulation
  !
  !
  !Moved this from wqm_modules.F
  !
      Logical :: BENTHIC_OUTPUT, STEADY_STATE_SED_G3, &
     & STEADY_STATE_SED_IC, STEADY_STATE_SED_TS, STEADY_STATE_SED_AVG, &
     & AGGREGATE_POM_FLUX ! If true, aggregate all POM reaching the sediments before distributing into 3 sediment labilities
  ! BENTHIC_OUTPUT,              & !flag for benthic model output
  ! STEADY_STATE_SED_G3,        & !flag for steady state simulation for G3
  ! STEADY_STATE_SED_IC,        & !flag for setting initial condition of sediments using steady state solution
  !                               !     with first time step forcing from water column to sediments and steady state formulation
  ! STEADY_STATE_SED_TS,        & ! quasi-steady state time series calculation of sediment flux
  ! STEADY_STATE_SED_AVG,       & !flag for calculating steady state using time averaged forcing
  !                               !    from water colum to sediments and using steady state formulation
  ! AGGREGATE_POM_FLUX            ! If true, aggregate all POM reaching the sediments before distributing into 3 sediment labilities
  !
  !                 BALGAE_CALC            !flag for activating benthic algae
  !                                        !calculation
  !
  ! variables that are local to sediment module
  !
      Real (SP) :: DLTS !sediment simulation time step (day)
  !==DLT/86400, where DLT is in seconds
!
      Real (SP), Allocatable, Dimension (:) :: HSED !sediment layer thickness (H2 (m))
      Real (SP), Allocatable, Dimension (:) :: HSED1 !sediment layer thickness (H1 (m))
  !
  ! input variables that control the sediment modules simulation
  !
      Integer :: NDTSED ! number of WC timesteps between sedim calculations (DTsed = NDTSED * DTwater)
      Integer :: SSTATEG3 !==1 for steady state calculation of G3 POC PON POP in sediments
  !    based on time average of the simulation at the end of the integration
  !
  !
      Integer :: SSTATEIC !==1 for setting initial condition of sediment using steady state solution
  !   based on overlying water column concentration and fluxes of the first time step
  !
      Integer :: QUASISS !==1 for quasi-steady state calculation, meaning that even the forcing flux (JPOC, JPON, JPOP, JPOS)
  !   and overlying water column is time varying
  !   the calculation of sediment module is using
  !   steady state.
  !
  !
  !
      Integer :: SSTATEAVG !==1 for calculating sediment model using time averaged overlying water column
  !   forcing and time averaged POM fluxes to sediments, and steady state formulation
  !   The results of this can be used for subsequent simulations as initial conditions
  !   when SSTATEIC /= 1
  !
      Integer :: SWITCH_LAB !==1 to aggregate all POM reaching the sediments before distributing into 3 sediment labilities
  !  if 0, then LPOMflux goes to G1 in the sediments and RPOM is divided in 65% to G2 and 35% to G3
  !
  !
      Integer :: ITVWCLF !==1 for forcing with time varying overlying water column condition input file
      Real (SP) :: WCLFDAYS !length of water column forcing for each WCLFN file (days)
  !
      Real (SP) :: HSEDALL !Depth of sediment layer (h1+h2 ~= h2) (cm)
      Real (SP) :: DIFFT !water column-sediment diffusion rate (m^2/sec)
  !
      Real (SP) :: SALTSW !Salinity concentration for determining whether methane or sulfide SOD formulation (ppt)
  !if salt water  (> SALTSW) use sulfide,
  !if freshwater (<=SALTSW) use methane
  !
      Real (SP) :: SALTND !Salinity for determing whether to use freshwater (<=SALTND) or salt water (>SALTND)
  !nitrification/denitrification (ppt)
  !
      Real (SP), Dimension (3) :: FRPALG1, FRPALG2, FRPALG3, FRNALG1, &
     & FRNALG2, FRNALG3, FRCALG1, FRCALG2, FRCALG3
  !FRPALG1,    & !Algal group 1 Phosphorus fraction for
  !                                         !G1,G2,G3 in phytoplankton
  !                         FRPALG2,    & !Algal group 2 Phosphorus fraction for.
  !                                         !G1,G2,G3 in phytoplankton
  !                         FRPALG3,    & !Algal group 3 Phosphorus fraction for.
  !                                         !G1,G2,G3 in phytoplankton
  !                         FRNALG1,    & !Algal group 1 Nitrogen fractions for
  !                                         !G1,G2,G3 in phytoplankton
  !                         FRNALG2,    & !Algal group 2 Nitrogen fractions for.
  !                                         !G1,G2,G3 in phytoplankton
  !                         FRNALG3,    & !Algal group 3 Nitrogen fractions for G1
  !                                         !G2 G3 in phytoplankton
  !                         FRCALG1,    & !Algal group 1 Carbon fractions for G1,
  !                                         !G2, G3 in phytoplankton
  !                         FRCALG2,    & !Algal group 2 Carbon fractions for G1,G2
  !                                         !G3 in phytoplankton
  !                         FRCALG3!,   & !Algal group 3 Carbon fractions for G1,
  !                                         !G2, G3 in phytoplankton
  !                         FRCPHB,    & !Fraction of POC generated for G1,G2,G3
  !                                         !due to predation on benthic algae
  !                         FRNPHB,    & !Fraction of PON genereated for G1,G2,G3
  !                                         !due to predation on benthic algae
  !                         FRPPHB       !Fraction of POP generated for G1,G2,G3
  !                                         !due to predation on benthic algae
  !
      Real (SP), Allocatable :: FRPOP (:, :), FRPON (:, :), FRPOC (:, &
     & :)!MTLOC x 3, Fraction of non-algal particulate organic carbon (dimensionless)
!
      Real (SP) :: FRPON1, FRPOP1, FRPOC1 ! Fraction of G1 in POC (detritus) flux from water column
  !
  !
      Real (SP) :: FRACL1, FRACL2, FRACL3 ! Fraction of total POM flux from water column sent to G2 in sediments
  !
  !Moved here from wqm_modules.F
  !sediment reaction kinetics.
  !
      Real (SP) :: KPDIAG (3), DPTHTA (3), KNDIAG (3), DNTHTA (3), &
     & KCDIAG (3), DCTHTA (3)!Temperature control theta for C
  !
      Real (SP) :: KSI, THTASI, THTASISAT !temperature control theta for Si saturation concentration
  !
  !
      Real (SP) :: M1, M2, THTADP, THTADD !temperature tehat for porewater diffusion coefficient
  !      M1,       & !layer 1 solids concentration  (kg/L)
  !      M2,       & !layer 2 solids concentration  (kg/L)
  !      THTADP,   & !temperature theta for particle diffusion coefficient
  !                  !(Dpp)
  !      THTADD      !temperature tehat for porewater diffusion coefficient
  !(Ddp)
  !
  !Ammonia nitrification kinetics.
      Real (SP) :: KAPPNH4F, KAPPNH4S, PIENH4, THTANH4, KMNH4, &
     & THTAKMNH4, KMNH4O2
  !     KAPPNH4F, & !Nitrification reaction velocity (m/day) (integrated over layer 1 thickness)
  !                 !for freshwater in layer 1
  !     KAPPNH4S, & !Nitrification reaction velocity (m/day) for salt water
  !                 !in layer 1 (integrated over layer 1 thickness)
  !     PIENH4,   & !Ammonia partition coefficient (L/kg)
  !     THTANH4,  & !temperature control theta for Theta for nitrification
  !                 !reaction velicities.
  !     KMNH4,    & !Nitrification half saturation constant for ammonia in layer 1
  !                 ! mg N/m^3 (for total NH4 in layer 1) (DMD book (3.30)
  !     THTAKMNH4,& !temperature contro theta for nitrification half saturation constant KMNH4
  !     KMNH4O2     !Nitrification half saturation constant for oxygen
  ! mgO2/L
  !
  !
  !NO3 denitrification kinetics
      Real (SP) :: KAPPNO3F, KAPPNO3S, K2NO3, THTANO3
  ! KAPPNO3F,  &  !Denitrification reaction velocity (m/day) for fresh
  !                   !water in layer 1 (integrated over layer 1 thickness)
  !     KAPPNO3S,  &  !Denitrification reaction velocity (m/day) for salt
  !                   !water in layer 1 (integrated over layer 1 thickness)
  !     K2NO3,     &  !Denitrification reaction velocity in layer 2 (m/day) (DMD p 105 table 4.1)
  !                   !(integrated over the thickness of layer 2)
  !     THTANO3       !temperature control theta for denitrification
  !
  !HS sulfide oxidation kinetics
      Real (SP) :: KAPP1HSD, KAPP1HSP, PIE1HS, PIE2HS, THTAH2S, KMHSO2
  ! KAPP1HSD,     & !Reaction velocity for dissolved sulfide oxidation in layer 1
  !                   !(m/day)
  !     KAPP1HSP,     & !Reaction velocity for particulate sulfide oxidation
  !                   !in layer 1 (m/day)
  !     PIE1HS,      & !Partition coefficient for sulfide in layer 1 (L/kg)
  !     PIE2HS,      & !Partition coefficient for sulfide in layer 2 (L/kg)
  !     THTAH2S,    & !temperature control theta for both dissolved and
  !                   !particulate sulfide
  !     KMHSO2        !Sulfide oxidation normalization constant
  !for oxygen (mg O2/L)
  !Silica kinetics
      Real (SP) :: CSISAT, CSISATT20, DPIE1SI, PIE2SI, KMPSI !Particulate biogenic silica half saturation for silica.
  !for dissolution (mgSi/m^3)
  !
      Real (SP) :: O2CRITSI, JSIDETR !Detrital biogenic silica source to
  ! O2CRITSI, & !Critical dissolved oxygen concentration for layer 1
  !                  ! incremental silica sorption  mg O2/L
  !         JSIDETR     !Detrital biogenic silica source to
  ! sediment mg Si/m2-d
!
  !PO4 sorption kinetics
      Real (SP) :: DPIE1PO4F, DPIE1PO4S, PIE2PO4, O2CRITPO4, KMO2DP !Particle mixing half saturation constant
  ! DPIE1PO4F,  &    !Incremental partition coefficient for
  !                      !phosphate in layer 1 (freshwater) (unitless)
  !      DPIE1PO4S, &    !Incremental partition coefficient for
  !                      !phosphate in layer 1 (saltwater)  (unitless)
  !      PIE2PO4,   &    !Partition coefficient for phosphate
  !                      !in layer 2                         (L/kg)
  !      O2CRITPO4,    &    !Critical dissolved oxygen concentration for
  !                      !layer 1 incremental phosphate sorption  (mg O2/L)
  !      KMO2DP          !Particle mixing half saturation constant
  !for oxygen (mgO2/L) (bioturbation)
  !
  !benthic stress and bioturbation etc mixing control
      Real (SP) :: TEMPBEN, KBENSTR, KLBNTH, DPMIN !Minimum particle diffusion coefficient (m2/day)
  !  TEMPBEN, &      !temperature at which benthic stress is reset to
  !                      !zero (degC)
  !      KBENSTR, &      !Decay constant for benthic stress (1/day)
  !      KLBNTH,  &      !Ratio of bio-irrigation to bioturbation
  !      DPMIN           !Minimum particle diffusion coefficient (m2/day)
  !
  !methane kinetics
      Real (SP) :: KAPPCH4, THTACH4, KMCH4O2, KMSO4 !saturation concentration for SO4 in layer 2 for methane oxidation by sulfate (mgO2/L)
  !
!
  !Moved from subroutine SED_READ to here
      Character (Len=8) :: SPVARS, PRINTS
!
  !Checked the units of these settling velocities
  !net settling rate is supposed to be smaller than normal settling rate
  !net settling rate = (normal settling - suspension rate)
  !so the suspension rate (m/d) = (normal settling - net settling ) (m/d)
      Real (SP), Allocatable, Dimension (:) :: WSSBNET, WSLBNET, &
     & WSRBNET, WS1BNET, WS2BNET, WS3BNET, WSUBNET !net settling rate of particulate biogenic silica (m/d)
  !
  !Moved from subroutine SED_READ
      Character (Len=8) :: SPVARB, PRINTB
  !
      Real (SP), Allocatable, Dimension (:) :: VSED, VPMIX, VDMIX !m^2/day, VSED in m/yr on input, model converts to m/day
  !
  !
  !Moved from subroutine SED_READ
      Character (Len=8) :: SPVARLR, PRINTLR
  !
  !other local variables that are used in the sediment module
  !
  !
  !Moved the following from wqm_modules.F
  !Renamed DP to DPP and DD to DDP here
      Real (SP) :: DPP, W2, DDP, H2
  ! DPP,       &   !Particulate diffusion coefficient (m2/day)
  !                           !                                   <==VPMIX
  !        W2,            &   !Burial velocity (m/day)            <==VSED
  !        DDP,           &   ! porewater diffusion rate (m2/day) <==VDMIX
  !        H2                 !sediment layer thickness (m).
  !
  !
      Real (SP) :: DPIE1PO4 !sediment 1st layer incremental partition coefficient for phosphate  (dimensionless)
  !
      Real (SP) :: KAPPNH4, KAPPNO3
  !KAPPNH4,   &   !sediment layer 1 nitrification velocity (m/day)
  !                        !not used, replaced by KAPPNH4F and KAPPNH4S
  !KAPPNO3        !sediment layer 1 denitrification velocity
  !(m/day), never used, replaced by KAPPNO3F,
  !KAPPNO3F (integrated over layer 1 thickness)
  !Commented KAPP1CH4, it is duplicated with KAPPCH4
  ! KAPP1CH4        !sediment layer 1 CH4 oxidation velocity (m/day)
!
      Integer :: KWC, IERR, ITEMP
  !
  !Moved from wqm_modules.F
      Character (Len=20) :: SSNAME (17)
  !
  !Mmoved from wqm_modules.F
      Real (SP) :: K0H1D, K0H1P, K1H1D, K1H1P, K2H2D, K2H2P, K3, J1, &
     & J2, KMC1, KL12, KL12NOM, KHS_1, KL12SO4 !(used for NH4, i.e. mgN/m^3, for other constituents, KMC1 is zero)
  !    K0H1D,   &        !dissolved constituent reaction velocity of first  layer (m/d) 
  !             K0H1P,   &        !particulate constituent reaction velocity of first layer (m/d) (always zero) 
  !             K1H1D,   &        !dissolved constituent reaction velocity in first layer (m/d) 
  !             K1H1P,   &        !particulate constituent reaction velocity of first layer (m/d) 
  !             K2H2D,   &
  !             K2H2P,   &
  !             K3,      &    !dissolution rate of particulate biogenic silicate  PSISED (m/d)
  !             J1,      &  !**** need to be careful with J1 here, better rename to JJ1, J1 is used else where
  !             J2,      &  !**** need to be careful with J2 here, better rename to JJ2, J2 is used else where
  !             KMC1,    &    !half reaction constant for total concentration in layer 1 (CT1)
  !                         !(used for NH4, i.e. mgN/m^3, for other constituents, KMC1 is zero)
  !             KL12,    &
  !             KL12NOM, &    !KL12 after temperature regulation (m/day)
  !             KHS_1,   &
  !             KL12SO4
  !
  !sediment POM concentration (all two layers ~= layer 2)
      Real (SP) :: PON1, PON1TM1, PON2, PON2TM1, PON3, PON3TM1, POC1, &
     & POC1TM1, POC2, POC2TM1, POC3, POC3TM1, POP1, POP1TM1, POP2, &
     & POP2TM1, POP3, POP3TM1 !sediment POP G3 (mgP/m^3) at previous time step
  !    PON1,          & !sediment PON G1 (mgN/m^3)    at new time step
  !             PON1TM1,       & !sediment PON G1 (mgN/m^3) at previous time step
  !
  !             PON2,          & !sediment PON G2 at new time step
  !             PON2TM1,       & !sediment PON G2 (mgN/m^3) at previous time step
  !
  !             PON3,          & !sediment PON G3 (mgN/m^3) at new time step
  !             PON3TM1,       & !sediment PON G3 (mgN/m^3) at previous time step
  !
  !             POC1,          & !sediment POC G1 (mgC/m^3) at new time step
  !             POC1TM1,       & !sediment POC G1 (mgC/m^3) at previous time step
  !
  !             POC2,          & !sediment POC G2 (mgC/m^3) at new time step
  !             POC2TM1,       & !sediment POC G2 (mgC/m^3) at previous time step
  !
  !             POC3,          & !sediment POC G3 (mgC/m^3) at new time step
  !             POC3TM1,       & !sediment POC G3 (mgC/m^3) at previous time step
  !
  !             POP1,          & !sediment POP G1 (mgP/m^3) at new time step
  !             POP1TM1,       & !sediment POP G1 (mgP/m^3) at previous time step
  !
  !             POP2,          & !sediment POP G2 (mgP/m^3) at new time step
  !             POP2TM1,       & !sediment POP G2 (mgP/m^3) at previous time step
  !
  !             POP3,          & !sediment POP G3 (mgP/m^3) at new time step
  !             POP3TM1             !sediment POP G3 (mgP/m^3) at previous time step
  !
      Real (SP) :: PSISED, PSISEDTM1 !sediment particulate biogenic Si (mgSi/m^3)
  !
      Real (SP) :: JNX, JCX, JPX !depth integrated POP diagenesis rate (mgP/m^2/day)
  !
      Real (SP) :: XJCNO3, XJCNO31, XJCNO32, XJCO2, XJC1 ! (gO2/m^2/day) cabon diagenesis flux converted to oxygen equivalents.
  !discounting organic matter used for denirification in layer
  !1 and layer 2
  !Moved from wqm_modules.F
  !
      Real (SP) :: NH40, NH41, NH41TM1, NH4T1, NH42, NH4T2, NH4T2TM1, &
     & JNH4 !mgN/m^2/day
  !
      Real (SP) :: NO30, NO31, NO31TM1, NO3T1, NO32, NO3T2, NO3T2TM1, &
     & JNO3 !mgN/m^2/day
  !
      Real (SP) :: PO40, PO41, PO42, PO4T1, PO4T2, PO41TM1, PO4T2TM1, &
     & JPO4 !(mgP/m^2/day)
  !
  !Moved these here from wqm_modules.F
  !
      Real (SP) :: SI0, SI1, SI2, SIT1, SIT2, SI1TM1, SIT2TM1, JSI !(mgSi/m^2/day)
  !
      Real (SP) :: HS0, HS1, HS2, HST1, HST2, HS1TM1, HST2TM1, JHS !(gO2/m^2/day)
  !
      Real (SP) :: CH40, CH41, CH41TM1, CH4T1, CH42, CH42AV, CH4T2, &
     & CH4T2AV, CH4T2TM1, JCH4G, JCH4AQ, JCH4 !(gO2/m^2/day)
  !
      Real (SP) :: SO40, SO40MG, SO41, SO41TM1, SO4T1, SO42, SO42AV, &
     & SO4T2, SO4T2AV, SO4T2TM1, JSO4, SO4T2SAV !never used
  !
      Real (SP) :: JO2NH4, CH4T2SAV, JCH4GASO2, JGAS !never used
  !
      Real (SP) :: KPOP1, KPOP2, KPOP3, KPON1, KPON2, KPON3, KPOC1, &
     & KPOC2, KPOC3 !G3 POC diagenesis rate (1/day)
  !
      Real (SP) :: ISEDMN, ISEDMP, ISEDMC !total amount of carbon in sediments (POC) (previous time step)
  !
      Real (SP) :: SEDMN, SEDMP, SEDMC !total sediment carbon  (kgC) (new time step)
  !
      Real (SP) :: WSSINETMMD, WLPOMNETMMD, WRPOMNETMMD, WB1NETMMD, &
     & WB2NETMMD, WB3NETMMD, WPSINETMMD !Particulate biogenic silica settling rate (mm/day)
  !
      Real (SP) :: STP20, DF, PO4AVL, BFOR, W12MIN !minimum mixing speed of layer 1 and 2 (cm/day?)
  !
      Real (SP) :: LOGICT
  !
      Real (SP) :: XKPOC1, XKPOC2, XKPOC3, XKPOP1, XKPOP2, XKPOP3, &
     & XKPON1, XKPON2, XKPON3 !diagenesis reaction velocity (layer integrated ~ kpon3*H2) (m/d) of PON3 with temperature control
  !
      Real (SP) :: TEMP_SED
  !
      Real (SP) :: SFEED_TMP_1, SFEED_TMP_2, SFEED_TMP_3 !3 types of suspension feeders
  !
      Real (SP) :: DLTS_H2 !sedimen time step (day) divided by second sediment layer thickness (m), final unit: 1/(m*d)
  !
      Real (SP) :: ERROR, A1, AO2N !oxygen to nitogren ratio for denitrification ( (4-45)*(4-46)) (gO2/gN)
  !
      Real (SP) :: FP1SO4, FP2SO4, HS2AV, XJ2, XJ2CH4, X1J2
  !
      Real (SP) :: PF, PPO4 !particulate iorganic ortho phosphorate
  !
      Real (SP) :: FLUXHS, FLUXHSCH4, VJCH4G
  !
  !sediment POM and particulate biogenic silica
      Real (SP), Allocatable, Dimension (:) :: PON1TM1S, PON2TM1S, &
     & PON3TM1S, POC1TM1S, POC2TM1S, POC3TM1S, POP1TM1S, POP2TM1S, &
     & POP3TM1S, PSISEDTM1S !sediment particulate organic silica  (mgSi/m^3)
  !
  !sediment inorganic chemicals.
      Real (SP), Allocatable, Dimension (:) :: NH41TM1S, NH42TM1S, &
     & NH4T2TM1S, NO31TM1S, NO32TM1S, NO3T2TM1S, PO41TM1S, PO42TM1S, &
     & PO4T2TM1S, SI1TM1S, SI2TM1S, SIT2TM1S, HS1TM1S, HS2TM1S, &
     & HST2TM1S, CH41TM1S, CH42TM1S, CH4T2TM1S, SO4T2TM1S !sediment total SO4 (2nd layer) (gO2/m^3)
  ! NH41TM1S,  &!sediment disslolved NH4 (1st layer) (mgN/m^3)
  !                                 NH42TM1S,  &!sediment dissolved NH4 (2nd layer) (mgN/m^3)
  !                                 NH4T2TM1S, &!sediment total NH4 (2nd layer) (mgN/m^3)
  !
  !                                 NO31TM1S,  &!sediment dissolved NO3 (1st layer) (mgN/m^3)
  !                                 NO32TM1S,  &!sediment dissolved NO3 (2nd layer) (mgN/m^3)
  !                                 NO3T2TM1S, &!sediment total NO3 (2nd layer) (mgN/m^3)
  !
  !
  !                                 PO41TM1S,  &!sediment dissolved PO4 (1st layer) (mgP/m^3)
  !                                 PO42TM1S,  &!sediment dissolved PO4 (2nd layer) (mgP/m^3)
  !                                 PO4T2TM1S, &!sediment total PO4 (2nd layer) (mgP/m^3)
  !
  !
  !                                 SI1TM1S,   &!sediment dissolved silica (1st layer) (mgSi/m^3)
  !                                 SI2TM1S,   &!sediment dissolved silica (2nd layer) (mgSi/m^3)
  !
  !                                 HS1TM1S,   &!sediment dissolved HS (1st layer) (gO2/m^3)
  !                                 HS2TM1S,   &!sediment dissolved HS (2nd layer) (gO2/m^3)
  !                                 HST2TM1S,  &!sediment total HS (2nd layer) (gO2/m^3)
  !
  !
  !                                 CH41TM1S,  &!sediment dissolved CH4 (1st layer) (gO2/m^3)
  !                                 CH42TM1S,  &!sediment dissolved CH4 (2nd layer) (gO2/m^3)
  !                                 CH4T2TM1S, &!sediment total CH4 (2nd layer) (gO2/m^3)
  !
  !                                 !SO41TM1S,  &!sediment dissolved SO4 (1st layer) (gO2/m^3)
  !                                 SO4T2TM1S   !sediment total SO4 (2nd layer) (gO2/m^3)
  !
  !Moved the following here from wqm_modules.F
  !benthic stress calculation
      Real (SP), Allocatable, Dimension (:) :: BENSTRTM1S, BFORMAXS, &
     & ISWBENS !Array indicating whether benthic stress accumulation is
  !set to zero
  !If sediment temperature > TEMPBEN, then
  !    yes, it is set to zero   and ISWBENS=1, BFORMAX=0
  !otherwise ISWBENS is 0,
  !      BFOR=BFORMAX
  !
      Real (SP), Allocatable, Dimension (:) :: SODTM1S, JNH4TM1S, &
     & JNO3TM1S, JPO4TM1S, JCH4TM1S, JCH4GTM1S, JHSTM1S, JSITM1S !dissolved Si flux to water column (gO2/m^2/d)
  !
  !Moved here from wqm_modules.F
  !sediment burial fluxes of C, N, P
      Real (SP), Allocatable, Dimension (:) :: BURIALC, BURIALN, &
     & BURIALP !Burial flux of phosphorus (mgP/m^2/day)
  !    BURIALC,    & !Burial flux of carbon    (mgC/m^2/day)
  !                                                     !(positive leaving sediment)
  !                                     BURIALN,    & !Burial flux of nitrogen    (mgN/m^2/day)
  !                                                     !(positive leaving sediment)
  !                                     BURIALP       !Burial flux of phosphorus (mgP/m^2/day)
  !(positive leaving sediments)
  !
  !Moved here from wqm_modules.F
      Real (SP), Allocatable, Dimension (:) :: DIAGENC !carbon diagenesis flux (gC/m^2/day)
  !
  !
  !Moved these to wqm_sed.F.
  !look up tables for temperathre control used for interpolation at given
  !temperature
      Real (SP), Dimension (350) :: ZHTADP, ZHTADD, ZHTANH4F, ZHTANH4S, &
     & ZHTANO3F, ZHTANO3S, ZHTAK2NO3, ZHTA1HSD, ZHTA1HSP, ZHTASI, &
     & ZL12NOM, ZW12NOM, ZHTAPON1, ZHTAPON2, ZHTAPON3, ZHTAPOC1, &
     & ZHTAPOC2, ZHTAPOC3, ZHTAPOP1, ZHTAPOP2, ZHTAPOP3, ZHTACH4 !sediment methane oxidation rate look up table (m/day) after applying temperature effect
  !
  !Moved this to wqm_sed.F.
  !POM diagenesis rate temperature theta.
  !used in setting the temperature control look up tables
  !
      Real (SP) :: THTAPOP1, THTAPOP2, THTAPOP3, THTAPON1, THTAPON2, &
     & THTAPON3, THTAPOC1, THTAPOC2, THTAPOC3 !(G3)
  !  THTAPOP1,  &  !sediment POP diagenesis reaction temperature theta
  !                    !(G1)
  !      THTAPOP2,  &  !(G2)
  !      THTAPOP3,  &  !(G3)
  !      THTAPON1,  &  !sediment PON diagenesis reaction temperature theta
  !                    !(G1)
  !      THTAPON2,  &  !(G2)
  !      THTAPON3,  &  !(G3)
  !      THTAPOC1,  &  !sediment POC diagenesis reaction temperature theta
  !                    !(G1)
  !      THTAPOC2,  &  !(G2)
  !      THTAPOC3      !(G3)
  !
  !Moved this here from wqm_modules.F
  !time averaged G3 particulate organic flux (mg/m^2/day)
      Real (SP), Allocatable, Dimension (:) :: AG3CFL, AG3NFL, AG3PFL !G3 sediment POP flux
  !									AG3CFL,  &     !G3 sediment POC flux
  !                                                  !(mgC/m^2/day)
  !                                 AG3NFL,  &       !G3 sediment PON flux
  !                                                  !(mgN/m^2/day)
  !                                 AG3PFL           !G3 sediment POP flux
  !                                                  !(mgP/m^2/day)
  !
   !time averaged sediment temperature (degC)
      Real (SP), Allocatable, Dimension (:) :: ASDTMP !time averaged sediment
  !temperature  (degC)
!
      Real (SP) :: PIE1, PIE2, W12, TEMPD, O20, CH4SAT, SAL !overlying water salinity  (ppt)
  !    	  PIE1,          & !sediment 1st layer partition coefficient (L/kg)
  !                             !temporary variable
  !         PIE2,          & !sediment 2nd layer partition coefficient (L/kg)
  !                             !temporary variable
  !         W12,           & !particle mixing velocity (cm/yr) between Layer 1
  !                             !and layer 2
  !         TEMPD,         & !sediment temperature (degC) temporary variable
  !         O20,           & !overlying water concentration (mgO2/L)
  !         CH4SAT,        & !saturation concentration of CH4 in sediment pore
  !                             !water (mgO2/L)
  !         SAL              !overlying water salinity  (ppt)
  !
  !
      Real (SP) :: XAPPNH4, XAPP1HSD, XAPP1HSP, XAPP1NO3, XK2NO3, XKSI, &
     & XAPPCH4, TEMP20, TEMP202, FD1, FP1, FD2, FP2, SOD, CSOD, S, &
     & W12NOM, HSO4, DDSO4, CSODHS, CSODCH4, BENSTR, BENSTRS, &
     & BENSTRTM1, SODTM1, ISWBEN, BFORMAX !benthic stress limitation temporary variable (dimensionless)
  !
      Real (SP) :: ZHTANH4, ZHTANO3 !loop up table value for denitrification rate after
  !     ZHTANH4,    & 	!look up table value for nitrification rate after applying temperature control
  !                       !not used, replaced by ZHTANH4S, ZHTANH4F
  !         ZHTANO3       !loop up table value for denitrification rate after
  !applying temperature control
  !not used, replaced by ZHTANO3F. ZHTANO3S
  !
  !Moved here from wqm_modules.F
      Real (SP) :: TINTIM !time (days) for steady state simulation of sediments
  !
  !global variables fror initial condition and also output
      Real (SP), Allocatable :: CTEMP_GL (:), PO4T2TM1S_GL (:), &
     & NH4T2TM1S_GL (:), NO3T2TM1S_GL (:), HST2TM1S_GL (:), &
     & CH4T2TM1S_GL (:), CH41TM1S_GL (:), SO4T2TM1S_GL (:), SIT2TM1S_GL &
     & (:), BENSTRTM1S_GL (:), SODTM1S_GL (:), CPOP_GL (:, :), CPON_GL &
     & (:, :), CPOC_GL (:, :), CPOS_GL (:), PO41TM1S_GL (:), &
     & NH41TM1S_GL (:), NO31TM1S_GL (:), HS1TM1S_GL (:), SI1TM1S_GL (:)
  !
  !used for outputing sediment flux variables into history and station files
  !
      Real (SP), Allocatable :: JPOC_GL (:, :), JPON_GL (:, :), JPOP_GL &
     & (:, :), JPOS_GL (:), CH40_GL (:), JNH4_GL (:), JNO3_GL (:), &
     & BENDEN_GL (:), JCH4_GL (:), JCH4G_GL (:), JHS_GL (:), JPO4_GL &
     & (:), JSI_GL (:), NH41_GL (:), NH42_GL (:), NO31_GL (:), NO32_GL &
     & (:), PO41_GL (:), PO42_GL (:), SI1_GL (:), SI2_GL (:), CH42_GL &
     & (:), HS1_GL (:), HS2_GL (:), HSED1_GL (:)!,    &
  !	JPOC_GL(:,:),    &
  !                             JPON_GL(:,:),    &
  !                             JPOP_GL(:,:),    &
  !                             JPOS_GL(:),    &
  !                             !O20_GL(:),        &!given by C2_GL(:,:,27)
  !                             !D_GL(:),        &!already defined
  !                             !T_GL(:,:),        &!already defined
  !                             !NH40_GL(:),    &!given by C2_GL
  !                             !NO30_GL(:)),    &!given by C2_GL
  !                             !SI0_GL(:)),    &!given by C2_GL
  !                             CH40_GL(:),        &
  !                             !SAL_GL(:),        &!given by C2_GL
  !                             !SOD_GL(:),        &!given by SODTM1S_GL
  !                             JNH4_GL(:),        &
  !                             JNO3_GL(:),        &
  !                             BENDEN_GL(:),    &
  !                             JCH4_GL(:),        &
  !                             JCH4G_GL(:),    &
  !                             JHS_GL(:),        &
  !                             JPO4_GL(:),        &
  !                             JSI_GL(:),        &
  !                             NH41_GL(:),        &
  !                             NH42_GL(:),        &
  !                             NO31_GL(:),        &
  !                             NO32_GL(:),        &
  !                             PO41_GL(:),        &
  !                             PO42_GL(:),        &
  !                             SI1_GL(:),        &
  !                             SI2_GL(:),        &
  !                             !CH41_GL(:),    &!given by CH41TM1S_GL
  !                             CH42_GL(:),        &
  !                             HS1_GL(:),        &
  !                             HS2_GL(:),        &
  !                             !POC1_GL(:),    &!given by CPOC_GL
  !                             !POC2_GL(:),    &
  !                             !POC3_GL(:),    &
  !                             !PON1_GL(:),    &!given by CPON_GL
  !                             !PON2_GL(:),    &
  !                             !PON3_GL(:),    &
  !                             !POP1_GL(:),    &!given by CPOP_GL
  !                             !POP2_GL(:),    &
  !                             !POP3_GL(:),    &
  !                             !PSISED_GL(:),    &!given by CPOS_GL
  !                             HSED1_GL(:)!,    &
  !BENSTR_GL(:)     !given by BENSTRTM1S_GL
  !
  !
  !diagnostic flux from P,N,C, Si for sediment module
      Real (SP), Allocatable, Dimension (:) :: DIAGP, DIAGN, DIAGC, &
     & DIAGS !!not used
  !
      Real (SP), Allocatable, Dimension (:) :: MTVEL !mass transfer velocity (m/sec) (=SOD/O20)
  !
  !The following needs to be moved to SED_INIT
      Real (SP), Allocatable :: CPOP (:, :), CPON (:, :), CPOC (:, :), &
     & CPOS (:)!Sediment particulate organic silica (mgSi/m^3)
  !
      Real (SP), Allocatable, Dimension (:) :: CPO4, CNO3, CNH4, CSI, &
     & CCH4, CSO4, CHS !Sediment HS concentration (gO2/m^3)
  !
  !Initial conditions of sediment constitutients (moved here from mod_wqminit.F) for individual cell
      Real (SP) :: CPOPI (3), CPONI (3), CPOCI (3), CPOSI, PO41TI, &
     & PO4T2I, NO31TI, NO3T2I, NH41TI, NH4T2I, CH4T2I, CH41TI, SO41TI, &
     & SO4T2I, HS1TI, HST2I, SI1TI, SIT2I, CTEMPI, BENSTI, SODI !Initial SOD mgO2/m^2/day
  !
  !
      Real (SP), Allocatable :: JPOP (:, :), JPON (:, :), JPOC (:, :), &
     & JPOS (:)!POSi settling flux(mgSi/m^2/day) positive into sediments
  !
      Real (SP), Allocatable :: JPOPaccum (:, :), JPONaccum (:, :), &
     & JPOCaccum (:, :), JPOSaccum (:)!accumulative POSi settling flux(mgSi/m^2/day) positive into sediments
  !
      Real (SP) POC1R !mgC/gSediments  reference POC concentration in sediments for calculating W12NOM (DMD page 277)
  !=0.1mgC/gSediments ~ 100000*M2 mgC/m^3 ( if M2 has unit kg/L for particle concentration of layer 2)
  !(DMD book page 278, table 13.1, eqn(13.1))
  !
      Integer :: IDEBUG_SED
  !
      Character (Len=72), Dimension (NFLP) :: WCLFN !sediment overlying water column forcing filenames
  !
  !variables related to records in WCLFN which provids overlying water column condition as forcing to
  !sediment flux model
  !
  !: will need to make this as field arrays for all TCE cells rather then a domain wide uniform value
  !
      Real (SP) :: NXJDAY
      Real (SP) :: NXWCL
      Real (SP) :: JDAY_R1 !JDAY of 1st record kept in memory                    (day)
      Real (SP) :: JDAY_R2 !JDAY of 2nd record kept in memory                    (day)
  !
      Real (SP) :: JCIN_R1, JCIN_R2 !1st and 2nd record of Jcin                         (gO2/m^2/d)
      Real (SP) :: JNIN_R1, JNIN_R2 !1st and 2nd record of Jnin                         (gN/m^2/d)
      Real (SP) :: JPIN_R1, JPIN_R2 !1st and 2nd record of Jpin                         (gP/m^2/d)
      Real (SP) :: JSIN_R1, JSIN_R2 !1st and 2nd record of Jsin                         (gSi/m^2/d)
  !
      Real (SP) :: O20_R1, O20_R2 !1st and 2nd record of O20                             (mgO2/L)
      Real (SP) :: D_R1, D_R2 !1st and 2nd record of depth                         (m)
      Real (SP) :: TW_R1, TW_R2 !1st and 2nd record of overlying water temperature     (degC)
  !
      Real (SP) :: NH30_R1, NH30_R2 !1st and 2nd record of overlying ammonia             (mgN/L)
      Real (SP) :: NO30_R1, NO30_R2 !1st and 2nd record of overlying nitrate             (mgN/L)
      Real (SP) :: PO40_R1, PO40_R2 !1st and 2nd record of overlying phosphate             (mgP/L)
      Real (SP) :: SIAT0_R1, SIAT0_R2 !1st and 2nd record of overlying total silicate     (mgSi/L)
      Real (SP) :: CH40_R1, CH40_R2 !1st and 2nd record of overlying water CH4             (mgO2/L)
      Real (SP) :: SALT0_R1, SALT0_R2 !1st and 2nd record of overlying water salinity     (psu)
  !
      Real (SP) :: SODMIN, SODMAX !min and max value of SOD
Contains
  !
  !subroutines:
  !subroutine SED_INIT()
  !subroutine SED_ALLOC()
  !subroutine SED_DEALLOC()
  !subroutine SED_READ()
  !subroutine POM_ACCUMUL
  !subroutine SED_CALC()
  !subroutine SED_DIAGENESIS_G3()
  !
  !subroutine SEDTSFNL()
  !subroutine SEDSSFNL()
  !
  !functions:
  !function SEDF()
  !function ZBRENT()
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   SED_INIT                          **
  !********************************************************************************
      Subroutine SED_INIT
         Use MOD_WQMINIT, Only: UNI_ICON_IN_SED_VAR
    !
    !Moved from wqm_modules.F
         SSNAME (1) = 'Sediment Temperature'
         SSNAME (2) = 'Sediment POP        '
         SSNAME (3) = 'Sediment PON        '
         SSNAME (4) = 'Sediment POC        '
         SSNAME (5) = 'Sediment PBS        '
         SSNAME (6) = 'Sediment PO4        '
         SSNAME (7) = 'Sediment NH4        '
         SSNAME (8) = 'Sediment NO3        '
         SSNAME (9) = 'Sediment HS         '
         SSNAME (10) = 'Sediment CH4        '
         SSNAME (11) = 'Sediment CH4        '
         SSNAME (12) = 'Sediment SO4        '
         SSNAME (13) = 'Sediment DSIL       '
         SSNAME (14) = 'Benthic Stress      '
         SSNAME (15) = 'Benthic Algae       '
         SSNAME (16) = 'Deposit Feeders     '
         SSNAME (17) = 'Suspension Feeders  '
    !
    !allocate variables related to sediment diagenesis
         Call SED_ALLOC
    !
    !initialize time varying overlying water colum input related variables
         ITVWCLF = 0 !Assume no input file is needed
         NTVWCLF = 0 !assume no time varying overlying water column condition input file
         WCLPTR = 0 !pointer in WCLFN() array
         WCLFDAYS = 0.0 !lenght of file set to zero
    !
         If ( .Not. UNI_ICON_IN_SED_VAR) Then !: added this if.  IF UNI_ICON_IN_SED_VAR is true (ie, ICIC=WQ_SEDV in wqm_con.npt), then sediment initialization occured in wqm_inputs.F
            Call SED_INIT_ICI !get the input initial conditions first, only if providing constant,single values for each sediment variable (ICIC=UNIFORM in wqm_con.npt)
         End If
    !
      End Subroutine SED_INIT
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   SED_ALLOC                         **
  !********************************************************************************
      Subroutine SED_ALLOC
    !
    !USE MOD_SED_DF_EXCHANGE_VARS, ONLY: &
    !    POC1TM1S_SED_DF,                &
    !    POC2TM1S_SED_DF
    !
    !Moved here from wqm_modules.F
         Allocate (FRPOP(MTLOC, 3))
         FRPOP = 0.0
         Allocate (FRPON(MTLOC, 3))
         FRPON = 0.0
         Allocate (FRPOC(MTLOC, 3))
         FRPOC = 0.0
    !
    !Moved here from wqm_modules.F
         Allocate (WSSBNET(MTLOC))
         WSSBNET = 0.0
         Allocate (WSLBNET(MTLOC))
         WSLBNET = 0.0
         Allocate (WSRBNET(MTLOC))
         WSRBNET = 0.0
         Allocate (WS1BNET(MTLOC))
         WS1BNET = 0.0
         Allocate (WS2BNET(MTLOC))
         WS2BNET = 0.0
         Allocate (WS3BNET(MTLOC))
         WS3BNET = 0.0
         Allocate (WSUBNET(MTLOC))
         WSUBNET = 0.0
    !
    !Moved here from wqm_modules.F
         Allocate (VSED(MTLOC))
         VSED = 0.0
         Allocate (VPMIX(MTLOC))
         VPMIX = 0.0
         Allocate (VDMIX(MTLOC))
         VDMIX = 0.0
    !
    !Moved here from wqm_modules.F
         Allocate (PON1TM1S(MTLOC))
         PON1TM1S = 0.0
         Allocate (PON2TM1S(MTLOC))
         PON2TM1S = 0.0
         Allocate (PON3TM1S(MTLOC))
         PON3TM1S = 0.0
         Allocate (POC1TM1S(MTLOC))
         POC1TM1S = 0.0
         Allocate (POC2TM1S(MTLOC))
         POC2TM1S = 0.0
         Allocate (POC3TM1S(MTLOC))
         POC3TM1S = 0.0
         Allocate (POP1TM1S(MTLOC))
         POP1TM1S = 0.0
         Allocate (POP2TM1S(MTLOC))
         POP2TM1S = 0.0
         Allocate (POP3TM1S(MTLOC))
         POP3TM1S = 0.0
         Allocate (PSISEDTM1S(MTLOC))
         PSISEDTM1S = 0.0
    !
    !Added POC1TM1S_SED_DF
    !     ALLOCATE(POC1TM1S_SED_DF(MTLOC));    POC1TM1S_SED_DF=0.0
    !     ALLOCATE(POC2TM1S_SED_DF(MTLOC));  POC2TM1S_SED_DF=0.0
    !
    !Moved here from wqm_modules.F
         Allocate (BENSTRTM1S(MTLOC))
         BENSTRTM1S = 0.0
         Allocate (BFORMAXS(MTLOC))
         BFORMAXS = 0.0
         Allocate (ISWBENS(MTLOC))
         ISWBENS = 0.0
    !
         Allocate (SODTM1S(MTLOC))
         SODTM1S = 0.0
         Allocate (JNH4TM1S(MTLOC))
         JNH4TM1S = 0.0
         Allocate (JNO3TM1S(MTLOC))
         JNO3TM1S = 0.0
         Allocate (JPO4TM1S(MTLOC))
         JPO4TM1S = 0.0
         Allocate (JSITM1S(MTLOC))
         JSITM1S = 0.0
         Allocate (JCH4TM1S(MTLOC))
         JCH4TM1S = 0.0
         Allocate (JCH4GTM1S(MTLOC))
         JCH4GTM1S = 0.0
         Allocate (JHSTM1S(MTLOC))
         JHSTM1S = 0.0
    !
    !Moved here from wqm_modules.F
         Allocate (NH41TM1S(MTLOC))
         NH41TM1S = 0.0
         Allocate (NH42TM1S(MTLOC))
         NH42TM1S = 0.0
         Allocate (NH4T2TM1S(MTLOC))
         NH4T2TM1S = 0.0
         Allocate (NO31TM1S(MTLOC))
         NO31TM1S = 0.0
         Allocate (NO32TM1S(MTLOC))
         NO32TM1S = 0.0
         Allocate (NO3T2TM1S(MTLOC))
         NO3T2TM1S = 0.0
         Allocate (HS1TM1S(MTLOC))
         HS1TM1S = 0.0
         Allocate (HS2TM1S(MTLOC))
         HS2TM1S = 0.0
         Allocate (HST2TM1S(MTLOC))
         HST2TM1S = 0.0
         Allocate (SI1TM1S(MTLOC))
         SI1TM1S = 0.0
         Allocate (SI2TM1S(MTLOC))
         SI2TM1S = 0.0
         Allocate (SIT2TM1S(MTLOC))
         SIT2TM1S = 0.0
         Allocate (PO41TM1S(MTLOC))
         PO41TM1S = 0.0
         Allocate (PO42TM1S(MTLOC))
         PO42TM1S = 0.0
         Allocate (PO4T2TM1S(MTLOC))
         PO4T2TM1S = 0.0
    !
    !Moved here from wqm_modules.F
         Allocate (CH41TM1S(MTLOC))
         CH41TM1S = 0.0
         Allocate (CH42TM1S(MTLOC))
         CH42TM1S = 0.0
         Allocate (CH4T2TM1S(MTLOC))
         CH4T2TM1S = 0.0
         Allocate (SO4T2TM1S(MTLOC))
         SO4T2TM1S = 0.0
    !
    !Moved here from wqm_modules.F
         Allocate (BURIALN(MTLOC))
         BURIALN = 0.0
         Allocate (BURIALP(MTLOC))
         BURIALP = 0.0
         Allocate (BURIALC(MTLOC))
         BURIALC = 0.0
    !
    !Moved here from wqm_modules.F
         Allocate (DIAGENC(MTLOC))
         DIAGENC = 0.0
    !
    !Moved this here from wqm_modules.F
         Allocate (AG3CFL(MTLOC))
         AG3CFL = 0.0
         Allocate (AG3NFL(MTLOC))
         AG3NFL = 0.0
         Allocate (AG3PFL(MTLOC))
         AG3PFL = 0.0
    !
    !Moved this here from wqm_modules.F
         Allocate (ASDTMP(MTLOC))
         ASDTMP = 0.0
    !
    !Moved this here from wqm_modules.F
    !   ALLOCATE(FIB(MTLOC));                FIB   = 0.0
    !   ALLOCATE(NLB(MTLOC));                NLB   = 0.0
    !   ALLOCATE(PLB(MTLOC));                PLB   = 0.0
    !   ALLOCATE(NPPB(MTLOC));               NPPB  = 0.0
    !   ALLOCATE(BBM(MTLOC));                BBM   = 0.0
    !   ALLOCATE(BLITE(MTLOC));              BLITE = 0.0
    !   ALLOCATE(BMB(MTLOC));                BMB   = 0.0
    !   ALLOCATE(PB(MTLOC));                 PB    = 0.0
    !   ALLOCATE(PRB(MTLOC));                PRB   = 0.0
    !   ALLOCATE(BANH4(MTLOC));              BANH4 = 0.0
    !   ALLOCATE(BANO3(MTLOC));              BANO3 = 0.0
    !   ALLOCATE(BAPO4(MTLOC));              BAPO4 = 0.0
    !   ALLOCATE(BADOC(MTLOC));              BADOC = 0.0
    !   ALLOCATE(BADO(MTLOC));               BADO  = 0.0
    !   ALLOCATE(BAPOC(MTLOC));              BAPOC = 0.0
    !   ALLOCATE(BAPON(MTLOC));              BAPON = 0.0
    !   ALLOCATE(BAPOP(MTLOC));              BAPOP = 0.0
    !
         Allocate (HSED(MTLOC))
         HSED = 0.0 !Moved from mod_wqm.F
         Allocate (HSED1(MTLOC))
         HSED1 = 0.0 !Moved from mod_wqm.F
    !
    !used for writing and reading restart file
         Allocate (CTEMP_GL(MGL))
         CTEMP_GL = 0.0
         Allocate (CPOP_GL(MGL, 3))
         CPOP_GL = 0.0
         Allocate (CPON_GL(MGL, 3))
         CPON_GL = 0.0
         Allocate (CPOC_GL(MGL, 3))
         CPOC_GL = 0.0
         Allocate (CPOS_GL(MGL))
         CPOS_GL = 0.0
         Allocate (PO4T2TM1S_GL(MGL))
         PO4T2TM1S_GL = 0.0
         Allocate (NH4T2TM1S_GL(MGL))
         NH4T2TM1S_GL = 0.0
         Allocate (NO3T2TM1S_GL(MGL))
         NO3T2TM1S_GL = 0.0
         Allocate (HST2TM1S_GL(MGL))
         HST2TM1S_GL = 0.0
         Allocate (CH4T2TM1S_GL(MGL))
         CH4T2TM1S_GL = 0.0
         Allocate (CH41TM1S_GL(MGL))
         CH41TM1S_GL = 0.0
         Allocate (SO4T2TM1S_GL(MGL))
         SO4T2TM1S_GL = 0.0
         Allocate (SIT2TM1S_GL(MGL))
         SIT2TM1S_GL = 0.0
    !
         Allocate (BENSTRTM1S_GL(MGL))
         BENSTRTM1S_GL = 0.0
         Allocate (SODTM1S_GL(MGL))
         SODTM1S_GL (:) = 0.d0
    !
    !: save other sediment variables needed for restarts
         Allocate (PO41TM1S_GL(MGL))
         PO41TM1S_GL = 0.0
         Allocate (NH41TM1S_GL(MGL))
         NH41TM1S_GL = 0.0
         Allocate (NO31TM1S_GL(MGL))
         NO31TM1S_GL = 0.0
         Allocate (HS1TM1S_GL(MGL))
         HS1TM1S_GL = 0.0
         Allocate (SI1TM1S_GL(MGL))
         SI1TM1S_GL = 0.0
    !used for outputing sediment flux variables into history and station files
    !
         Allocate (JPOC_GL(MGL, 3))
         JPOC_GL = 0.d0
         Allocate (JPON_GL(MGL, 3))
         JPON_GL = 0.d0
         Allocate (JPOP_GL(MGL, 3))
         JPOP_GL = 0.d0
         Allocate (JPOS_GL(MGL))
         JPOS_GL = 0.d0
    !ALLOCATE(O20_GL(MGL));            O20_GL        =0.d0    !given by C2_GL(:,:,27)
    !ALLOCATE(D_GL(MGL));                D_GL        =0.d0    !already defined
    !ALLOCATE(T_GL(MGL,KBM1));           T_GL        =0.d0    !already defined
    !ALLOCATE(NH40_GL(MGL));          NH40_GL        =0.d0    !given by C2_GL
    !ALLOCATE(NO30_GL(MGL));          NO30_GL        =0.d0    !given by C2_GL
    !ALLOCATE(SI0_GL(MGL));            SI0_GL        =0.d0    !given by C2_GL
         Allocate (CH40_GL(MGL))
         CH40_GL = 0.d0
    !ALLOCATE(SAL_GL(MGL));            SAL_GL        =0.d0    !given by C2_GL
    !ALLOCATE(SOD_GL(MGL));            SOD_GL        =0.d0    !given by SODTM1S_GL
         Allocate (JNH4_GL(MGL))
         JNH4_GL = 0.d0 !
         Allocate (JNO3_GL(MGL))
         JNO3_GL = 0.d0
         Allocate (BENDEN_GL(MGL))
         BENDEN_GL = 0.d0
         Allocate (JCH4_GL(MGL))
         JCH4_GL = 0.d0
         Allocate (JCH4G_GL(MGL))
         JCH4G_GL = 0.d0
         Allocate (JHS_GL(MGL))
         JHS_GL = 0.d0
         Allocate (JPO4_GL(MGL))
         JPO4_GL = 0.d0
         Allocate (JSI_GL(MGL))
         JSI_GL = 0.d0
         Allocate (NH41_GL(MGL))
         NH41_GL = 0.d0
         Allocate (NH42_GL(MGL))
         NH42_GL = 0.d0
    !
         Allocate (NO31_GL(MGL))
         NO31_GL = 0.d0
         Allocate (NO32_GL(MGL))
         NO32_GL = 0.d0
         Allocate (PO41_GL(MGL))
         PO41_GL = 0.d0
         Allocate (PO42_GL(MGL))
         PO42_GL = 0.d0
         Allocate (SI1_GL(MGL))
         SI1_GL = 0.d0
         Allocate (SI2_GL(MGL))
         SI2_GL = 0.d0
    !ALLOCATE(CH41_GL(MGL));            CH41_GL        =0.d0    !given by CH41TM1S_GL
         Allocate (CH42_GL(MGL))
         CH42_GL = 0.d0
         Allocate (HS1_GL(MGL))
         HS1_GL = 0.d0
         Allocate (HS2_GL(MGL))
         HS2_GL = 0.d0
    !
    !ALLOCATE(POC1_GL(MGL));            POC1_GL        =0.d0    !given by CPOC_GL
    !ALLOCATE(POC2_GL(MGL));            POC2_GL        =0.d0
    !ALLOCATE(POC3_GL(MGL));            POC3_GL        =0.d0
    !
    !ALLOCATE(PON1_GL(MGL));            PON1_GL        =0.d0    !given by CPON_GL
    !ALLOCATE(PON2_GL(MGL));            PON2_GL        =0.d0
    !ALLOCATE(PON3_GL(MGL));            PON3_GL        =0.d0
    !
    !ALLOCATE(POP1_GL(MGL));            POP1_GL        =0.d0    !given by CPOP_GL
    !ALLOCATE(POP2_GL(MGL));            POP2_GL        =0.d0
    !ALLOCATE(POP3_GL(MGL));            POP3_GL        =0.d0
    !
    !ALLOCATE(PSISED_GL(MGL));            PSISED_GL    =0.d0    !given by CPOS_GL
         Allocate (HSED1_GL(MGL))
         HSED1_GL = 0.d0 !
    !ALLOCATE(BENSTR_GL(MGL));            BENSTR_GL    =0.d0    !given by BENSTRTM1S_GL
    !
    !these are the outputs of the excel version (Greg Pelletier)
    !
    !!',',(JPOC(I,1)+JPOC(I,2)+JPOC(I,3))/1000.0*2.667,    &    !'Jcin'                mgC/m^2/d            gO2/m^2/d
    !!',',(JPON(I,1)+JPON(I,2)+JPON(I,3))/1000.0,    &            !'Jnin'                *mgN/m^2/d            gN/m^2/d
         !!',',(JPOP(I,1)+JPOP(I,2)+JPOP(I,3))/1000.0,    &            !'Jpin'                *mgP/m^2/d            gP/m^2/d
    !!',',JPOS(I)/1000.0,    &                                    !'Jsin                 *mgSi/m^2/d            gSi/m^2/d
         !!',',O20,            &                                    !'O20'                mgO2/L                mgO2/L
    !!',',D(I),            &                                    !'Depth'            m                    m
         !!',',T(I,KBM1),        &                                    !'Tw'                degC                degC
    !!',',NH40/1000.0,    &                                    !'NH30'                *mgN/m^3            mgN/L
         !!',',NO30/1000.0,    &                                    !'NO30'                *mgN/m^3            mgN/L
    !!',',PO40/1000.0,    &                                    !'PO40'                *mgP/m^3            mgP/L
         !!',',SI0/1000.0,        &                                    !'SI0'                *mgSi/m^3             mgSi/L
    !!',',CH40,            &                                    !'CH40'                gO2/m^3                mgO2/L
         !!',',SAL,            &                                    !'SALw'              ppt                    ppt
    !!',',SOD,            &                                    !'SOD'                gO2/m^2/day            gO2/m^2/d
         !!',',JNH4/1000.0,    &                                    !'Jnh4                *mgN/m^2/day        gN/m^2/d
    !!',',JNO3/1000.0,    &                                    !'Jno3                *mgN/m^2/day        gN/m^2/d
         !!',',BENDEN(1),        &                                    !'JDenitT'            gN/m^2/day            gN/m^2/d
    !!                        !JDenit(1) = Denit(1) * NO3(1)
    !!                        !JDenit(2) = Denit(2) * NO3(2)
    !!                        !JDenitT = JDenit(1) + JDenit(2)
    !!',',JCH4,            &                                    !'Jch4'                gO2/m^2/day            gO2/m^2/d
    !!',',JCH4G,            &                                    !'Jch4g'            gO2/m^2/day            gO2/m^2/d
         !!',',JHS,            &                                    !'Jhs'                gO2/m^2/day            gO2/m^2/d
    !!',',JPO4/1000.0,    &                                    !'Jpo4'                *mgP/m^2/day        gP/m^2/d
         !!',',JSI/1000.0,        &                                    !'Jsi'                *mgSi/m^2/day        gSi/m^2/d
    !
    !!',',NH41/1000.0,    &                    !'NH3(1)'            *mgN/m^3             mgN/L    (dissolved)
    !!',',NH42/1000.0,    &                    !'NH3(2)'            *mgN/m^3             mgN/L    (dissolved)
         !!',',NO31/1000.0,    &                    !'NO3(1)'            *mgN/m^3             mgN/L    (dissolved)
    !!',',NO32/1000.0,    &                    !'NO3(2)'            *mgN/m^3             mgN/L    (dissolved)
         !!',',PO41/1000.0,    &                    !'PO4(1)'            *mgP/m^3             mgP/L    (dissolved)
    !!',',PO42/1000.0,    &                    !'PO4(2)'            *mgP/m^3             mgP/L    (dissolved)
         !!',',SI1/1000.0,        &                    !'Si(1)'            *mgSi/m^3             mgSi/L    (dissolved)
    !!',',SI2/1000.0,        &                    !'Si(2)'            *mgSi/m^3             mgSi/L    (dissolved)
         !!',',CH41,            &                    !'CH4(1)'            gO2/m^3             mgO2/L    (dissolved)
    !!',',CH42,            &                    !'CH4(2)'            gO2/m^3             mgO2/L    (dissolved)
         !!',',HS1,            &                    !'HS(1)'            gO2/m^3             mgO2/L    (dissolved)
    !!',',HS2,            &                    !'HS(2)'            gO2/m^3             mgO2/L    (dissolved)
         !!',',POC1/1000.0*2.667,&                    !'POC2(1)'            *mgC/m^3             gO2/m^3    (G1 of POC in layer 2)
    !!',',POC2/1000.0*2.667,&                    !'POC2(2)'             *mgC/m^3             gO2/m^3    (G2 of POC in layer 2)
         !!',',POC3/1000.0*2.667,&                    !'POC2(3)'             *mgC/m^3             gO2/m^3    (G3 of POC in layer 2)
    !!',',PON1/1000.0,    &                    !'PON2(1)'            *mgN/m^3             gN/m^3    (G1 of PON in layer 2)
         !!',',PON2/1000.0,    &                    !'PON2(2)'             *mgN/m^3             gN/m^3    (G2 of PON in layer 2)
    !!',',PON3/1000.0,    &                    !'PON2(3)'             *mgN/m^3             gN/m^3    (G3 of PON in layer 2)
         !!',',POP1/1000.0,    &                    !'POP2(1)'            *mgP/m^3             gP/m^3    (G1 of POP in layer 2)
    !!',',POP2/1000.0,    &                    !'POP2(2)'            *mgP/m^3             gP/m^3    (G2 of POP in layer 2)
         !!',',POP3/1000.0,    &                    !'POP2(3)'             *mgP/m^3             gP/m^3    (G3 of POP in layer 2)
    !!',',PSISED/1000.0,    &                    !'POS2'                *mgSi/m^3             gSi/m^3    (particulate biogenic Silic
         !!',',HSED1(I)*1000.0,&                    !'H1'                m                    mm
    !!',',BENSTR                                 !'BEN_STR'             dimensionless        dimensionless
    !
         Allocate (MTVEL(MTLOC))
         MTVEL = 0.0 !Moved from mod_wqm.F
         Allocate (CPOC(MTLOC, 3))
         CPOC = 0.0
         Allocate (CPOP(MTLOC, 3))
         CPOP = 0.0
         Allocate (CPON(MTLOC, 3))
         CPON = 0.0
         Allocate (CPOS(MTLOC))
         CPOS = 0.0
         Allocate (CPO4(MTLOC))
         CPO4 = 0.0
         Allocate (CNO3(MTLOC))
         CNO3 = 0.0
         Allocate (CNH4(MTLOC))
         CNH4 = 0.0
    !
    !
         Allocate (CCH4(MTLOC))
         CCH4 = 0.0
         Allocate (CSO4(MTLOC))
         CSO4 = 0.0
         Allocate (CHS(MTLOC))
         CHS = 0.0
         Allocate (CSI(MTLOC))
         CSI = 0.0
    !
         Allocate (DIAGP(MTLOC))
         DIAGP = 0.0
         Allocate (DIAGN(MTLOC))
         DIAGN = 0.0
         Allocate (DIAGC(MTLOC))
         DIAGC = 0.0
         Allocate (DIAGS(MTLOC))
         DIAGS = 0.0
    !
         Allocate (JPOC(MTLOC, 3))
         JPOC = 0.0
         Allocate (JPOP(MTLOC, 3))
         JPOP = 0.0
         Allocate (JPON(MTLOC, 3))
         JPON = 0.0
         Allocate (JPOS(MTLOC))
         JPOS = 0.0
    !
         Allocate (JPOCaccum(MTLOC, 3))
         JPOCaccum = 0.0
         Allocate (JPOPaccum(MTLOC, 3))
         JPOPaccum = 0.0
         Allocate (JPONaccum(MTLOC, 3))
         JPONaccum = 0.0
         Allocate (JPOSaccum(MTLOC))
         JPOSaccum = 0.0
    !
      End Subroutine SED_ALLOC
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   SED_DEALLOC                       **
  !********************************************************************************
  !
      Subroutine SED_DEALLOC
    !
    !Moved here from wqm_main.F
    !
         If (ALLOCATED(FRPOP)) DEALLOCATE (FRPOP)
         If (ALLOCATED(FRPON)) DEALLOCATE (FRPON)
         If (ALLOCATED(FRPOC)) DEALLOCATE (FRPOC)
    !
    !Moved here from wqm_main.F
         If (ALLOCATED(WSSBNET)) DEALLOCATE (WSSBNET)
         If (ALLOCATED(WSLBNET)) DEALLOCATE (WSLBNET)
         If (ALLOCATED(WSRBNET)) DEALLOCATE (WSRBNET)
         If (ALLOCATED(WS1BNET)) DEALLOCATE (WS1BNET)
         If (ALLOCATED(WS2BNET)) DEALLOCATE (WS2BNET)
         If (ALLOCATED(WS3BNET)) DEALLOCATE (WS3BNET)
         If (ALLOCATED(WSUBNET)) DEALLOCATE (WSUBNET)
    !
    !Moved here from wqm_main.F
         If (ALLOCATED(VSED)) DEALLOCATE (VSED)
         If (ALLOCATED(VPMIX)) DEALLOCATE (VPMIX)
         If (ALLOCATED(VDMIX)) DEALLOCATE (VDMIX)
    !
    !Moved here from wqm_main.F
         If (ALLOCATED(PON1TM1S)) DEALLOCATE (PON1TM1S)
         If (ALLOCATED(PON2TM1S)) DEALLOCATE (PON2TM1S)
         If (ALLOCATED(PON3TM1S)) DEALLOCATE (PON3TM1S)
         If (ALLOCATED(POC1TM1S)) DEALLOCATE (POC1TM1S)
         If (ALLOCATED(POC2TM1S)) DEALLOCATE (POC2TM1S)
         If (ALLOCATED(POC3TM1S)) DEALLOCATE (POC3TM1S)
         If (ALLOCATED(POP1TM1S)) DEALLOCATE (POP1TM1S)
         If (ALLOCATED(POP2TM1S)) DEALLOCATE (POP2TM1S)
         If (ALLOCATED(POP3TM1S)) DEALLOCATE (POP3TM1S)
         If (ALLOCATED(PSISEDTM1S)) DEALLOCATE (PSISEDTM1S)
    !
    !Added POC1TM1S_SED_DF POC2TM1S_SED_DF
    !
    !     IF(ALLOCATED(POC1TM1S_SED_DF))DEALLOCATE(POC1TM1S_SED_DF)
    !     IF(ALLOCATED(POC2TM1S_SED_DF))DEALLOCATE(POC2TM1S_SED_DF)
    !
    !
    !Moved here from wqm_main.F
         If (ALLOCATED(BENSTRTM1S)) DEALLOCATE (BENSTRTM1S)
         If (ALLOCATED(BFORMAXS)) DEALLOCATE (BFORMAXS)
         If (ALLOCATED(ISWBENS)) DEALLOCATE (ISWBENS)
    !
         If (ALLOCATED(SODTM1S)) DEALLOCATE (SODTM1S)
         If (ALLOCATED(JNH4TM1S)) DEALLOCATE (JNH4TM1S)
         If (ALLOCATED(JNO3TM1S)) DEALLOCATE (JNO3TM1S)
         If (ALLOCATED(JPO4TM1S)) DEALLOCATE (JPO4TM1S)
         If (ALLOCATED(JSITM1S)) DEALLOCATE (JSITM1S)
         If (ALLOCATED(JCH4TM1S)) DEALLOCATE (JCH4TM1S)
         If (ALLOCATED(JCH4GTM1S)) DEALLOCATE (JCH4GTM1S)
         If (ALLOCATED(JHSTM1S)) DEALLOCATE (JHSTM1S)
    !
    !Moved here from wqm_main.F
         If (ALLOCATED(NH41TM1S)) DEALLOCATE (NH41TM1S)
         If (ALLOCATED(NH42TM1S)) DEALLOCATE (NH42TM1S)
         If (ALLOCATED(NH4T2TM1S)) DEALLOCATE (NH4T2TM1S)
         If (ALLOCATED(NO31TM1S)) DEALLOCATE (NO31TM1S)
         If (ALLOCATED(NO32TM1S)) DEALLOCATE (NO32TM1S)
         If (ALLOCATED(NO3T2TM1S)) DEALLOCATE (NO3T2TM1S)
         If (ALLOCATED(HS1TM1S)) DEALLOCATE (HS1TM1S)
         If (ALLOCATED(HS2TM1S)) DEALLOCATE (HS2TM1S)
         If (ALLOCATED(HST2TM1S)) DEALLOCATE (HST2TM1S)
         If (ALLOCATED(SI1TM1S)) DEALLOCATE (SI1TM1S)
         If (ALLOCATED(SI2TM1S)) DEALLOCATE (SI2TM1S)
         If (ALLOCATED(SIT2TM1S)) DEALLOCATE (SIT2TM1S)
         If (ALLOCATED(PO41TM1S)) DEALLOCATE (PO41TM1S)
         If (ALLOCATED(PO42TM1S)) DEALLOCATE (PO42TM1S)
         If (ALLOCATED(PO4T2TM1S)) DEALLOCATE (PO4T2TM1S)
    !
    !Moved these here   wqm_main.F
         If (ALLOCATED(CH41TM1S)) DEALLOCATE (CH41TM1S)
         If (ALLOCATED(CH42TM1S)) DEALLOCATE (CH42TM1S)
         If (ALLOCATED(CH4T2TM1S)) DEALLOCATE (CH4T2TM1S)
         If (ALLOCATED(SO4T2TM1S)) DEALLOCATE (SO4T2TM1S)
    !Moved these here   wqm_main.F
    !
         If (ALLOCATED(BURIALN)) DEALLOCATE (BURIALN)
         If (ALLOCATED(BURIALP)) DEALLOCATE (BURIALP)
         If (ALLOCATED(BURIALC)) DEALLOCATE (BURIALC)
    !
    !Moved here from wqm_main.F
         If (ALLOCATED(DIAGENC)) DEALLOCATE (DIAGENC)
    !
    !Moved here from wqm_main.F
         If (ALLOCATED(AG3CFL)) DEALLOCATE (AG3CFL)
         If (ALLOCATED(AG3NFL)) DEALLOCATE (AG3NFL)
         If (ALLOCATED(AG3PFL)) DEALLOCATE (AG3PFL)
         If (ALLOCATED(ASDTMP)) DEALLOCATE (ASDTMP)
    !
         If (ALLOCATED(HSED)) DEALLOCATE (HSED)
         If (ALLOCATED(HSED)) DEALLOCATE (HSED1)
    !
    !global arrays used as input initial conditions and output
         If (ALLOCATED(CTEMP_GL)) DEALLOCATE (CTEMP_GL)
         If (ALLOCATED(CPOP_GL)) DEALLOCATE (CPOP_GL)
         If (ALLOCATED(CPON_GL)) DEALLOCATE (CPON_GL)
         If (ALLOCATED(CPOC_GL)) DEALLOCATE (CPOC_GL)
         If (ALLOCATED(CPOS_GL)) DEALLOCATE (CPOS_GL)
         If (ALLOCATED(PO4T2TM1S_GL)) DEALLOCATE (PO4T2TM1S_GL)
         If (ALLOCATED(NH4T2TM1S_GL)) DEALLOCATE (NH4T2TM1S_GL)
         If (ALLOCATED(NO3T2TM1S_GL)) DEALLOCATE (NO3T2TM1S_GL)
         If (ALLOCATED(HST2TM1S_GL)) DEALLOCATE (HST2TM1S_GL)
         If (ALLOCATED(CH4T2TM1S_GL)) DEALLOCATE (CH4T2TM1S_GL)
         If (ALLOCATED(CH41TM1S_GL)) DEALLOCATE (CH41TM1S_GL)
         If (ALLOCATED(SO4T2TM1S_GL)) DEALLOCATE (SO4T2TM1S_GL)
         If (ALLOCATED(SIT2TM1S_GL)) DEALLOCATE (SIT2TM1S_GL)
         If (ALLOCATED(BENSTRTM1S_GL)) DEALLOCATE (BENSTRTM1S_GL)
         If (ALLOCATED(SODTM1S_GL)) DEALLOCATE (SODTM1S_GL)
    !
    !: save other sediment variables needed for restarts
         If (ALLOCATED(PO41TM1S_GL)) DEALLOCATE (PO41TM1S_GL)
         If (ALLOCATED(NH41TM1S_GL)) DEALLOCATE (NH41TM1S_GL)
         If (ALLOCATED(NO31TM1S_GL)) DEALLOCATE (NO31TM1S_GL)
         If (ALLOCATED(HS1TM1S_GL)) DEALLOCATE (HS1TM1S_GL)
         If (ALLOCATED(SI1TM1S_GL)) DEALLOCATE (SI1TM1S_GL)
    !
    !
    !global arrays used for output to history and station outputs
         If (ALLOCATED(JPOC_GL)) DEALLOCATE (JPOC_GL)
         If (ALLOCATED(JPON_GL)) DEALLOCATE (JPON_GL)
         If (ALLOCATED(JPOP_GL)) DEALLOCATE (JPOP_GL)
         If (ALLOCATED(JPOS_GL)) DEALLOCATE (JPOS_GL)
    !IF(ALLOCATED(O20_GL))        DEALLOCATE(O20_GL)        !given by C2_GL(:,:,27)
    !IF(ALLOCATED(D_GL))        DEALLOCATE(D_GL)        !already defined
    !IF(ALLOCATED(T_GL))        DEALLOCATE(T_GL)        !already defined
    !IF(ALLOCATED(NH40_GL))    DEALLOCATE(NH40_GL)        !given by C2_GL
    !IF(ALLOCATED(NO30_GL))    DEALLOCATE(NO30_GL)        !given by C2_GL
    !IF(ALLOCATED(SI0_GL))        DEALLOCATE(SI0_GL)         !given by C2_GL
         If (ALLOCATED(CH40_GL)) DEALLOCATE (CH40_GL)
    !IF(ALLOCATED(SAL_GL))        DEALLOCATE(SAL_GL)        !given by C2_GL
    !IF(ALLOCATED(SOD_GL))        DEALLOCATE(SOD_GL)        !given by SODTM1S_GL
         If (ALLOCATED(JNH4_GL)) DEALLOCATE (JNH4_GL)!
         If (ALLOCATED(JNO3_GL)) DEALLOCATE (JNO3_GL)
         If (ALLOCATED(BENDEN_GL)) DEALLOCATE (BENDEN_GL)
         If (ALLOCATED(JCH4_GL)) DEALLOCATE (JCH4_GL)
         If (ALLOCATED(JCH4G_GL)) DEALLOCATE (JCH4G_GL)
         If (ALLOCATED(JHS_GL)) DEALLOCATE (JHS_GL)
         If (ALLOCATED(JPO4_GL)) DEALLOCATE (JPO4_GL)
         If (ALLOCATED(JSI_GL)) DEALLOCATE (JSI_GL)
         If (ALLOCATED(NH41_GL)) DEALLOCATE (NH41_GL)
         If (ALLOCATED(NH42_GL)) DEALLOCATE (NH42_GL)
    !
         If (ALLOCATED(NO31_GL)) DEALLOCATE (NO31_GL)
         If (ALLOCATED(NO32_GL)) DEALLOCATE (NO32_GL)
         If (ALLOCATED(PO41_GL)) DEALLOCATE (PO41_GL)
         If (ALLOCATED(PO42_GL)) DEALLOCATE (PO42_GL)
         If (ALLOCATED(SI1_GL)) DEALLOCATE (SI1_GL)
         If (ALLOCATED(SI2_GL)) DEALLOCATE (SI2_GL)
    !IF(ALLOCATED(CH41_GL))    DEALLOCATE(CH41_GL)    !given by CH41TM1S_GL
         If (ALLOCATED(CH42_GL)) DEALLOCATE (CH42_GL)
         If (ALLOCATED(HS1_GL)) DEALLOCATE (HS1_GL)
         If (ALLOCATED(HS2_GL)) DEALLOCATE (HS2_GL)
    !
    !IF(ALLOCATED(POC1_GL))    DEALLOCATE(POC1_GL)
    !IF(ALLOCATED(POC2_GL))    DEALLOCATE(POC2_GL)
    !IF(ALLOCATED(POC3_GL))    DEALLOCATE(POC3_GL)
    !
    !IF(ALLOCATED(PON1_GL))    DEALLOCATE(PON1_GL)    !given by CPON_GL
    !IF(ALLOCATED(PON2_GL))    DEALLOCATE(PON2_GL)
    !IF(ALLOCATED(PON3_GL))    DEALLOCATE(PON3_GL)
    !
    !IF(ALLOCATED(POP1_GL))    DEALLOCATE(POP1_GL)    !given by CPOP_GL
    !IF(ALLOCATED(POP2_GL))    DEALLOCATE(POP2_GL)
    !IF(ALLOCATED(POP3_GL))    DEALLOCATE(POP3_GL)
    !
    !IF(ALLOCATED(PSISED_GL))  DEALLOCATE(PSISED_GL)    !given by CPOS_GL
         If (ALLOCATED(HSED1_GL)) DEALLOCATE (HSED1_GL)
    !IF(ALLOCATED(BENSTR_GL))    DEALLOCATE(BENSTR_GL)    !given by BENSTRTM1S_GL
    !
    !used as local calculations
         If (ALLOCATED(MTVEL)) DEALLOCATE (MTVEL)!Moved to mod_sed.F
    !
         If (ALLOCATED(CPOC)) DEALLOCATE (CPOC)
         If (ALLOCATED(CPON)) DEALLOCATE (CPON)
         If (ALLOCATED(CPOP)) DEALLOCATE (CPOP)
         If (ALLOCATED(CPOS)) DEALLOCATE (CPOS)
    !
         If (ALLOCATED(CPO4)) DEALLOCATE (CPO4)
         If (ALLOCATED(CNO3)) DEALLOCATE (CNO3)
         If (ALLOCATED(CNH4)) DEALLOCATE (CNH4)
    !
         If (ALLOCATED(CCH4)) DEALLOCATE (CCH4)
         If (ALLOCATED(CSO4)) DEALLOCATE (CSO4)
         If (ALLOCATED(CHS)) DEALLOCATE (CHS)
         If (ALLOCATED(CSI)) DEALLOCATE (CSI)
    !
         If (ALLOCATED(DIAGC)) DEALLOCATE (DIAGC)
         If (ALLOCATED(DIAGN)) DEALLOCATE (DIAGN)
         If (ALLOCATED(DIAGP)) DEALLOCATE (DIAGP)
         If (ALLOCATED(DIAGS)) DEALLOCATE (DIAGS)
    !
         If (ALLOCATED(JPOC)) DEALLOCATE (JPOC)
         If (ALLOCATED(JPON)) DEALLOCATE (JPON)
         If (ALLOCATED(JPOP)) DEALLOCATE (JPOP)
         If (ALLOCATED(JPOS)) DEALLOCATE (JPOS)
    !
         If (ALLOCATED(JPOCaccum)) DEALLOCATE (JPOCaccum)
         If (ALLOCATED(JPONaccum)) DEALLOCATE (JPONaccum)
         If (ALLOCATED(JPOPaccum)) DEALLOCATE (JPOPaccum)
         If (ALLOCATED(JPOSaccum)) DEALLOCATE (JPOSaccum)
    !
      End Subroutine SED_DEALLOC
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   S E D _ R E A D                   **
  !********************************************************************************
      Subroutine SED_READ
    !
         Use MOD_SED_DF_EXCHANGE_VARS, Only: POC1TM1S_SED_DF, &
        & POC2TM1S_SED_DF, M1_SED_DF, M2_SED_DF
    !
#  if defined (MULTIPROCESSOR)
         Use MOD_PAR, Only: NGID, NHN, HN_LST
#  endif
    !
         Use MOD_WQMINIT, Only: UNI_ICON_IN_SED_VAR
    !
         Implicit None
         Save
    !
    !***** Variable declarations
    !
         Character (Len=24) :: FRNAME (14)
    !
         Integer :: I, J, JG, JT
         Real (SP), Allocatable :: RTMP11 (:), RTMP12 (:), RTMP13 (:), &
        & RTMP14 (:), RTMP15 (:), RTMP16 (:), RTMP17 (:)
    !
    !***** Data declarations
!
         Data FRNAME / 'Group 1 algal phosphorus', 'Group 2 algal phosp&
        &horus', 'Group 3 algal phosphorus', 'Detrital org phosphorus '&
        & , 'Group 1 algal nitrogen  ', 'Group 2 algal nitrogen  ', 'Gr&
        &oup 3 algal nitrogen  ', 'Detrital org nitrogen   ', 'Group 1 &
        &algal carbon    ', 'Group 2 algal carbon    ', 'Group 3 algal &
        &carbon    ', 'Benthic algal carbon    ', 'Benthic algal nitrog&
        &en  ', 'Benthic algal phosphorus' /!14   Benthic algae P
!
         If (BENTHIC_OUTPUT) Then
            If (MSR) Then
               Write (*,*) 'BFOFN=', BFOFN
               Open (BFO, File=BFOFN)
            End If
         End If
    !
    !********************************************************************************
    !**                                  Inputs                                    **
    !********************************************************************************
    !
    !
999      Format (:///8 X, F8.0, I8, I8, I8, I8, I8)! added ,I8 for NDTSED
9999     Format (/)
1020     Format (://(8 X, 7 F8.1))
1022     Format (://(8 X, 3 F8.1))
1040     Format (://(8 X, 6 F8.1))
1060     Format (/ / 8 X, 2 A8, F8.0)
1070     Format (://8 X, I8, I8, F8.0)
1080     Format (://8 X, 9 F8.0)
1090     Format (/ / (8 X, A72))
1099     Format (://(8 X, I8, 3 F8.1))
    !
         Read (BFI, 9999, Err=10100)!move to third line with format '(/)',
    !READ empty moves to second line, / moves another line
    !
         Read (BFI, 999, Err=10100) HSEDALL, SSTATEG3, SSTATEIC, &
        & SSTATEAVG, NDTSED, QUASISS ! depth of sediments, steady state flag for G3, number of WCtimesteps per sedim timestep (DTsed=NDTSED*DTwc) -added by LB
    !
         If (MSR) WRITE (*,*), 'HSEDALL=', HSEDALL, 'SSTATEG3=', &
        & SSTATEG3, 'SSTATEIC=', SSTATEIC, 'SSTATEAVG=', SSTATEAVG, 'QU&
        &ASISS=', QUASISS
    !
    !
         Read (BFI, 1070, Err=10100) ITVWCLF, NTVWCLF, WCLFDAYS !flag for focing with time varying water column information, length in days
    !in the file (all files sholuld have same length)
    !
         If (MSR) WRITE (*,*) 'ITVWCLF=', ITVWCLF !==1 to have overyling water forcing input files
         If (MSR) WRITE (*,*) 'NTVWCLF=', NTVWCLF !number of files
         If (MSR) WRITE (*,*) 'WCLFDAYS=', WCLFDAYS !lenght for each file (days)
    !
         If (ITVWCLF == 1 .And. WCLFDAYS <= 0.0) Then
            Write (*,*) 'Error, WCLFDAYS should be greater than zero'
            Go To 10100
         End If
         If (ITVWCLF == 1 .And. NTVWCLF < 1) Then
            Write (*,*) 'Error, NTVWCLF should be >=1'
            Go To 10100
         End If
    !
         Read (BFI, 1080, Err=10100) DIFFT !diffusion rate at sediment water interface (m^2/sec)
         !DIFFT = 0.0001 * DIFFT !convert from cm^2/sec to m^2/sec !
		 !Canceled conversion,so input value should be prepared in m^2/sec
    !
         If (MSR) WRITE (*,*) 'DIFFT=', DIFFT
    !
         Read (BFI, 1080, Err=10100) SALTSW, SALTND !freshwater salt water formulation for SOD  and nitrification/denitrification
         If (MSR) WRITE (*,*) 'SALTSW=', SALTSW, 'SALTND=', SALTND
    !
         Read (BFI, 1080, Err=10100) FRPALG1 !Fractions of P in water column algae group 1 for G1,G2,G3
         Read (BFI, 1080, Err=10100) FRPALG2 !Fractions of P in water column algae group 2 for G1,G2,G3
         Read (BFI, 1080, Err=10100) FRPALG3 !Fractions of P in water column algae group 3 for G1,G2,G3
         If (MSR) WRITE (*,*) 'FRPALG1=', FRPALG1
         If (MSR) WRITE (*,*) 'FRPALG2=', FRPALG2
         If (MSR) WRITE (*,*) 'FRPALG3=', FRPALG3
    !
         Read (BFI, 1080, Err=10100) FRNALG1 !Fractions of N in water column algae group 1 for G1,G2,G3
         Read (BFI, 1080, Err=10100) FRNALG2 !Fractions of N in water column algae group 2 for G1,G2,G3
         Read (BFI, 1080, Err=10100) FRNALG3 !Fractions of N in water column algae group 3 for G1,G2,G3
         If (MSR) WRITE (*,*) 'FRNALG1=', FRNALG1
         If (MSR) WRITE (*,*) 'FRNALG2=', FRNALG2
         If (MSR) WRITE (*,*) 'FRNALG3=', FRNALG3
    !
         Read (BFI, 1080, Err=10100) FRCALG1 !Fractions of C in water column algae group 1 for G1, G2, G3
         Read (BFI, 1080, Err=10100) FRCALG2 !Fractions of C in water column algae group 2 for G1, G2, G3
         Read (BFI, 1080, Err=10100) FRCALG3 !Fractions of C in water column algae group 3 for G1, G2, G3
         If (MSR) WRITE (*,*) 'FRCALG1=', FRCALG1
         If (MSR) WRITE (*,*) 'FRCALG2=', FRCALG2
         If (MSR) WRITE (*,*) 'FRCALG3=', FRCALG3
    !
         Read (BFI, 1080, Err=10100) (KPDIAG(JG), DPTHTA(JG), JG=1, 3)
    !
         KPOP1 = KPDIAG (1)
         KPOP2 = KPDIAG (2)
         KPOP3 = KPDIAG (3)
    !
         THTAPOP1 = DPTHTA (1)
         THTAPOP2 = DPTHTA (2)
         THTAPOP3 = DPTHTA (3)
    !
         If (MSR) WRITE (*,*) 'KPOP1=', KPOP1
         If (MSR) WRITE (*,*) 'KPOP2=', KPOP2
         If (MSR) WRITE (*,*) 'KPOP3=', KPOP3
    !
         If (MSR) WRITE (*,*) 'THTAPOP1=', THTAPOP1
         If (MSR) WRITE (*,*) 'THTAPOP2=', THTAPOP2
         If (MSR) WRITE (*,*) 'THTAPOP3=', THTAPOP3
    !
         Read (BFI, 1080, Err=10100) (KNDIAG(JG), DNTHTA(JG), JG=1, 3)
    !
    !Added the following and removed the equivalence declaration in wqm_modules.F
         KPON1 = KNDIAG (1)
         KPON2 = KNDIAG (2)
         KPON3 = KNDIAG (3)
         THTAPON1 = DNTHTA (1)
         THTAPON2 = DNTHTA (2)
         THTAPON3 = DNTHTA (3)
    !
         If (MSR) WRITE (*,*) 'KPON1=', KPON1
         If (MSR) WRITE (*,*) 'KPON2=', KPON2
         If (MSR) WRITE (*,*) 'KPON3=', KPON3
    !
         If (MSR) WRITE (*,*) 'THTAPON1=', THTAPON1
         If (MSR) WRITE (*,*) 'THTAPON2=', THTAPON2
         If (MSR) WRITE (*,*) 'THTAPON3=', THTAPON3
    !
         Read (BFI, 1080, Err=10100) (KCDIAG(JG), DCTHTA(JG), JG=1, 3)
    !
    !Added the following and removed the equivalence declaration in wqm_modules.F
         KPOC1 = KCDIAG (1)
         KPOC2 = KCDIAG (2)
         KPOC3 = KCDIAG (3)
         THTAPOC1 = DCTHTA (1)
         THTAPOC2 = DCTHTA (2)
         THTAPOC3 = DCTHTA (3)
    !
         If (MSR) WRITE (*,*) 'KPOC1=', KPOC1
         If (MSR) WRITE (*,*) 'KPOC2=', KPOC2
         If (MSR) WRITE (*,*) 'KPOC3=', KPOC3
    !
         If (MSR) WRITE (*,*) 'THTAPOC1=', THTAPOC1
         If (MSR) WRITE (*,*) 'THTAPOC2=', THTAPOC2
         If (MSR) WRITE (*,*) 'THTAPOC3=', THTAPOC3
    !
         Read (BFI, 1080, Err=10100) KSI, THTASI, THTASISAT
         If (MSR) WRITE (*,*) 'KSI=', KSI
         If (MSR) WRITE (*,*) 'THTASI=', THTASI
         If (MSR) WRITE (*,*) 'THTASISAT=', THTASISAT
    !
    !
         Read (BFI, 1080, Err=10100) M1, M2, THTADP, THTADD
         If (MSR) WRITE (*,*) 'M1=', M1
         If (MSR) WRITE (*,*) 'M2=', M2
         If (MSR) WRITE (*,*) 'THTADP=', THTADP
         If (MSR) WRITE (*,*) 'THTADD=', THTADD
    !
    !
         Read (BFI, 1080, Err=10100) KAPPNH4F, KAPPNH4S, PIENH4, &
        & THTANH4, KMNH4, KMNH4O2, THTAKMNH4
    !
         If (MSR) WRITE (*,*) 'KAPPNH4F=', KAPPNH4F
         If (MSR) WRITE (*,*) 'KAPPNH4S=', KAPPNH4S
         If (MSR) WRITE (*,*) 'PIENH4=', PIENH4
         If (MSR) WRITE (*,*) 'THTANH4=', THTANH4
         If (MSR) WRITE (*,*) 'KMNH4=', KMNH4
         If (MSR) WRITE (*,*) 'KMNH4O2=', KMNH4O2
         If (MSR) WRITE (*,*) 'THTAKMNH4=', THTAKMNH4
    !
         Read (BFI, 1080, Err=10100) KAPPNO3F, KAPPNO3S, K2NO3, THTANO3
    !
         If (MSR) WRITE (*,*) 'KAPPNO3F=', KAPPNO3F
         If (MSR) WRITE (*,*) 'KAPPNO3S=', KAPPNO3S
         If (MSR) WRITE (*,*) 'K2NO3=', K2NO3
         If (MSR) WRITE (*,*) 'THTANO3=', THTANO3
    !
         Read (BFI, 1080, Err=10100) KAPP1HSD, KAPP1HSP, PIE1HS, &
        & PIE2HS, THTAH2S, KMHSO2
    !
         If (MSR) WRITE (*,*) 'KAPP1HSD=', KAPP1HSD
         If (MSR) WRITE (*,*) 'KAPP1HSP=', KAPP1HSP
         If (MSR) WRITE (*,*) 'PIE1HS=', PIE1HS
         If (MSR) WRITE (*,*) 'PIE2HS=', PIE2HS
         If (MSR) WRITE (*,*) 'THTAH2S=', THTAH2S
         If (MSR) WRITE (*,*) 'KMHSO2=', KMHSO2
    !
         Read (BFI, 1080, Err=10100) CSISATT20, DPIE1SI, PIE2SI, KMPSI !KMPSI is gSi/m^3 on input!!!!
         KMPSI = KMPSI * 1000.d0 !convert from gSi/m^3 to mgSi/m^3
         If (MSR) WRITE (*,*) 'CSISATT20=', CSISATT20
         If (MSR) WRITE (*,*) 'DPIE1SI=', DPIE1SI
         If (MSR) WRITE (*,*) 'PIE2SI=', PIE2SI
         If (MSR) WRITE (*,*) 'KMPSI=', KMPSI
    !
         Read (BFI, 1080, Err=10100) O2CRITSI, JSIDETR
         If (MSR) WRITE (*,*) 'O2CRITSI=', O2CRITSI
         If (MSR) WRITE (*,*) 'JSIDETR=', JSIDETR
    !
         Read (BFI, 1080, Err=10100) DPIE1PO4F, DPIE1PO4S, PIE2PO4, &
        & O2CRITPO4, KMO2DP
         If (MSR) WRITE (*,*) 'DPIE1PO4F=', DPIE1PO4F
         If (MSR) WRITE (*,*) 'DPIE1PO4S=', DPIE1PO4S
         If (MSR) WRITE (*,*) 'PIE2PO4=', PIE2PO4
         If (MSR) WRITE (*,*) 'O2CRITPO4=', O2CRITPO4
         If (MSR) WRITE (*,*) 'KMO2DP=', KMO2DP
    !
         Read (BFI, 1080, Err=10100) TEMPBEN, KBENSTR, KLBNTH, DPMIN
    !
         If (MSR) WRITE (*,*) 'TEMPBEN=', TEMPBEN
         If (MSR) WRITE (*,*) 'KBENSTR=', KBENSTR
         If (MSR) WRITE (*,*) 'KLBNTH=', KLBNTH
         If (MSR) WRITE (*,*) 'DPMIN=', DPMIN
    !
    !
         Read (BFI, 1080, Err=10100) KAPPCH4, THTACH4, KMCH4O2, KMSO4
    !
         If (MSR) WRITE (*,*) 'KAPPCH4=', KAPPCH4
         If (MSR) WRITE (*,*) 'THTACH4=', THTACH4
         If (MSR) WRITE (*,*) 'KMCH4O2=', KMCH4O2
         If (MSR) WRITE (*,*) 'KMSO4=', KMSO4
    !
    !read information of water colum time-varying conditions (forcing with water column input file)
    !
         If (ITVWCLF >= 1) Then
            Read (BFI, 1090, Err=10100) (WCLFN(J), J=1, NTVWCLF)
       !
            If (MSR) Then
               Do J = 1, NTVWCLF
                  Write (*,*) 'WCLFN(', J, ')=', WCLFN (J)
               End Do
            End If
         End If
    !
         M1_SED_DF = M1
         M2_SED_DF = M2 !will be used by mod_df
    !
    ! net settling rates
    !
         Read (BFI, 1060, Err=10100) SPVARS, PRINTS
         If (MSR) WRITE (*,*) 'SPVARS=', SPVARS
         If (MSR) WRITE (*,*) 'PRINTS=', PRINTS
    !
         If (SPVARS == 'CONSTANT') Then
       !
       !
            Read (BFI, 1020, Err=10100) WSSBNET (1), WSLBNET (1), &
           & WSRBNET (1), WS1BNET (1), WS2BNET (1), WS3BNET (1), &
           & WSUBNET (1)!settling rate of particulate biogenic silica (m/d)
       !
            If (MSR) WRITE (*,*) 'WSSBNET=', WSSBNET (1)
            If (MSR) WRITE (*,*) 'WSLBNET=', WSLBNET (1)
            If (MSR) WRITE (*,*) 'WSRBNET=', WSRBNET (1)
            If (MSR) WRITE (*,*) 'WS1BNET=', WS1BNET (1)
            If (MSR) WRITE (*,*) 'WS2BNET=', WS2BNET (1)
            If (MSR) WRITE (*,*) 'WS3BNET=', WS3BNET (1)
            If (MSR) WRITE (*,*) 'WSUBNET=', WSUBNET (1)
       !
            Do I = 2, MTLOC
               WSSBNET (I) = WSSBNET (1)
               WSLBNET (I) = WSLBNET (1)
               WSRBNET (I) = WSRBNET (1)
               WS1BNET (I) = WS1BNET (1)
               WS2BNET (I) = WS2BNET (1)
               WS3BNET (I) = WS3BNET (1)
               WSUBNET (I) = WSUBNET (1)
            End Do
         Else
            Allocate (RTMP11(MGL))
            RTMP11 = 0.0
            Allocate (RTMP12(MGL))
            RTMP12 = 0.0
            Allocate (RTMP13(MGL))
            RTMP13 = 0.0
            Allocate (RTMP14(MGL))
            RTMP14 = 0.0
            Allocate (RTMP15(MGL))
            RTMP15 = 0.0
            Allocate (RTMP16(MGL))
            RTMP16 = 0.0
            Allocate (RTMP17(MGL))
            RTMP17 = 0.0
            Do I = 1, MGL
               If (MSR) WRITE (*,*) 'I=', I
               Read (BFI, 1020, Err=10100) RTMP11 (I), RTMP12 (I), &
              & RTMP13 (I), RTMP14 (I), RTMP15 (I), RTMP16 (I), RTMP17 &
              & (I)
            End Do
            If (SERIAL) Then
               WSSBNET = RTMP11
               WSLBNET = RTMP12
               WSRBNET = RTMP13
               WS1BNET = RTMP14
               WS2BNET = RTMP15
               WS3BNET = RTMP16
               WSUBNET = RTMP17
          !
            End If
#    if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do I = 1, MLOC
                  WSSBNET (I) = RTMP11 (NGID(I))
                  WSLBNET (I) = RTMP12 (NGID(I))
                  WSRBNET (I) = RTMP13 (NGID(I))
                  WS1BNET (I) = RTMP14 (NGID(I))
                  WS2BNET (I) = RTMP15 (NGID(I))
                  WS3BNET (I) = RTMP16 (NGID(I))
                  WSUBNET (I) = RTMP17 (NGID(I))
               End Do
               Do I = 1, NHN
                  WSSBNET (I+MLOC) = RTMP11 (HN_LST(I))
                  WSLBNET (I+MLOC) = RTMP12 (HN_LST(I))
                  WSRBNET (I+MLOC) = RTMP13 (HN_LST(I))
                  WS1BNET (I+MLOC) = RTMP14 (HN_LST(I))
                  WS2BNET (I+MLOC) = RTMP15 (HN_LST(I))
                  WS3BNET (I+MLOC) = RTMP16 (HN_LST(I))
                  WSUBNET (I+MLOC) = RTMP17 (HN_LST(I))
               End Do
            End If
#    endif
       !
            If (MSR) WRITE (*,*) 'WSSBNET=', WSSBNET (1)
            If (MSR) WRITE (*,*) 'WSLBNET=', WSLBNET (1)
            If (MSR) WRITE (*,*) 'WSRBNET=', WSRBNET (1)
            If (MSR) WRITE (*,*) 'WS1BNET=', WS1BNET (1)
            If (MSR) WRITE (*,*) 'WS2BNET=', WS2BNET (1)
            If (MSR) WRITE (*,*) 'WS3BNET=', WS3BNET (1)
            If (MSR) WRITE (*,*) 'WSUBNET=', WSUBNET (1)
       !
            Deallocate (RTMP11, RTMP12, RTMP13, RTMP14, RTMP15, RTMP16, &
           & RTMP17)
         End If
    !
    ! burial and mixing rates
    !
         Read (BFI, 1080, Err=10100) POC1R !mgC/gSediment
         If (MSR) WRITE (*,*) 'POC1R=', POC1R
    !
         Read (BFI, 1060, Err=10100) SPVARB, PRINTB
         If (MSR) WRITE (*,*) 'SPVARB=', SPVARB
         If (MSR) WRITE (*,*) 'PRINTB=', PRINTB
    !
    !
         If (SPVARB == 'CONSTANT') Then
            Read (BFI, 1022, Err=10100) VSED (1), VPMIX (1), VDMIX (1)
            Do I = 2, MTLOC
               VSED (I) = VSED (1)
               VPMIX (I) = VPMIX (1)
               VDMIX (I) = VDMIX (1)
            End Do
       !
            If (MSR) WRITE (*,*) 'VSED(1)=', VSED (1)
            If (MSR) WRITE (*,*) 'VPMIX(1)=', VPMIX (1)
            If (MSR) WRITE (*,*) 'VDMIX(1)=', VDMIX (1)
       !
       !
         Else
            Allocate (RTMP11(MGL))
            RTMP11 = 0.0
            Allocate (RTMP12(MGL))
            RTMP12 = 0.0
            Allocate (RTMP13(MGL))
            RTMP13 = 0.0
            Do I = 1, MGL
               Read (BFI, 1022, Err=10100) RTMP11 (I), RTMP12 (I), &
              & RTMP13 (I)
            End Do
            If (SERIAL) Then
               VSED = RTMP11
               VPMIX = RTMP12
               VDMIX = RTMP13
            End If
#    if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do I = 1, MLOC
                  VSED (I) = RTMP11 (NGID(I))
                  VPMIX (I) = RTMP12 (NGID(I))
                  VDMIX (I) = RTMP13 (NGID(I))
               End Do
               Do I = 1, NHN
                  VSED (I+MLOC) = RTMP11 (HN_LST(I))
                  VPMIX (I+MLOC) = RTMP12 (HN_LST(I))
                  VDMIX (I+MLOC) = RTMP13 (HN_LST(I))
               End Do
            End If
#    endif
       !
            If (MSR) WRITE (*,*) 'VSED(1)=', VSED (1)
            If (MSR) WRITE (*,*) 'VPMIX(1)=', VPMIX (1)
            If (MSR) WRITE (*,*) 'VDMIX(1)=', VDMIX (1)
       !
            Deallocate (RTMP11, RTMP12, RTMP13)
         End If !end of SPVARB switch
    !
    ! splits of refractory water column into G2, G3 sediments
    !
         Read (BFI, 1060, Err=10100) SPVARLR, PRINTLR
    !
         If (MSR) WRITE (*,*) 'SPVARLR=', SPVARLR
         If (MSR) WRITE (*,*) 'PRINTLR=', PRINTLR
    !
         If (SPVARLR == 'CONSTANT') Then
            Read (BFI, 1040, Err=10100) FRPOP (1, 2), FRPOP (1, 3), &
           & FRPON (1, 2), FRPON (1, 3), FRPOC (1, 2), FRPOC (1, 3)
       !
            If (MSR) WRITE (*,*) ' FRPOP(1,2)=', FRPOP (1, 2)
            If (MSR) WRITE (*,*) ' FRPOP(1,3)=', FRPOP (1, 3)
       !
            If (MSR) WRITE (*,*) ' FRPON(1,2)=', FRPON (1, 2)
            If (MSR) WRITE (*,*) ' FRPON(1,3)=', FRPON (1, 3)
       !
            If (MSR) WRITE (*,*) ' FRPOC(1,2)=', FRPOC (1, 2)
            If (MSR) WRITE (*,*) ' FRPOC(1,3)=', FRPOC (1, 3)
       !
            Do I = 2, MTLOC
               FRPOP (I, 2) = FRPOP (1, 2)
               FRPOP (I, 3) = FRPOP (1, 3)
               FRPON (I, 2) = FRPON (1, 2)
               FRPON (I, 3) = FRPON (1, 3)
               FRPOC (I, 2) = FRPOC (1, 2)
               FRPOC (I, 3) = FRPOC (1, 3)
            End Do
         Else
            Allocate (RTMP11(MGL))
            RTMP11 = 0.0
            Allocate (RTMP12(MGL))
            RTMP12 = 0.0
            Allocate (RTMP13(MGL))
            RTMP13 = 0.0
            Allocate (RTMP14(MGL))
            RTMP14 = 0.0
            Allocate (RTMP15(MGL))
            RTMP15 = 0.0
            Allocate (RTMP16(MGL))
            RTMP16 = 0.0
            Do I = 1, MGL
               Read (BFI, 1040, Err=10100) RTMP11 (I), RTMP12 (I), &
              & RTMP13 (I), RTMP14 (I), RTMP15 (I), RTMP16 (I)
            End Do
       !
            If (SERIAL) Then
               FRPOP (:, 2) = RTMP11
               FRPOP (:, 3) = RTMP12
               FRPON (:, 2) = RTMP13
               FRPON (:, 3) = RTMP14
               FRPOC (:, 2) = RTMP15
               FRPOC (:, 3) = RTMP16
            End If
#    if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do I = 1, MLOC
                  FRPOP (I, 2) = RTMP11 (NGID(I))
                  FRPOP (I, 3) = RTMP12 (NGID(I))
                  FRPON (I, 2) = RTMP13 (NGID(I))
                  FRPON (I, 3) = RTMP14 (NGID(I))
                  FRPOC (I, 2) = RTMP15 (NGID(I))
                  FRPOC (I, 3) = RTMP16 (NGID(I))
               End Do
               Do I = 1, NHN
                  FRPOP (I+MLOC, 2) = RTMP11 (HN_LST(I))
                  FRPOP (I+MLOC, 3) = RTMP12 (HN_LST(I))
                  FRPON (I+MLOC, 2) = RTMP13 (HN_LST(I))
                  FRPON (I+MLOC, 3) = RTMP14 (HN_LST(I))
                  FRPOC (I+MLOC, 2) = RTMP15 (HN_LST(I))
                  FRPOC (I+MLOC, 3) = RTMP16 (HN_LST(I))
               End Do
            End If
#    endif
       !
            If (MSR) WRITE (*,*) ' FRPOP(1,2)=', FRPOP (1, 2)
            If (MSR) WRITE (*,*) ' FRPOP(1,3)=', FRPOP (1, 3)
            If (MSR) WRITE (*,*) ' FRPON(1,2)=', FRPON (1, 2)
            If (MSR) WRITE (*,*) ' FRPON(1,3)=', FRPON (1, 3)
            If (MSR) WRITE (*,*) ' FRPOC(1,2)=', FRPOC (1, 2)
            If (MSR) WRITE (*,*) ' FRPOC(1,3)=', FRPOC (1, 3)
       !
       !
            Deallocate (RTMP11, RTMP12, RTMP13, RTMP14, RTMP15, RTMP16)
         End If !end of SPVARLR switch
    !
    !
         Read (BFI, 1099, Err=10100) SWITCH_LAB, FRACL1, FRACL2, FRACL3
         If (MSR) WRITE (*,*) 'SWITCH_LAB=', SWITCH_LAB
         If (MSR) WRITE (*,*) 'FRACTION Lab1=', FRACL1
         If (MSR) WRITE (*,*) 'FRACTION Lab2=', FRACL2
         If (MSR) WRITE (*,*) 'FRACTION Lab3=', FRACL3
    !
    !***** Define logical variables
    !
         STEADY_STATE_SED_G3 = SSTATEG3 == 1
         STEADY_STATE_SED_IC = ( .Not. UNI_ICON_IN_SED_VAR .And. &
        & SSTATEIC == 1)!: if reading initial conds from ini file, automatically assume not steady state for ini conds
         STEADY_STATE_SED_TS = QUASISS == 1
         STEADY_STATE_SED_AVG = SSTATEAVG == 1
         AGGREGATE_POM_FLUX = SWITCH_LAB == 1
    !
    !WRITE(*,*)'STEADY_STATE_SED_IC=,',STEADY_STATE_SED_IC,' UNI_ICON_IN_SED_VAR=',UNI_ICON_IN_SED_VAR,' SSTATEIC=',SSTATEIC
    !********************************************************************************
    !**                                 Outputs                                    **
    !********************************************************************************
    !
         If (BENTHIC_OUTPUT) Then
            If (MSR) Then
               Write (BFO, 2000)
2000           Format (/ // 34 X, 'Sediment-water column linkages and s&
              &ediment ', 'depths and volumes' /)!and add one more line
               Write (BFO, 2020) HSEDALL
2020           Format (/ ' ACTIVE LAYER DEPTH ', F8.3, ' CM')!print on 7th line
          !
               If (STEADY_STATE_SED_G3) Then
                  Write (BFO, 2022)
2022              Format (/ ' STEADY-STATE VALUES OF G3 COMPUTED' /)!print on 9th line and add one line
               Else
                  Write (BFO, 2025)
2025              Format (/ ' NO STEADY-STATE VALUES OF G3 COMPUTED' /)
               End If
               Write (BFO, 2030)
2030           Format (/ // / 33 X, 'S E D I M E N T   I N I T I A L   &
              &C O N D I T ', 'I O N S' /)
          !
               Write (BFO, 2040) SSNAME (1)!temperature
2040           Format (/ / 25 X, 'Sediment initial conditions for ', &
              & A20 /)
          !
               Write (BFO, 2050) (CTEMP(I), I=1, MLOC)!: need to use global array
2050           Format (13 X, 3(7 X, 1 PE11.4))
          !
               Write (BFO, 2060) SSNAME (2)!POP
2060           Format (/ / 25 X, 'Sediment initial conditions for ', &
              & A20 / 37 X, 'G1', 22 X, 'G2', 22 X, 'G3' /)
               Write (BFO, 2070) ((CPOP(I, JG), JG=1, 3), I=1, MLOC)!:Need to use global array
2070           Format (18 X, 3(2 X, 1 PE11.4))
          !
               Write (BFO, 2060) SSNAME (3)!PON
               Write (BFO, 2070) ((CPON(I, JG), JG=1, 3), I=1, MLOC)!:Need to use global array
          !
               Write (BFO, 2060) SSNAME (4)!POC
               Write (BFO, 2070) ((CPOC(I, JG), JG=1, 3), I=1, MLOC)!:Need to use global array
          !
               Write (BFO, 2060) SSNAME (5)!PBS    particulate biogenic silica
               Write (BFO, 2070) (CPOS(I), I=1, MLOC)!:Need to use global array
          !
               Write (BFO, 2040) SSNAME (6)!total PO4 in layer 2
               Write (BFO, 2050) (PO4T2TM1S(I), I=1, MLOC)!:Need to use global array
          !
               Write (BFO, 2040) SSNAME (7)!total NH4 in layer 2
               Write (BFO, 2050) (NH4T2TM1S(I), I=1, MLOC)!:Need to use global array
          !
               Write (BFO, 2040) SSNAME (8)!total NO3 in layer 2
               Write (BFO, 2050) (NO3T2TM1S(I), I=1, MLOC)!:Need to use global array
          !
               Write (BFO, 2040) SSNAME (9)!total H2S in layer 2
               Write (BFO, 2050) (HST2TM1S(I), I=1, MLOC)!:Need to use global array
          !
               Write (BFO, 2040) SSNAME (10)!total DSIL (dissolved silicate) in layer 2
               Write (BFO, 2050) (SIT2TM1S(I), I=1, MLOC)!:Need to use global array
          !
               Write (BFO, 2040) SSNAME (11)!benthic stress
               Write (BFO, 2050) (BENSTRTM1S(I), I=1, MLOC)!:Need to use global array
          !
               Write (BFO, 2040) SSNAME (13)!total CH4 in layer 2
               Write (BFO, 2050) (CH4T2TM1S(I), I=1, MLOC)!:Need to use global array
          !
               Write (BFO, 2040) SSNAME (14)!total SO4 in layer 2
               Write (BFO, 2050) (SO4T2TM1S(I), I=1, MLOC)!:Need to use global array
          !
               Write (BFO, 2080) DIFFT !m^2/s
2080           Format (/ / 30 X, 'Temperature diffusion coefficient ', &
              & E10.3, ' m**2/sec')
          !
               Write (BFO, 2090) SALTSW, SALTND
2090           Format (/ / 31 X, 'If salinity < ', F10.3, ' ppt, methan&
              &e formed', / 31 X, 'If salinity < ', F10.3, ' ppt, high &
              &nit/denit used')
               Write (BFO, 2100)
2100           Format (/ / 30 X, 'Particulate organic matter G-model sp&
              &lits' / 10 X, 'fraction of....', 5 X, 'recycled to', 5 &
              & X, 'G1', 5 X, 'G2', 5 X, 'G3')
               Write (BFO, 2110) FRNAME (1), FRPALG1
               Write (BFO, 2110) FRNAME (2), FRPALG2
               Write (BFO, 2110) FRNAME (3), FRPALG3
               Write (BFO, 2110) FRNAME (5), FRNALG1
               Write (BFO, 2110) FRNAME (6), FRNALG2
               Write (BFO, 2110) FRNAME (7), FRNALG3
               Write (BFO, 2110) FRNAME (9), FRCALG1
               Write (BFO, 2110) FRNAME (10), FRCALG2
               Write (BFO, 2110) FRNAME (11), FRCALG3
2110           Format (6 X, A24, 11 X, 3 F7.2)
          !
               Write (BFO, 2250) (I, FRPOP(I, 2), FRPOP(I, 3), FRPON(I, &
              & 2), FRPON(I, 3), FRPOC(I, 2), FRPOC(I, 3), I=1, MLOC)
2250           Format (/ / 31 X, 'G2 - G3 splits for Refractory Particu&
              &lates' / '  NODE   FRG2P   FRG3P   FRG2N   FRG3N   FRG2C&
              &   FRG3C' / (I5, 6 F8.3))
          !
               Write (BFO, 2120) (KPDIAG(JG), DPTHTA(JG), JG=1, 3), &
              & (KNDIAG(JG), DNTHTA(JG), JG=1, 3), (KCDIAG(JG), &
              & DCTHTA(JG), JG=1, 3), KSI, THTASI
2120           Format (/ / 30 X, 'Diagenesis rates (1/day) | Temp corr &
              &factor' / 30 X, 'Phosphorus' / 39 X, 'G1', E11.3, 5 X, &
              & F7.3 / 39 X, 'G2', E11.3, 5 X, F7.3 / 39 X, 'G3', &
              & E11.3, 5 X, F7.3 / 30 X, 'Nitrogen' / 39 X, 'G1', &
              & E11.3, 5 X, F7.3 / 39 X, 'G2', E11.3, 5 X, F7.3 / 39 X, &
              & 'G3', E11.3, 5 X, F7.3 / 30 X, 'Carbon' / 39 X, 'G1', &
              & E11.3, 5 X, F7.3 / 39 X, 'G2', E11.3, 5 X, F7.3 / 39 X, &
              & 'G3', E11.3, 5 X, F7.3 / 30 X, 'Silica' / 41 X, E11.3, &
              & 5 X, F7.3)
          !
               Write (BFO, 2170) M1, M2, THTADP, THTADD
2170           Format (/ / 35 X, 'Additional constants' / 30 X, 'M1....&
              &....', F8.2, ' kg/L' / 30 X, 'M2........', F8.2, ' kg/L' &
              & / 30 X, 'THTADP....', F8.3, / 30 X, 'THTADD....', F8.3)
               Write (BFO, 2180) KAPPNH4F, KAPPNH4S, PIENH4, THTANH4, &
              & KMNH4, KMNH4O2
2180           Format (30 X, 'KAPPNH4F..', F8.3, ' m/day'/30 X, 'KAPPNH4S..', F8.3, ' m/day'/30 X, 'PIENH4....', F8.3, ' L/kg'/30 &
              & X, 'THTANH4...', F8.3,/30 X, 'KMNH4.....', F8.3, ' mg n/m**3'/30 X, 'KMNH4O2...', F8.3, ' mg o2/L')
               Write (BFO, 2190) KAPPNO3F, KAPPNO3S, K2NO3, THTANO3
2190           Format (30 X, 'KAPPNO3F..', F8.3, ' m/day'/30 X, 'KAPPNO&
              &3S..', F8.3, ' m/day'/30 X, 'K2NO3.....', F8.3, ' /day'/&
              & 30 X, 'THTANO3...', F8.3)
               Write (BFO, 2200) KAPP1HSD, KAPP1HSP, PIE1HS, PIE2HS, &
              & THTAH2S, KMHSO2
2200           Format (30 X, 'KAPP1HSD....', F8.3, ' m/day'/30 X, 'KAPP&
              &1HSP....', F8.3, ' m/day'/30 X, 'PIE1HS......', F8.3, ' &
              &L/kg'/30 X, 'PIE2HS......', F8.3, ' L/kg'/30 X, 'THTAH2S&
              &...', F8.3,/30 X, 'KMHSO2....', F8.3, ' mg o2/L')
               Write (BFO, 2210) CSISATT20, DPIE1SI, PIE2SI, KMPSI, &
              & O2CRITSI, JSIDETR
2210           Format (30 X, 'CSISATT20..', F8.1, ' mg si/m**3'/30 X, '&
              &DPIE1SI....', F8.3, ' L/kg'/30 X, 'PIE2SI....', F8.3, ' &
              &L/kg'/30 X, 'KMPSI.....', E8.2, ' mg si/m**3'/30 X, 'O2C&
              &RITSI..', E8.2, ' mg O2/L'/30 X, 'JSIDETR...', E8.2, ' m&
              &g Si/m2-d')
               Write (BFO, 2220) DPIE1PO4F, DPIE1PO4S, PIE2PO4, &
              & O2CRITPO4, KMO2DP
2220           Format (30 X, 'DPIE1PO4F..', F8.3, ' l/kg'/30 X, 'DPIE1P&
              &O4S..', F8.3, ' l/kg'/30 X, 'PIE2PO4..', F8.3, ' l/kg'/3&
              & 0 X, 'O2CRITPO4.', F8.3, ' mg o2/l'/30 X, 'KMO2DP....', &
              & F8.3, ' mg o2/l')
               Write (BFO, 2230) TEMPBEN, KBENSTR, KLBNTH, DPMIN
2230           Format (30 X, 'TEMPBEN...', F8.3, ' deg c'/30 X, 'KBENST&
              &R...', F8.3, ' /day'/30 X, 'KLBNTH....', F8.3, '---'/30 &
              & X, 'DPMIN.....', F8.3, ' m2/d')
               Write (BFO, 2240) KAPPCH4, THTACH4
2240           Format (30 X, 'KAPPCH4...', F8.3, ' m/day'/30 X, 'THTACH&
              &4...', F8.3)
          !
               Write (BFO, 2130) (I, WSSBNET(I), WSLBNET(I), &
              & WSRBNET(I), WS1BNET(I), WS2BNET(I), WS3BNET(I), &
              & WSUBNET(I), I=1, MLOC)!:Need to use global array
2130           Format (/ / 3 X, 'NODE', 6 X, 'WSSNET', 4 X, 'WSLNET', 4 &
              & X, 'WSRNET', 4 X, 'WS1NET', 4 X, 'WS2NET', 4 X, 'WS3NET&
              &', 4 X, 'WSUNET' / (I7, 7 F10.3))
               Write (BFO, 2140) (I, VSED(I), I=1, MLOC)!:Need to use global array
          !
2140           Format (/ / 31 X, 'Sedimentation (burial) rates (cm/yr)' &
              & / 10 X, 8 (I5, F6.2))
               Write (BFO, 2150) (I, VPMIX(I), I=1, MLOC)!:Need to use global array
2150           Format (/ / 30 X, 'Sediment solid-phase mixing rates (m*&
              &*2/day)' / 10 X, 8 (I5, F6.2))
               Write (BFO, 2160) (I, VDMIX(I), I=1, MLOC)!:Need to use global array
2160           Format (/ / 30 X, 'Sediment dissolved-phase mixing rates&
              & (m**2/day)' / 10 X, 8 (I5, F6.2))
               Close (BFO)
            End If
         End If
    !
    !***** Convert cell heights and burial velocities to sediment units
    !
    !MTLOC    --include the halo nodes
    !MLOC    --does not include halo
    !
         Do I = 1, MLOC
            HSED (I) = HSEDALL * 0.01 !Convert from cm to m
            VSED (I) = VSED (I) * 2.73791E-5 !cm/yr to m/day = 0.01/365.0 = 2.73791E-5
         End Do
    !
    !
    !Open the first time variable water column condition file if need be
    !and the read past the header lines
         If (ITVWCLF == 1) Then
            WCLPTR = 1
            Open (WCL, File=WCLFN(WCLPTR), Status='OLD')
            Read (WCL, 1000)!read past the 3 header lines
       !
1000        Format (/ /)
1010        Format (14(F12.5, 1 X))
         End If
    !
    !CALL SED_INIT2
    !
         Return
    !
    !***** Error traps
    !
10100    Continue
    !
         If (BENTHIC_OUTPUT) Then
            Write (BFO, 3010)
            Close (BFO)
         End If
3010     Format (/ ' Read error in sediment input deck')
         Stop 'SED_READ'
    !
         Return
      End Subroutine SED_READ
  !
  !********************************************************************************
  !**                             Initializations                                **
  !********************************************************************************
  !
  !set the initial condition to sediment if it is read from ICI file
      Subroutine SED_INIT_ICI
    !
         Use MOD_WQM, Only: C2
    !
         Implicit None
         Save
         Integer :: I, JG
    !
    !if(MSR)write(*,*)'Init temp sedim = ',CTEMPI,'*******************************************'
         Do I = 1, MLOC
       !
       !initialize temperature of sediments from water column instead
            CTEMP (I) = T (I, KBM1)!take bottom layer of water column instead
       !
       !#if defined (test2a) || defined (test2b) ||defined(test2c) ||defined(test2d) ||defined (test2e)
       !            !initial condition of concentrations in sediments
       !            CTEMP(I)= CTEMPI
       !#endif
       !
            Do JG = 1, 3
               CPOP (I, JG) = CPOPI (JG)
               CPON (I, JG) = CPONI (JG)
               CPOC (I, JG) = CPOCI (JG)
            End Do
       !
            CPOS (I) = CPOSI
       !
            PO41TM1S (I) = PO41TI
            PO4T2TM1S (I) = PO4T2I
       !
            NH41TM1S (I) = NH41TI
            NH4T2TM1S (I) = NH4T2I
       !
            NO31TM1S (I) = NO31TI
            NO3T2TM1S (I) = NO3T2I
       !
            HS1TM1S (I) = HS1TI
            HST2TM1S (I) = HST2I
       !
            CH41TM1S (I) = CH41TI
            CH4T2TM1S (I) = CH4T2I
       !
       !SO41TM1S(I)   = SO41TI
            SO4T2TM1S (I) = SO4T2I
       !
            SI1TM1S (I) = SI1TI
            SIT2TM1S (I) = SIT2I
       !
            BENSTRTM1S (I) = BENSTI
            SODTM1S (I) = SODI
       !

       !
         End Do
    !
    !
      End Subroutine SED_INIT_ICI
  !
  !
      Subroutine SED_INIT2
    !
         Use MOD_SED_DF_EXCHANGE_VARS, Only: POC1TM1S_SED_DF, &
        & POC2TM1S_SED_DF, M1_SED_DF, M2_SED_DF
	!
		 Use MOD_SED_SAV_EXCHANGE_VARS, Only: NH4T2TM1S_SHARE, &
        & PO4T2TM1S_SHARE !AN
    !
#  if defined (MULTIPROCESSOR)
         Use MOD_PAR, Only: NGID, NHN, HN_LST
#  endif
    !
         Implicit None
         Save
         Integer :: I, JG, JT
    !
!
    !***** Set sediment concentrations to initial concentrations
    !
!
         Do I = 1, MLOC
            POP1TM1S (I) = CPOP (I, 1)
            POP2TM1S (I) = CPOP (I, 2)
            POP3TM1S (I) = CPOP (I, 3)
            PON1TM1S (I) = CPON (I, 1)
            PON2TM1S (I) = CPON (I, 2)
            PON3TM1S (I) = CPON (I, 3)
            POC1TM1S (I) = CPOC (I, 1)
            POC2TM1S (I) = CPOC (I, 2)
            POC3TM1S (I) = CPOC (I, 3)
            PSISEDTM1S (I) = CPOS (I)
       !
       !Added POC1TM1S_SED_DF and POC2TM1S_SED_DF
       !
            If (DFEEDER) Then
               POC1TM1S_SED_DF (I) = CPOC (I, 1)
               POC2TM1S_SED_DF (I) = CPOC (I, 2)
            End If
       !
         End Do
    !
    !***** Initialize mass balance variables (kgC, kgN, kgP)
         ISEDMC = 0.
         ISEDMN = 0.
         ISEDMP = 0.
    !
    !
         Do I = 1, MLOC
            CPO4 (I) = PO4T2TM1S (I)
            CNO3 (I) = NO3T2TM1S (I)
            CNH4 (I) = NH4T2TM1S (I)
       !
       !initial total C, N, P in sediments (kg)
            ISEDMC = ISEDMC + (CPOC(I, 1)+CPOC(I, 2)+CPOC(I, 3)) * ART1 &
           & (I) * HSED (I) / 1.E6
       !unit of ISEDMC is mgC/m^3 * m^2 * m /1000/1000
       !==> mgC/1000/1000 ==>gC/1000 ==> kgC
       !
            ISEDMN = ISEDMN + (CPON(I, 1)+CPON(I, 2)+CPON(I, &
           & 3)+CNH4(I)+CNO3(I)) * ART1 (I) * HSED (I) / 1.E6
       !unit of ISEDMN is mgN/m^3 * m^2 * m /1000/1000
       !==> mgN/1000/1000 ==>gN/1000 ==> kgN
            ISEDMP = ISEDMP + (CPOP(I, 1)+CPOP(I, 2)+CPOP(I, &
           & 3)+CPO4(I)) * ART1 (I) * HSED (I) / 1.E6
       !unit of ISEDMP is mgP/m^3 * m^2 * m /1000/1000
       !==> mgP/1000/1000 ==>gP/1000 ==> kgP

         End Do
    !
    !***** Set up reaction rates in table look-up form
    !
         Do JT = 1, 350
            TEMP_SED = REAL (JT-1) / 10. + 0.05
            TEMP20 = TEMP_SED - 20.
            TEMP202 = TEMP20 / 2.
       !
            ZHTANH4F (JT) = KAPPNH4F * THTANH4 ** TEMP202
            ZHTANH4S (JT) = KAPPNH4S * THTANH4 ** TEMP202
            ZHTA1HSD (JT) = KAPP1HSD * THTAH2S ** TEMP202
            ZHTA1HSP (JT) = KAPP1HSP * THTAH2S ** TEMP202
            ZHTANO3F (JT) = KAPPNO3F * THTANO3 ** TEMP202
            ZHTANO3S (JT) = KAPPNO3S * THTANO3 ** TEMP202
       !
            ZHTAK2NO3 (JT) = K2NO3 * THTANO3 ** TEMP20 !m/d
       !
            ZL12NOM (JT) = THTADD ** TEMP20
            ZW12NOM (JT) = THTADP ** TEMP20 !NO UNIT
       !
            ZHTAPON1 (JT) = KPON1 * THTAPON1 ** TEMP20
            ZHTAPON2 (JT) = KPON2 * THTAPON2 ** TEMP20
            ZHTAPON3 (JT) = KPON3 * THTAPON3 ** TEMP20
            ZHTAPOC1 (JT) = KPOC1 * THTAPOC1 ** TEMP20
            ZHTAPOC2 (JT) = KPOC2 * THTAPOC2 ** TEMP20
       !
            IDEBUG_SED = 0
       !
            ZHTAPOC3 (JT) = KPOC3 * THTAPOC3 ** TEMP20
            ZHTAPOP1 (JT) = KPOP1 * THTAPOP1 ** TEMP20
            ZHTAPOP2 (JT) = KPOP2 * THTAPOP2 ** TEMP20
            ZHTAPOP3 (JT) = KPOP3 * THTAPOP3 ** TEMP20
            ZHTASI (JT) = KSI * THTASI ** TEMP20 !
            ZHTACH4 (JT) = KAPPCH4 * THTACH4 ** TEMP202
       !
         End Do
    !
    !***** Turn off settling
    !
         If ( .Not. SETTLING) Then
            Do I = 1, MLOC
               WSSBNET (I) = 0.
               WSLBNET (I) = 0.
               WSRBNET (I) = 0.
               WS1BNET (I) = 0.
               WS2BNET (I) = 0.
               WS3BNET (I) = 0.
               WSUBNET (I) = 0.
            End Do
         End If
    !
    !***** Initialize accumulators for steady-state computations
    !
         If (STEADY_STATE_SED_G3) Then
            TINTIM = 0.
            Do I = 1, MLOC
               AG3CFL (I) = 0.
               AG3NFL (I) = 0.
               AG3PFL (I) = 0.
               ASDTMP (I) = 0.
            End Do
         End If
    !
         Return
    !
      End Subroutine SED_INIT2
  !
  !
  !
  !********************************************************************************
  !**              S U B R O U T I N E   P O M _ A C C U M U L                   **
  !**      Accumulate POM from water column, before reaching the sedims          **
  !********************************************************************************
  !
      Subroutine POM_ACCUMUL
    !
         Use MOD_SAV, Only: NSAVCELL, SAVCELL, NSAVSPC, LEAF, STEM, &
        & SAVFRAC, WSSSAV, WSLSAV, WSRSAV, WS1SAV, WS2SAV, WS3SAV, &
        & WSUSAV, SEDPOCSAV, FRPOCSAV, SEDPONSAV, FRPONSAV, SEDPOPSAV, &
        & FRPOPSAV, SEDNH4SAV, SEDNO3SAV, SEDDOSAV, SEDPO4SAV
# if defined (WET_DRY)
	USE  MOD_WD, ONLY: ISWETN
# endif
         Implicit None
    !
         Integer :: I, N,K
         Real (SP) :: SAVEFCT
    !******* Assign base net settling rates
!
         Do I = 1, MLOC
            WSSNET (I) = WSSBNET (I)
            WSLNET (I) = WSLBNET (I)
            WSRNET (I) = WSRBNET (I)
            WS1NET (I) = WS1BNET (I)
            WS2NET (I) = WS2BNET (I)
            WS3NET (I) = WS3BNET (I)
            WSUNET (I) = WSUBNET (I)
         End Do
    !
    !******* Adjust net settling for SAV effect
    !
         If (SAV_CALC) Then
            Do I = 1, NSAVCELL
               B = SAVCELL (I)
               Do N = 1, NSAVSPC (B)
                  SAVEFCT = (LEAF(B, N)+STEM(B, N)) * SAVFRAC (B, N)
                  WSSNET (B) = WSSNET (B) + WSSSAV * SAVEFCT !increase of settling rate (m/d) to suspended solids
                  WSLNET (B) = WSLNET (B) + WSLSAV * SAVEFCT !increase of settling rate (m/d) to LPOM
                  WSRNET (B) = WSRNET (B) + WSRSAV * SAVEFCT !increase of settling rate (m/d) to RPOM
                  WS1NET (B) = WS1NET (B) + WS1SAV * SAVEFCT !increase of settling rate (m/d) to alg 1
                  WS2NET (B) = WS2NET (B) + WS2SAV * SAVEFCT !increase of settling rate (m/d) to alg 2
                  WS3NET (B) = WS3NET (B) + WS3SAV * SAVEFCT !increase of settling rate (m/d) to alg 3
                  WSUNET (B) = WSUNET (B) + WSUSAV * SAVEFCT !increase of settling rate (m/d) to particulate biogenic silicate (unavaiable)
               End Do
            End Do
         End If
    !
    !***** Calculate fluxes
    !
         Do I = 1, MLOC
            KWC = KBM1
       !
       !******* Settling rate (LHS has unit mm/day) ****
       !
            WSSINETMMD = 1000. * WSSNET (I)
            WLPOMNETMMD = 1000. * WSLNET (I)
            WRPOMNETMMD = 1000. * WSRNET (I)
            WB1NETMMD = 1000. * WS1NET (I)
            WB2NETMMD = 1000. * WS2NET (I)
            WB3NETMMD = 1000. * WS3NET (I)
            WPSINETMMD = 1000. * WSUNET (I)
       !
       !*** Fluxes***(mg/m^2/d) (positive into sediments)***
       !
            If (ITVWCLF /= 1) Then !use settling calculated from overlying water settling
          !
          !G1 flux of POC, Unit : [(mm/d)*(1)*gC/m^3 ] ==> (mgC/m^2/day)
          !     JPOCaccum (I, 1) = WB1NETMMD * FRCALG1 (1) * B1 (I, KWC) &
          !    & + WB2NETMMD * FRCALG2 (1) * B2 (I, KWC) + WB3NETMMD * &
          !    & FRCALG3 (1) * B3 (I, KWC) + WLPOMNETMMD * LPOC (I, KWC) &
          !    & + JPOCaccum (I, 1)!full channel has LPOC meaningful
			   JPOCaccum (I, 1) = WB1NETMMD * FRCALG1 (1) * B1 (I, KWC) &
              & + WB2NETMMD * FRCALG2 (1) * B2 (I, KWC) &
              & + WLPOMNETMMD * LPOC (I, KWC) &
              & + JPOCaccum (I, 1)!full channel has LPOC meaningful
          !
               If (SFEEDER) Then
                  JPOCaccum (I, 1) = JPOCaccum (I, 1) + JLPOC_SED_SF &
                 & (I) * FRCALG3 (1)!suspension feeders
               End If
          !
          !
          !G2 flux of POC, Unit : [(mm/d)*(1)*gC/m^3 ] ==> (mgC/m^2/day)
          !     JPOCaccum (I, 2) = WB1NETMMD * FRCALG1 (2) * B1 (I, KWC) &
          !    & + WB2NETMMD * FRCALG2 (2) * B2 (I, KWC) + WB3NETMMD * &
          !    & FRCALG3 (2) * B3 (I, KWC) + WRPOMNETMMD * RPOC (I, KWC) &
          !    & * FRPOC (I, 2) / (FRPOC(I, 2)+FRPOC(I, 3)) + JPOCaccum &
          !    & (I, 2)
			   JPOCaccum (I, 2) = WB1NETMMD * FRCALG1 (2) * B1 (I, KWC) &
              & + WB2NETMMD * FRCALG2 (2) * B2 (I, KWC) &
              & + WRPOMNETMMD * RPOC (I, KWC) &
              & * FRPOC (I, 2) / (FRPOC(I, 2)+FRPOC(I, 3)) + JPOCaccum &
              & (I, 2)
          !
               If (SFEEDER) Then
                  JPOCaccum (I, 2) = JPOCaccum (I, 2) + JLPOC_SED_SF &
                 & (I) * FRCALG3 (2)!suspension feeders
               End If
          !
          !
          !G3 flux of POC, Unit : [(mm/d)*(1)*gC/m^3 ] ==> (mgC/m^2/day)
          !     JPOCaccum (I, 3) = WB1NETMMD * FRCALG1 (3) * B1 (I, KWC) &
          !    & + WB2NETMMD * FRCALG2 (3) * B2 (I, KWC) + WB3NETMMD * &
          !    & FRCALG3 (3) * B3 (I, KWC) + WRPOMNETMMD * RPOC (I, KWC) &
          !    & * FRPOC (I, 3) / (FRPOC(I, 2)+FRPOC(I, 3)) + JPOCaccum &
          !    & (I, 3)
			   JPOCaccum (I, 3) = WB1NETMMD * FRCALG1 (3) * B1 (I, KWC) &
              & + WB2NETMMD * FRCALG2 (3) * B2 (I, KWC) &
              & + WRPOMNETMMD * RPOC (I, KWC) &
              & * FRPOC (I, 3) / (FRPOC(I, 2)+FRPOC(I, 3)) + JPOCaccum &
              & (I, 3)
          !
               If (SFEEDER) Then
                  JPOCaccum (I, 3) = JPOCaccum (I, 3) + JLPOC_SED_SF &
                 & (I) * FRCALG3 (3)!suspension feeders
               End If
          !
          !G1 flux of PON, Unit : [(mm/d)*(gN/gC)*(1)*gC/m^3 ] ==> (mgN/m^2/day)
          !     JPONaccum (I, 1) = WB1NETMMD * ANC1 * FRNALG1 (1) * B1 &
          !    & (I, KWC) + WB2NETMMD * ANC2 * FRNALG2 (1) * B2 (I, KWC) &
          !    & + WB3NETMMD * ANC3 * FRNALG3 (1) * B3 (I, KWC) + &
          !    & WLPOMNETMMD * LPON (I, KWC) + JPONaccum (I, 1)
			   JPONaccum (I, 1) = WB1NETMMD * ANC1 * FRNALG1 (1) * B1 &
              & (I, KWC) + WB2NETMMD * ANC2 * FRNALG2 (1) * B2 (I, KWC) &
              & + WLPOMNETMMD * LPON (I, KWC) + JPONaccum (I, 1)
          !
               If (SFEEDER) Then
                  JPONaccum (I, 1) = JPONaccum (I, 1) + JLPON_SED_SF &
                 & (I) * FRNALG3 (1)!suspension feeders
               End If
          !
          !
          !G2 flux of PON, Unit : [(mm/d)*(gN/gC)*(1)*gC/m^3 ] ==> (mgN/m^2/day)
          !     JPONaccum (I, 2) = WB1NETMMD * ANC1 * FRNALG1 (2) * B1 &
          !    & (I, KWC) + WB2NETMMD * ANC2 * FRNALG2 (2) * B2 (I, KWC) &
          !    & + WB3NETMMD * ANC3 * FRNALG3 (2) * B3 (I, KWC) + &
          !    & WRPOMNETMMD * RPON (I, KWC) * FRPON (I, 2) / (FRPON(I, &
          !    & 2)+FRPON(I, 3)) + JPONaccum (I, 2)
			   JPONaccum (I, 2) = WB1NETMMD * ANC1 * FRNALG1 (2) * B1 &
              & (I, KWC) + WB2NETMMD * ANC2 * FRNALG2 (2) * B2 (I, KWC) &
              & + WRPOMNETMMD * RPON (I, KWC) * FRPON (I, 2) / (FRPON(I, &
              & 2)+FRPON(I, 3)) + JPONaccum (I, 2)
          !
               If (SFEEDER) Then
                  JPONaccum (I, 2) = JPONaccum (I, 2) + JLPON_SED_SF &
                 & (I) * FRNALG3 (2)!suspension feeders
               End If
          !
          !G3 flux of PON, Unit : [(mm/d)*(gN/gC)*(1)*gC/m^3 ] ==> (mgN/m^2/day)
          !     JPONaccum (I, 3) = WB1NETMMD * ANC1 * FRNALG1 (3) * B1 &
          !    & (I, KWC) + WB2NETMMD * ANC2 * FRNALG2 (3) * B2 (I, KWC) &
          !    & + WB3NETMMD * ANC3 * FRNALG3 (3) * B3 (I, KWC) + &
          !    & WRPOMNETMMD * RPON (I, KWC) * FRPON (I, 3) / (FRPON(I, &
          !    & 2)+FRPON(I, 3)) + JPONaccum (I, 3)
			   JPONaccum (I, 3) = WB1NETMMD * ANC1 * FRNALG1 (3) * B1 &
              & (I, KWC) + WB2NETMMD * ANC2 * FRNALG2 (3) * B2 (I, KWC) &
              & + WRPOMNETMMD * RPON (I, KWC) * FRPON (I, 3) / (FRPON(I, &
              & 2)+FRPON(I, 3)) + JPONaccum (I, 3)
          !
               If (SFEEDER) Then
                  JPONaccum (I, 3) = JPONaccum (I, 3) + JLPON_SED_SF &
                 & (I) * FRNALG3 (3)!suspension feeders
               End If
          !
          !G1 flux of POP, Unit : [(mm/d)*(gP/gC)*(1)*gC/m^3 ] ==> mgP/m^2/day)
          !     JPOPaccum (I, 1) = WB1NETMMD * Q1 (I, KWC) * FRPALG1 (1) &
          !    & * B1 (I, KWC) + WB2NETMMD * Q2 (I, KWC) * FRPALG2 (1) * &
          !    & B2 (I, KWC) + WB3NETMMD * Q3 (I, KWC) * FRPALG3 (1) * &
          !    & B3 (I, KWC) + WLPOMNETMMD * LPOP (I, KWC) + JPOPaccum &
          !    & (I, 1)
               JPOPaccum (I, 1) = WB1NETMMD * Q1 (I, KWC) * FRPALG1 (1) &
              & * B1 (I, KWC) + WB2NETMMD * Q2 (I, KWC) * FRPALG2 (1) * &
              & B2 (I, KWC) &
              & + WLPOMNETMMD * LPOP (I, KWC) + JPOPaccum &
              & (I, 1)
          !
               If (SFEEDER) Then
                  JPOPaccum (I, 1) = JPOPaccum (I, 1) + JLPOP_SED_SF &
                 & (I) * FRPALG3 (1)!suspension feeders
               End If
          !
          !G2 flux of POP, Unit : [(mm/d)*(gP/gC)*(1)*gC/m^3 ] ==> (mgP/m^2/day)
          !     JPOPaccum (I, 2) = WB1NETMMD * Q1 (I, KWC) * FRPALG1 (2) &
          !    & * B1 (I, KWC) + WB2NETMMD * Q2 (I, KWC) * FRPALG2 (2) * &
          !    & B2 (I, KWC) + WB3NETMMD * Q3 (I, KWC) * FRPALG3 (2) * &
          !    & B3 (I, KWC) + WRPOMNETMMD * RPOP (I, KWC) * FRPOP (I, &
          !    & 2) / (FRPOP(I, 2)+FRPOP(I, 3)) + JPOPaccum (I, 2)
               JPOPaccum (I, 2) = WB1NETMMD * Q1 (I, KWC) * FRPALG1 (2) &
              & * B1 (I, KWC) + WB2NETMMD * Q2 (I, KWC) * FRPALG2 (2) * &
              & B2 (I, KWC) &
              & + WRPOMNETMMD * RPOP (I, KWC) * FRPOP (I, &
              & 2) / (FRPOP(I, 2)+FRPOP(I, 3)) + JPOPaccum (I, 2)
          !
               If (SFEEDER) Then
                  JPOPaccum (I, 2) = JPOPaccum (I, 2) + JLPOP_SED_SF &
                 & (I) * FRPALG3 (2)!suspension feeders
               End If
          !
          !G3 flux of POP, Unit : [(mm/d)*(gP/gC)*(1)*gC/m^3 ] ==> (mgP/m^2/day)
          !     JPOPaccum (I, 3) = WB1NETMMD * Q1 (I, KWC) * FRPALG1 (3) &
          !    & * B1 (I, KWC) + WB2NETMMD * Q2 (I, KWC) * FRPALG2 (3) * &
          !    & B2 (I, KWC) + WB3NETMMD * Q3 (I, KWC) * FRPALG3 (3) * &
          !    & B3 (I, KWC) + WRPOMNETMMD * RPOP (I, KWC) * FRPOP (I, &
          !    & 3) / (FRPOP(I, 2)+FRPOP(I, 3)) + JPOPaccum (I, 3)
               JPOPaccum (I, 3) = WB1NETMMD * Q1 (I, KWC) * FRPALG1 (3) &
              & * B1 (I, KWC) + WB2NETMMD * Q2 (I, KWC) * FRPALG2 (3) * &
              & B2 (I, KWC) &
              & + WRPOMNETMMD * RPOP (I, KWC) * FRPOP (I, &
              & 3) / (FRPOP(I, 2)+FRPOP(I, 3)) + JPOPaccum (I, 3)
          !
               If (SFEEDER) Then
                  JPOPaccum (I, 3) = JPOPaccum (I, 3) + JLPOP_SED_SF &
                 & (I) * FRPALG3 (3)!suspension feeders
               End If
          !
          !Particulate Biogeonic Silica flux to sediments (positive to sediments)
               JPOSaccum (I) = WB1NETMMD * ASC1 * B1 (I, KWC) + &
              & WB2NETMMD * ASC2 * B2 (I, KWC) &
              & + WPSINETMMD * SIUPB (I, KWC) + JPOSaccum (I)
          !
          !Should add RPOM contribution to Silicate based Si/C ratio in detritus
          !and JSIDETR is supposed to caculate that instead of being given from input
          !
               If (SFEEDER) Then
                  JPOSaccum (I) = JPOSaccum (I) + JSU_SED_SF (I)!+    SU_SED_SF(I)    !suspension feeder
               End If
          !
            Else
          !
          !
               JPOCaccum (I, 1) = JCIN_R1 * FRCALG1 (1) * 1000.0 / &
              & 2.667 !G1 POC flux converted from gO2/m^2/day to mgC/m^2/day
               JPOCaccum (I, 2) = JCIN_R1 * FRCALG1 (2) * 1000.0 / &
              & 2.667 !G2 POC flux converted from gO2/m^2/day to mgC/m^2/day
               JPOCaccum (I, 3) = JCIN_R1 * FRCALG1 (3) * 1000.0 / &
              & 2.667 !G3 POC flux converted from gO2/m^2/day to mgC/m^2/day
          !
               JPONaccum (I, 1) = JNIN_R1 * FRNALG1 (1) * 1000.0 !G1 PON flux converted from gN/m^2/day to mgN/m^2/day
               JPONaccum (I, 2) = JNIN_R1 * FRNALG1 (2) * 1000.0 !G2 PON flux converted from gN/m^2/day to mgN/m^2/day
               JPONaccum (I, 3) = JNIN_R1 * FRNALG1 (3) * 1000.0 !G3 PON flux converted from gN/m^2/day to mgN/m^2/day
          !
               JPOPaccum (I, 1) = JPIN_R1 * FRPALG1 (1) * 1000.0 !G1 POP flux converted from gP/m^2/day to mgP/m^2/day
               JPOPaccum (I, 2) = JPIN_R1 * FRPALG1 (2) * 1000.0 !G2 POP flux converted from gP/m^2/day to mgP/m^2/day
               JPOPaccum (I, 3) = JPIN_R1 * FRPALG1 (3) * 1000.0 !G3 POP flux converted from gP/m^2/day to mgP/m^2/day
          !
          !
          !
               JPOSaccum (I) = JSIN_R1 * 1000.0 !Silicate flux from gSi/m^2/day to mgSi/m^2/day
          !
          !
            End If
         End Do

		 DO I=1, MLOC
#  if !defined (WET_DRY)
			IF(D(I) <= 0.0_SP)THEN
#  else
			IF(ISWETN(I) ==0)THEN
#  endif
				DO K=1,3
					 JPOCaccum (I, K)=0.0_SP
					 JPONaccum (I, K)=0.0_SP
					 JPOPaccum (I, K)=0.0_SP
				ENDDO
			    JPOSaccum (I)=0.0_SP
			ENDIF
		 ENDDO


    !
      End Subroutine POM_ACCUMUL
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   S E D _ C A L C                   **
  !**                           Sediment Calculations                            **
  !********************************************************************************
  !
      Subroutine SED_CALC (TF_SSTATE)
    !
         Use MOD_SED_DF_EXCHANGE_VARS, Only: POC1TM1S_SED_DF, &
        & POC2TM1S_SED_DF
    !
		 Use MOD_SED_SAV_EXCHANGE_VARS, Only: NH4T2TM1S_SHARE, &
        & PO4T2TM1S_SHARE, HST2TM1S_SHARE !AN
	!
         Use MOD_SAV, Only: NSAVCELL, SAVCELL, NSAVSPC, LEAF, STEM, &
        & SAVFRAC, WSSSAV, WSLSAV, WSRSAV, WS1SAV, WS2SAV, WS3SAV, &
        & WSUSAV, SEDPOCSAV, FRPOCSAV, SEDPONSAV, FRPONSAV, SEDPOPSAV, &
        & FRPOPSAV, SEDNH4SAV, SEDNO3SAV, SEDDOSAV, SEDPO4SAV
#if defined (MULTIPROCESSOR)
         Use MOD_PAR, Only: NGID
#endif
# if defined (WET_DRY)
	USE  MOD_WD, ONLY: ISWETN
# endif
    !
          Implicit None
    !
         Save
         Integer :: I, J, JSF, N, IERR, IDEBUG

         Real (SP) :: SAVEFCT
         Real (SP) :: SOD !,ZBOUT
    !
         Real (SP) :: DFEED_TMP, DF_SOD_TMP
    !
         Real (SP) :: SAV_SOD_TMP
    !
         Logical :: TF_SSTATE !Type of integration  FALSE - time intergration from current time step to next time step
    !                        TRUE  - no time integration, calculate steady steate solution
!
         Real (SP) :: Ibreak 
!
         DLTS = NDTSED * DLT / 86400. !time step in days for sediment flux model - sediment timestep is NDTSED times larger than water column time step (DLT) -added by LB
    !
         If (STEADY_STATE_SED_G3) TINTIM = TINTIM + DLTS
    !
    !***** Initialize sediment nutrient masses
    !
    !
         SEDMC = 0.
         SEDMN = 0.
         SEDMP = 0.
    !
    !
    
    !
         Do I = 1, MLOC ! moved DO here because of moved calculations of JPOM
            KWC = KBM1
       !
       ! *** Sum particulate fluxes to water column, negative into sediments (postitive into water)
       ! used for detecting the mass balance from the water column's point of view
       ! (LHS is g/m^2/day)
       !
       !Sum of particuplate P, N, C, Particulate Biogenic Si, Suspended Solids to water column (g/m^2/day)
            PPFWS (I) = - 0.001 * (JPOP(I, 1)+JPOP(I, 2)+JPOP(I, 3))!converted from mgP/m^2/d --> gP/m^2/day
            PNFWS (I) = - 0.001 * (JPON(I, 1)+JPON(I, 2)+JPON(I, 3))!converted from mgN/m^2/d --> gN/m^2/day
            PCFWS (I) = - 0.001 * (JPOC(I, 1)+JPOC(I, 2)+JPOC(I, 3))!converted from mgC/m^2/d --> gC/m^2/day
       !
            PSFWS (I) = - 0.001 * JPOS (I)!converted from mgSi/m^2/d --> gSi/m^2/day
       !
       !Unit: (m/day)*(g/m^3) --> gSolids/m^2/day
            SSFWS (I) = - WSSNET (I) * SSI (I, KWC)
            If (SFEEDER) Then
               SSFWS (I) = SSFWS (I) - 0.001 * SSI_SED_SF (I)!converted from mgSolids/m^2/day -->gSolids/m^2/day
            End If
       !
         End Do
    !
    ! ADD IN THE FLUX FROM ROOT MORTALITY OF SAV
    !
         If (SAV_CALC) Then
            Do I = 1, NSAVCELL
               B = SAVCELL (I)
               Do J = 1, 3
                  JPOC (B, J) = JPOC (B, J) + 1000. * SEDPOCSAV (B) * &
                 & FRPOCSAV (J)
                  JPON (B, J) = JPON (B, J) + 1000. * SEDPONSAV (B) * &
                 & FRPONSAV (J)
                  JPOP (B, J) = JPOP (B, J) + 1000. * SEDPOPSAV (B) * &
                 & FRPOPSAV (J)
               End Do
            End Do
         End If
    !
    !    *** Accumulate fluxes for steady-state computation for G3
    !
         If (STEADY_STATE_SED_G3) Then
            Do I = 1, MLOC
               AG3CFL (I) = AG3CFL (I) + JPOC (I, 3) * DLTS
               AG3NFL (I) = AG3NFL (I) + JPON (I, 3) * DLTS
               AG3PFL (I) = AG3PFL (I) + JPOP (I, 3) * DLTS
            End Do
         End If
!
    !
    !***** Assign previous timestep concentrations to particulate organics
!
         Do I = 1, MLOC
            CPOP (I, 1) = POP1TM1S (I)
            CPOP (I, 2) = POP2TM1S (I)
            CPOP (I, 3) = POP3TM1S (I)
            CPON (I, 1) = PON1TM1S (I)
            CPON (I, 2) = PON2TM1S (I)
            CPON (I, 3) = PON3TM1S (I)
            CPOC (I, 1) = POC1TM1S (I)
            CPOC (I, 2) = POC2TM1S (I)
            CPOC (I, 3) = POC3TM1S (I)
            CPOS (I) = PSISEDTM1S (I)
         End Do
!
    !
    !***** Update sediment concentrations
    !
    !calculation deposition feeder and flux terms
    !due to deposition feeder
         If (DFEEDER) Then
            Call DF_CALC
         End If
	   !
         Do I = 1, MLOC

       !
       !get previous time step concentrations of POM and inorganics
       !
            POC1TM1 = POC1TM1S (I)
            POC2TM1 = POC2TM1S (I)
            POC3TM1 = POC3TM1S (I)
       !
            PON1TM1 = PON1TM1S (I)
            PON2TM1 = PON2TM1S (I)
            PON3TM1 = PON3TM1S (I)
       !
            POP1TM1 = POP1TM1S (I)
            POP2TM1 = POP2TM1S (I)
            POP3TM1 = POP3TM1S (I)
       !
            PSISEDTM1 = PSISEDTM1S (I)
       !
            BENSTRTM1 = BENSTRTM1S (I)
       !
            SODTM1 = SODTM1S (I)
       !
            NH41TM1 = NH41TM1S (I)
            NH4T2TM1 = NH4T2TM1S (I)
       !
            NO31TM1 = NO31TM1S (I)
            NO3T2TM1 = NO3T2TM1S (I)
       !
            PO41TM1 = PO41TM1S (I)
            PO4T2TM1 = PO4T2TM1S (I)
       !
            SI1TM1 = SI1TM1S (I)
            SIT2TM1 = SIT2TM1S (I)
       !
       !
            HS1TM1 = HS1TM1S (I)
            HST2TM1 = HST2TM1S (I)
       !
       !
            CH41TM1 = CH41TM1S (I)! CH4 of layer 1
            CH4T2TM1 = CH4T2TM1S (I)! total CH4 of layer 2
       !
            SO4T2TM1 = SO4T2TM1S (I)! total SO4 of layer 2
       !
       ! dfeedm1  = dfeedm1s(I)  !deposition feeder befor integrating to next time step
       !moved to mod_df.F
       !
       ! ACCOUNT FOR SAV NUTRIENT UPTAKE
       ! NOTE THIS IS OVER ALL CELLS, NOT JUST SAV CELLS SO SEDNH4SAV
       ! MUST BE ZEROED OUT FOR ALL CELLS
       !
            If (SAV_CALC) Then
               NH4T2TM1 = NH4T2TM1 - 1000. * SEDNH4SAV (I) * DLTS / &
              & HSED (I)
               NO3T2TM1 = NO3T2TM1 - 1000. * SEDNO3SAV (I) * DLTS / &
              & HSED (I)!: Added SEDNO3SAV here
               PO4T2TM1 = PO4T2TM1 - 1000. * SEDPO4SAV (I) * DLTS / &
              & HSED (I)
            End If
!
       !
            BFORMAX = BFORMAXS (I)
            ISWBEN = ISWBENS (I)
            H2 = HSED (I)
       !
       !******* Sedimentation, mixing rates, and sediment temperature
       !
            W2 = VSED (I)!m/d
            DPP = VPMIX (I)!Particle mixing coefficient m**2/day
       !Renamed DD to DDP so that is is not duplicated with DD elsewhere out of sediment module
       !DD    = VDMIX(I)
            DDP = VDMIX (I)!dissolved phase diffusion (mixing) m**2/day
            TEMPD = CTEMP (I)!temperature of sediments
            STP20 = TEMPD - 20.
       !
       !
            CSISAT = CSISATT20 * THTASISAT ** STP20 ! calculation saturation concentration based on temperature control
            KWC = KBM1
       !
       !PO40 (mgP/m^3)
       !
       !steady state estimate (isotherm) of dissolved PO4 in overlying water column
            DF = 1. / (1.+KADPO4*SSI(I, KWC))
            PO4AVL = DF * PO4 (I, KWC)!avaiable PO4 for algae growth (dissolved)
            PO40 = PO4AVL * 1000. !convert from gP/m^3 to mgP/m^3
       !
       !NH40 (mgN/m^3)
            NH40 = NH4 (I, KWC) * 1000. !convert from gN/m^3 to mgN/m^3
            NO30 = NO3 (I, KWC) * 1000. !convert from gN/m^3 to mgN/m^3
       !
       !SI0 (mgSi/m^3)
            DF = 1. / (1.+KADSA*SSI(I, KWC))
            SI0 = DF * SIAT (I, KWC) * 1000. !convert from gSi/m^3 to mgSi/m^3
       !
       !O20 (mgO2/L) = (gO2/m^3)
            O20 = AMAX1 (DOXG(I, KWC), 0.010)!mgO2/L
       !
       !HS0 (gO2/m^3)
            HS0 = COD (I, KWC)!mgO2/L (in oxygen equivalents)
       !
       !salinity (ppt)
            SAL = SALT (I, KWC)!overlying water salinity (ppt)
       !
       !
       !SO40 (gO2/m^3)
       !
       ! lets flag DOXG if it starts near zero as it may cause root finder problems
       ! with suspension and deposit feeders on.
       !
       !
       !       Regression function to get SO4 concentration from SAL
       !       [SO4] = 20 mgO2/L        for        [Cl] < 6 mgO2/L
       !             = (10/3)[Cl]       for        [Cl] > 6 mgO2/L
       !       1 ppt = 607.445 mg/L Cl
       !
            If (SAL > 0.0099) Then
               SO40MG = 20.0 + (27./190.) * 607.445 * SAL !mgSO4/L
            Else
               SO40MG = 20.0 !mgSO4/L
            End If
       !**** units: so4 in o2 equivalents
       !     SO4 (mg so4/L)* 1 mmol SO4 /98 mg SO4 * 2 mmol O2/ 1 mmol SO4
       !     * 32 mg O2 / mmol O2= 0.65306122
       !
       !==> 1 mgSO4/L ~ 0.65306122 mgO2/L
       !
            SO40 = SO40MG * 0.65306122 !mgO2/L, i.e gO2/m^3
       !
       !CH40 (gO2/m^3)
            If (ITVWCLF /= 1) Then
               CH40 = 0.
            Else
               CH40 = CH40_R1
            End If
       !
       !Methane saturation (gO2/m^3) in oxygen equivalents
       !
       !   CH4SAT = 99.*(1.+(ZD(I,KWC)+D(I)*DZ2D(I,KWC)+HSED(I))/10.)*0.9759**STP20
       !     CH4SAT = 99.*(1.+(D(I)+HSED(I))/10.)*0.9759**STP20
            CH4SAT = 100. * (1.+(D(I)+HSED(I))/10.) * 1.024 ** (-STP20)!, 0.9759 should be 1.024 according to QAPP
       !Also exponent is 20-T rather than T-20
       !And 1/1.024 = 0.9765625
       !     1/0.9759=1.0246951531919255
       !which means they are close!
!
       !
       !******* Evaluate the temperature dependent coefficients
       !
            ITEMP = 10. * TEMPD + 1
       !
       !******* Salinity dependence of nitrification and denitrification
       !
            If (SAL < SALTND) Then
               XAPPNH4 = GET_ZHTANH4F (TEMPD)
               XAPP1NO3 = GET_ZHTANO3F (TEMPD)
            Else
               XAPPNH4 = GET_ZHTANH4S (TEMPD)
               XAPP1NO3 = GET_ZHTANO3S (TEMPD)
            End If
       !
            XAPP1HSD = GET_ZHTA1HSD (TEMPD)
            XAPP1HSP = GET_ZHTA1HSP (TEMPD)
       !
            XK2NO3 = GET_ZHTAK2NO3 (TEMPD)!*H2        !m/d     !
            XKSI = GET_ZHTASI (TEMPD) * H2 !m/d     !DMD book eqn(7.7)
       !
            XAPPCH4 = GET_ZHTACH4 (TEMPD)
       !
       !    KL12NOM  = DDP/H2*GET_ZL12NOM(TEMPD)         !OLD code uses H2 (this depends on definition of DD)
       !
            KL12NOM = DDP / (H2/2) * GET_ZL12NOM (TEMPD)!normalized porewater diffusion rate across layer1,2 (m/day)

       !
       !    W12NOM   = DPP/H2*GET_ZW12NOM(TEMPD)*POC1/1.0E5     !Old code uses H2
       !
            If (TF_SSTATE) Then
          !
               W12NOM = DPP / (H2/2) * GET_ZW12NOM (TEMPD)!
          !
            Else
          !
               W12NOM = DPP / (H2/2) * GET_ZW12NOM (TEMPD) * POC1TM1 / &
              & (POC1R*1.0E6*M2)!New code uses H2/2
          !we see new DPP is half of old DPP to keep  WL12NOM same
          !
          !Here "*1.0E6*M2 converts  POC1R from mgC/gSediment to mgC/m^3
          !
          !normalized particle diffusion rate across layer1,2 !changed POC1 to POC1TM1
          !: POC1R = (1.0E5*M2) (mgC/m^3) if POC1R= 0.1mgC/gSediment
          ! eqn(3-26) (DMD book page 278, table 13.1, eqn(13.1))
          !
          !M2 has unit kg/L = 1000kg/m^3 =1000,000 g/m^3
          !
          !POC1R = 0.1mgC/gSediment *M2 ~ 0.1mgC *1000,000 gSediment/m^3  = 1.0E5 *M2 (mgC/m^3)
          !
          !
          !
          !In excel program: w12base = Dp * ThtaDp ^ (Tw - 20) / (H2 / 2) * (POC2(1) / 1000# / m2) / POC1R
          !where H2 is repaced by H2/2
          !POC2(1) is POC in layer 2 for G1 class in gO2/m^3Sediment unit (~ 1 mgO2/Lsediment)
          !m2 is sediment concentration in kg/L
          !POC1R is reference POC2(1) for bioturbation , ~= 0.2667 mgO2/gSediment
          !POC2(1)/1000/m2 = POC2(1)/(1000*m2) converts to mgO2/gSediment
          !==> POC2(1)/1000/m2 /POC1R finally has unit of (mgO2/gSediment)/(mgO2/gSediment) ==> unitness
          !
          ! Derivation:
          !POC1R in DMD eqn (13.1) is 0.1mgC/g sediment ~ 0.1gC/kg
          !m2 has unit  kg/L            gC      kg
          !==> POC1R = 0.1gC/kg ~0.1 ----- *m2---- = 0.1*m2 gC/L  ~ 0.1*m2 * 1000gC/m^3 = 100*m2 (gC/m^3) =100000*m2 (mgC/m^3)
          !                             kg      L
          !            ~ 0.1*10E6*m2 (mgC/m^3)    !i.e. POC1R*1.0E6*M2 used in W12NOM calculation above
          !
          !
            End If
       !
       !---------------
       !
            BFOR = KMO2DP / (KMO2DP+O20)
       !
            If ( .Not. TF_SSTATE) Then
          !
          !predict benthic pressure at new time setp based on previous value
               BENSTR = (BENSTRTM1+DLTS*BFOR) / (1.+KBENSTR*DLTS)! ()
          !
            Else !steady state
          !
               BENSTR = BFOR / KBENSTR ! ()
          !
               BENSTR = 0.d0 !
          !
            End If
       !
            BFORMAX = AMAX1 (BENSTR, BFORMAX)
       !
       !     W12MIN = DPMIN/H2 !Old Code
            W12MIN = DPMIN / (H2/2.0)!New Code ()
       !We see that new DPMIN is half of old DPMIN to keep WL12MIN the same
            W12 = W12NOM * (1.-KBENSTR*BFORMAX)!+W12MIN    ! removed minimum mixing
       !
       !## -- klbnth is ratio of bio-irrigation to bio-particle mixing
            KL12 = KL12NOM + KLBNTH * W12 !dissolved phase mixing (KL12) affected by bio-particle mixing W12
       !
       !******* Lookup reaction rates (m/d)
       !
            ITEMP = 10. * TEMPD + 1
       !
            XKPOC1 = GET_ZHTAPOC1 (TEMPD) * H2
            XKPOC2 = GET_ZHTAPOC2 (TEMPD) * H2
            XKPOC3 = GET_ZHTAPOC3 (TEMPD) * H2
            XKPOP1 = GET_ZHTAPOP1 (TEMPD) * H2
            XKPOP2 = GET_ZHTAPOP2 (TEMPD) * H2
            XKPOP3 = GET_ZHTAPOP3 (TEMPD) * H2
            XKPON1 = GET_ZHTAPON1 (TEMPD) * H2
            XKPON2 = GET_ZHTAPON2 (TEMPD) * H2
            XKPON3 = GET_ZHTAPON3 (TEMPD) * H2
!
       !******* Calculate sediment concentrations
       !
            DLTS_H2 = DLTS / H2 !multipler dt/H2 in front of sediment diagensis equations 
       !
            FRPON1 = 1. - (FRPON(I, 2)+FRPON(I, 3))
            FRPOP1 = 1. - (FRPOP(I, 2)+FRPOP(I, 3))
            FRPOC1 = 1. - (FRPOC(I, 2)+FRPOC(I, 3))
       !
       !     IDEBUG=8
       !     IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
       !
       !sediment silicate
            FD2 = 1. / (1.+M2*PIE2SI)! (4-6-7)(3-7)
            K3 = XKSI * (CSISAT-FD2*SIT2TM1) / (PSISEDTM1+KMPSI)!Dissolution rate (m/day) ()
       !
            If ( .Not. TF_SSTATE) Then !non-steady state calculation, integration one time step DLTS
          !
               POC1 = (JPOC(I, 1)*DLTS_H2+POC1TM1) / &
              & (1.+(XKPOC1+W2)*DLTS_H2)
               POC2 = (JPOC(I, 2)*DLTS_H2+POC2TM1) / &
              & (1.+(XKPOC2+W2)*DLTS_H2)
          !
               POC3 = (JPOC(I, 3)*DLTS_H2+POC3TM1) / &
              & (1.+(XKPOC3+W2)*DLTS_H2)
          !
               PON1 = (JPON(I, 1)*DLTS_H2+PON1TM1) / &
              & (1.+(XKPON1+W2)*DLTS_H2)
          !
               PON2 = (JPON(I, 2)*DLTS_H2+PON2TM1) / &
              & (1.+(XKPON2+W2)*DLTS_H2)
          !
               PON3 = (JPON(I, 3)*DLTS_H2+PON3TM1) / &
              & (1.+(XKPON3+W2)*DLTS_H2)
          !
               POP1 = (JPOP(I, 1)*DLTS_H2+POP1TM1) / &
              & (1.+(XKPOP1+W2)*DLTS_H2)
          !
               POP2 = (JPOP(I, 2)*DLTS_H2+POP2TM1) / &
              & (1.+(XKPOP2+W2)*DLTS_H2)
          !
               POP3 = (JPOP(I, 3)*DLTS_H2+POP3TM1) / &
              & (1.+(XKPOP3+W2)*DLTS_H2)
          !
          !## -- modification for detrital Si input to sediment
          !##  PSISED  = (JPOS(I)*DLTS/H2+PSISEDTM1)/(1.+(K3+W2)*DLTS/H2)
          !
               PSISED = ((JPOS(I)+JSIDETR)*DLTS_H2+PSISEDTM1) / &
              & (1.+(K3+W2)*DLTS_H2)!(DMD book 13.31, 13.33)
          !
               PSISED = (JPOS(I)*DLTS_H2+PSISEDTM1) / &
              & (1.+(XKSI+W2)*DLTS_H2)!
          !
          !JSIDETR here should be included in JPOS calculation
          !see DMD p391 JdetrSi (16A.21)
          !Modification to PON, POC, POP due to deposition feeder fluxes
               If (DFEEDER) Then
             !
                  POC1 = POC1 + (-DF_GROW_POC1(I)*DLTS/H2+&
                 & FRPOC1*(DF_MORT(I)+DF_PRED(I))*DLTS/H2) / (1.+&
                 & (XKPOC1+W2)*DLTS_H2)!contribute to POC1 at all?
             !
                  POC2 = POC2 + (-DF_GROW_POC2(I)*DLTS/H2+FRPOC(I, &
                 & 2)*(DF_MORT(I)+DF_PRED(I))*DLTS/H2) / &
                 & (1.+(XKPOC2+W2)*DLTS_H2)
             !
                  PON1 = PON1 + (-DF_GROW_POC1(I)*DLTS/H2/AMCN+&
                 & FRPON1*(DF_MORT(I)+DF_PRED(I))*DLTS/H2/AMCN) / (1.+&
                 & (XKPON1+W2)*DLTS_H2)!contribute to PON1 at all? Yes
             !
                  PON2 = PON2 + (-DF_GROW_POC2(I)*DLTS/H2/AMCN+FRPON(I, &
                 & 2)*(DF_MORT(I)+DF_PRED(I))*DLTS/H2/AMCN) / &
                 & (1.+(XKPON2+W2)*DLTS_H2)
             !
                  POP1 = POP1 + (-DF_GROW_POC1(I)*DLTS/H2/AMCP+&
                 & FRPOP1*(DF_MORT(I)+DF_PRED(I))*DLTS/H2/AMCP) / (1.+&
                 & (XKPOP1+W2)*DLTS_H2)!contribute to POP1 at all? Yes
             !
                  POP2 = POP2 + (-DF_GROW_POC2(I)*DLTS/H2/AMCP+FRPOP(I, &
                 & 2)*(DF_MORT(I)+DF_PRED(I))*DLTS/H2/AMCP) / &
                 & (1.+(XKPOP2+W2)*DLTS_H2)
             !
               End If
          !
            Else !STEADY STATE CALCULATION
!
               POC1 = JPOC (I, 1) / (XKPOC1+W2)!
               POC2 = JPOC (I, 2) / (XKPOC2+W2)
               POC3 = JPOC (I, 3) / (XKPOC3+W2)
!
               PON1 = JPON (I, 1) / (XKPON1+W2)!
               PON2 = JPON (I, 2) / (XKPON2+W2)
               PON3 = JPON (I, 3) / (XKPON3+W2)
          !
               POP1 = JPOP (I, 1) / (XKPOP1+W2)!
               POP2 = JPOP (I, 2) / (XKPOP2+W2)
               POP3 = JPOP (I, 3) / (XKPOP3+W2)
          !
          !## -- modification for detrital Si input to sediment
               PSISED = (JPOS(I)+JSIDETR) / (K3+W2)! (DMD book 13.31, 13.33)
               PSISED = JPOS (I) / (XKSI+W2)!
          !
               If (DFEEDER) Then
             !
                  POC1 = (JPOC(I, &
                 & 1)-DF_GROW_POC1(I)+FRPOC1*(DF_MORT(I)+DF_PRED(I))) / &
                 & (XKPOC1+W2)!
             !
             !
                  POC2 = (JPOC(I, 2)-DF_GROW_POC2(I)+FRPOC(I, &
                 & 2)*(DF_MORT(I)+DF_PRED(I))) / (XKPOC2+W2)!
             !
                  PON1 = (JPON(I, 1)-DF_GROW_POC1(I)/AMCN+&
                 & FRPON1*(DF_MORT(I)+DF_PRED(I))/AMCN) / (XKPON1+W2)! 
             !
             !
                  PON2 = (JPON(I, 2)-DF_GROW_POC2(I)/AMCN+FRPON(I, &
                 & 2)*(DF_MORT(I)+DF_PRED(I))/AMCN) / (XKPON2+W2)! 
             !
             !
                  POP1 = (JPOP(I, 1)-DF_GROW_POC1(I)/AMCP+&
                 & FRPOP1*(DF_MORT(I)+DF_PRED(I))/AMCP) / (XKPOP1+W2)! 
             !
             !
                  POP2 = (JPOP(I, 2)-DF_GROW_POC2(I)/AMCP+FRPOP(I, &
                 & 2)*(DF_MORT(I)+DF_PRED(I))/AMCP) / (XKPOP2+W2)! 
             !
             !
               End If
          !
            End If
       !
            If (PON1 < 0.0) PON1 = 0.0
            If (PON2 < 0.0) PON2 = 0.0
            If (POC1 < 0.0) POC1 = 0.0
            If (POC2 < 0.0) POC2 = 0.0
            If (POP1 < 0.0) POP1 = 0.0
            If (POP2 < 0.0) POP2 = 0.0
            If (PSISED < 0.0) PSISED = 0.0 !Make sure PSISED >=0
       !
       !******* Assign diagenesis values for sediment model
       !basically the decomposition of PON, POC, POP in sediments
       !aided by respiration of deposition feeder
       !
            JCX = XKPOC1 * POC1 + XKPOC2 * POC2 + XKPOC3 * POC3 !             (m/d * mgC/m^3  ==> mgC/m^2/day)
            JNX = XKPON1 * PON1 + XKPON2 * PON2 + XKPON3 * PON3 ! JN term (mgN/m^2/day)
            JPX = XKPOP1 * POP1 + XKPOP2 * POP2 + XKPOP3 * POP3 ! JP term (mgP/m^2/day)
       !
!
            If (DFEEDER) Then
               JCX = JCX + DF_RESP (I)!deposition feeder contribution
               JNX = JNX + DF_RESP (I) * (1./AMCN)!deposition feeder contribution converted from C to N
               JPX = JPX + DF_RESP (I) * (1./AMCP)!deposition feeder contribution converted from C to P
            End If
       !
       !: The sediment flux are solved in the following procedures
       ! (1) take initial guess of SOD and H1
       ! (2) solve for CT1 and CT2 at t+dt time for NH4,NO3, HS and CH4
       !       (a) solve for NH4 at t+dt in layer 1 and layer 2
       !       (b) compute oxygen demand consumed by nitrification, called NSOD
       !       (c) solve for NO3 at t+dt in layer 1 and layer 2
       !       (d) compute CH4 (freshwater) or HS( salt water) oxidation
       !             <i> for salt water, compute HS reaction terms and compute SOD due to HS oxidation
       !             <ii> for fresh water, compute CH4 flux by establishing the chemical specific conditions
       !                    [1] compare computed and saturation concentrations of CH4 and correct it
       !                     [2] calculate CSOD due to methane oxidation
       !             <1>compute CSOD  due to HS or CH4
       !             <2>compute flux of JCH4 or JHS
       !             <3>compute total SOD = NSOD + CSOD
       !       (e) refine estimate of SOD, make new estimate of SOD
       !  (3)go to (2) if no covergence of SOD and H1
       !  (4)solve Si and PO4 using steady state equations
       !
       !
       !****** Evaluate the NH4, NO3, SO4, HS, CH4 and SOD equations *****
       !
       !
       !initial guess of SOD (including SOD due to deposition feeder)
            If (DFEEDER) Then
               DF_SOD_TMP = DF_SOD (I)
            Else
               DF_SOD_TMP = 0.0
            End If
       !
            If (SAV_CALC) Then
               SAV_SOD_TMP = - SEDDOSAV (I)!SAV_SOD(I)  !SEDDOSAV
            Else
               SAV_SOD_TMP = 0.0
            End If
!
       ! REAL(SP) ::     DFSOD_SED !SOD due to deposition feeder
       !REAL(SP) ::     SAVSOD_SED !SOD due to SAV root respiration
       !
       !
       !solve NH4, NO3, SO4, HS, CH4 using ZBRENT function which
       !calls SEDF, where SEDF solves NH4, NO3, SO4, HS, CH4
       !and finally find out what SOD is until ZBRENT converges
       !
            If (TF_SSTATE) Then
               SODTM1 = 2.667 * 1.0E-3 * JCX + 1.714 * 1.0E-3 * JNX !mgN/m^2/day converted to mgO2/m^2/day then to gO2/m^2/day
            Else
			   !SODTM1=SODTM1  !no change
               !
            End If
       !
            SODMIN = 1.E-4
            SODMAX = 100.

#  if !defined (WET_DRY)
		IF(D(I) > 0.0_SP)THEN
#  else
		IF(ISWETN(I) == 1)THEN
#  endif

            SOD = ZBRENT (IERR, DF_SOD_TMP, SAV_SOD_TMP, TF_SSTATE, &
           & SODTM1)!returns IERR and SOD
       !
       !ZBRENT returns IERR and value ZBRENT as SOD
       !if IERR = 0   good, find soultion and SOD converges
       !   IERR = 1   did not find solution of SOD within [SODMIN, SODMAX] range
       !   IERR = 2   did not find the solution within the IMAX interations (did not converge)
       !
            If (IERR == 1 .Or. IERR == 2) Then
          !
               If (DFEEDER) Then
                  DF_SOD_TMP = DF_SOD (I)
                  DFEED_TMP = DFEEDM1S (I)
               Else
                  DF_SOD_TMP = 0.0
                  DFEED_TMP = 0.0
               End If
          !
               If (SAV_CALC) Then
                  SAV_SOD_TMP = SEDDOSAV (I)!SAV_SOD(I)
               Else
                  SAV_SOD_TMP = 0.0
               End If
          !
          !this is to use SODMIN and SODMAX to estimate the bounds of error
               SODMIN = 0.0001
               SODMAX = 100.
          !
               If (SFEEDER) Then !three species of suspension feeder and their contribution to sediment diagenesis
             !    SFEED_TMP_1= SFEED(I,1)    !this needs to be treated as SFEED_SED_SF in mod_sed_sf_exchange_vars.F
             !    SFEED_TMP_2= SFEED(I,2)
             !    SFEED_TMP_3= SFEED(I,3)
                  SFEED_TMP_1 = 0. !disabled for now
                  SFEED_TMP_1 = 0.
                  SFEED_TMP_1 = 0.
               Else
                  SFEED_TMP_1 = 0.
                  SFEED_TMP_2 = 0.
                  SFEED_TMP_3 = 0.
               End If
          !
               If (MSR) WRITE (6, 9000) JDAY, IERR, I, SAL, SO40MG, &
              & DFEED_TMP, SFEED_TMP_1, SFEED_TMP_2, SFEED_TMP_3, &
              & SODMIN, SODMAX
          !
               If (MSR) WRITE (6, 9911) CSODHS, CSODCH4, CSOD
9911           Format (/ 1 X, ' CSODHS, CSODCH4, CSOD' / 3E10.3)
               If (MSR) WRITE (6, 9910) CH41, CH42, HST1, HS1, HS2
9910           Format (/ 1 X, ' CH41   CH42   HST1   HS1   HS2' / &
              & 5E10.3)
          !
               If (IERR == 2) Then !SOD=ZBRENT(SOD) did not converge
             !
                  If (MSR) WRITE (6, 9900) JDAY, CTEMP (I), POP1, POP2, &
                 & POP3
                  If (MSR) WRITE (6, 9901) PON1, PON2, PON3, POC1, &
                 & POC2, POC3
                  If (MSR) WRITE (6, 9902) PO4T2, HST2, SIT2, PSISED
                  If (MSR) WRITE (6, 9903) (JPOP(I, 1)+JPOP(I, &
                 & 2)+JPOP(I, 3)), (JPON(I, 1)+JPON(I, 2)+JPON(I, 3)), &
                 & (JPOC(I, 1)+JPOC(I, 2)+JPOC(I, 3))
                  If (MSR) WRITE (6, 9904) O20, CSOD, DF_SOD_TMP, SOD, &
                 & S, H2, HSED (I), VSED (I)
                  If (MSR) WRITE (6, 9905) JPX, JNX, JCX, JO2NH4, XJC1
                  If (MSR) WRITE (6, 9906) JPO4, JNH4, JNO3, JHS, JSI, &
                 & JCH4AQ, JCH4G, BENSTR
                  If (MSR) WRITE (6, 9907) PO40, PO41, PO42, PO4T2, &
                 & NH40, NH41, NH42, NH4T2
                  If (MSR) WRITE (6, 9908) NO30, NO31, NO32, NO3T2, &
                 & HS1, HS2, HST2
                  If (MSR) WRITE (6, 9909) SI0, SI1, SI2, SIT2

             !
             !stop the calculation
             !
                  Ibreak = I
#       if defined (MULTIPROCESSOR)
                  If (PAR) Then
                     Ibreak = NGID (I)
                  End If
#       endif
                  Write (*,*) 'Oops, zbrent did not converge for box nu&
                 &mber I=', Ibreak
                  Stop

             !
               Else !(IERR==1) !did not find solution within SODMIN and SODMAX
             !
             !If IERR /=0 and /=2 then redo the calculation with minimum SOD and find out the size of error
             !                    also redo the calculation with maximum SOD and find out the size of error
             !
                  If (DFEEDER) Then
                     DF_SOD_TMP = DF_SOD (I)
                  Else
                     DF_SOD_TMP = 0.0
                  End If
             !
                  If (SAV_CALC) Then
                     SAV_SOD_TMP = SEDDOSAV (I)!SAV_SOD(I)
                  Else
                     SAV_SOD_TMP = 0.0
                  End If
             !calculate NO3, NH4, SO4, HS etc using SEDF function with a)minimum SOD, b)maximum SOD
             !
             !
             !calculate ERROR using minimum SOD, basically report the lower bound of error
             !
             !
                  ERROR = SEDF (SODMIN, DF_SOD_TMP, SAV_SOD_TMP, &
                 & TF_SSTATE)
             !
                  If (MSR) WRITE (6, 9889) JDAY, SODMIN, ERROR
                  If (MSR) WRITE (6, 9900) JDAY, CTEMP (I), POP1, POP2, &
                 & POP3
                  If (MSR) WRITE (6, 9901) PON1, PON2, PON3, POC1, &
                 & POC2, POC3
                  If (MSR) WRITE (6, 9902) PO4T2, HST2, SIT2, PSISED
                  If (MSR) WRITE (6, 9903) (JPOP(I, 1)+JPOP(I, &
                 & 2)+JPOP(I, 3)), (JPON(I, 1)+JPON(I, 2)+JPON(I, 3)), &
                 & (JPOC(I, 1)+JPOC(I, 2)+JPOC(I, 3))
             !
                  If (MSR) WRITE (6, 9904) O20, CSOD, DF_SOD_TMP, SOD, &
                 & S, H2, HSED (I), VSED (I)
             !
                  If (MSR) WRITE (6, 9911) CSODHS, CSODCH4, CSOD
                  If (MSR) WRITE (6, 9905) JPX, JNX, JCX, JO2NH4, XJC1
                  If (MSR) WRITE (6, 9906) JPO4, JNH4, JNO3, JHS, JSI, &
                 & JCH4AQ, JCH4G, BENSTR
                  If (MSR) WRITE (6, 9907) PO40, PO41, PO42, PO4T2, &
                 & NH40, NH41, NH42, NH4T2
                  If (MSR) WRITE (6, 9908) NO30, NO31, NO32, NO3T2, &
                 & HS1, HS2, HST2
                  If (MSR) WRITE (6, 9909) SI0, SI1, SI2, SIT2

             !
             !calculate ERROR using maximum SOD, basically report the upper bound of error
             !
                  If (DFEEDER) Then
                     DF_SOD_TMP = DF_SOD (I)
                  Else
                     DF_SOD_TMP = 0.0
                  End If
             !
                  If (SAV_CALC) Then
                     SAV_SOD_TMP = SEDDOSAV (I)!SAV_SOD(I)
                  Else
                     SAV_SOD_TMP = 0.0
                  End If
             !
             !calculate ERROR using maximum SOD
                  ERROR = SEDF (SODMAX, DF_SOD_TMP, SAV_SOD_TMP, &
                 & TF_SSTATE)
             !
                  If (MSR) WRITE (6, 9889) JDAY, SODMAX, ERROR
                  If (MSR) WRITE (6, 9900) JDAY, CTEMP (I), POP1, POP2, &
                 & POP3
                  If (MSR) WRITE (6, 9901) PON1, PON2, PON3, POC1, &
                 & POC2, POC3
                  If (MSR) WRITE (6, 9902) PO4T2, HST2, SIT2, PSISED
             !
                  If (MSR) WRITE (6, 9903) (JPOP(I, 1)+JPOP(I, &
                 & 2)+JPOP(I, 3)), (JPON(I, 1)+JPON(I, 2)+JPON(I, 3)), &
                 & (JPOC(I, 1)+JPOC(I, 2)+JPOC(I, 3))
             !
                  If (MSR) WRITE (6, 9904) O20, CSOD, DF_SOD_TMP, SOD, &
                 & S, H2, HSED (I), VSED (I)
                  If (MSR) WRITE (6, 9911) CSODHS, CSODCH4, CSOD
                  If (MSR) WRITE (6, 9905) JPX, JNX, JCX, JO2NH4, XJC1
                  If (MSR) WRITE (6, 9906) JPO4, JNH4, JNO3, JHS, JSI, &
                 & JCH4AQ, JCH4G, BENSTR
                  If (MSR) WRITE (6, 9907) PO40, PO41, PO42, PO4T2, &
                 & NH40, NH41, NH42, NH4T2
                  If (MSR) WRITE (6, 9908) NO30, NO31, NO32, NO3T2, &
                 & HS1, HS2, HST2
                  If (MSR) WRITE (6, 9909) SI0, SI1, SI2, SIT2


				  WRITE (*, 9889) JDAY, SODMAX, ERROR
                  WRITE (*, 9900) JDAY, CTEMP (I), POP1, POP2, &
                 & POP3
                  WRITE (*, 9901) PON1, PON2, PON3, POC1, &
                 & POC2, POC3
                  WRITE (*, 9902) PO4T2, HST2, SIT2, PSISED
             !
                  WRITE (*, 9903) (JPOP(I, 1)+JPOP(I, &
                 & 2)+JPOP(I, 3)), (JPON(I, 1)+JPON(I, 2)+JPON(I, 3)), &
                 & (JPOC(I, 1)+JPOC(I, 2)+JPOC(I, 3))
             !
                  WRITE (*, 9904) O20, CSOD, DF_SOD_TMP, SOD, &
                 & S, H2, HSED (I), VSED (I)
                  WRITE (*, 9911) CSODHS, CSODCH4, CSOD
                  WRITE (*, 9905) JPX, JNX, JCX, JO2NH4, XJC1
                  WRITE (*, 9906) JPO4, JNH4, JNO3, JHS, JSI, &
                 & JCH4AQ, JCH4G, BENSTR
                  WRITE (*, 9907) PO40, PO41, PO42, PO4T2, &
                 & NH40, NH41, NH42, NH4T2
                  WRITE (*, 9908) NO30, NO31, NO32, NO3T2, &
                 & HS1, HS2, HST2
                  WRITE (*, 9909) SI0, SI1, SI2, SIT2
             !
             !stop the calculation
             !
                  Ibreak = I
#       if defined (MULTIPROCESSOR)
                  If (PAR) Then
                     Ibreak = NGID (I)
                  End If
#       endif
                  Write (*,*) 'Oops, did not find SOD soultion with SOD&
                 &MIN<SOD<SODMAX, SODMIN=', SODMIN, 'SODMAX=', SODMAX, &
                 & 'for box I=', Ibreak

				  Write (*,*) 'SAV_SOD_TMP',SAV_SOD_TMP

                  Stop
             !
               End If
          !
9889           Format (/ 5 X, 'ZBRENT DIAGNOSTICS AT TIME =', F8.3, '(DAYS) F&
              &OR SOD =', F8.4, ' ERROR =', E12.3 /)
9000           Format (/ 5 X, 'ZBRENT FAILURE AT TIME =', F8.3, ' WITH &
              &IERR=', I2 / 5 X, 'IN SEDIMENT SEGMENT IR=', I5 / 5 X, '&
              &WITH SALT, SO40MG     =', 2E10.3 / 5 X, 'DFEED=', F10.3, &
              & ' SFEED=', 3 F11.3 / 5 X, '(SODMIN,SODMAX=', F6.4, &
              & F6.1, ')' / 5 X, 'PROGRAM TERMINATION FOLLOWS DIAGNOSTI&
              &C DUMPS')
9900           Format (/ 1 X, ' TIME,CTEMP,POP1,POP2,POP3' / 8E10.3)
9901           Format (/ 1 X, ' PON1,PON2,PON3,POC1,POC2,POC3' / &
              & 8E10.3)
9902           Format (/ 1 X, ' PO4T2,HST2,SIT2,PSISED' / 8E10.3)
9903           Format (/ 1 X, ' JPOP,JPON,JPOC' / 8E10.3)
9904           Format (/ 1 X, ' O20,CSOD,DFSOD,SOD,S,H2', ',HSED,VSED' &
              & / 10E10.3)
9905           Format (/ 1 X, ' JP,JN,JC,JO2NH4,XJC1' / 8E10.3)
9906           Format (/ 1 X, ' JPO4,JNH4,JNO3,JHS,JSI,JCH4AQ,JCH4G,BEN&
              &STR' / 8E10.3)
9907           Format (/ 1 X, ' PO40,PO41,PO42,PO4T2,NH40,NH41,NH42,NH4&
              &T2' / 8E10.3)
9908           Format (/ 1 X, ' NO30,NO31,NO32,NO3T2,HS1,HS2,HST2' / &
              & 8E10.3)
9909           Format (/ 1 X, ' SI0,SI1,SI2,SIT2' / 8E10.3)
          !
          !
            End If !Done with error reporting for node I


	   ELSE		!Dry places do not do sediment calculation for SOD, simply use the SOD of last time step

			SOD=max(SODTM1S(I),SODMIN)
			S=SOD/max(O20,0.001)

	   ENDIF
	   !
       !******* accumulate temperature for steady state calculation of POC, PON, POP for G3
       !
       !
            If (STEADY_STATE_SED_G3) Then
               ASDTMP (I) = ASDTMP (I) + TEMPD * DLTS
            End If
       !
       !Calculate Silica and PO4 for they do not affect SOD
       !
       !calculate Si in sediments
       !
            K0H1D = 0.
            K0H1P = 0.
            KMC1 = 0.
            K1H1D = S
            K1H1P = 0.
            K2H2D = 0.
            K2H2P = 0.
            J1 = S * SI0
       !
       ! Oxygen dependency of pie1 (L/kg) (DMD book page 156, 391, 571)
       !
            If (O20 < O2CRITSI) Then
               PIE1 = PIE2SI * DPIE1SI ** (O20/O2CRITSI)! DMD  (7.19)
            Else
               PIE1 = PIE2SI * DPIE1SI !DMD (7.18)
            End If
       !
            PIE2 = PIE2SI
       !
       !******* Silica dissolution kinetics
       !
            FD2 = 1. / (1.+M2*PIE2)
       !
            K3 = XKSI * PSISED / (PSISED+KMPSI) * FD2 !
       !
       !particle fraction in overlying water column
            PF = KADSA * SSI (I, KWC) / (1.+KADSA*SSI(I, KWC))
       !
            J2 = XKSI * PSISED / (PSISED+KMPSI) * CSISAT ! &        !
       !
            If (SFEEDER) Then
               J2 = J2 + SA_SED_SF (I)! Suspension Feeders
            End If
!
            Call SEDTSFNL (SI1, SI2, SIT1, SIT2, SI1TM1, SIT2TM1, &
           & TF_SSTATE)
       !
            JSI = 0.0 !S*(SI1-SI0)  !
       !
!
       !
       !calculate PO4 in sediments
       !
            K0H1D = 0.
            K0H1P = 0.
            KMC1 = 0.
            K1H1D = S
            K1H1P = 0.
            K2H2D = 0.
            K2H2P = 0.
       !
       !
            J1 = S * PO40
            K3 = 0.
       !
            PF = KADPO4 * SSI (I, KWC) / (1.+KADPO4*SSI(I, KWC))
       !
            PPO4 = PF * PO4 (I, KWC)!PO4 adsorbed to SSI in overlying water (gP/m^3)
       !and that will get into sediments through SSI settling
       !
            J2 = JPX !                                &!diagenesis of PO4 from POP decay (mgP/m^2/day)
       !+ WSSINETMMD*PPO4                     &!settling of dissolved PO4 adsorbed on to Suspended solids (mgP/m^2/day)
       !+ WSSINETMMD*PIP(I,KWC)                 !settling of particulate inorganic phosphorus (mgP/m^2/day)
       !
            If (SFEEDER) Then
               J2 = J2 + PIP_SED_SF (I)! Suspension Feeders
            End If
       !
       !Salinity dependence of DPIE1PO4 which is  (delta pi PO4,1 ) in eqn (16A.18) of DMD page 191
       !
            If (SAL < SALTSW) Then
               DPIE1PO4 = DPIE1PO4F !use fresh water incremental partitioning coefficient
            Else
               DPIE1PO4 = DPIE1PO4S
            End If
       !
       ! Oxygen dependency of pie1
       !
            If (O20 < O2CRITPO4) Then
               PIE1 = PIE2PO4 * DPIE1PO4 ** (O20/O2CRITPO4)!DMD (16A.19)
            Else
               PIE1 = PIE2PO4 * DPIE1PO4 !DMD (16A.18)
            End If
       !
!
       !
            PIE2 = PIE2PO4 !(L/kg), partitioning coef of PO4
       !
            Call SEDTSFNL (PO41, PO42, PO4T1, PO4T2, PO41TM1, PO4T2TM1, &
           & TF_SSTATE)
       !
       !
!
            JPO4 = S * (PO41-PO40)
       !

       !
       !******* Assign flux-flux results to mod_wqm arrays
       !
            BENDO (I) = - SOD !gO2/m^2/day
            If (SFEEDER) Then
               BENDO (I) = BENDO (I) - SOD_SED_SF (I)!suspension feeders
            End If
       !
       !
            MTVEL (I) = SOD / O20 !mass transfer velocity (m/day)
       !
       !calculate layer 1 thickness (m)
            HSED1 (I) = KL12 * H2 / S !Estimated thickness of layer 1 (H1) (Di Toro Appendix B) Book page 576
       ![m/d ]* [m] /[m/d] ==> [m]
       !
       !BENNH4 is in gN/m^2/day
       !where JNH4 is in mgN/m^2/day
       !
!
            BENNH4 (I) = JNH4 / 1000.
            If (SFEEDER) Then
               BENNH4 (I) = BENNH4 (I) + JNH4_SED_SF (I) / 1000. ! suspension feeders
            End If
       !
       !IDEBUG=162
       !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
       !
            BENNO3 (I) = JNO3 / 1000.
            BENPO4 (I) = JPO4 / 1000.
            If (SFEEDER) Then
               BENPO4 (I) = BENPO4 (I) + JPO4_SED_SF (I) / 1000. ! suspension feeders
            End If
       !
       !IDEBUG=163
       !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
       !
            BENDOC (I) = 0.0_SP
            BENCOD (I) = (JHS+JCH4AQ)		!: gO2/m^2/day
            BENCH4G (I) = JCH4G 			!: gO2/m^2/day
            BENCH4A (I) = JCH4AQ 			!: gO2/m^2/day
       !
       !benthic avaiable silica flux
       !converted from mgSi/m^2/day to gSi/m^2/day (positive into water column)
            BENSA (I) = JSI / 1000.
            If (SFEEDER) Then
               BENSA (I) = BENSA (I) + JSA_SED_SF (I) / 1000. ! suspension feeders
            End If
       !
#if defined(CSOD_CALC)
            CSODflux (I) = CSOD + XJCNO3 !
#endif
!
       !denitrification flux
            ITEMP = 10 * TEMPD + 1
            If (SAL < SALTND) Then
               XAPP1NO3 = GET_ZHTANO3F (TEMPD)
            Else
               XAPP1NO3 = GET_ZHTANO3S (TEMPD)
            End If
       !
            XK2NO3 = GET_ZHTAK2NO3 (TEMPD)!m/day    !: no need to multiply by H2
       !
       !convert from mgN/m^2/day to gN/m^2/day
            BENDEN (I) = (XAPP1NO3*XAPP1NO3*NO31/S+XK2NO3*NO32) / 1000. !gN/m^2/day
       !
       !: This is same as the calculation in Excel code:
       !Denit(1) = (KappaNO3_1 ^ 2 * ThtaNO3 ^ (Tw - 20) / s)
       !Denit(2) = KappaNO3_2 * ThtaNO3 ^ (Tw - 20)
       !JdenitT =Denit(1)*NO31 + Denit(2)*NO32
       !
       !Check units above: XK2NO3*NO32/1000 ~ gN/m^2/d
       !                         *(mgN/m^3)/1000
       !==> XK2NO3 must have unit m/d
       !
       !and KappaNO3_2 also has unit m/d
       !
       !in Excel: KappaNO3_2 is
       !
       !
       !Fluxes due to burial of particulates (mg/m^2/day)
       !
            BURIALN (I) = (PON1+PON2+PON3+NO3T2+NH4T2) * W2 !mgN/m^3*m/d ==> mgN/m^2/d
       !
            BURIALP (I) = (POP1+POP2+POP3+PO4T2) * W2
!
       !
            BURIALC (I) = (POC1+POC2+POC3) * W2
       !
       ! Diagenesis of carbon forms (gC/m^2/day)
       !
            DIAGENC (I) = JCX / 1000. !gC/m^2/day
!
            If (BALGAE_CALC) Then
          !
               Call BA_CALC
          !
          !******* Phosphorus available for benthic algal growth
          !
          !: adjust BENNH4 BENNO3, BENPO4, BENDOC, BENDO etc flux based on benthic algae
               BENNH4 (I) = BENNH4 (I) + BANH4 (I)!(gN/m2/day)
               BENNO3 (I) = BENNO3 (I) + BANO3 (I)
          !
               BENPO4 (I) = BENPO4 (I) + BAPO4 (I)!(gP/m^2/day)
          !
               BENDOC (I) = BENDOC (I) + BADOC (I)!(gC/m^2/day)
               BENDO (I) = BENDO (I) + BADO (I)!(gO2/m^2/day)
          !
          ! COMPUTE EFFECTS OF ALGAL ACTIVITY ON ORGANIC PARTICULATES (MG/M**3)
          !
          !
          ! : acumulate POCi using DTPOCi_BA returned from MOD_BA
               POC1 = POC1 + DTPOC1_BA (I)!(mgC/m^3)
               POC2 = POC2 + DTPOC2_BA (I)!(mgC/m^3)
               POC3 = POC3 + DTPOC3_BA (I)!(mgC/m^3)
          !
          !: accumulate PONi using DTPONi_BA returned from MOD_BA
               PON1 = PON1 + DTPON1_BA (I)!mgN/m^3
               PON2 = PON2 + DTPON2_BA (I)!mgN/m^3
               PON3 = PON3 + DTPON3_BA (I)!mgN/m^3
          !
          !: accumulate POPi using DTPOPi_BA returned from MOD_BA
               POP1 = POP1 + DTPOP1_BA (I)
               POP2 = POP2 + DTPOP2_BA (I)
               POP3 = POP3 + DTPOP3_BA (I)
          !
          !******* Accumulate fluxes for steady-state computation for G3
          !
               If (STEADY_STATE_SED_G3) Then
             !
                  AG3CFL (I) = AG3CFL (I) + DTAG3CFL_BA (I)!1000.*PRB(I)*FRCPHB(3)*BBM(I)*DLTS         !(mgC/m^2/day)
                  AG3NFL (I) = AG3NFL (I) + DTAG3NFL_BA (I)!1000.*PRB(I)*FRNPHB(3)*ANCB*BBM(I)*DLTS    !(mgN/m^2/day
                  AG3PFL (I) = AG3PFL (I) + DTAG3PFL_BA (I)!1000.*PRB(I)*FRPPHB(3)*APCB*BBM(I)*DLTS  !(mgP/m^2/day)
             !
               End If
          !
            End If


#  if !defined (WET_DRY)
			IF(D(I) <=0.0_SP)THEN
#  else
			IF(ISWETN(I) == 0)THEN
#  endif
			   BENNH4 (I) = 0.0_SP
               BENNO3 (I) = 0.0_SP
               BENPO4 (I) = 0.0_SP
               BENDOC (I) = 0.0_SP
               BENDO (I)  = 0.0_SP
			   BENCH4G(I) = 0.0_SP
			   BENCH4A(I) = 0.0_SP
#if defined(CSOD_CALC)
				CSODflux(I) = 0.0_SP  !also set CSOD flux to zero
#endif
			ENDIF

       !******* Total sediment nutrient mass
       !
       !!     SEDMN = SEDMN+(PON1+PON2+PON3+NH4T2+NO3T2)*SFA(I)*H2/1.E6
       !!     SEDMP = SEDMP+(POP1+POP2+POP3+PO4T2)*SFA(I)*H2/1.E6
       !!     SEDMC = SEDMC+(POC1+POC2+POC3)*SFA(I)*H2/1.E6
       !
            SEDMN = SEDMN + (PON1+PON2+PON3+NH4T2+NO3T2) * ART1 (I) * &
           & H2 / 1.E6 !(kgN) ! replaced SFA by ART1
            SEDMP = SEDMP + (POP1+POP2+POP3+PO4T2) * ART1 (I) * H2 / &
           & 1.E6 !(kgP) ! replaced SFA by ART1
            SEDMC = SEDMC + (POC1+POC2+POC3) * ART1 (I) * H2 / 1.E6 !(kgC) ! replaced SFA by ART1
       !
       !                (mgC/m^3) *(m^2)*m /1000/1000 = mgC/1000/1000 = gC/1000 = kgC
       !
       !
       !cycle the new clauclations to n-1 time step for next step marching
       !
       !
!
!#  if !defined (WET_DRY)
!			IF(D(I) > 0.0_SP)THEN
!#  else
!			IF(ISWETN(I) == 1)THEN
!#  endif
				BENSTRTM1S (I) = BENSTR
       !
				SODTM1S (I) = SOD
				JNH4TM1S (I) = JNH4
				JNO3TM1S (I) = JNO3

				JPO4TM1S (I) = JPO4
				JCH4TM1S (I) = JCH4
				JCH4GTM1S (I) = JCH4G
				JHSTM1S (I) = JHS
				JSITM1S (I) = JSI
       !
				PON1TM1S (I) = PON1
				PON2TM1S (I) = PON2
				PON3TM1S (I) = PON3
!
				POC1TM1S (I) = POC1
				POC2TM1S (I) = POC2
				POC3TM1S (I) = POC3
!
       !
				POP1TM1S (I) = POP1
				POP2TM1S (I) = POP2
				POP3TM1S (I) = POP3
       !
				PSISEDTM1S (I) = PSISED !sediment particulate organic silicate
       !
				NH41TM1S (I) = NH41 !dissolved ammonia in layer 1
				NH42TM1S (I) = NH42 !dissolved ammonia in layer 2
				NH4T2TM1S (I) = NH4T2 !total ammonia in layer 2
       !
				NO31TM1S (I) = NO31 !dissolved nitrate in layer 1
				NO32TM1S (I) = NO32 !dissolved nitrate in layer 2
				NO3T2TM1S (I) = NO3T2 !total nitrate in layer 2
       !
				PO41TM1S (I) = PO41 !dissolved PO4 in layer 1
				PO42TM1S (I) = PO42 !dissolved PO4 in layer 2
				PO4T2TM1S (I) = PO4T2 !total PO4 in layer 2
       !
				HS1TM1S (I) = HS1 !dissolved HS in layer 1
				HS2TM1S (I) = HS2 !dissolved HS in layer 2
				HST2TM1S (I) = HST2 !total HS in layer 2
       !
				SI1TM1S (I) = SI1 !dissolved Si in layer 1
				SI2TM1S (I) = SI2 !dissolved Si in layer 2
				SIT2TM1S (I) = SIT2 !total Si in layer 2
       !
				!Added POC1TM1S_SED_DF
				If (DFEEDER) Then
					POC1TM1S_SED_DF (I) = POC1
					POC2TM1S_SED_DF (I) = POC2
				End If
       !
				BFORMAXS (I) = BFORMAX
				ISWBENS (I) = ISWBEN
       !
				!     DFEEDM1S(I)  = DFEED   !moved to mod_df.F
				CH41TM1S (I) = CH41 ! dissolved CH4 in layer 1
				CH42TM1S (I) = CH42 ! dissolved CH4 in layer 2
				CH4T2TM1S (I) = CH4T2 ! total CH4 in sediment layer 2
				SO4T2TM1S (I) = SO4T2 ! SO4
       !
				DIAGN (I) = JNX ! CFC

            If (SAV_CALC) Then
               NH4T2TM1S_SHARE (I) = NH4T2
               PO4T2TM1S_SHARE (I) = PO4T2
		 	   HST2TM1S_SHARE (I) = HST2
            End If


!		 	ENDIF
       !
         End Do
    !
    !***** Assign concentrations to plot variables !
    !IF(MSR)write(*,*)'LBnote vi: POC2TM1S(38)=',POC2TM1S(38),'  CPOC(38,2)=',CPOC(38,2)
         Do I = 1, MLOC

            CPOC (I, 1) = POC1TM1S (I)
            CPOC (I, 2) = POC2TM1S (I)
            CPOC (I, 3) = POC3TM1S (I)
            CPON (I, 1) = PON1TM1S (I)
            CPON (I, 2) = PON2TM1S (I)
            CPON (I, 3) = PON3TM1S (I)
            CPOP (I, 1) = POP1TM1S (I)
            CPOP (I, 2) = POP2TM1S (I)
            CPOP (I, 3) = POP3TM1S (I)
       !
            CNH4 (I) = NH4T2TM1S (I)
            CNO3 (I) = NO3T2TM1S (I)
            CPO4 (I) = PO4T2TM1S (I)
       !
            CPOS (I) = PSISEDTM1S (I)
            CCH4 (I) = CH4T2TM1S (I)
            CSO4 (I) = SO4T2TM1S (I)
            CHS (I) = HST2TM1S (I)
            CSI (I) = SIT2TM1S (I)

         End Do
!
    !
    !***** Take temperature integration step
    !

		IF(.NOT.TF_SSTATE) THEN  !non steady state calculation of sediment temperature
			Do I = 1, MLOC

#if defined (SED_TEMP)

			IF(SEDTEMP_CALC)THEN
# if !defined (WET_DRY)
				IF(D(I) > 0.0_SP)THEN
# else
				IF(ISWETN(I) == 1)THEN
# endif

				!Explicit scheme for calculating sediment temperature by diffusion with water

					!IF(DLT*DIFFT/HSED(I)/HSED(I)>0.5)THEN

					!   WRITE(*,*)'Oops time step does not meet stability criterion for solving sediment temperature'
					!   WRITE(IREPORT,*)'statbility criterion: DLT*DIFFT/HSED(I)/HSED(I)=',DLT*DIFFT/HSED(I)/HSED(I)
					!   STOP 'time step too large for sediment temperature calculation'
					!ENDIF

					!CTEMP(I) = CTEMP(I)+DLTS*86400.0_SP*DIFFT/HSED(I)/HSED(I)          &
					!            *(T(I,KWC)-CTEMP(I))

					!Note here DIFFT has unit m^2/s so convert DLTS to seconds

				!Implicit scheme for calculating sediment temperature by diffusionw with water

					CTEMP(I) = (CTEMP(I)+DLTS*86400.0_SP*DIFFT/HSED(I)/HSED(I)*T(I,KWC))	&
								/(1+DLTS*86400.0_SP*DIFFT/HSED(I)/HSED(I))
				!Note here DIFFT has unit m^2/s so use 86400.0 to convert to sec

				ELSE

					CTEMP(I)=CTEMP(I) !hold it not changing??
					!should use air temperature here

					!
					!USE IATBOT(I,KBM1) as input heat flux
					!How about latent heat loss?
					!


				ENDIF
			ELSE
				CTEMP(I) = T(I,KWC) !no calculation, simply use water colcumn temperature
			ENDIF
#else

# if !defined (WET_DRY)
				IF(D(I) > 0.0_SP)THEN
# else
				IF(ISWETN(I) == 1)THEN
# endif
				!Wet use water temperature
					CTEMP(I) = T(I,KWC)  !simply use water temperature for sediment temperature
				ELSE
					!Dry, or should use air temperature here
					CTEMP(I)= CTEMP(I)   !hold it with no change
				ENDIF
#endif
			End Do

		ELSE !Quasi-Steady state calculation of sedimen temperature
			DO I=1,MLOC
				CTEMP(I)=T(I,KWC)		!use quasi-steady state calculation for sediment temperature
										!i.e. sediment temperature is always same as water temperature
			ENDDO
		ENDIF

        Return
      End Subroutine SED_CALC
  !
      Subroutine SED_DIAGENESIS_G3 ()
    !
         Integer :: I

    !
         Do I = 1, MLOC
            CPOC (I, 3) = AG3CFL (I) / &
           & (KCDIAG(3)*DCTHTA(3)**(ASDTMP(I)-20.)*HSED(I)+VSED(I))
            CPON (I, 3) = AG3NFL (I) / &
           & (KNDIAG(3)*DNTHTA(3)**(ASDTMP(I)-20.)*HSED(I)+VSED(I))
            CPOP (I, 3) = AG3PFL (I) / &
           & (KPDIAG(3)*DPTHTA(3)**(ASDTMP(I)-20.)*HSED(I)+VSED(I))
         End Do
    !
         Return
      End Subroutine SED_DIAGENESIS_G3
  !
  !
  !********************************************************************************
  !**                          F U N C T I O N   S E D F                         **
  !********************************************************************************
  !
      Function SEDF (SOD1, DF_SOD_TMP, SAV_SOD_TMP, TF_SSTATE) Result &
     & (SEDFOUT)
    !
    !
    !Function that calculates the imbalance of equation
    !
    ! SOD=CSOD(SOD1)+NSOD(SOD1) + DFSOD(SOD1) + SAVSOD(SOD1)
    !
    !where SOD1 is given SOD, and CSOD is calculated using CSOD_HS , CSOD_CH4
    !NSOD is calculated as NSOD_NH4
    !
    !     CSOD = CSOD_HS  for salt > SALTSW (salt water system)
    !or
    !     CSOD = CSOD_CH4 for salt < SALTSW (fresh water system)
    !
    !     NSOD = oxygen consumed during nitrification of NH4 to NO3 (DMD book eqn (1.25) reaction)
    !
    !Imbalance = SOD-SOD1
    !
    !Aim of the program is to find SOD=SOD1 so that Imbalance is zero, this is done through zbrent() function
    !
    !
    !
    !First of all,give any SOD1 value (gO2/m^2/day), calculate NO3, NH4, SO4, CH4, HS based on XJC1
    !
    !Then caluclate CSOD and NSOD, and put together CSOD+NSOD+DFSOD+SAVSOD
    !
    !Finally calculate the Imbalance SEDF=SOD-SOD1=(CSOD+NSOD+DFSOD+SAVSOD) - SOD1
    !
    !
         Implicit None
    !
         Save
         Logical :: TF_SSTATE
         Real (SP) :: SOD1, DF_SOD_TMP, SAV_SOD_TMP, SEDFOUT !returned error SOD-SOD1 for checking convergence (gO2/m^3/d)
    !
         Real (SP) :: DCH4T2, DHST2 !rate of HS generation in layer 2 (gO2/m^3/d * m/d = gO2/m^2/d)
    !
         Real (SP) :: SAVSOD_SED, DFSOD_SED
    !
    ! KURT GLAESEMANN Changed H to HH to avoid name conflict with module
         Real (SP) AD (4, 4), BX (4), G (2), HH (2, 2)
    !
         Real (SP) DBLSO41, DBLSO42, RA0, RA1, RA2, R1, R2, DISC, SN1
    !
         Real (SP) CSODMAX, SECH_ARG !argument of hyperbolic secant function (no unit)
    !
    !***** Compute the NH4, NO3, and SOD fluxes
    !
    !surface diffusion velocity (m/day)
    !
         S = SOD1 / O20 !(mgO2/m^2/day) / (mgO2/m^3) ==> (m/day)
    !
    !***** Ammonia flux
    !
    !layer integrated reaction rates (m/d)
         K0H1P = 0.
         K1H1P = 0.
    !
         K2H2D = 0.
         K2H2P = 0.
    !
         If (KMNH4 /= 0.) Then !KMNH4 is K_{NH4} )
            K0H1D = XAPPNH4 ** 2 / S * KMNH4 * THTAKMNH4 ** (STP20) * &
           & (O20/(2*KMNH4O2+O20))!O20/(KMNH4O2+O20) is F_{O2} in 
       ! (3-32-8)
            K1H1D = S ! (3-32-6)
         Else
       !
            K0H1D = 0. ! (3-32-7)
            K1H1D = XAPPNH4 ** 2 / S * (O20/(2*KMNH4O2+O20)) + S ! (3-32-5)
         End If
    !
         J1 = S * NH40 !(m/d) * (mgN/m^3) ==> mgN/m^2/day 
         K3 = 0. ! (4-8)
         J2 = JNX ! (4-9)
         PIE1 = PIENH4
         PIE2 = PIENH4
         KMC1 = KMNH4 * THTAKMNH4 ** (STP20)!Michael=Mention half constant for layer 1 reaction, DMD eqn (3.30)
    !
         Call SEDTSFNL (NH41, NH42, NH4T1, NH4T2, NH41TM1, NH4T2TM1, &
        & TF_SSTATE)
    !
         If (TF_SSTATE) Then !for steady state make sure NH41TM1=NH41 just calculated
            NH41TM1 = NH41 !: this is important for steady state
         End If
    !
         JNH4 = S * (NH41-NH40)!(m/d)*(mgN/m^3) ==> mgN/m^2/d
    !
    !***** Oxygen consumed by nitrification (NSOD)
    !
    !
         A1 = 0.0045714 ! notes (4-14-1) see DMD book eqn (1.25) and p206,p 189
    !
    !calculate NSOD (gO2/m^2/d)
    !
         If (KMNH4 /= 0.) Then
            JO2NH4 = A1 * K0H1D * NH41 / &
           & (KMNH4*THTAKMNH4**(STP20)+NH41TM1)
       ! (3-32-14) (4-14)
       !Note oxygen regulation is already in K0H1D, DMD book eqn (9.12)(9.14)
         Else
            JO2NH4 = A1 * (K1H1D-S) * NH41 ! notes (4-14-2) (3-32-15)
         End If
    !(mgO2/mgN)* (m/d) *(mgN/m^3)/1000 ==> mgO2/m^2/d/1000 ==> gO2/m^2/d
    !
    !
    !***** Denitrification
    !
         K0H1D = 0. ! (3-32-7)
         K0H1P = 0. ! (3-32-9)
         KMC1 = 0.
         K1H1D = XAPP1NO3 ** 2 / S + S ! (3-32-5)
         K1H1P = 0. ! (3-32-10, should not be zero but equivalent to zero since FP1 = 0 for NO3)
    !    it is multiplied by FP1 anyway
         K2H2D = XK2NO3 ! (3-36-7)
         K2H2P = 0. ! (3-36-8)
         If (KMNH4 /= 0.) Then
            J1 = S * NO30 + XAPPNH4 ** 2 / S * KMNH4 * THTAKMNH4 ** &
           & (STP20) / (KMNH4*THTAKMNH4**(STP20)+NH41TM1) * &
           & (O20/(2*KMNH4O2+O20)) * NH41 ! (4-40-1) with KMNH4 !DMD(3.34)
       !
         Else
            J1 = S * NO30 + XAPPNH4 ** 2 / S * (O20/(2*KMNH4O2+O20)) * &
           & NH41 ! (4-40-1) with no KMNH4
         End If
         K3 = 0.
         J2 = 0.
         PIE1 = 0.
         PIE2 = 0.
    !
         Call SEDTSFNL (NO31, NO32, NO3T1, NO3T2, NO31TM1, NO3T2TM1, &
        & TF_SSTATE)
         JNO3 = S * (NO31-NO30)! (4-44)   !fixed to -0.011 in experiment 'PSMdiag_fixedJno3' by LB
!
    !
    !***** Sulfide/methane oxidation
    !
         AO2N = 2.857142857142857 !AO2N= 2.85714 is alpha_O2,C and alpha_C,N in 
    !alpha_O2C = 2.667gO2/gC 
    !alpha_C,N = 1.071 gC/gN 
    !
    != 2.0*16/12.0*10/8.0*12/14.0
    !
    !
    !==> AO2N=2.667*1.071 = 2.85714 gO2/gN
    !
    !
         XJCNO31 = AO2N * XAPP1NO3 ** 2 / S * NO31 / 1000.0 !layer 1 denitrification flux (gO2/m^2/d) 
         XJCNO32 = AO2N * XK2NO3 * NO32 / 1000.0 !layer 2 denitrification flux (gO2/m^2/d) 
    !mgO2/mgN  * m/d * mgN/m^3 /1000= mgO2/m^2/d /1000 = gO2/m^2/d
    !
    !
    !
    !***** Add the aerobic and first anaerobic layer to keep mass balance
    !
         XJCNO3 = XJCNO31 + XJCNO32 !JO2NO3T in Greg's excel version
    !
    !***** Convert carbon diagenesis flux to O2 units
    !
         XJC1 = AMAX1 (2.667E-3*JCX-XJCNO3, 1.0e-10)!, note here XJC1 is in gO2/m^2/day as 2.667 converts from C to O2
    !and 1.0E-3 converts from mgC to gC
    ! notes eqn(4-57) + eqn (4-56) for both layer 1 and layer 2
    !
    !XJC1 is JH2S,1+JH2S,2 
    !
    !
#if(1) /* James L. Martins method of calculating HS and CH4 : solve HS if SAL<SALTSW using DMD chapter 10 */
         If (SAL >= SALTSW) Then
       !
       !dissolved and particle fractions for HS
       !
            PIE1 = PIE1HS
            PIE2 = PIE2HS
       !
            KMC1 = 0.0
            K0H1D = 0. ! (3-32-7)
            K0H1P = 0. ! (3-32-9)
       !
       !   K1H1D=XAPP1HSD**2/S*(O20/KMHSO2) + S
       !   K1H1P=XAPP1HSP**2/S*(O20/KMHSO2)
       !
            K1H1D = XAPP1HSD ** 2 / S * (O20/(KMHSO2)) + S ! notes (4-58)
       ! (3-32-5)
            K1H1P = XAPP1HSP ** 2 / S * (O20/(KMHSO2))! (3-32-10)
       !
            K2H2D = 0. ! (3-36-8)
            K2H2P = 0. ! (3-36-8-1)
       !
            J1 = 0. ! (4-52) J1 term with J1=S*fd0*HST0 + JH2ST,1, where HST0=0, JH2ST,1 is put in J2
       !instead
       !
            K3 = 0.0 ! (3-36-6)
       !
            J2 = XJC1 !gO2/m2/day, carbon diagenesis in layer 2 in oxygen units after discounting diaggenesis
       ! (4-57) + (4-56) with JH2ST,1 and JH2ST,2 are put in JH2ST,2 together and ignoring
       !layer 1
       !
            Call SEDTSFNL (HS1, HS2, HST1, HST2, HS1TM1, HST2TM1, &
           & TF_SSTATE)
       !
            JHS = S * (HS1-HS0)!gO2/m^2/day (flux of dissolved HS to overlying water column)
       ! (4-61)
       !
            FD1 = 1. / (1.+M1*PIE1)
            FP1 = M1 * PIE1 / (1.+M1*PIE1)
            FD2 = 1. / (1.+M2*PIE2)
            FP2 = M2 * PIE2 / (1.+M2*PIE2)
       !
       !
       !CSOD due to HS oxidation in layer 1 limited by avaiable O2 in overlying water
       ! (m/d)*(gO2/m^3)*(1) ==> unit is gO2/m^2/d
       !
       !
       ! fix slight error, should use XAPP1HSD here
       !            CSODHS= (  XAPP1HSP**2/S*FD1*HST1                 &!
       !
            CSODHS = (XAPP1HSD**2/S*FD1*HST1+XAPP1HSP**2/S*FP1*HST1) * &
           & (O20/KMHSO2)!limitation by avaiable O20
       ! (4-62)
       !
       !No methane generation, no CSOD due to methane oxidation to CO2
       !
            CSODCH4 = 0.d0
            JCH4AQ = 0.d0
            JCH4G = 0.d0
            JCH4 = 0.0d0 !
       !
            CSOD = CSODCH4 + CSODHS
       !
       !
       !dummy CH4 solution as no CH4 will be generated??: should we retrain them to be
       !same as previous time step? consider SAL<SALTSW is altenatingly satisfied.
       !: may need to set CH40 to zero as well?
       !
       !
            CH4T2 = 0
            CH4T1 = 0
            CH41 = 0
            CH42 = 0
       !
         Else !methane generation with salt is < saltsw
       !
       !
            CSODMAX = Min (Sqrt(2.d0*KL12*CH4SAT*XJC1), XJC1)! 4-65)
       !
            SECH_ARG = XAPPCH4 / S ! (4-70), DMD book (10.39)
       !
            CSODCH4 = 0.d0
       !
       !
            If (.False.) Then !--the analytical solution
               If (SECH_ARG < 400) Then
                  CSODCH4 = CSODMAX * &
                 & (1.0-(2.0/(Exp(SECH_ARG)+Exp(-SECH_ARG))))!CSODMAX(1-sech(sech_arg))
             !where hyperbolic secant is : sech(x)=2/[exp(x)+exp(-x)]
               Else
                  CSODCH4 = CSODMAX
               End If
          !
            Else !--the numerical solution
          !
               CH41 = (CSODMAX+S*CH40) / (XAPPCH4**2/S+S)!analytical solution of CH4(1)
               CSODCH4 = (XAPPCH4**2/S) * CH41
            End If
       !
            JCH4AQ = CSODMAX - CSODCH4 !DMD book (10.35a)
            JCH4G = XJC1 - JCH4AQ - CSODCH4 !DMD book (10.33)
       !
            CH41 = (CSODMAX+S*CH40) / (XAPPCH4**2/S+S)!analytical solution of CH4(1)
       !
       !
            JCH4 = JCH4AQ !S*(CH41-CH40)
       !
       !
       !dummy HS solution as no HS will be generated??: should we retrain them to be
       !same value as previous time step? Consider SAL>SALTSW condition is alternatingly satisfied.
       !May need to set HS0 to zero as well?
       !
       !
            HST1 = 0
            HST2 = 0
            HS1 = 0
            HS2 = 0
       !
            CSODHS = 0.d0
            JHS = 0.d0
       !
            CSOD = CSODHS + CSODCH4
       !
         End If
    !
    !do not solve for SO4 at all (make all SO4 terms zero)
    !Need to set SO40 to zero as well ?
         SO4T1 = 0
         SO4T2 = 0
         SO41 = 0
         SO42 = 0
    !
#else /*solve HS, SO4 and CH4 together following DMD book chapter 11*/
    !
         Call GET_SO4_HS_CH4_CSOD_FLUXES (TF_SSTATE)
    !
#endif
    !
    !**** Volumetric methane and total gas flux (L/m2-d)
         VJCH4G = 22.4 / 64. * JCH4G !JCH4G is mgO2/m^2/day
         SOD = CSOD + JO2NH4 ! (4-9-14)

         DFSOD_SED = DF_SOD_TMP ! deposit feeders
         SOD = SOD + DFSOD_SED ! gO2/m^2/day

    !include SOD due to SAV root respiratioon
    !IF(SAVC_ALC)THEN
         SAVSOD_SED = SAV_SOD_TMP
         SOD = SOD + SAVSOD_SED !gO2/m^2/day
    !ENDIF
    !
         SEDFOUT = SOD - SOD1 !calculate error in SOD and SOD1
    !
         Return
      End Function SEDF
  !
  !
  !********************************************************************************
  !**                        F U N C T I O N   Z B R E N T                       **
  !********************************************************************************
  !
      Function ZBRENT (IERR, DF_SOD_TMP, SAV_SOD_TMP, TF_SSTATE, &
     & SODold) Result (ZBOUT)
    !
    !Trying to solve the equation : SEDF[SOD]== SOD-(CSOD+NSOD)[SOD] =0
    !by guessing SOD and then call SEDF, here [] means argument to function
    !and and eventually converge to a value within the range of SODMIN and SODMAX
    !
    !In order to solve for f(SOD)=0, i.e. SEDF[SOD]=0, with f(SOD)=SOD-CSOD(SOD)+NSOD(SOD)
    !   where SOD is unknown, CSOD and NSOD are functions of SOD as well, this program ZBRENT
    !takes guess of what SOD is by using SODMIN and SODMAX initially, then narrow down the
    !search successively by moving B point, [SODMIN, A, B, SODMAX],
    !until A, B gets closer and closer with f(A),f(B) ~=0
    !
    !zbrent finds the root of a function without knowing the derivative
    !http://en.wikipedia.org/wiki/Brent%27s_method
    !
    !
    !ZBRENT returns IERR and value of SOD as ZBOUT
    !
    !if IERR = 0   good, find soultion and SOD converges and the solution is in ZBOUT=B
    !   IERR = 1    did not find solution of SOD within [SODMIN, SODMAX] range
    !   IERR = 2    did not find the solution within the IMAX interations, last point B is saved in ZBOUT
    !
    !
         Implicit None
         Logical :: TF_SSTATE
         Integer, Parameter :: IMAX = 500
         Real (SP), Parameter :: EPS = 3.E-10, TOL = 1.E-6
    !
         Integer :: IERR, I !                by ZBRENT

         Real (SP) :: ZBOUT, TOL1 !                            SODMAX,       &!maximum SOD (gO2/m^2/day)
    !     ZBOUT,      &
    !                      SODMIN,       &!minimum SOD (gO2/m^2/day) !moved to module definition
    !                      SODMAX,       &!maximum SOD (gO2/m^2/day)
    !                     TOL1
    !
         Real (SP) :: A, B, C, D, E, P, Q, R, XM
         Real (SP) :: FA, FB, FC
         Real (SP) :: DF_SOD_TMP, SAV_SOD_TMP
         Real (SP) :: SODold !previous time step SOD
    !

         ZBOUT = 0.0
    !***** Initialize upper and lower limits for solution
    !
         IERR = 0
!
         A = SODMIN
         B = SODMAX
    !
         FA = SEDF (A, DF_SOD_TMP, SAV_SOD_TMP, TF_SSTATE)!sediment flux calculation (SOD) using minimum SOD
         FB = SEDF (B, DF_SOD_TMP, SAV_SOD_TMP, TF_SSTATE)!sediment flux calculation (SOD) using maximum SOD
    !
    !***** Root must bracket ZBRENT
    !
         If (FB*FA > 0.) Then !if f(A)*f(B)>0, then no solution within [A B] bracket
            IERR = 1
            Return
         End If

         C = B ! followed JLM code and initialize C as B instead of (A+B)/2
    !
         E = B - A !
    ! following JLM, no need to initiate E here
    !
    !start with FC = FB, i.e. searching from SODMAX end towards SODMIN
         FC = FB
    !
    !iterate IMAX times until found solution or failed
         Do I = 1, IMAX
       !
       !
       !     IF(MSR)WRITE(*,*) 'ZBRENT=3, IMAX=',IMAX,'I=',I
            If (FB*FC > 0.) Then !If f(B)f(C) >0 then start lookfing for A side instead
               C = A
               FC = FA
               D = B - A !D keeps the search window width = B-A
               E = D !E also keeps a copy of D
            End If
       !
       !     IF(MSR)WRITE(*,*) 'ZBRENT=4, IMAX=',IMAX,'I=',I
            If (Abs(FC) < Abs(FB)) Then !if f(C) is closer to zero than f(B) then
          !replace search widnow [A B] by new window [B C]
          !and again start guessing root from the C side
          !
          !
          !                  ~
          !                /  \
          !before ------A(C)-------B-----
          !               /     \
          !         .- ~ *          `_
          !              f(C)        \
          !                        * f(B)
          !
          !               |<---------->|
          !                    E
          !
          !
          !                  ~
          !                /  \
          !after  -------(B)----(C)A-----
          !               _^     \
          !         .- ~ *          `_
          !              f(B)        \
          !                        * f(A),f(C)
          !
          !               |<------->|
          !                  XM*2        !new window width
               A = B
          !
          !       IF(MSR)WRITE(*,*) 'ZBRENT=71, IMAX=',IMAX,'I=',I
          !
               B = C
          !       IF(MSR)WRITE(*,*) 'ZBRENT=72, IMAX=',IMAX,'I=',I
          !
               C = A
          !       IF(MSR)WRITE(*,*) 'ZBRENT=73, IMAX=',IMAX,'I=',I
          !
               FA = FB
          !       IF(MSR)WRITE(*,*) 'ZBRENT=74, IMAX=',IMAX,'I=',I
          !
               FB = FC
          !       IF(MSR)WRITE(*,*) 'ZBRENT=75, IMAX=',IMAX,'I=',I
          !
               FC = FA
          !       IF(MSR)WRITE(*,*) 'ZBRENT=76, IMAX=',IMAX,'I=',I
          !
            End If
       !
       !
       !now we have |f(A)| greater than |f(B)|, which is either of the following two cases
       !
       !case 1:
       !
       !                  ~
       !                /  \
       !       ------A(C)-------B--------------->x
       !                /      \
       !               /       `\* f(B)
       !            ~ *              \
       !          /                `
       !
       !               |<---------->|
       !                    E
       !
       !             |<--------->|
       !                  -XM*2
       !
       !
       !case 2:
       !
       !                  ~
       !                /  \
       !      -------(B)----(C)A---------------->x
       !               _^     \
       !         .- ~ *          `_
       !              f(B)        \
       !                        * f(A),f(C)
       !
       !               |<------->|
       !                  XM*2        !new window width
       !
       !
       !
       !     IF(MSR)WRITE(*,*) 'ZBRENT=5, IMAX=',IMAX,'I=',I
            TOL1 = 2. * EPS * Abs (B) + 0.5 * TOL !update tolerance value based on size of B
            XM = 0.5 * (C-B)!XM records half width of new window [B C]
       !
            If (Abs(XM) <= TOL1 .Or. FB == 0.) Then !if widnow is very samall or FB is already zero
          !then B must be the solution to f(x)=0
               ZBOUT = B
          !
          !ZBOUT = (ZBOUT+SODold)/2.d0
          !
          !       IF(MSR)WRITE(*,*) 'ZBRENT=6, IMAX=',IMAX,'I=',I
               Return
            End If
       !
       !
            If (Abs(E) >= TOL1 .And. Abs(FA) > Abs(FB)) Then !If window XM is still too large
          !then check if previoius widnow E is greater
          !than tolreance and if f(A) is farther away
          !from zero than f(B)
          !
          !if yes, then calculate the slope of ratio FB/FA
          !which should have |f(B)/f(A)|<1
          !
               S = FB / FA
          !
               If (A == C) Then !if A and C are the same then set P to scaled down size of window
                  P = 2. * XM * S !by S, and have Q as ratio shortened
                  Q = 1. - S
               Else
                  Q = FA / FC !if A and C are not the same, e.g. C is with B, and f(C)*f(B)
             !
             !            ----A-------(C)B-----
             !
             !then have Q=f(A)/f(C), R = f(B)/f(C)
             !
             !
             !
                  R = FB / FC
                  P = S * (2.*XM*Q*(Q-R)-(B-A)*(R-1.))
                  Q = (Q-1.) * (R-1.) * (S-1.)
               End If
          !
          !recalculate P and Q
          !
               If (P > 0.) Q = - Q
               P = Abs (P)
          !
               If (2.*P < Min(3.*XM*Q-Abs(TOL1*Q), Abs(E*Q))) Then
                  E = D
                  D = P / Q
               Else
                  D = XM
                  E = D
               End If
            Else
               D = XM
               E = D
            End If
       !
       !     IF(MSR)WRITE(*,*) 'ZBRENT=8, IMAX=',IMAX,'I=',I
       !
       !set A and f(A) to B point
       !
            A = B
            FA = FB
       !
       !
       !move B by half window size XM
       !and E still keeps a copy of windows size D
       !
       !
            If (Abs(D) > TOL1) Then
               B = B + D
            Else
               B = B + SIGN (TOL1, XM)
            End If
       !     IF(MSR)WRITE(*,*) 'ZBRENT=9, IMAX=',IMAX,'I=',I
       !
       !recalculate f(B)
       !
            FB = SEDF (B, DF_SOD_TMP, SAV_SOD_TMP, TF_SSTATE)
       !
       !     IF(MSR)WRITE(*,*) 'ZBRENT=10, IMAX=',IMAX,'I=',I
       !
         End Do
    !
    ! do loop of IMAX interation finished and still not converging
    !save the current results in ZBOUT
         IERR = 2
         ZBOUT = B
    !ZBOUT = (ZBOUT+SODold)/2.d0   !, use gp method to average with previous step
    !
    !   IF(MSR)WRITE(*,*) 'ZBRENT=11, IERR=',IERR,'B=',B
    !
         Return
      End Function ZBRENT
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   S E D T S F N L                   **
  !********************************************************************************
  !
  !Time series solution of sediment diagenesis
  !
  !
      Subroutine SEDTSFNL (C1S, C2S, CT1S, CT2S, C1TM1S, CT2TM1S, &
     & STEADYSTATE)
    !inputs: C1TM1S         !total concentration at last time step
    !        CT2TM1S      !total concentration at last time step
    !
    !outputs: C1S        !dissolved concentration in layer 1
    !         C2S        !dissolved concentration in layer 2
    !         CT1S        !total concentration in layer 1
    !         CT2S        !total concentration in layer 2
    !
         Implicit None
         Real (SP) :: C1TM1S, CT2TM1S
         Logical :: STEADYSTATE
    !
         Real (SP) :: C1S, C2S, CT1S, CT2S
         Real (SP) :: A11, A12, A21, A22, B_1, B_2
         Real (SP) :: F12, F21, XK0, XK1, XK2, DELTA
    !
    !***** Initialize constants
    !
         FD1 = 1. / (1.+M1*PIE1)
         FP1 = M1 * PIE1 / (1.+M1*PIE1)
         FD2 = 1. / (1.+M2*PIE2)
         FP2 = M2 * PIE2 / (1.+M2*PIE2)
    !
         F12 = W12 * FP1 + KL12 * FD1 ! (3-32-1-0)
    !
         F21 = W12 * FP2 + KL12 * FD2 ! (3-36-2)
    !
    !***** Evaluate the MM term at time level t-1
    !
         If (KMC1 /= 0.) Then
       !
            XK0 = (K0H1D*FD1+K0H1P*FP1) / (KMC1+C1TM1S)! (3-32-4)
       !Note C1TM1S is dissolved fraction in layer 1
         Else
            XK0 = 0. ! (3-32-3)
         End If
    !
         XK1 = XK0 + K1H1D * FD1 + K1H1P * FP1
         XK2 = K2H2D * FD2 + K2H2P * FP2
    !
    !Excel
    !If steadystate Then
    !  a22 = -fd2 * KL12 - fp2 * w12 - w2
    !  b2 = -Jn
    !
    !Else  !unsteady state
    !  a22 = -fd2 * KL12 - fp2 * w12 - w2 -H2 / tc
    !  b2 = -Jn - H2 / tc * NH3Tp(2)
    !End If
    !
    !
         A11 = - F12 - XK1 - W2 !eqn(3-32)
         A21 = F12 + W2 !eqn(3-35)
         A12 = F21 !eqn(3-33)
         B_1 = - J1 !eqn(3-34)
    !
         If (STEADYSTATE) Then
            A22 = - F21 - XK2 - W2 - K3 !eqn(3-36)
            B_2 = - J2 !eqn(3-37)
         Else
            A22 = - F21 - XK2 - W2 - K3 - H2 / DLTS !eqn(3-36)
            B_2 = - J2 - H2 / DLTS * CT2TM1S !eqn(3-37)
         End If
    !
    !***** Solve the 2x2 set of linear equations
    !
         DELTA = A11 * A22 - A12 * A21
         If (DELTA == 0.) Then
            Print *, 'DELTA is singular: A11,A12,A21,A22'
            Write (*,*) F12, XK1, W2
            Write (*,*) W12, FP1, KL12, FD1
            Write (*,*) XK0, K1H1D, FD1, K1H1P, FP1
            Print *, A11, A12, A21, A22
            Stop
         End If
    !
    !
    !***** Assign results
    !
         CT1S = (B_1*A22-B_2*A12) / DELTA !Total concentration in layer 1
         CT2S = (B_2*A11-B_1*A21) / DELTA !Total concentration in layer 2
    !
         C1S = FD1 * CT1S !dissolved concentration in layer 1
         C2S = FD2 * CT2S !dissolved concentration in layer 2
    !
    !
         Return
      End Subroutine SEDTSFNL
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   S E D S S F N L                   **
  !********************************************************************************
  !steady state solution of sediment diagenesis
      Subroutine SEDSSFNL (C1, C2, C2AV, CT1, CT2, CT2AV, C1TM1, &
     & CT2TM1, ITYPE)!TM1 ~ t-1
         Implicit None
         Integer ITYPE
         Real (SP) :: C1, C2, C2AV, CT1, CT2, CT2AV, C1TM1, CT2TM1
    !
    !  This subroutine translates between SEDTSFNL and SEDSSFNL
    !
         If (ITYPE .Eq. 1) Then
            Call SEDTSFNL (C1, C2, CT1, CT2, C1TM1, CT2TM1, .True.)!call time series solution with steadystate flag
            C2AV = C2
            CT2AV = CT2
         Else
            Stop 'Oops, you are not supposed to call this function with&
           & ITYPE/=1'
         End If
      !
         Return
      End Subroutine SEDSSFNL
    !
      Subroutine GET_SO4_HS_CH4_CSOD_FLUXES (TF_SSTATE)
      !**** *******************************************************************************
      !**** New code for methane formation.  CH4 starts forming                           *
      !**** once all sulfate is used up. (Chapter 11 of DMD book instead of Chapter 10)   *
      !**** *******************************************************************************
      !
         Logical :: TF_SSTATE !Type of integration FALSE - time intergration from current step to next step
      !                    TRUE  - steady state solution
      !**** Sulfide and sulfate in O2 equivalents
      !
      !
         Real (SP) :: DCH4T2, DHST2 !rate of HS generation in layer 2 (gO2/m^3/d * m/d = gO2/m^2/d)
      !
      !   REAL(SP) :: SAVSOD_SED,        &
      !               DFSOD_SED
         Real (SP) AD (4, 4), BX (4), G (2), HH (2, 2)
      !
         Real (SP) DBLSO41, DBLSO42, RA0, RA1, RA2, R1, R2, DISC, SN1
      !
         KMC1 = 0.0
         K0H1D = 0. ! (3-32-7)
         K0H1P = 0. ! (3-32-9)
      !
      !
         K1H1D = XAPP1HSD ** 2 / S * (O20/(KMHSO2)) + S ! notes 4-58 with O20 downstairs
      ! (3-32-5)
         K1H1P = XAPP1HSP ** 2 / S * (O20/(KMHSO2))! (3-32-10)
      !
         K2H2D = 0. ! (3-36-8)
         K2H2P = 0. ! (3-36-8-1)
      !
         J1 = 0. ! (4-52) J1 term with J1=S*fd0*HST0 + JH2ST,1, where HST0=0, JH2ST,1 is put in J2
      !instead
      !
         K3 = 0.0 ! (3-36-6)
      !
         J2 = XJC1 !gO2/m2/day, carbon diagenesis in layer 2 in oxygen units after discounting diaggenesis
      !
         PIE1 = PIE1HS
         PIE2 = PIE2HS
      !
      !**** Set KL12 using HSO4
      !
         ITEMP = 10. * TEMPD + 1
      !
         DDSO4 = KL12NOM * H2 / 2.0 

      !
#if(0)  /* !: not very clear on how the quadratic equation for SO42 is derived, based on WASP cand ICM code */
      !
      !
      !*
      !* overall it is a combination of HS and SO4 equations, with HS elimited to arrive two linear equations
      !* for SO41 and SO42, and then according to DMD book pp.232 eqn(11.27c) a third equation is needed for solving
      !* the third unknown HSO4
      !*
      !* However the following statement calculates HSO4 directly based on DMD book pp.229 eqn (11.17)
      !* since HSO4 already calculated, then it is hard to understand why we need the third equation (11.27c) anymore
      !* and why we need to derive a quadratic equation for SO42
      !*
      !
      !   HSO4  =SQRT(2.*DDSO4*SO40*H2/XJC1/)                !
      !wihich is 1/2 *********
      !Well, it depends on how DDSO4 is defined
      !
         HSO4 = Sqrt (2.*DDSO4*SO40*H2/(XJC1*(1./2.0)))! included 1/2.0 here which is alpha_SO4,C in eqn (4-142)
      !where XJC1 already included alpha_O2,C = 2.667gO2/gC
      !i.e. 
      !
      !
      !**** No deeper than H2
         If (HSO4 > H2) HSO4 = H2
      !
         KL12SO4 = KL12 * H2 / HSO4 !KL12=DD/H2 ==> DD= KL12*H2
      !KL12SO4 = DD/HSO4
      !==> KL12SO4 = KL12*H2/HSO4  (see)
      !
         IDEBUG_SED = IDEBUG_SED + 1
      !IF(MSR)THEN
      !        WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
      !        WRITE(*,*)'K0H1D=',K0H1D
      !        WRITE(*,*)'K0H1P=',K0H1P
      !        WRITE(*,*)'K1H1D=',K1H1D
      !        WRITE(*,*)'K2H2D=',K2H2D
      !        WRITE(*,*)'K1H1P=',K1H1P
      !        WRITE(*,*)'K2H2P=',K2H2P
      !        WRITE(*,*)'J1=',J1
      !        WRITE(*,*)'K3=',K3
      !        WRITE(*,*)'J2=',J2
      !        WRITE(*,*)'PIE1=',PIE1
      !        WRITE(*,*)'PIE2=',PIE2
      !        WRITE(*,*)'DDSO4=',DDSO4
      !        WRITE(*,*)'HSO4=',HSO4
      !        WRITE(*,*)'KL12SO4=',KL12SO4
      !        READ(*,*)
      !ENDIF
      !
      !**** Fractions and overall decay reaction velocity
         FD1 = 1. / (1.+M1*PIE1)
         FP1 = M1 * PIE1 / (1.+M1*PIE1)
         FD2 = 1. / (1.+M2*PIE2)
         FP2 = M2 * PIE2 / (1.+M2*PIE2)
         FP1SO4 = FP1 !: 
         FP2SO4 = FP2 !: 
      !
      !
      !DMD book page 244 to 246
      ! notes (4-133), (4-144)
      !
      !
      !equations for HS and SO4
      !
      !i.e.  AD(1,1)*SO41+ AD(1,2)*SO42+AD(1,3)*HST1 +BX(1) = 0  !
      !                                                          !
      !      AD(2,1)*SO41+ AD(2,2)*SO42              +BX(2) = 0  !
      !
      !       AD(3,3)*HST1+ AD(3,4)*HST2 +BX(3) = 0               
      !      AD(4,3)*HST1+ AD(4,4)*HST2 +BX(4) = 0               !
      !
      !
      !
      !
      !
      !equations for SO4
      !
      !
         AD (1, 1) = - S - KL12SO4 ! eqn (4-8-7) assuming H1~=0, W2~=0
         AD (1, 2) = KL12SO4 ! eqn (4-8-8) assuming H1~=0
         AD (2, 1) = KL12SO4 ! eqn (4-8-10) H1~= constant (i.e. H1 ~=0) and W2=0
      !
         If ( .Not. TF_SSTATE) Then
            AD (2, 2) = - (DLTS*KL12SO4+H2) / DLTS ! eqn (4-8-11) assuming H1~=0, W2=0 and H2,SO4 ~= H2
         !   AD(2,2) = -(DLTS*KL12SO4+HSO4)/DLTS        ! 
         Else
            AD (2, 2) = - KL12SO4 !  assuming H1~=0, W2=0 and H2,SO4 ~= H2
         End If
      !
      !
         BX (1) = S * SO40 ! (4-8-15) assuming H1~=0
      !
         If ( .Not. TF_SSTATE) Then
            BX (2) = H2 * SO4T2TM1 / DLTS ! (4-8-12) assuming JSO4,2 = 0 and H2,SO4 ~= H2
         !BX(2) = HSO4*SO4T2TM1/DLTS                    ! ( assuming JSO4,2 = 0 and H2,SO4 ~= H2
         Else
            BX (2) = 0.d0 ! (4-8-12) assuming JSO4,2 = 0
         End If
      !
         KHS_1 = FP1 * XAPP1HSP ** 2 / S * (O20/(KMHSO2)) + FD1 * &
        & XAPP1HSD ** 2 / S * (O20/(KMHSO2))
      !
         AD (1, 3) = KHS_1 !Basically diagenesis of SO4 (oxidation of HS) in layer 1, i.e. JSO4,1 at t+dt time
      !
      !
      !
      !equations for HS
      !
      !
      !
      !(after assuming H1~=0, JH2ST,2=0, JH2ST,1=0 (***: ok to assume JH2ST,1=0, but
      !should not assume JH2ST,2=0???), i.e. BX(4) should  have -J2 in it. However, here it is added to RA1
      !term in the solution for SO42)
      !
      !
         AD (3, 3) = - (W2) - (FP1) * (W12) - (FD1) * (S) - (FD1) * &
        & (KL12SO4) - (KHS_1)
      !             -(FD1)*(KL12)-(KHS_1)                            !, *****should be KL12 *****
      ! notes eqn (4-8-19)
      !
         AD (3, 4) = (FP2) * (W12) + (FD2) * (KL12SO4)
      !   AD(3,4) = (FP2)*(W12)+(FD2)*(KL12)                        !, *****should be KL12 ****
      ! notes eqn (4-8-20)
      !
         AD (4, 3) = (W2) + (FP1) * (W12) + (FD1) * (KL12SO4)
      !   AD(4,3) = (W2)+(FP1)*(W12)+(FD1)*(KL12)                    !, **** should be KL12*****
      ! eqn (4-8-22)
      !
         If ( .Not. TF_SSTATE) Then
            AD (4, 4) = - ((DLTS)*(FP2)*(W12)+(DLTS)*(FD2)*(KL12SO4)+&
           & (DLTS)*(W2)+(H2)) / (DLTS)!              +(DLTS)*(FD2)*(KL12)                          &!!, ****should be KL12 ****
         ! AD(4,4) = -((DLTS)*(FP2)*(W12)                        &!! eqn (4-8-23)
         !             +(DLTS)*(FD2)*(KL12SO4)                      &!
         !            +(DLTS)*(FD2)*(KL12)                          &!!, ****should be KL12 ****
         !             +(DLTS)*(W2)                               &!
         !             +(H2)                                      &!
         !             )/(DLTS)
         Else
            AD (4, 4) = - ((FP2)*(W12)+(FD2)*(KL12SO4)+(W2))!               +(FD2)*(KL12)                      &!, ****should be KL12 ****
         ! AD(4,4) = -((FP2)*(W12)                       &! eqn (4-8-23)
         !             +(FD2)*(KL12SO4)                  &!
         !             +(FD2)*(KL12)                      &!, ****should be KL12 ****
         !             +(W2)                               &!
         !             )
         End If
      !
         BX (3) = (HS0) * (S)! eqn (4-8-21)
      !
         If ( .Not. TF_SSTATE) Then
            BX (4) = (H2) * (HST2TM1) / (DLTS)!&!            ! eqn (4-8-24)
         !           +J2                                
         !
         Else
            BX (4) = 0.d0 !&!! eqn (4-8-24)
         !           +J2        
         End If
      !
         G (1) = (BX(1)*AD(3, 3)*AD(4, 4)-BX(3)*AD(1, 3)*AD(4, &
        & 4)-BX(1)*AD(3, 4)*AD(4, 3)+BX(4)*AD(1, 3)*AD(3, &
        & 4)+BX(2)*AD(1, 3)*AD(3, 4)) / (AD(1, 3)*AD(3, 4))

         G (2) = (BX(1)*AD(3, 3)*AD(4, 4)-BX(3)*AD(1, 3)*AD(4, &
        & 4)-BX(1)*AD(3, 4)*AD(4, 3)+BX(4)*AD(1, 3)*AD(3, 4)) / (AD(1, &
        & 3)*AD(3, 4))
      !
         HH (1, 1) = (AD(1, 1)*AD(3, 3)*AD(4, 4)-AD(1, 1)*AD(3, &
        & 4)*AD(4, 3)+AD(1, 3)*AD(2, 1)*AD(3, 4)) / (AD(1, 3)*AD(3, 4))
      !
         HH (1, 2) = (AD(1, 2)*AD(3, 3)*AD(4, 4)-AD(1, 2)*AD(3, &
        & 4)*AD(4, 3)+AD(1, 3)*AD(2, 2)*AD(3, 4)) / (AD(1, 3)*AD(3, 4))
      !
         HH (2, 1) = (AD(1, 1)*AD(3, 3)*AD(4, 4)-AD(1, 1)*AD(3, &
        & 4)*AD(4, 3)) / (AD(1, 3)*AD(3, 4))
      !
         HH (2, 2) = (AD(1, 2)*AD(3, 3)*AD(4, 4)-AD(1, 2)*AD(3, &
        & 4)*AD(4, 3)) / (AD(1, 3)*AD(3, 4))
      !

      !Find the solution for SO42 based on quadratic equation
      !
         RA2 = HH (1, 1) * HH (2, 2) - HH (1, 2) * HH (2, 1)
      !
         RA0 = (HH(1, 1)*G(2)-HH(2, 1)*G(1)) * (KMSO4)!Here KMSO4 is saturation concentration of SO4 in oxygen equivalents
      !(mgO2/L) for methane generation, DMD book does not have this term
      !
      !
         RA1 = (HH(1, 1)*G(2)-HH(2, 1)*G(1)) + (HH(1, 1)*HH(2, 2)-HH(1, &
        & 2)*HH(2, 1)) * (KMSO4) + HH (1, 1) * J2 !, why add J2 here ???? This shouuld be in the JHST,2 in BX(4) above
      !
      !Sign infront of sqrt
         SN1 = 1. !solution of a2*x^2+a1*x+a0=0
      !
         If (RA1 <= 0.0) SN1 = - 1. !see Num Rec p178
      !
         DISC = - (RA1+SN1*DSQRT(RA1**2-4.0D0*RA2*RA0)) / 2.0D0
      !
         If (Abs(DISC) /= 0.0D0) Then
         !
            R1 = DISC / RA2 !first root
         !==> DISC= R1*RA2
            R2 = RA0 / DISC !==> basically x2 = c/(a*x1)
         !
         !
         Else ! vjp 11/17/2005 added logic for case (a2*a0) = 0
            If (Abs(RA2) == 0.0D0) Then !  a2 = 0, a1/=0
            !
            !RA1*x + RA0 = 0 ==> x = - RA0/RA1
            !
               R1 = - RA0 / RA1 !first root
               R2 = R1 !2nd root
            !
            Else !  a0 = 0, a2 /= 0
            !
            !(RA2*x+RA1)x = 0
            !==> x = 0 or x= -RA1/RA2
               R1 = - RA1 / RA2 !first root
               R2 = 0.0D0 !2nd root
            End If
         End If
      !
         DBLSO42 = R1
      !
         If (DBLSO42 < 0.) DBLSO42 = R2
      !
      !
      !solve for SO41 from
      !
      !  SO41*HH(1,1)+HH(1,2)*SO42 + G(1) = 0
      !
      !
         DBLSO41 = - (HH(1, 2)*DBLSO42+G(1)) / HH (1, 1)
      !
      !   HST1 = -( AD(1,2)*DBLSO42                &        !James Martin code
      !             +AD(1,1)*DBLSO41+BX(1)        &
      !             )/AD(1,3)
      !
         HST1 = - (AD(1, 2)*DBLSO42+AD(1, 1)*DBLSO41+BX(1)) / AD (1, 3)
      !
      !i.e.  AD(1,1)*SO41+ AD(1,2)*SO42+AD(1,3)*HST1 +BX(1) = 0  !
      !                                                          !
      !      AD(2,1)*SO41+ AD(2,2)*SO42+             +BX(2) = 0  !
      !
      !
         HST2 = (AD(1, 2)*AD(3, 3)*DBLSO42+AD(1, 1)*AD(3, &
        & 3)*DBLSO41+BX(1)*AD(3, 3)-AD(1, 3)*BX(3)) / (AD(1, 3)*AD(3, &
        & 4))
     
#else /* ! included the following code for calculating HSO4, SO42,SO41, HS1, HS2, Chapter 11 of DMD */
      
      !**** Fractions and overall decay reaction velocity
         FD1 = 1. / (1.+M1*PIE1)
         FP1 = M1 * PIE1 / (1.+M1*PIE1)
         FD2 = 1. / (1.+M2*PIE2)
         FP2 = M2 * PIE2 / (1.+M2*PIE2)
         FP1SO4 = FP1
         FP2SO4 = FP2
      !
         RA2 = 2 * 0.5 * S * J2 !here J2 is in gO2/m^2/d and has carbon diagenesis discounted by denitrification
         RA1 = DDSO4 * 0.5 * J2
      !
         RA0 = - S * H2 * DDSO4 * SO40

      !Sign in front of sqrt
         SN1 = 1. !solution of a2*x^2+a1*x+a0=0
      !
         If (RA1 <= 0.0) SN1 = - 1. !see Num Rec p178
      !
         DISC = - (RA1+SN1*DSQRT(RA1**2-4.0D0*RA2*RA0)) / 2.0D0
      !
         If (Abs(DISC) /= 0.0D0) Then
         !
            R1 = DISC / RA2 !first root
         !==> DISC= R1*RA2
            R2 = RA0 / DISC !==> basically x2 = c/(a*x1)
         !
         !
         Else !logic for case (a2*a0) = 0
            If (Abs(RA2) == 0.0D0) Then !  a2 = 0, a1/=0
            !
            !RA1*x + RA0 = 0 ==> x = - RA0/RA1
            !
               R1 = - RA0 / RA1 !first root
               R2 = R1 !2nd root
            !
            Else !  a0 = 0, a2 /= 0
            !
            !(RA2*x+RA1)x = 0
            !==> x = 0 or x= -RA1/RA2
               R1 = - RA1 / RA2 !first root
               R2 = 0.0D0 !2nd root
            End If
         End If
      !
      !take the positive root as HSO4
      !
         HSO4 = R1
         If (HSO4 < 0.0) HSO4 = R2
      !
         If (HSO4 < 0.0) Then
            Write (*,*) 'oops, HSO4 is less than zero!!'
            HSO4 = 0.0000001 !make a very small number to avoid singularity
         End If
      !
         If (HSO4 > H2) HSO4 = H2
         KL12SO4 = DDSO4 / HSO4 !KL12=DD/H2 ==> DD= KL12*H2
      !
         AD (1, 1) = - (S) - (KL12SO4)!(4-8-7) assuming H1~=0, W2~=0
         AD (1, 2) = (KL12SO4)!(4-8-8) assuming H1~=0
         BX (1) = (S) * (SO40)!(4-8-15) assuming H1~=0
      !
         AD (2, 1) = (KL12SO4)!(4-8-10) H1~= constant (i.e. H1 ~=0) and W2=0
      !
         AD (2, 2) = - (KL12SO4)!-(HSO4))/(DLTS)                    !(4-8-11) assuming H1~=0, W2=0
      !
      !
      !  BX(2)     = (H2)*(SO4T2TM1)/(DLTS)                        !(4-8-12) assuming JSO4,2 = 0
         BX (2) = - J2 * 0.5 !+(HSO4)*(SO4T2TM1)/(DLTS)            !(4-8-12) assuming steady state
      !
      !here J2 is JH2ST in gO2/m^2/d and 0.5 is alpha_SO4C in gO2/gO2
      !
      !
      !solve linear equation for SO42 and SO41 
      !
         DBLSO41 = (-AD(2, 2)*BX(1)+AD(1, 2)*BX(2)) / (AD(1, 1)*AD(2, &
        & 2)-AD(1, 2)*AD(2, 1))
      !
      !
      !
         DBLSO42 = (-AD(1, 1)*BX(2)+AD(2, 1)*BX(1)) / (AD(1, 1)*AD(2, &
        & 2)-AD(1, 2)*AD(2, 1))
      !
      !
         If (DBLSO41 < 0.0) DBLSO41 = 0.0
         If (DBLSO42 < 0.0) DBLSO42 = 0.0
      !
      !check if SO42= SO41/2.0 is ture
      !
         If (Abs(DBLSO41-DBLSO42*2.0) > 0.00001) Then
            Write (*,*) 'Oops, SO41 /= SO42*2 '
            Write (*,*) 'SO41=', DBLSO41, 'SO42*2=', DBLSO42 * 2
         End If
      !
         KHS_1 = FP1 * XAPP1HSP ** 2 / S * (O20/(KMHSO2)) + FD1 * &
        & XAPP1HSD ** 2 / S * (O20/(KMHSO2))
      !
      !Assuming H1~=0, JH2ST,1=0.0
      !
         AD (3, 3) = - (W2) - (FP1) * (W12) - (FD1) * (S) - (FD1) * &
        & (KL12) - (KHS_1)!(4-8-19)

         AD (3, 4) = (FP2) * (W12) + (FD2) * (KL12)!(4-8-20)
      !
      !  AD(4,3) = (W2)+(FP1)*(W12)+(FD1)*(KL12SO4)
         AD (4, 3) = (W2) + (FP1) * (W12) + (FD1) * (KL12)!(4-8-22)
      !
         If ( .Not. TF_SSTATE) Then
            AD (4, 4) = - &
           & ((DLTS)*(FP2)*(W12)+(DLTS)*(FD2)*(KL12)+(DLTS)*(W2)+(H2)) &
           & / (DLTS)   !                   +(DLTS)*(FD2)*(KL12SO4)                 &!
         Else
            AD (4, 4) = - ((FP2)*(W12)+(FD2)*(KL12)+(W2)) !                   +(FD2)*(KL12SO4)                     &!

         End If
      !
         BX (3) = (HS0) * (S)!(4-8-21)
      !
         If ( .Not. TF_SSTATE) Then
            BX (4) = (H2) * (HST2TM1) / (DLTS) + J2 * 0.5 !(4-8-24)
         !here J2 is JH2ST,2 in gO2/m^2/d, 0.5 is to convert to SO4 unit also in gO2 (4-57)
         Else
            BX (4) = + J2 * 0.5 !(4-8-24)
         End If
      !
      !
      !solve for HST1, HST2 
      !
      !
         HST1 = (-AD(4, 4)*BX(3)+AD(3, 4)*BX(4)) / (AD(3, 3)*AD(4, &
        & 4)-AD(3, 4)*AD(4, 3))
      !
         HST2 = (-AD(3, 3)*BX(4)+AD(4, 3)*BX(3)) / (AD(3, 3)*AD(4, &
        & 4)-AD(3, 4)*AD(4, 3))
      !
      !check to make sure both HST1 HST2 >0
         If (HST1 < 0.0) HST1 = 0
         If (HST2 < 0.0) HST2 = 0
      !
#endif    /*end solving SO4 and HSO4 and HS */
      !
         HS1 = FD1 * HST1
         HS2 = FD2 * HST2
         HS2AV = FD2 * HST2
      !
         IDEBUG_SED = IDEBUG_SED + 1
      !
         SO41 = DBLSO41
         SO42 = DBLSO42
         SO42AV = SO42
         SO4T2 = SO42
      !
         JHS = S * (HS1-HS0)!gO2/m^2/day (flux of dissolved HS to overlying water column)
      ! (4-61)
      !
      !CSOD due to HS oxidation in layer 1 limited by avaiable O2 in overlying water
      ! (m/d)*(gO2/m^3)*(1) ==> unit is gO2/m^2/d
      !
         CSODHS = (XAPP1HSP**2/S*FD1*HST1+XAPP1HSP**2/S*FP1*HST1) * &
        & (O20/KMHSO2)!limitation by avaiable O20
      ! (4-62)
      !
         IDEBUG_SED = IDEBUG_SED + 1
      !solve for methane equation
      !
      !
         XJ2 = J2 * KMSO4 / (SO42+KMSO4)! (4-70) (4-157-2)
      !Where J2 = XJC1 = JCX-XJCNO3
      !
      !gO2/m^2/day !generation of CH4 in layer 2 regulated by abundance level of SO4 in layer 2
         XJ2CH4 = XJ2 ! (4-70) flux of methane generation (gO2/m^2/d)
      !
         X1J2 = J2 * DBLSO42 / (SO42+KMSO4)!never used = J2-XJ2
      !e.g. the amount of J2 oxidized by SO42
      !
         IDEBUG_SED = IDEBUG_SED + 1
      !**** Methane solved assuming steady state
      !
         KMC1 = 0.0
         PIE1 = 0.0
         PIE2 = 0.0
      !
         K0H1D = 0. ! (3-32-7)
         K0H1P = 0. ! (3-32-9)
      !
         K1H1P = 0. ! (3-32-10, which should be non-zero but equivalent to zero with fp1 =0 in (3-32-11)
         K1H1D = XAPPCH4 ** 2 / S * (O20/(KMCH4O2+O20)) + S ! (4-162-1)(3-32-11 with fp1=0) (3-32-5)
      !
         K2H2D = 0. ! (3-36-8)
         K2H2P = 0. ! (3-36-8-1)
      !
         J1 = S * CH40 ! (4-164)
      !
         K3 = 0.0 ! (3-36-6)
      !
         J2 = XJ2 !gO2/m^2/day  !carbon diagenesis flux converted from mgC/m^3/day to gO2/m^3/day
      !call steady state calculations
         Call SEDSSFNL (CH41, CH42, CH42AV, CH4T1, CH4T2, CH4T2AV, &
        & CH41TM1, CH4T2TM1, 1)
      !
         If (CH42 > CH4SAT) Then
            CH42 = CH4SAT
            CH41 = (CH40*S**2+CH42*KL12*S) / &
           & (S**2+KL12*S+XAPPCH4**2*(O20/(KMCH4O2+O20)))! (4-169) with W2~=0
         !
         End If
      !
      !**** Calculate changes in CH4 and HS stored in the sediment
      !
         If ( .Not. TF_SSTATE) Then
            DCH4T2 = (CH4T2-CH4T2TM1) * H2 / DLTS
            DHST2 = (HST2-HST2TM1) * H2 / DLTS !
         !
         Else
            DCH4T2 = 0.d0
            DHST2 = 0.d0
         End If
      !
      !**** Calculate CSOD (gO2/M2/day)
      !
         CSODCH4 = XAPPCH4 ** 2 / S * (O20/(KMCH4O2+O20)) * CH41 ! (4-170)
      !
         CSOD = CSODCH4 + CSODHS
      !
      !**** Calculate Fluxes
         JCH4 = S * (CH41-CH40)
         JCH4AQ = S * (CH41-CH40)!same as JCH4
      !
      !   FLUXHS    = S*FD1*HS1               ! (this should be HST1)
         FLUXHS = S * FD1 * HST1
         FLUXHSCH4 = JCH4AQ + FLUXHS !Never used
      !
      !**** If not flux or SOD or stored then it must escape as gas flux
         JCH4G = 0.
         If (CH42 == CH4SAT) Then
            JCH4G = XJ2CH4 - DCH4T2 - CSODCH4 - JCH4AQ
            JCH4G = XJC1 - DCH4T2 - DHST2 - CSOD - FLUXHSCH4 ! used Jame Martin's method
         !
         !
         End If
         Return
      End Subroutine GET_SO4_HS_CH4_CSOD_FLUXES
    !
      Function GET_ZHTANH4F (TEMPVAL)!m/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTANH4F
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTANH4F = KAPPNH4F * THTANH4 ** TEMP202
         Return
      End Function GET_ZHTANH4F
    !
      Function GET_ZHTANH4S (TEMPVAL)!m/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTANH4S
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTANH4S = KAPPNH4S * THTANH4 ** TEMP202
         Return
      End Function GET_ZHTANH4S
    !
      Function GET_ZHTA1HSD (TEMPVAL)!m/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTA1HSD
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTA1HSD = KAPP1HSD * THTAH2S ** TEMP202
         Return
      End Function GET_ZHTA1HSD
    !
      Function GET_ZHTA1HSP (TEMPVAL)!m/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTA1HSP
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTA1HSP = KAPP1HSP * THTAH2S ** TEMP202
         Return
      End Function GET_ZHTA1HSP
    !
      Function GET_ZHTANO3F (TEMPVAL)!m/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTANO3F
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTANO3F = KAPPNO3F * THTANO3 ** TEMP202
         Return
      End Function GET_ZHTANO3F
    !
      Function GET_ZHTANO3S (TEMPVAL)!m/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTANO3S
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTANO3S = KAPPNO3S * THTANO3 ** TEMP202
         Return
      End Function GET_ZHTANO3S
    !
      Function GET_ZHTAK2NO3 (TEMPVAL)!DMD p. 105 eqn (4.9) (m/d) (anerobic layer denitrification velocity)
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTAK2NO3
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTAK2NO3 = K2NO3 * THTANO3 ** TEMP20 !m/d
         Return
      End Function GET_ZHTAK2NO3
    !
      Function GET_ZL12NOM (TEMPVAL)!temperature contol on diffusion rate of layer 1, 2 (unitless)
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZL12NOM
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZL12NOM = THTADD ** TEMP20 !unitless
         Return
      End Function GET_ZL12NOM
    !
      Function GET_ZW12NOM (TEMPVAL)!temperature control on particle mixing rate of layer 1,2 (unitless)
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZW12NOM
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZW12NOM = THTADP ** TEMP20 !unitless
         Return
      End Function GET_ZW12NOM
    !
      Function GET_ZHTAPON1 (TEMPVAL)!1/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTAPON1
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTAPON1 = KPON1 * THTAPON1 ** TEMP20
         Return
      End Function GET_ZHTAPON1
    !
      Function GET_ZHTAPON2 (TEMPVAL)!1/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTAPON2
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTAPON2 = KPON2 * THTAPON2 ** TEMP20
         Return
      End Function GET_ZHTAPON2
    !
      Function GET_ZHTAPON3 (TEMPVAL)!1/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTAPON3
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTAPON3 = KPON3 * THTAPON3 ** TEMP20
         Return
      End Function GET_ZHTAPON3
    !
      Function GET_ZHTAPOC1 (TEMPVAL)!1/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTAPOC1
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTAPOC1 = KPOC1 * THTAPOC1 ** TEMP20
         Return
      End Function GET_ZHTAPOC1
    !
      Function GET_ZHTAPOC2 (TEMPVAL)!1/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTAPOC2
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTAPOC2 = KPOC2 * THTAPOC2 ** TEMP20
         Return
      End Function GET_ZHTAPOC2
    !
      Function GET_ZHTAPOC3 (TEMPVAL)!1/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTAPOC3
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTAPOC3 = KPOC3 * THTAPOC3 ** TEMP20
         Return
      End Function GET_ZHTAPOC3
    !
      Function GET_ZHTAPOP1 (TEMPVAL)!1/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTAPOP1
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTAPOP1 = KPOP1 * THTAPOP1 ** TEMP20
         Return
      End Function GET_ZHTAPOP1
    !
      Function GET_ZHTAPOP2 (TEMPVAL)!1/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTAPOP2
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTAPOP2 = KPOP2 * THTAPOP2 ** TEMP20
         Return
      End Function GET_ZHTAPOP2
    !
      Function GET_ZHTAPOP3 (TEMPVAL)!1/d
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTAPOP3
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTAPOP3 = KPOP3 * THTAPOP3 ** TEMP20
         Return
      End Function GET_ZHTAPOP3
    !
      Function GET_ZHTASI (TEMPVAL)!Silicate dissolution rate (1/day), DMD book page P 151
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTASI
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTASI = KSI * THTASI ** TEMP20
         Return
      End Function GET_ZHTASI
    !
      Function GET_ZHTACH4 (TEMPVAL)!m/d ! (4-70), DMD (40.39)
         Real (SP) :: TEMPVAL !temperature (degC)
         Real (SP) :: GET_ZHTACH4
         TEMP20 = TEMPVAL - 20.d0
         TEMP202 = TEMP20 / 2.d0
         GET_ZHTACH4 = KAPPCH4 * THTACH4 ** TEMP202 ! (4-70)
         Return
      End Function GET_ZHTACH4
    !
    !
      Subroutine SEDTEST_UPDATE_WATERCOLUMN (NXJDAY)
      !
      !read the overlying water column data information to make sure overlying water column concentration of current and next time
      !steps are given correctly
      !
         Implicit None
         Real (SP) :: NXJDAY
      !
      !Input file looks like the following (e.g. for test2c copied from SedFlux_ver17b12_test2c.xlsm
      !
      !======================
      !Sediment overyling water column condition forcing file
      !        Segment ID    Date+time    Jcin        Jnin        Jpin        Jsin        O20            Depth        Tw        NH30
      !                                  (gO2/m^2/d)    (gN/m^2/d)    (gP/m^2/d)    (gSi/m^2/d)    (mgO2/L)    (m)        (deg C)
      !        1    1/1/00 12:00 AM        0.093450    0.051800    0.000750    0.097125    12.0        4.30    0.0            0.100
      !        1    1/1/00 6:00 AM        0.092830    0.051066    0.000743    0.095748    12.1        4.32    0.3            0.100
      !        1    1/1/00 12:00 PM        0.092670    0.050807    0.000740    0.095264    12.1        4.34    0.3            0.100
      !        1    1/1/00 6:00 PM        0.092371    0.050443    0.000736    0.094580    12.2        4.35    0.4            0.100
      !        1    1/2/00 12:00 AM        0.092201    0.050215    0.000732    0.094153    12.2        4.35    0.4            0.100
      !        1    1/2/00 6:00 AM        0.092018    0.049986    0.000729    0.093723    12.2        4.35    0.4            0.100
      !        1    1/2/00 12:00 PM        0.091872    0.049799    0.000726    0.093373    12.2        4.36    0.4            0.100
      !        1    1/2/00 6:00 PM        0.091734    0.049624    0.000723    0.093045    12.2        4.36    0.4            0.100
      !        1    1/3/00 12:00 AM        0.091610    0.049466    0.000721    0.092749    12.3        4.36    0.4            0.100
      !        1    1/3/00 6:00 AM        0.091494    0.049317    0.000718    0.092469    12.3        4.36    0.4            0.100
      !        1    1/3/00 12:00 PM        0.091386    0.049176    0.000716    0.092205    12.3        4.36    0.4            0.100
      !....
      !        1    12/30/03 6:00 AM    0.139803    0.059920    0.001033    0.112350    12.2        4.43    0.6            0.073
      !        1    12/30/03 12:00 PM    0.139610    0.059949    0.001032    0.112405    12.2        4.43    0.6            0.073
      !        1    12/30/03 6:00 PM    0.139419    0.059978    0.001031    0.112460    12.2        4.43    0.6            0.072
      !======================
      !
133      Continue
         Do WHILE (JDAY >= NXWCL)
         !
         !cycle back the record
            JCIN_R1 = JCIN_R2
            JNIN_R1 = JNIN_R2
            JPIN_R1 = JPIN_R2
            JSIN_R1 = JSIN_R2
            O20_R1 = O20_R2
            D_R1 = D_R2
            TW_R1 = TW_R2
            NH30_R1 = NH30_R2
            NO30_R1 = NO30_R2
            PO40_R1 = PO40_R2
            SIAT0_R1 = SIAT0_R2
            CH40_R1 = CH40_R2
            SALT0_R1 = SALT0_R2
         !
         !read in a new record, if end of file, jump to 134 and open next file
         !
            Read (WCL, 1010, END=134) NXWCL, JCIN_R2, JNIN_R2, JPIN_R2, &
           & JSIN_R2, O20_R2, D_R2, TW_R2, NH30_R2, NO30_R2, PO40_R2, &
           & SIAT0_R2, CH40_R2, SALT0_R2 !salinity
         !
            If (MSR) Then
               Write (*, 1010) NXWCL, JCIN_R2, JNIN_R2, JPIN_R2, &
              & JSIN_R2, O20_R2, D_R2, TW_R2, NH30_R2, NO30_R2, &
              & PO40_R2, SIAT0_R2, CH40_R2, SALT0_R2 !salinity
            !
            End If
         !
            NXWCL = (WCLPTR-1) * WCLFDAYS + NXWCL !global record time in all forcing of WCL
         !
         End Do
      !
         Go To 135
      !
      !
      !Open next data file and read one record if necesary
      !
      !
134      Continue
      !
         WCLPTR = WCLPTR + 1
         Write (*,*) 'Opening water column condition ', 'file ', WCLFN &
        & (WCLPTR), ' at ', 'day ', JDAY
         Close (WCL)
      !
         Open (WCL, File=WCLFN(WCLPTR), Status='OLD')
      !
         Read (WCL, 1000)!Read past the header lines
         Read (WCL, 1010, END=134) NXWCL, JCIN_R2, JNIN_R2, JPIN_R2, &
        & JSIN_R2, O20_R2, D_R2, TW_R2, NH30_R2, NO30_R2, PO40_R2, &
        & SIAT0_R2, CH40_R2, SALT0_R2 !salinity
         If (MSR) Then
            Write (*, 1010) NXWCL, JCIN_R2, JNIN_R2, JPIN_R2, JSIN_R2, &
           & O20_R2, D_R2, TW_R2, NH30_R2, NO30_R2, PO40_R2, SIAT0_R2, &
           & CH40_R2, SALT0_R2 !salinity
         !
         !WRITE(*,*)'here I ''m doing nothing 2'
         !READ(*,*)
         End If
      !
         NXWCL = (WCLPTR-1) * WCLFDAYS + NXWCL !global record time in all forcing of WCL
      !
         If (JDAY >= NXWCL) Go To 133!keep reading if JDAY still larger than NXWCL
      !
      !
135      Continue
      !
         NXJDAY = Max (NXJDAY, NXWCL)!global record time in all forcing of WCL
         Write (*,*) 'NXJDAY=', NXJDAY
1000     Format (/ /)
!!!---1010 FORMAT(13F8.0,:/(:8X,14F8.0))
1010     Format (14(F12.5, 1 X))
      !
      End Subroutine SEDTEST_UPDATE_WATERCOLUMN
    !
End Module MOD_SED
  !
  !
  !
