!mod_sf.F
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
!
!************************************************************************
!**              MODULE MOD_SF	: Suspension feeder module             **
!************************************************************************
!**                                                                    **
!** Initial development based on Chesapeake Bay Suspension Feeder Model**
!**                                                                    **
!**      -- HydroQual, Inc.                                            **
!**         Version: 18 Apr 97                                         **
!**                                                                    **
!**         hypoxic mortality - SAL 4.16.97                            **
!**                                                                    **
!**                                                                    **
!**      Latest change: form a separate module, input, output          **
!**   and calculation are independent of main program                  **
!************************************************************************
!
Module MOD_SF
  !
      Use MOD_SIZES, Only: NSSFP, MGL
      Use MOD_LIMS, Only: MTLOC, MLOC, KBM1
      Use MOD_PREC, Only: SP
      Use MOD_SED_SF_EXCHANGE_VARS, Only: JLPOC_SED_SF, JLPOC_SED_SF, &
     & JLPON_SED_SF, JLPOP_SED_SF, JRPOC_SED_SF, JRPON_SED_SF, &
     & JRPOP_SED_SF, JNH4_SED_SF, JPO4_SED_SF, SOD_SED_SF, JSA_SED_SF, &
     & JSU_SED_SF, SSI_SED_SF, SU_SED_SF, SA_SED_SF, PIP_SED_SF
  !
      Implicit None
      Save
  !
  !Moved the following from mod_wqm.F
      Logical :: SFEEDER, HYPOXFX_SF
      Real (SP) :: FRSASF
      Integer :: NSPECIES
  !
      Integer, Allocatable :: SEDTYPE (:, :)
  !
  !
      Real (SP), Dimension (NSSFP) :: FILT, SFA1, SFA2, SFA3, SFA4, &
     & SFA5, MAXING, SFCN, SFCP, SFRESP, SFPRED, SFTMN, THTAFILT, &
     & THTARESP, THTAPRED, XKPO2, SFTD, FILTFACT, RESPFACT, SFDOH, &
     & SFDOQ, SFATURB, SFBTURB !Maximum potential ASSIMILATED carbon
  !	FILT,		&   !filtration rate of suspension feeders
  !!Maximum potential ASSIMILATED carbon
  !SFA1,		&	! assimilation rate cyanobacteria
  !SFA2,		&	! assimilation rate diatoms
  !SFA3,    	&   ! assimilation rate greens
  !SFA4,		&   ! assimilation rate lpoc
  !SFA5,       &   ! assimilation rate rpoc
  !MAXING,		&
  !SFCN,		&	!C:N ratio of suspension feeder
  !SFCP,       &   !C:P ratio of suspension feeder
  !SFRESP,		&   !
  !SFPRED,		&
  !SFTMN, 		&
  !THTAFILT,   &
  !THTARESP, 	&
  !THTAPRED, 	&
  !XKPO2,      &
  !SFTD,  		&
  !FILTFACT, 	&
  !RESPFACT,   &
  !SFDOH, 		&
  !SFDOQ, 		&
  !SFATURB, 	&
  !SFBTURB
  !
  !MBM 971021 added turbidity sensitivity parameters
      Real (SP), Allocatable :: SF_JLPOC (:), SF_JLPON (:), SF_JLPOP &
     & (:), SF_JRPOC (:), SF_JRPON (:), SF_JRPOP (:), JNH4SF (:), &
     & JPO4SF (:), SODSF (:), JSASF (:), JSUSF (:), SF_SSI (:), SF_SU &
     & (:), SF_SA (:), SF_PIP (:)
  !
      Real (SP), Allocatable :: SFEED (:, :)
  !
      Real (SP), Allocatable :: SFFILTCT (:, :), SFRESPCT (:, :), &
     & SFPREDCT (:, :), SFRATEI (:, :), SFGCMAX (:, :), SFGMIN (:, :), &
     & SFGCING (:, :), SFCFECES (:, :), SFNFECES (:, :), SFPFECES (:, &
     & :), SFCPSFEC (:, :), SFNPSFEC (:, :), SFPPSFEC (:, :)
  !
      Real (SP), Allocatable :: SF_SFGC (:, :), SF_RESP (:, :), SF_PRED &
     & (:, :), SF_RMORT (:, :), ASF_SFGC (:, :), ASF_RESP (:, :), &
     & ASF_PRED (:, :), ASF_RMORT (:, :)
  !
      Real (SP), Allocatable :: SFEED_GL (:, :)
      Real (SP), Allocatable :: SEDTYPE_GL (:, :)
  !
  !
Contains
  !
  !subroutines:
  !	subroutine 	SUSPFEED()
  !	subroutine 	SF_ALLOC()
  !	subroutine 	SF_DEALLOC()
  !	subroutine	SF_READ()
  !	subroutine	SF_INIT()
  !
      Subroutine SUSPFEED (N, DTB1, DTB2, DTB3, DTPIB1, DTPIB2, DTPIB3, &
     & DTLPOC, DTLPOP, DTLPON, DTRPOC, DTRPOP, DTRPON, DTSSI, DTSIUPB, &
     & DTSIAT, DTPO4)
!
         Use MOD_WQM, Only: ANC1, ANC2, ANC3, B1, B2, B3, CTEMP, DOXG, &
        & LPOC, LPON, LPOP, SIUPB, Q1, Q2, Q3, RPOC, RPON, RPOP, SSI, &
        & ASC1, ASC2, ASC3, DLT, SIAT, KADPO4, PO4, WSS, WSSNET !,         &!
!
         Use MOD_HYDROVARS, Only: DZ, DZ2D, D
         Implicit None
    !dummy variables
    !
         Integer :: N !suspension feeder species index
         Real (SP), Dimension (0:MTLOC, KBM1) :: DTB1, DTB2, DTB3, &
        & DTPIB1, DTPIB2, DTPIB3, DTLPOC, DTLPOP, DTLPON, DTRPOC, &
        & DTRPOP, DTRPON, DTSSI, DTSIUPB, DTSIAT, DTPO4 !rate of change of PO4 (available) gP/m^3/sec)
    !
    !
    ! THESE ARE LOCAL WORK ARRAYS (end in X)
         Real (SP), Dimension (MLOC) :: SF, JNH4SFX, JPO4SFX, SODSFX, &
        & JSASFX, JSUSFX, SEDTYPEX, SF_LPOCX, SF_LPONX, SF_LPOPX, &
        & SF_RPOCX, SF_RPONX, SF_RPOPX, SF_SUX, SF_SASX, SF_PIPX, &
        & SF_SSIX !  positive adding PO4 to water column
    !
         Real (SP) :: MAXINGX
    !
         Real (SP) :: LOGICT, CLOGI, LPOCASF, LPONASF, LPOPASF
    !
         Real (SP) :: NFECES, NFILT, NPSFEC
    !
         Integer :: I
         Real (SP) :: FILTX, SFA1X, SFA2X, SFA3X, SFA4X, SFA5X, SFCNX, &
        & SFCPX, SFRESPX, SFPREDX, SFTMNX, THTAFILTX, THTARESPX, &
        & THTAPREDX, XKPO2X, TDX, RESPFACTX, FILTFACTX, DOhx, DOqx, &
        & ATURB, BTURB, RD, RMORT, B1ASF, B2ASF, B3ASF, RPOCASF, &
        & RPONASF, RPOPASF
    !
         Real (SP) :: SSIASF, SUASF, FILTCT, RESPCT, PREDCT, RATEI, &
        & GCMAX, SFGC, GMIN, CFECES, RCFECES, RNFECES, PFECES, RPFECES, &
        & CFILT, PFILT, CPSFEC, RCPSFEC, RNPSFEC, PPSFEC, RPPSFEC, DOX, &
        & SFD, TSOLID, TURBRED, RELING, GN, GP, SFAC1, SFAC2, SFAC3, &
        & SFAC4, SFAC5, SFAN1, SFAN2, SFAN3, SFAN4, SFAN5, SFAP1, &
        & SFAP2, SFAP3, SFAP4, SFAP5, PF, DTDAY
    !
    !
         Do I = 1, MLOC
            SF (I) = SFEED (I, N)
            SEDTYPEX (I) = SEDTYPE (I, N)
       !
       !MBM 961226 diagnostic arrays reset everywhere, every time step
            SFFILTCT (I, N) = 0.
            SFRESPCT (I, N) = 0.
            SFPREDCT (I, N) = 0.
            SFRATEI (I, N) = 0.
            SFGCMAX (I, N) = 0.
            SFGMIN (I, N) = 0.
            SFCFECES (I, N) = 0.
            SFNFECES (I, N) = 0.
            SFPFECES (I, N) = 0.
            SFCPSFEC (I, N) = 0.
            SFNPSFEC (I, N) = 0.
            SFPPSFEC (I, N) = 0.
            SF_SFGC (I, N) = 0.
            SF_RESP (I, N) = 0.
            SF_PRED (I, N) = 0.
            SF_RMORT (I, N) = 0.
       !
         End Do
    !
         FILTX = FILT (N)
         SFA1X = SFA1 (N)
         SFA2X = SFA2 (N)
         SFA3X = SFA3 (N)
         SFA4X = SFA4 (N)
         SFA5X = SFA5 (N)
         MAXINGX = MAXING (N)
         SFCNX = SFCN (N)
         SFCPX = SFCP (N)
         SFRESPX = SFRESP (N)
         SFPREDX = SFPRED (N)
         SFTMNX = SFTMN (N)
         THTAFILTX = THTAFILT (N)
         THTARESPX = THTARESP (N)
         THTAPREDX = THTAPRED (N)
         XKPO2X = XKPO2 (N)
         TDX = SFTD (N)
         RESPFACTX = RESPFACT (N)
         FILTFACTX = FILTFACT (N)
         DOhx = SFDOH (N)
         DOqx = SFDOQ (N)
         ATURB = SFATURB (N)
         BTURB = SFBTURB (N)
    !
    ! do some initial calcs for low o2 respiration stress
    !
         RD = 4.605 / TDX ! ln(1/100)=4.6=99% mortality
    !
         RMORT = 0.0
    !
    ! reset local flux contribution to zero at start of each timestep as
    ! each species contribution is added to total at bottom
    !
         Do I = 1, MLOC
            JNH4SFX (I) = 0.0
            JPO4SFX (I) = 0.0
            SODSFX (I) = 0.0
            JSASFX (I) = 0.0
            JSUSFX (I) = 0.0
            SF_LPOCX (I) = 0.0
            SF_LPONX (I) = 0.0
            SF_LPOPX (I) = 0.0
            SF_RPOCX (I) = 0.0
            SF_RPONX (I) = 0.0
            SF_RPOPX (I) = 0.0
            SF_SSIX (I) = 0.0
            SF_SUX (I) = 0.0
            SF_SASX (I) = 0.0
            SF_PIPX (I) = 0.0
         End Do
    !
    ! Re-init global particulate fluxes
    !  labile (SFLUX{C,N,P}), refract, & unavail Si
    ! Only re-initialize prior to call for first species (N=1)...
         If (N == 1) Then
            Do I = 1, MLOC
               JNH4SF (I) = 0.0
               JPO4SF (I) = 0.0
               SODSF (I) = 0.0
               JSASF (I) = 0.0
               JSUSF (I) = 0.0
               SF_JLPOC (I) = 0.0
               SF_JLPON (I) = 0.0
               SF_JLPOP (I) = 0.0
               SF_JRPOC (I) = 0.0
               SF_JRPON (I) = 0.0
               SF_JRPOP (I) = 0.0
               SF_SSI (I) = 0.0
               SF_SU (I) = 0.0
               SF_SA (I) = 0.0
               SF_PIP (I) = 0.0
            End Do
         End If
    !
    ! main loop
    !
         Do I = 1, MLOC
       !     KWC=BBN(I)
       !
       ! COMPUTE FOOD,PARTICULATES AVAILABLE TO SUSPENSION FEEDERS
       !
            B1ASF = Max (1.E-6, B1(I, KBM1))
            B2ASF = Max (1.E-6, B2(I, KBM1))
            B3ASF = Max (1.E-6, B3(I, KBM1))
            LPOCASF = Max (1.E-6, LPOC(I, KBM1))
            RPOCASF = Max (1.E-6, RPOC(I, KBM1))
            LPONASF = Max (1.E-6, LPON(I, KBM1))
            RPONASF = Max (1.E-6, RPON(I, KBM1))
            LPOPASF = Max (1.E-6, LPOP(I, KBM1))
            RPOPASF = Max (1.E-6, RPOP(I, KBM1))
            SSIASF = Max (1.E-6, SSI(I, KBM1))
            SUASF = Max (1.E-6, SIUPB(I, KBM1))
       !
            FILTCT = 0.0
            RESPCT = 0.0
            PREDCT = 0.0
            RMORT = 0.0 !MBM 972830
            RATEI = 0.0
            GCMAX = 0.0
            SFGC = 0.0
            GMIN = 0.0
            CFECES = 0.0
            RCFECES = 0.0
            NFECES = 0.0
            RNFECES = 0.0
            PFECES = 0.0
            RPFECES = 0.0
            CFILT = 0.0
            NFILT = 0.0
            PFILT = 0.0
            CPSFEC = 0.0
            RCPSFEC = 0.0
            NPSFEC = 0.0
            RNPSFEC = 0.0
            PPSFEC = 0.0
            RPPSFEC = 0.0
       !
            If (CTEMP(I) < SFTMNX) Go To 255
       !
            If (DOXG(I, KBM1) < 0.0) Then
               DOX = 0.0
            Else
               DOX = DOXG (I, KBM1)
            End If
       !
       ! take into account habitat type
       !
            If (SEDTYPEX(I) == 0) Then
               SF (I) = 0.1 !  10 mgC/m2 "refuge" population
               Go To 255
            End If
       !
       ! calculate biomass dependent respiration and filtration
       !
            SFD = SF (I) / 1000. ! convert to g C /m^2 for this calc.
            RESPCT = SFRESPX / SFD ** RESPFACTX
            FILTCT = FILTX / SFD ** FILTFACTX
       !
            FILTCT = FILTCT * THTAFILTX ** (CTEMP(I)-20.0)
            RESPCT = RESPCT * THTARESPX ** (CTEMP(I)-20.0)
            PREDCT = SFPREDX * THTAPREDX ** (CTEMP(I)-20.0)
       !
       ! MBM 970928 Temperature dependence on Time-to-Death, DOXG-mortality rate
       !     RD = RD*(1.08**(CTEMP(I)-20.0))
       !
       ! O2 limitations on filtration, predation and respiration
       !
       ! MBM 961224 added logistical control on filtr, hypoxic resp fcts
       !   logict ranges 0.0-1.0 with DOX; use 1-logict for increased mort/resp.
       !   MAX() prevents underflows in the event of superhigh btm DOXG.
       !
       !  One HARDWIRED constant!
            CLOGI = 1.1
       !
       !MBM 970109  hypoxfx_SF switch added
            If (HYPOXFX_SF) Then
          !
               LOGICT = 1.0 / &
              & (1.0+Exp(Max(CLOGI*(DOhx-DOX)/(DOhx-DOqx),-25.)))
          !
          ! MBM 961224 added logistical control on filtr, hypoxic resp fcts
               FILTCT = FILTCT * LOGICT
          ! MBM 970612 added logistical control on resp to simulate
          ! MBM        greatly reduced (anaerobic) respiration rates as FILT -> 0
          ! MBM        Relieves SOD from SODSF component when SF are not using oxygen.
               RESPCT = RESPCT * LOGICT
               PREDCT = PREDCT * DOX / (DOX+XKPO2X)
          !
          ! MBM 961224 added logistical control on filtr, hypoxic resp fcts
          ! SAL 970416 now use RMORT instead of respiration
          ! MBM 970928 Re-establishing logistic function for hypoxic mortality
               RMORT = RD * (1.0-LOGICT)
          ! MBM 970429 New mortality formula from analysis of benthic anoxia LD50 data
          !     RMORT = 0.04773*exp( 0.05319*CTEMP(B) - 0.69538*DOX )
          !
            End If
       !
            TSOLID = 2.5 * (B1ASF+B2ASF+B3ASF+LPOCASF+RPOCASF) + SSIASF
       !
       ! MBM 971021
       ! Reduction in filtration from total particulate load from
       !   PHK oyster model; solids as ** gm/L **; red as % reduction of FILT
       !   Same function, values as in PHK, but has been algebraically modified
       !   to be more comprehensible (originally = [log(turb) + 3.38] / 0.0418).
       !
       !     TURBRED = 81. + 24. * LOG10( 0.001*TSOLID )
            TURBRED = Max (ATURB+BTURB*Log10(0.001*TSOLID), 0.0)
            FILTCT = Max (FILTCT*(1.0-0.01*TURBRED), 0.0)
       !
       !
       ! CARBON INGESTION RATE
       !   FILTCT - REALIZED FILTRATION RATE LITERS/D/MG C BIOMASS
       !
            RATEI = FILTCT * (B1ASF+B2ASF+B3ASF+LPOCASF+RPOCASF)
       !
       !     COMPUTE SUSPENSION FEEDER GROWTH RATE, COMPARING TO MAX
       !     INGESTION RATE - MAX AMOUNT OF FOOD GUT CAN TAKE IN
       !        if ratei < rateimax  then use normal growth rate
       !        if ratei > rateimax  then limit to max rate
       !
            RELING = Min (1.0, MAXINGX/(RATEI+1.0E-10))
       !
       ! Total Filtered from Organic Matter Pools    (MBM 971010 moved up from below)
            CFILT = SF (I) * RATEI
       !
            NFILT = SF (I) * FILTCT * &
           & (B1ASF*ANC1+B2ASF*ANC2+B3ASF*ANC3+LPONASF+RPONASF)
            PFILT = SF (I) * FILTCT * (B1ASF*Q1(I, KBM1)+B2ASF*Q2(I, &
           & KBM1)+B3ASF*Q3(I, KBM1)+LPOPASF+RPOPASF)
       !
       ! Maximum potential ASSIMILATED carbon
       !    sfa1 - assimilation rate cyanobacteria
       !    sfa2 - assimilation rate diatoms
       !    sfa3 - assimilation rate greens
       !    sfa4 - assimilation rate lpoc
       !    sfa5 - assimilation rate rpoc
       !
            GCMAX = SF (I) * FILTCT * (SFA1X*B1ASF+SFA2X*B2ASF+&
           & SFA3X*B3ASF+SFA4X*LPOCASF+SFA5X*RPOCASF)
       !
       ! Limit carbon growth to a max. theoret. value (mg C/m2.d)
            SFGC = RELING * GCMAX
       !
       !     NOW CHECK TO SEE IF THERE IS A NUTRIENT LIMITATION ON GROWTH
       !
       !        COMPUTE GROWTH RATE FOR NITROGEN
       !
            GN = SF (I) * FILTCT * (SFA1X*B1ASF*ANC1+SFA2X*B2ASF*ANC2+&
           & SFA3X*B3ASF*ANC3+SFA4X*LPONASF+SFA5X*RPONASF)
       !
       !        COMPUTE GROWTH RATE FOR PHOSPHORUS
       !
            GP = SF (I) * FILTCT * (SFA1X*B1ASF*Q1(I, &
           & KBM1)+SFA2X*B2ASF*Q2(I, KBM1)+SFA3X*B3ASF*Q3(I, &
           & KBM1)+SFA4X*LPOPASF+SFA5X*RPOPASF)
       !
       !        REQUIREMENTS FOR C, N AND P ARE DIFFERENT, SEE WHICH ONE IS
       !        LIMITING BY USING ANIMAL STOICHIOMETRY
       !          SFCN - C/N RATIO   SFCP - C/P RATIO
       !        CONVERT TO EQUIVALENT CARBON UNITS TO PERMIT DIRECT COMPARISON
       !
            GN = GN * SFCNX
            GP = GP * SFCPX
            GMIN = Min (GN, GP)
       !
       !        COMPARE TO CARBON GROWTH RATE - IF NUTRIENT LIMITED
       !        THEN RECOMPUTE CARBON ASSIMILATION EFFICIENCIES USING
       !        LIMITATION.  THIS HAS THE EFFECT OF LIMITING THE UPTAKE
       !        OF CARBON TO MATCH ANIMAL STOICHIOMETRY
       !
            If (GMIN < SFGC) Then
               SFAC1 = SFA1X * GMIN / SFGC
               SFAC2 = SFA2X * GMIN / SFGC
               SFAC3 = SFA3X * GMIN / SFGC
               SFAC4 = SFA4X * GMIN / SFGC
               SFAC5 = SFA5X * GMIN / SFGC
               SFGC = GMIN
            Else
               SFAC1 = SFA1X
               SFAC2 = SFA2X
               SFAC3 = SFA3X
               SFAC4 = SFA4X
               SFAC5 = SFA5X
            End If
       !
       !        KNOWING HOW MUCH CARBON WAS TAKEN UP, COMPUTE N AND P UPTAKE
       !        KNOWING HOW MUCH N AND P WAS FILTERED, COMPUTE HOW MUCH GOES
       !        TO FECES AND PSEUDOFECES
       !
       !         NITROGEN ASSIMILATION RATES
            If (GN > SFGC) Then
               SFAN1 = SFA1X * SFGC / GN
               SFAN2 = SFA2X * SFGC / GN
               SFAN3 = SFA3X * SFGC / GN
               SFAN4 = SFA4X * SFGC / GN
               SFAN5 = SFA5X * SFGC / GN
            Else
               SFAN1 = SFA1X
               SFAN2 = SFA2X
               SFAN3 = SFA3X
               SFAN4 = SFA4X
               SFAN5 = SFA5X
            End If
       !
       !        PHOSPHORUS ASSIMILATION RATES
            If (GP > SFGC) Then
               SFAP1 = SFA1X * SFGC / GP
               SFAP2 = SFA2X * SFGC / GP
               SFAP3 = SFA3X * SFGC / GP
               SFAP4 = SFA4X * SFGC / GP
               SFAP5 = SFA5X * SFGC / GP
            Else
               SFAP1 = SFA1X
               SFAP2 = SFA2X
               SFAP3 = SFA3X
               SFAP4 = SFA4X
               SFAP5 = SFA5X
            End If
       !
       !        C FECES - C NOT ASSILMIATED     mg / m2 /d
       ! MBM 970930 Refract feces added
            CFECES = SF (I) * FILTCT * RELING * ((1.-SFAC1)*B1ASF+(1.-&
           & SFAC2)*B2ASF+(1.-SFAC3)*B3ASF+(1.-SFAC4)*LPOCASF)
       !MBM .       +  (1.-SFAC5)*RPOCASF )
       !
            RCFECES = SF (I) * FILTCT * RELING * (1.-SFAC5) * RPOCASF
       !
       !        N FECES - N NOT ASSILMIATED
            NFECES = SF (I) * FILTCT * RELING * ((1.-SFAN1)*B1ASF*ANC1+&
           & (1.-SFAN2)*B2ASF*ANC2+(1.-SFAN3)*B3ASF*ANC3+(1.-&
           & SFAN4)*LPONASF)
       !MBM .       + (1.-SFAN5)*RPONASF )
       !
            RNFECES = SF (I) * FILTCT * RELING * (1.-SFAN5) * RPONASF
       !
       !        P FECES - P NOT ASSILMIATED
            PFECES = SF (I) * FILTCT * RELING * ((1.-SFAP1)*B1ASF*Q1(I, &
           & KBM1)+(1.-SFAP2)*B2ASF*Q2(I, KBM1)+(1.-SFAP3)*B3ASF*Q3(I, &
           & KBM1)+(1.-SFAP4)*LPOPASF)
       !MBM .       + (1.-SFAP5)*RPOPASF )
       !
            RPFECES = SF (I) * FILTCT * RELING * (1.-SFAP5) * RPOPASF
       !
       !
            CPSFEC = Max (0.0, CFILT-SFGC-CFECES-RCFECES)
            RCPSFEC = CPSFEC * RPOCASF / &
           & (B1ASF+B2ASF+B3ASF+LPOCASF+RPOCASF)
            CPSFEC = CPSFEC - RCPSFEC
       !
            NPSFEC = Max (0.0, NFILT-SFGC/SFCNX-NFECES-RNFECES)
            RNPSFEC = NPSFEC * RPONASF / &
           & (ANC1*B1ASF+ANC2*B2ASF+ANC3*B3ASF+LPONASF+RPONASF)
            NPSFEC = NPSFEC - RNPSFEC
       !
            PPSFEC = Max (0.0, PFILT-SFGC/SFCPX-PFECES-RPFECES)
            RPPSFEC = PPSFEC * RPOPASF / (Q1(I, KBM1)*B1ASF+Q2(I, &
           & KBM1)*B2ASF+Q3(I, KBM1)*B3ASF+LPOPASF+RPOPASF)
            PPSFEC = PPSFEC - RPPSFEC
       !
       !        APPORTION FECES TO DEPOSITIONAL FLUXES.  FECES ARE PRODUCED
       !        FROM FOOD TAKEN IN BUT NOT ASSIMILATED.  (NOTE: HERE THE
       !        NUTRIENT LIMITED ASSIMILATION RATES (IF APPROPRIATE) ARE USED
       !        ALSO ASSUME THAT PREDATION LOSSES END UP IN LABILE POC,PON,POP POOLS
       !
       ! Susp fdr particulate fluxes (local) MG/M2/DAY
       !
       ! MBM 970930 --> Error in RPO{C,N,P} logic below since RPOM is CONSUMED with
       !            Assim Effy SFA5. Flux ALREADY accounted for in Feces,
       !            Pseudofeces. New logic for RFeces, RPSFec added.
       !            Prev. error resulted in 1) double-acctg for RPOM (SF) flux,
       !            2) Converting Jpom_ref --> labile splits (90%) in mod_sed module
       !        NOW: Labile and refractory conponents of feces are separated.
       !
       ! LPOC
            SF_LPOCX (I) = CFECES + CPSFEC + (PREDCT*SF(I)+RMORT) * SF &
           & (I)
       !
       ! LPON
            SF_LPONX (I) = NFECES + NPSFEC + (PREDCT*SF(I)+RMORT) * SF &
           & (I) / SFCNX
       !
       ! LPOP
            SF_LPOPX (I) = PFECES + PPSFEC + (PREDCT*SF(I)+RMORT) * SF &
           & (I) / SFCPX
       !
       ! RPOC
            SF_RPOCX (I) = RCFECES + RCPSFEC
       ! RPON
            SF_RPONX (I) = RNFECES + RNPSFEC
       ! RPOP
            SF_RPOPX (I) = RPFECES + RPPSFEC
       !
       ! Particulate, unavailable Si uptaken by suspension feeder (positive leaving water column)
            SF_SUX (I) = FILTCT * SF (I) * SUASF
       ! Inorganic solids
            SF_SSIX (I) = FILTCT * SF (I) * SSIASF * WSSNET (I) / WSS &
           & (I, 1)
       ! P sorbed onto SSI (PIP)
            PF = KADPO4 * SSI (I, KBM1) / (1.+KADPO4*SSI(I, KBM1))
            SF_PIPX (I) = FILTCT * SF (I) * WSSNET (I) / WSS (I, 1) * &
           & PF * Max (0., PO4(I, KBM1))!generation of PO4 due to filtration
       !of SSI
       ! Si sorbed onto SSI (Si removed from water column)
            PF = KADPO4 * SSI (I, KBM1) / (1.+KADPO4*SSI(I, KBM1))
            SF_SASX (I) = FILTCT * SF (I) * PF * SIAT (I, KBM1) * &
           & WSSNET (I) / WSS (I, 1)
       !
       !        CORRECT DERIVATIVES FOR LOSSES
       !        NOTE: TOTAL LOSS IS TOTAL FILTERING RATE
       ! note: water quality derivative units are G/M3/S
       !       therefore divide by 1000. to get g from mg
       !       and divide by 86400. to get from days to seconds
       !       BL(KWC,3) is the deltaZ of the bottom layer wc box in meters !replaced by D(I)*DZ2D(I,KBM1)
       !
            DTB1 (I, KBM1) = DTB1 (I, KBM1) - (FILTCT*SF(I)*B1ASF) / &
           & (D(I)*DZ2D(I,KBM1)) / 1000. / 86400. ! units: g/m3/s
       !
            DTB2 (I, KBM1) = DTB2 (I, KBM1) - (FILTCT*SF(I)*B2ASF) / &
           & (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
            DTB3 (I, KBM1) = DTB3 (I, KBM1) - (FILTCT*SF(I)*B3ASF) / &
           & (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
            DTPIB1 (I, KBM1) = DTPIB1 (I, KBM1) - Q1 (I, KBM1) * &
           & (FILTCT*SF(I)*B1ASF) / (D(I)*DZ2D(I,KBM1)) / 1000. / 86400. ! units: g/m3/s
       !
            DTPIB2 (I, KBM1) = DTPIB2 (I, KBM1) - Q2 (I, KBM1) * &
           & (FILTCT*SF(I)*B2ASF) / (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
            DTPIB3 (I, KBM1) = DTPIB3 (I, KBM1) - Q3 (I, KBM1) * &
           & (FILTCT*SF(I)*B3ASF) / (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
            DTLPOC (I, KBM1) = DTLPOC (I, KBM1) - &
           & (FILTCT*SF(I)*LPOCASF) / (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
            DTLPOP (I, KBM1) = DTLPOP (I, KBM1) - &
           & (FILTCT*SF(I)*LPOPASF) / (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
            DTLPON (I, KBM1) = DTLPON (I, KBM1) - &
           & (FILTCT*SF(I)*LPONASF) / (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
            DTRPOC (I, KBM1) = DTRPOC (I, KBM1) - &
           & (FILTCT*SF(I)*RPOCASF) / (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
            DTRPOP (I, KBM1) = DTRPOP (I, KBM1) - &
           & (FILTCT*SF(I)*RPOPASF) / (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
            DTRPON (I, KBM1) = DTRPON (I, KBM1) - &
           & (FILTCT*SF(I)*RPONASF) / (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
            DTSSI (I, KBM1) = DTSSI (I, KBM1) - SF_SSIX (I) / &
           & (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
       !Unit: (mgSi/m^2/day) /(m) /1000 /86400 ==> gSi/m^3/sec
            DTSIUPB (I, KBM1) = DTSIUPB (I, KBM1) - SF_SUX (I) / &
           & (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
            DTSIAT (I, KBM1) = DTSIAT (I, KBM1) - SF_SASX (I) / &
           & (D(I)*DZ2D(I,KBM1)) / 1000. / 86400.
       !
            DTPO4 (I, KBM1) = DTPO4 (I, KBM1) + SF_PIPX (I) / &
           & (D(I)*DZ2D(I,KBM1)) / 1000. / 86400. !Wen Long: Sign should be negative?
       !Wen Long: Should this be removing PO4 from water column?
       !
       ! calculate respiration fluxes
       !
            JNH4SFX (I) = RESPCT * SF (I) / SFCNX ! units: mg/m2/d
            JPO4SFX (I) = RESPCT * SF (I) / SFCPX
       !
       !Oxygen Demand due to SF respiration, Unit: (1/d)*(mgC/m^2)*(mgO2/mgC)/1000 = mgO2/m^2/1000 ==> gO2/m^2/day
            SODSFX (I) = RESPCT * SF (I) * 2.67 / 1000. !SOD due to respiration of suspension feeder (gO2/m^2/day)
       !
       ! Si recycled from diatom consumption             [mg Si / m2 / day]
            JSASFX (I) = FILTCT * SF (I) * &
           & (B1ASF*ASC1+B2ASF*ASC2+B3ASF*ASC3) * FRSASF !  CFC
            JSUSFX (I) = FILTCT * SF (I) * &
           & (B1ASF*ASC1+B2ASF*ASC2+B3ASF*ASC3) * (1.-FRSASF)!  CFC
       !
       !        TAKE INTEGRATION STEP FOR SUSPENSION FEEDERS
       !
            DTDAY = DLT / 86400.
       !
            SF (I) = SF (I) + DTDAY * &
           & (SFGC-RESPCT*SF(I)-RMORT*SF(I)-PREDCT*SF(I)*SF(I))
       !
       ! DONT LET GO NEGATIVE
       !
            SF (I) = Max (SF(I), 0.1)
       !
       !MBM 961226 diagnostic arrays
            SFFILTCT (I, N) = FILTCT
            SFRESPCT (I, N) = RESPCT
            SFPREDCT (I, N) = PREDCT
            SFRATEI (I, N) = RATEI
            SFGCMAX (I, N) = GCMAX
            SFGMIN (I, N) = GMIN
            SFCFECES (I, N) = CFECES
            SFNFECES (I, N) = NFECES
            SFPFECES (I, N) = PFECES
            SFCPSFEC (I, N) = CPSFEC
            SFNPSFEC (I, N) = NPSFEC
            SFPPSFEC (I, N) = PPSFEC
            SF_SFGC (I, N) = SFGC
            SF_RESP (I, N) = RESPCT * SF (I)
            SF_PRED (I, N) = PREDCT * SF (I) * SF (I)
            SF_RMORT (I, N) = RMORT * SF (I)
       !
         End Do
255      Continue
    !
    ! CONVERT BACK INTO SAVED ARRAYS
    !
         Do I = 1, MLOC
            SFEED (I, N) = SF (I)
         End Do
    !
    ! ADD UP EACH SPECIES FLUX CONTRIBUTION
    ! Partic. fluxes (labile, refrac. C,N,P,Inorg solids,unavail Si,sorbed P,Si)
    !
         Do I = 1, MLOC
       !
       !nutrient flux into water column
            JNH4SF (I) = JNH4SF (I) + JNH4SFX (I)!NH4 to water column
            JPO4SF (I) = JPO4SF (I) + JPO4SFX (I)!PO4 to water column
            SODSF (I) = SODSF (I) + SODSFX (I)!positive removing DOXG from water column
            JSASF (I) = JSASF (I) + JSASFX (I)!avaiable dissolved silicate to water column
       !
       !particulate flux into sediments
            JSUSF (I) = JSUSF (I) + JSUSFX (I)!Particulate biogenic silica
       !
            SF_JLPOC (I) = SF_JLPOC (I) + SF_LPOCX (I)!LPOC
            SF_JLPON (I) = SF_JLPON (I) + SF_LPONX (I)!LPON
            SF_JLPOP (I) = SF_JLPOP (I) + SF_LPOPX (I)!LPOP
            SF_JRPOC (I) = SF_JRPOC (I) + SF_RPOCX (I)!RPOC
            SF_JRPON (I) = SF_JRPON (I) + SF_RPONX (I)!RPON
            SF_JRPOP (I) = SF_JRPOP (I) + SF_RPOPX (I)!RPOP
       !
            SF_SSI (I) = SF_SSI (I) + SF_SSIX (I)!solids (mgSolids/m^2/day)
            SF_SU (I) = SF_SU (I) + SF_SUX (I)!Particulate biogenic silica leaving water colum into sediment (suspension feeder) (mgSi/m^2/day)
            SF_SA (I) = SF_SA (I) + SF_SASX (I)
            SF_PIP (I) = SF_PIP (I) + SF_PIPX (I)
       !
       !
       !make them accessible by mod_sed.F
       !
            JLPOC_SED_SF (I) = SF_JLPOC (I)!sediment module will pick up JLPOC_SED_SF
            JLPON_SED_SF (I) = SF_JLPON (I)!sediment module will pick up JLPON_SED_SF
            JLPOP_SED_SF (I) = SF_JLPOP (I)!sediment module will pick up JLPOP_SED_SF
       !
            JRPOC_SED_SF (I) = SF_JRPOC (I)!sediment module will pick up JRPOC_SED_SF
            JRPON_SED_SF (I) = SF_JRPON (I)!sediment module will pick up JRPON_SED_SF
            JRPOP_SED_SF (I) = SF_JRPOP (I)!sediment module will pick up JRPOP_SED_SF
       !
            JNH4_SED_SF (I) = JNH4SF (I)!sediment module will pick up JNH4_SED_SF
            JPO4_SED_SF (I) = JPO4SF (I)!sediment module will pick up JPO4_SED_SF
            SOD_SED_SF (I) = SODSF (I)!sediment module will pick up SOD_SED_SF
            JSA_SED_SF (I) = JSASF (I)!sediment module will pick up JSA_SED_SF
            JSU_SED_SF (I) = JSUSF (I)!sediment module will pick up JSU_SED_SF
       !
            SSI_SED_SF (I) = SF_SSI (I)!sediment module will pick up SSI_SED_SF
            SU_SED_SF (I) = SF_SU (I)!sediment module will pick up SU_SED_SF
            SA_SED_SF (I) = SF_SA (I)!sediment module will pick up SA_SED_SF
            PIP_SED_SF (I) = SF_PIP (I)!sediment module will pick up PIP_SED_SF
       !
         End Do
    !
         Return
      End Subroutine SUSPFEED
  !
      Subroutine SF_ALLOC
    !
         Implicit None
    !
         Allocate (SFEED(MTLOC, NSSFP))
         SFEED = 0.0
    !
         Allocate (SF_JLPOC(MTLOC))
         SF_JLPOC = 0.0
         Allocate (SF_JLPON(MTLOC))
         SF_JLPON = 0.0
         Allocate (SF_JLPOP(MTLOC))
         SF_JLPOP = 0.0
         Allocate (SF_JRPOC(MTLOC))
         SF_JRPOC = 0.0
         Allocate (SF_JRPON(MTLOC))
         SF_JRPON = 0.0
         Allocate (SF_JRPOP(MTLOC))
         SF_JRPOP = 0.0
         Allocate (JNH4SF(MTLOC))
         JNH4SF = 0.0
         Allocate (JPO4SF(MTLOC))
         JPO4SF = 0.0
         Allocate (SODSF(MTLOC))
         SODSF = 0.0
         Allocate (JSASF(MTLOC))
         JSASF = 0.0
         Allocate (JSUSF(MTLOC))
         JSUSF = 0.0
         Allocate (SF_SSI(MTLOC))
         SF_SSI = 0.0
         Allocate (SF_SU(MTLOC))
         SF_SU = 0.0
         Allocate (SF_SA(MTLOC))
         SF_SA = 0.0
         Allocate (SF_PIP(MTLOC))
         SF_PIP = 0.0
    !
         Allocate (SFFILTCT(MTLOC, NSSFP))
         SFFILTCT = 0.0
         Allocate (SFRESPCT(MTLOC, NSSFP))
         SFRESPCT = 0.0
         Allocate (SFPREDCT(MTLOC, NSSFP))
         SFPREDCT = 0.0
         Allocate (SFRATEI(MTLOC, NSSFP))
         SFRATEI = 0.0
         Allocate (SFGCMAX(MTLOC, NSSFP))
         SFGCMAX = 0.0
         Allocate (SFGMIN(MTLOC, NSSFP))
         SFGMIN = 0.0
         Allocate (SFGCING(MTLOC, NSSFP))
         SFGCING = 0.0
         Allocate (SFCFECES(MTLOC, NSSFP))
         SFCFECES = 0.0
         Allocate (SFNFECES(MTLOC, NSSFP))
         SFNFECES = 0.0
         Allocate (SFPFECES(MTLOC, NSSFP))
         SFPFECES = 0.0
         Allocate (SFCPSFEC(MTLOC, NSSFP))
         SFCPSFEC = 0.0
         Allocate (SFNPSFEC(MTLOC, NSSFP))
         SFNPSFEC = 0.0
         Allocate (SFPPSFEC(MTLOC, NSSFP))
         SFPPSFEC = 0.0
    !
         Allocate (SF_SFGC(MTLOC, NSSFP))
         SF_SFGC = 0.0
         Allocate (SF_RESP(MTLOC, NSSFP))
         SF_RESP = 0.0
         Allocate (SF_PRED(MTLOC, NSSFP))
         SF_PRED = 0.0
         Allocate (SF_RMORT(MTLOC, NSSFP))
         SF_RMORT = 0.0
         Allocate (ASF_SFGC(MTLOC, NSSFP))
         ASF_SFGC = 0.0
         Allocate (ASF_RESP(MTLOC, NSSFP))
         ASF_RESP = 0.0
         Allocate (ASF_PRED(MTLOC, NSSFP))
         ASF_PRED = 0.0
         Allocate (ASF_RMORT(MTLOC, NSSFP))
         ASF_RMORT = 0.0
    !
         Allocate (SEDTYPE(MTLOC, NSSFP))
         SEDTYPE = 0.0
    !
         Allocate (SFEED_GL(1:MGL, 1:NSSFP))
         SFEED_GL = 0.0
         Allocate (SEDTYPE_GL(1:MGL, 1:NSSFP))
         SEDTYPE_GL = 0.0
    !
         Return
      End Subroutine SF_ALLOC
  !
      Subroutine SF_DEALLOC
    !
         Implicit None
    !
         Deallocate (SFEED)
         Deallocate (SF_JLPOC)
         Deallocate (SF_JLPON)
         Deallocate (SF_JLPOP)
         Deallocate (SF_JRPOC)
         Deallocate (SF_JRPON)
         Deallocate (SF_JRPOP)
         Deallocate (JNH4SF)
         Deallocate (JPO4SF)
         Deallocate (SODSF)
         Deallocate (JSASF)
         Deallocate (JSUSF)
         Deallocate (SF_SSI)
         Deallocate (SF_SU)
         Deallocate (SF_SA)
         Deallocate (SF_PIP)
    !
         Deallocate (SFFILTCT)
         Deallocate (SFRESPCT)
         Deallocate (SFPREDCT)
         Deallocate (SFRATEI)
         Deallocate (SFGCMAX)
         Deallocate (SFGMIN)
         Deallocate (SFGCING)
         Deallocate (SFCFECES)
         Deallocate (SFNFECES)
         Deallocate (SFPFECES)
         Deallocate (SFCPSFEC)
         Deallocate (SFNPSFEC)
         Deallocate (SFPPSFEC)
    !
         Deallocate (SF_SFGC)
         Deallocate (SF_RESP)
         Deallocate (SF_PRED)
         Deallocate (SF_RMORT)
         Deallocate (ASF_SFGC)
         Deallocate (ASF_RESP)
         Deallocate (ASF_PRED)
         Deallocate (ASF_RMORT)
    !
         Deallocate (SEDTYPE)
         Deallocate (SFEED_GL)
         Deallocate (SEDTYPE_GL)
    !
         Return
    !
      End Subroutine SF_DEALLOC
  !
      Subroutine SF_READ
    !
         Use MOD_FILEINFO, Only: SFI
         Use MOD_WQMINIT, Only: SFIFN
         Use MOD_CONTROL, Only: SERIAL, PAR
# if defined (MULTIPROCESSOR)
         Use MOD_PAR, Only: NGID, NHN, HN_LST
# endif
         Implicit None
         Integer :: N
         Integer :: B
         Integer :: NUMBOX
    !
         Open (SFI, File=SFIFN, Status='OLD')
         Read (SFI, 1032)
         Read (SFI, 1170) NSPECIES ! # of susp. feeders species to be modeled
         Read (SFI, 1005) FRSASF ! Recycle fraction of grazed diatom Si
         Do N = 1, NSPECIES
            Read (SFI, 1032)
            Read (SFI, 1003) FILT (N), SFA1 (N), SFA2 (N), SFA3 (N), &
           & SFA4 (N), SFA5 (N), MAXING (N), SFCN (N)
            Read (SFI, 1002) SFCP (N), SFRESP (N), SFPRED (N), SFTMN &
           & (N), THTAFILT (N), THTARESP (N), THTAPRED (N)
            Read (SFI, 1003) XKPO2 (N), SFTD (N), FILTFACT (N), &
           & RESPFACT (N), SFDOH (N), SFDOQ (N), SFATURB (N), SFBTURB &
           & (N)
       !
       !MBM 971021 added turbidity sensitivity parameters ^^^^
       !
         End Do
    !
         Read (SFI, 1032)
         Do B = 1, MGL
            Read (SFI, 1001) NUMBOX, (SEDTYPE_GL(B, N), SFEED_GL(B, N), &
           & N=1, NSPECIES)!Check here
         End Do
    !
         Close (SFI)
1001     Format (I6, 10(I4, F10.1))
1002     Format (7 F10.3)
1003     Format (8 F10.3)
1005     Format (/ / F10.1)
1032     Format (/)
1170     Format (://(2 I10))
    !
    !
         If (SERIAL) Then
            SEDTYPE = SEDTYPE_GL
            SFEED = SFEED_GL
         End If
    !
# if defined (MULTIPROCESSOR)
         If (PAR) Then
            Do N = 1, NSPECIES
               Do B = 1, MLOC
                  SEDTYPE (B, N) = SEDTYPE_GL (NGID(B), N)
                  SFEED (B, N) = SFEED_GL (NGID(B), N)
               End Do
               Do B = 1, NHN
                  SEDTYPE (B+MLOC, N) = SEDTYPE_GL (HN_LST(B), N)
                  SFEED (B+MLOC, N) = SFEED_GL (HN_LST(B), N)
               End Do
            End Do
         End If
# endif

    !
         Return
      End Subroutine SF_READ
  !
      Subroutine SF_INIT
         Implicit None
    !: moved these here wqm_inputs.F
         SF_JLPOC = 0.0
         SF_JLPON = 0.0
         SF_JLPOP = 0.0
    !
         SF_JRPOC = 0.0
         SF_JRPON = 0.0
         SF_JRPOP = 0.0
    !
         JNH4SF = 0.0
         JPO4SF = 0.0
         SODSF = 0.0
         JSUSF = 0.0
         JSASF = 0.0
         SF_SSI = 0.0
         SF_SU = 0.0
         SF_SA = 0.0
         SF_PIP = 0.0
    !
      End Subroutine SF_INIT
  !
End Module MOD_SF
