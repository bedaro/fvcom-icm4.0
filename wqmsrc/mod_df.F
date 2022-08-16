!mod_df.f
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
!** Deposition Feeder Module                                           **
!**                                                                    **
!**    Wen Long,10/05/2014                                             **
!**                                                                    **
!************************************************************************
!**                                                                    **
!** Inputs:                                                            **
!**                                                                    **
!**            Required inputs for benthic algae model                 **
!**                                                                    **
!** Outpus:                                                            **
!**                                                                    **
!**            Required inputs for benthic algae model                 **
!**                                                                    **
!************************************************************************
!
Module MOD_DF
  !
      Use MOD_WQM, Only: CTEMP, DFIFN, DFOFN, DLT, DOXG, MGL, MTLOC !,         &!
  !
  !
      Use MOD_FILEINFO, Only: DFI, DFO
      Use MOD_SIZES, Only: MGL
  !
      Use MOD_LIMS, Only: MLOC, KBM1, MTLOC
  !
  !Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: SERIAL, MSR, PAR
  !
      Use MOD_PREC, Only: SP
  !
      Implicit None
      Save
  !
      Real (SP) :: XKMI0, ING0, THTAI0, R, THTAR, BETA, THBETA, AMCN, &
     & AMCP, AA1, AA2, XKMG1, XKMG2, XKBO2, TDD
  !
  !
      Real (SP) :: DFDOH, DFDOQ, DOLOW !Never used
  !
  !
  !Moved the following from mod_wqm.F
  !
      Logical :: DFEEDER, HYPOXFX_DF !Flag for hypoxia effect on deposition feeder
  !
      Logical :: DF_INITIALIZED !flag to check if deposition feeder is initialized
  !
      Real (SP) :: XKR, XKI0, XKBETA
  !
      Real (SP) :: RDD, RMORT, XPOC1LIM, TEMP, XPOC2LIM, DFSOD !SOD due to deposition feeder respiration
  !
      Real (SP) :: DFEED, DFEEDM1, ZHTAI0 (350), ZHTAR (350), ZHTABETA &
     & (350)
  !
      Real (SP), Allocatable :: DFEEDM1S (:), DF_GROW (:), DF_GROW_POC1 &
     & (:), DF_GROW_POC2 (:), DF_RESP (:), DF_PRED (:), DF_MORT (:), &
     & DF_SOD (:), ADF_GROW (:), ADF_RESP (:), ADF_PRED (:), ADF_MORT &
     & (:)
  !
      Real (SP), Allocatable :: DFEEDM1S_GL (:)
  !
      Logical :: DFEEDER_OUTPUT !flag for deposition feeder model output
  !
      Real (SP) :: DLT_DF !time step (days) for deposition feeder
  !
Contains
  !
  !subroutines:
  !subroutine DF_INIT()
  !subroutine DF_ALLOC()
  !subroutine DF_DEALLOC()
  !subroutine DF_READ()
  !subroutine DF_CALC()
  !subroutine DF_INT()
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   DF_INIT                           **
  !********************************************************************************
      Subroutine DF_INIT
    !
    !Integer :: I  !
    !
         If ( .Not. DFEEDER) Then
            ING0 = 0.0
            R = 0.0
            BETA = 0.0
         End If
    !
         DF_INITIALIZED = .True.
    !
      End Subroutine DF_INIT
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   DF_ALLOC                          **
  !********************************************************************************
      Subroutine DF_ALLOC
    !
    !
         Allocate (DFEEDM1S(MTLOC))
         DFEEDM1S = 0.0
         Allocate (DF_GROW(MTLOC))
         DF_GROW = 0.0
         Allocate (DF_GROW_POC1(MTLOC))
         DF_GROW_POC1 = 0.0
         Allocate (DF_GROW_POC2(MTLOC))
         DF_GROW_POC2 = 0.0
         Allocate (DF_RESP(MTLOC))
         DF_RESP = 0.0
         Allocate (DF_PRED(MTLOC))
         DF_PRED = 0.0
         Allocate (DF_MORT(MTLOC))
         DF_MORT = 0.0
         Allocate (DF_SOD(MTLOC))
         DF_SOD = 0.0
         Allocate (ADF_GROW(MTLOC))
         ADF_GROW = 0.0
         Allocate (ADF_RESP(MTLOC))
         ADF_RESP = 0.0
         Allocate (ADF_PRED(MTLOC))
         ADF_PRED = 0.0
         Allocate (ADF_MORT(MTLOC))
         ADF_MORT = 0.0
         Allocate (DFEEDM1S_GL(MGL))
         DFEEDM1S_GL = 0.0
    !
    !
      End Subroutine DF_ALLOC
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   DF_DEALLOC                        **
  !********************************************************************************
  !
      Subroutine DF_DEALLOC
    !
    !Moved here from wqm_main.F
    !
    !
         If (ALLOCATED(DFEEDM1S)) DEALLOCATE (DFEEDM1S)
         If (ALLOCATED(DF_GROW)) DEALLOCATE (DF_GROW)
         If (ALLOCATED(DF_GROW_POC1)) DEALLOCATE (DF_GROW_POC1)
         If (ALLOCATED(DF_GROW_POC2)) DEALLOCATE (DF_GROW_POC2)
         If (ALLOCATED(DF_RESP)) DEALLOCATE (DF_RESP)
         If (ALLOCATED(DF_PRED)) DEALLOCATE (DF_PRED)
         If (ALLOCATED(DF_MORT)) DEALLOCATE (DF_MORT)
         If (ALLOCATED(DF_SOD)) DEALLOCATE (DF_SOD)
         If (ALLOCATED(ADF_GROW)) DEALLOCATE (ADF_GROW)
         If (ALLOCATED(ADF_RESP)) DEALLOCATE (ADF_RESP)
         If (ALLOCATED(ADF_PRED)) DEALLOCATE (ADF_PRED)
         If (ALLOCATED(ADF_MORT)) DEALLOCATE (ADF_MORT)
    !
         If (ALLOCATED(DFEEDM1S_GL)) DEALLOCATE (DFEEDM1S_GL)
    !
      End Subroutine DF_DEALLOC
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   D F _ R E A D                     **
  !********************************************************************************
      Subroutine DF_READ
    !
         Implicit None
         Save
    !
    !***** Variable declarations
    !
         If (MSR) Then
            If (DFEEDER_OUTPUT) Then
               Open (DFO, File=DFOFN)
            End If
         End If
    !
         Open (DFI, File=DFIFN, Status='OLD')
    !
    !Read the deposition feeder input file
         Read (DFI, 1015, Err=10100) XKMI0, ING0, THTAI0, R, THTAR, &
        & BETA, THBETA
         Read (DFI, 1015, Err=10100) AMCN, AMCP, AA1, AA2, XKMG1, XKMG2
         Read (DFI, 1015, Err=10100) XKBO2, TDD, DOLOW, DFDOH, DFDOQ
    !
1015     Format (/ / 8 X, 8 F8.1)
    !
         Close (DFI)
    !
         DF_INITIALIZED = .False.
    !
    !***** Error output FORMAT's
    !
    !
    !
         Return
    !
    !***** Error traps
    !
10100    Continue
    !
         If (MSR) Then
            If (DFEEDER_OUTPUT) Then
               Write (DFO, 3010)
               Close (DFO)
            End If
         End If
3010     Format (/ ' Read error in Deposition Feeder input deck')
         Stop 'DF_READ'
    !
         Return
      End Subroutine DF_READ
  !
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   D F _ C A L C                     **
  !**                       Deposition Feeder Calculations                       **
  !********************************************************************************
  !
      Subroutine DF_CALC
         Use MOD_SED_DF_EXCHANGE_VARS, Only: POC1TM1S_SED_DF, &
        & POC2TM1S_SED_DF, M2_SED_DF
    !
         Implicit None
    !
         Save
         Integer :: I, JT, ITEMP, KWC
         Real (SP) :: TEMP20, TEMP202
         Real (SP) :: TEMPDAY !temperature day of year
         Real (SP) :: LOGICT
         Real (SP) :: O20 !oxygen concentration of overylying water (mgO2/L)
         Real (SP) :: POC1TM1S_DF, POC2TM1S_DF !(mgC/m^3)
    !
    !Look up table for tempreature control of deposition feeder rates
    !
         Do JT = 1, 350
            TEMP = REAL (JT-1) / 10. + 0.05
            TEMP20 = TEMP - 20.
            TEMP202 = TEMP20 / 2.
       !
            ZHTAI0 (JT) = ING0 * THTAI0 ** TEMP20 ! deposit feeders
            ZHTAR (JT) = R * THTAR ** TEMP20 ! deposit feeders
            ZHTABETA (JT) = BETA * THBETA ** TEMP20 ! deposit feeders
         End Do
    !
    !***** Update Depisition feedr concentration
    !
    !
         KWC = KBM1
         Do I = 1, MLOC
       !
            O20 = AMAX1 (DOXG(I, KWC), 0.010)
            TEMPDAY = CTEMP (I)!Sediment temperature (degC) temporary variable
            ITEMP = 10. * TEMPDAY + 1 !calculate temperature days (using temperature to figure out the lookup
       !index in the temperature dependent rates
       !
            POC1TM1S_DF = POC1TM1S_SED_DF (I)!take POC G1 concentration of sediment
            POC2TM1S_DF = POC2TM1S_SED_DF (I)!take POC G2 concentration of sediment
       !
            If (ITEMP < 1) ITEMP = 1 !make sure index >=1
            If (ITEMP > 350) ITEMP = 350 !mkake sure index <=350
       !
            DFEEDM1 = DFEEDM1S (I)!deposition feeder befor integrating to next time step
       !
       ! Deposit feeding ingestion rate
       !
            XKI0 = ZHTAI0 (ITEMP)
       !
       ! respiration rate
       !
            XKR = ZHTAR (ITEMP)! deposit feeders
       !
       ! quadratic predation
       !
            XKBETA = ZHTABETA (ITEMP)! deposit feeders
       !
       !
       !MBM 970109 control switch for hypoxic effects
            If (HYPOXFX_DF) Then
          !
               LOGICT = 1.0 / &
              & (1.0+Exp(Max(1.1*(DFDOH-O20)/(DFDOH-DFDOQ),-25.)))
          !
          !  Reduce ingestion when o2 is low
          !       xki0=xki0*(o20/(o20+xkmi0))
               XKI0 = XKI0 * LOGICT
          !
          ! Mortality due to hypoxia (adds to sediment POM pools)
               RDD = 4.6 / TDD ! ln(1/100) for 99% kill in time tdd
               RMORT = RDD * (1.0-LOGICT)
          !
          ! Reduce predation when o2 low
               XKBETA = XKBETA * O20 / (O20+XKBO2)
            End If
       !
       !
       ! Growth rate limitation by available POC1
       !
            XPOC1LIM = XKMG1 / (POC1TM1S_DF+XKMG1)! deposit feeders
            XPOC2LIM = XKMG2 / (POC2TM1S_DF+XKMG2)! deposit feeders
       !
       ! calculate deposit feeders biomass (in layer 2 of sediments, layer 1 ignored)
       !
       !
            DF_GROW_POC1 (I) = AA1 * XKI0 / (M2_SED_DF*1e+09) * &
           & POC1TM1S_DF * XPOC1LIM * DFEEDM1 !growth based on grazing POC1
            DF_GROW_POC2 (I) = AA2 * XKI0 / (M2_SED_DF*1e+09) * &
           & POC2TM1S_DF * XPOC2LIM * DFEEDM1 !growth based on grazing POC2
       !
            DF_GROW (I) = DF_GROW_POC1 (I) + DF_GROW_POC2 (I)
       !
       !
            DF_RESP (I) = XKR * DFEEDM1
            DF_PRED (I) = XKBETA * DFEEDM1 * DFEEDM1
            DF_MORT (I) = RMORT * DFEEDM1
            DLT_DF = DLT / 86400.0 !day
            DFEED = DFEEDM1 + DLT_DF * &
           & (DF_GROW(I)-DF_RESP(I)-DF_PRED(I)-DF_MORT(I))
       !
       !
       ! dont let it go negative
            DFEED = Max (DFEED, 0.1)
       !
       !calculation of SOD due to deposition feeder
       !
            DFSOD = XKR * DFEEDM1 * 2.667E-3 ! DEPOSIT FEEDERS RESP.
       !
            DF_SOD (I) = DFSOD !This will be returned to sediment module
       !
            DFEEDM1S (I) = DFEED
       !
         End Do
    !
         Return
      End Subroutine DF_CALC
  !
      Subroutine DF_INT ()
    !
         Integer :: I
    !
         Return
      End Subroutine DF_INT
  !
  !
End Module MOD_DF
