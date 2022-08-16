!mod_pH_pCO2.F
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
Module MOD_PH_PCO2
!
      Use MOD_PREC, Only: SP
!
      Use MOD_LIMS, Only: MLOC, KBM1, MTLOC
!
      Use MOD_CONTROL, Only: MSR
!
      Use MOD_WQM, Only: pH, pCO2, TALK, TDIC
!
      Use MOD_CO2SYS, Only: KW, K0, K1, K2, TB, KBb, TF, KF, TS, KS, &
     & IonS, FugFac, CO2star_surf, CO2SYSCONSTANTS       !KSi,    &!
!f90pprUSE MOD_CO2SYS, ONLY :    &
!f90pprKW,    &
!f90pprK0,    &
!f90pprK1,    &
!f90pprK2,    &
!f90pprTB,    &
!f90pprKBb,   &  !KBb instead of KB because KB is an integer = number of sigma levels (mod_lims.F)
!f90pprTF,    &
!f90pprKF,    &
!f90pprTS,    &
!f90pprKS,    &
!f90pprIonS,  &
!f90ppr!TP,     & ! Not considering protonation of phosphate or sulfate in pH calculation
!f90ppr!KP1,    &!
!f90ppr!KP2,    &!
!f90ppr!KP3,    &!
!f90ppr!TSi,    &!
!f90ppr!KSi,    &!
!f90pprFugFac,  &
!f90pprCO2star_surf, &
!f90pprCO2SYSCONSTANTS
!
      Use MOD_HYDROVARS, Only: RHO
!
!
!
      Implicit None
      Save
!
Contains
  !CALC_PH_PCO2
  !==============================================================================|
  !    Subroutine CALC_PH_PCO2 is taken from:                                             !
  !   CalculatepHfCO2fromTATC from CO2SYS.m   (v1.1 sept 2011)                   !
  !                                                                               !
  !    Subroutine PH_INIT initializes pH as pHGuess (set in mod_CO2SYSconstants.F)!
  !                                                                               !
  !                           !
  !==============================================================================|
!
  !!************************************************************************
  !!**                 S U B R O U T I N E   P H _ I N I T                **
  !!************************************************************************
  !   SUBROUTINE PH_INIT(pH)
  !     USE WQM, ONLY : pH !or whatever module defines pH
  !
  !     pH = 8.0
  !
  !   RETURN
  !   END SUBROUTINE PH_INIT
!
!
  !************************************************************************
  !**              S U B R O U T I N E   C A L C _ P H _ P C O 2         **
  !************************************************************************
!
      Subroutine CALC_PH_PCO2 !(pHlocal,pCO2local)
!
         Implicit None
         Save
!
         Real (SP) :: H, Denom, CAlk, BAlk, OH, FREEtoTOT, Hfree, HSO4, &
        & HF, dTalk, Slope, deltapH, fCO2, H2, TA, TC
         Real (SP) :: auxLB
    !REAL(SP) ::  PhosTop, PhosBot, PAlk, SiAlk
!
         Real (SP), Dimension (0:MTLOC, KBM1) :: CO2star
    !REAL(SP),DIMENSION(0:MTLOC,KBM1),INTENT(INOUT) :: pHlocal
    !REAL(SP),DIMENSION(0:MTLOC,KBM1),INTENT(OUT)   :: pCO2local
!
         Integer :: I, K
!
         Real (SP) :: pHTol = 0.0001 ! tolerance for iterations end
         Real (SP) :: ln10 = Log (10.0)
!
         deltapH = pHTol + 1.
         CO2star = 0.0
!
    !!Constants for pH calculations
         Call CO2SYSCONSTANTS
!
!
    !in the case that WHICHK1K2=7 (Peng et al), alkalinity is corrected
    !! by PO4 in CO2SYSCONSTANT:
    !! TALK = TALK - total phosphate.
    !! Protonation of silicate and phosphate
    !!              is ignored in the pH calculation, assumed negligible 
!
!
    !!**TALK and TDIC must be in mol/kgSW but in the model they are in mmol/m3
    !!  Used to Assume density = 1000 kg/m^3 (such that mmol/m3=umol/kg)
    !!            and divide TALK and TDIC by 10^6 to get mol/kg
    !!  NOW, use RHO to change mmol/m3 to umol/kg-SW (TDIC/rho*10^3), and then divide by 10^6  (8 jan 2016)
    !
    !!start calculation
         Do K = 1, KBM1
            Do I = 1, MLOC

               TC = TDIC (I, K) / RHO (I, K) * 1.E-3_SP !conversion of mmol/m3 to mol/kg,
               TA = TALK (I, K) / RHO (I, K) * 1.E-3_SP !using RHO in kg/m3
!
!
               Do WHILE (Abs(deltapH) > pHTol)
                  H = 10.0_SP ** (-pH(I, K))
                  Denom = H * H + K1 (I, K) * H + K1 (I, K) * K2 (I, K)
                  CAlk = TC * K1 (I, K) * (H+2*K2(I, K)) / Denom
                  BAlk = TB (I, K) * KBb (I, K) / (KBb(I, K)+H)
                  OH = KW (I, K) / H
             !            PhosTop   = KP1(I,K)*KP2(I,K)*H + 2.0_SP*KP1(I,K)*KP2(I,K)*KP3(I,K) &
             !                           - H*H*H
             !            PhosBot   = H*H*H + KP1(I,K)*H*H + KP1(I,K)*KP2(I,K)*H         &
             !                           + KP1(I,K)*KP2(I,K)*KP3(I,K)
             !            PAlk      = TP(I,K)*PhosTop/PhosBot
             !            SiAlk     = TSi(I,K)*KSi(I,K)/(KSi(I,K) + H)
                  FREEtoTOT = (1.0_SP+TS(I, K)/KS(I, K))! pH scale conversion factor
                  Hfree = H / FREEtoTOT ! for H on the total scale
                  HSO4 = TS (I, K) / (1.0_SP+KS(I, K)/Hfree)! since KS is on the free scale
                  HF = TF (I, K) / (1.0_SP+KF(I, K)/Hfree)! since KF is on the free scale
                  dTalk = TA - CAlk - BAlk - OH + Hfree + HSO4 + HF !&
             !- PAlk - SiAlk                              !ignoring PAlk and SiAlk
!
             !! find Slope dTA/dpH
             !! (this is not exact, but keeps all important terms)
                  Slope = ln10 * (TC*K1(I, K)*H*(H*H+K1(I, K)*K2(I, &
                 & K)+4.0_SP*H*K2(I, K))/Denom/Denom+BAlk*H/(KBb(I, &
                 & K)+H)+OH+H)
                  deltapH = dTalk / Slope ! this is Newton's method
             !! to keep the jump from being too big
                  Do WHILE (Abs(deltapH) > 1.0_SP)
                     deltapH = deltapH / 2.0_SP
                  End Do
!
             !!---Get pH
                  pH (I, K) = pH (I, K) + deltapH ! It is on the same scale as K1 and K2 were calculated...
!
               End Do
!
               deltapH = pHTol + 1._SP !
!
          !!---Calculate pCO2
               H = 10._SP ** (-pH(I, K))
               H2 = H * H
               CO2star (I, K) = TC * H2 / (H2+K1(I, K)*H+K1(I, K)*K2(I, &
              & K))!mole/kg-SW
               fCO2 = CO2star (I, K) / K0 (I, K)!fugacity (untis: atm)
               pCO2 (I, K) = fCO2 / FugFac (I, K) * 1.e6_SP !uatm  !FugFac assumes that the pressure is at one atmosphere, or close to it.
          !but approximation holds at P of natural environments (Weiss, Mar.Chem., 1974)
            End Do
         End Do
!
    !!--Keep CO2star at surface (to compute air-sea CO2 flux in calc_DIC_TALK.F)
         Do I = 1, MLOC
       !CO2star_surf(I) = CO2star(I,1) * 1.e6_SP  !mmol/m3 (assume density=1000 kg-SW/m3) !
            auxLB = CO2star (I, 1) * RHO (I, 1) * 1.e3_SP !!
            CO2star_surf (I) = auxLB !mmol/m3!
         End Do

         Return
      End Subroutine CALC_PH_PCO2
!
End Module MOD_PH_PCO2
