!mod_CO2SYSconstants.F
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
Module MOD_CO2SYS 
!
      Use MOD_PREC, Only: SP
!
      Use MOD_LIMS, Only: MTLOC, KBM1, MLOC
      Use MOD_SIZES, Only: MGL
      Use MOD_CONTROL, Only: MSR !LB: needed only for debugging (writing text on screen)
!
      Implicit None
      Save
!
      Real (SP), Allocatable, Dimension (:, :) :: KW, K0, K1, K2, TB, &
     & KBb, TF, KF, TS, KS, IonS, FugFac, VPFac       !KSi,    &!
!f90pprREAL(SP), ALLOCATABLE, DIMENSION (:,:)  :: &
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
!f90ppr!TP,     &! Not considering protonation of phosphate or sulfate in pH calculation
!f90ppr!KP1,    &!
!f90ppr!KP2,    &!
!f90ppr!KP3,    &!
!f90ppr!TSi,    &!
!f90ppr!KSi,    &!
!f90pprFugFac, &
!f90pprVPFac
!
      Real (SP), Allocatable, Dimension (:) :: pCO2atm, CO2star_sat, &
     & CO2star_surf
!
      Real (SP), Allocatable, Dimension (:, :) :: DICUPT, DICBMP, &
     & DICPRD, DICMNL, DICDEN, DICGAS, DICSED, ALKNH4, ALKNO3, ALKNIT, &
     & ALKDEN, ALKREM, ALKNH4SED, ALKNO3SED
  !DICNIT
!
      Real (SP), Allocatable, Dimension (:, :) :: DICUPT_GL, DICBMP_GL, &
     & DICPRD_GL, DICMNL_GL, DICDEN_GL, DICGAS_GL, DICSED_GL, &
     & ALKNH4_GL, ALKNO3_GL, ALKNIT_GL, ALKDEN_GL, ALKREM_GL, &
     & ALKNH4SED_GL, ALKNO3SED_GL
  !DICNIT_GL
!
      Integer :: AIRSEA_OPTION !choice of parameterization to calculate O2 and CO2 air-sea exchange
!
      Real (SP) :: NXPCO2, pco2atmNX
      Character (Len=7) :: MAPCO2ATM
      Character (Len=4) :: UNITCO2ATM
!
	  Real (SP), Dimension (2) :: JDstart1 !AN
Contains
  !subroutine CO2SYSCONSTANTS
  !subroutine CO2SYSCONST_ALLOC
  !subroutine CO2SYSCONST_DEALLOC
!
  !************************************************************************
  !**         S U B R O U T I N E   CO2SYSCONST_ALLOC                    **
  !************************************************************************
      Subroutine CO2SYSCONST_ALLOC
!
         Allocate (KW(0:MTLOC, KBM1))
         KW = 0.0
         Allocate (K0(0:MTLOC, KBM1))
         K0 = 0.0
         Allocate (K1(0:MTLOC, KBM1))
         K1 = 0.0
         Allocate (K2(0:MTLOC, KBM1))
         K2 = 0.0
         Allocate (TB(0:MTLOC, KBM1))
         TB = 0.0
         Allocate (KBb(0:MTLOC, KBM1))
         KBb = 0.0
         Allocate (TF(0:MTLOC, KBM1))
         TF = 0.0
         Allocate (KF(0:MTLOC, KBM1))
         KF = 0.0
         Allocate (TS(0:MTLOC, KBM1))
         TS = 0.0
         Allocate (KS(0:MTLOC, KBM1))
         KS = 0.0
         Allocate (IonS(0:MTLOC, KBM1))
         IonS = 0.0
    !ALLOCATE(TP(0:MTLOC,KBM1));   TP  = 0.0
    !ALLOCATE(KP1(0:MTLOC,KBM1));  KP1  = 0.0
    !ALLOCATE(KP2(0:MTLOC,KBM1));  KP2  = 0.0
    !ALLOCATE(KP3(0:MTLOC,KBM1));  KP3  = 0.0
    !ALLOCATE(TSi(0:MTLOC,KBM1));  TSi  = 0.0
    !ALLOCATE(KSi(0:MTLOC,KBM1));  KSi  = 0.0
         Allocate (FugFac(0:MTLOC, KBM1))
         FugFac = 0.0
         Allocate (VPFac(0:MTLOC, KBM1))
         VPFac = 0.0
!
         Allocate (pCO2atm(0:MTLOC))
         pCO2atm = 0.0
         Allocate (CO2star_sat(0:MTLOC))
         CO2star_sat = 0.0
         Allocate (CO2star_surf(0:MTLOC))
         CO2star_surf = 0.0
!
    !allocation of variables to save DTTALK and DTTIC outputs
         Allocate (DICUPT(0:MTLOC, KBM1))
         DICUPT = 0.0 !DICUPT=-CP1-CP2  TDIC uptake
         Allocate (DICBMP(0:MTLOC, KBM1))
         DICBMP = 0.0 !DICBMP=DOR1+DOR2 TDIC source by basal met and photoresp
         Allocate (DICPRD(0:MTLOC, KBM1))
         DICPRD = 0.0 !DICPRD=DOP1+DOP2  "   "     by predation
         Allocate (DICMNL(0:MTLOC, KBM1))
         DICMNL = 0.0 !DICMNL=MNLLDOC+MNLRDOC
         Allocate (DICDEN(0:MTLOC, KBM1))
         DICDEN = 0.0 !DICDEN=DENIT
    !ALLOCATE(DICNIT(0:MTLOC,KBM1));   DICNIT  = 0.0  !DICNIT=-NT
         Allocate (DICGAS(0:MTLOC, KBM1))
         DICGAS = 0.0 !DICGAS=FLUXCO2  !only surface value is different to zero !Wen Long: clould have been just one layer
         Allocate (DICSED(0:MTLOC, KBM1))
         DICSED = 0.0 !DICSED=BENDIC   !only bottom value is different to zero  !Wen Long: could have just been one layer
!
         Allocate (ALKNH4(0:MTLOC, KBM1))
         ALKNH4 = 0.0 !
         Allocate (ALKNO3(0:MTLOC, KBM1))
         ALKNO3 = 0.0 !
         Allocate (ALKNIT(0:MTLOC, KBM1))
         ALKNIT = 0.0 !
         Allocate (ALKDEN(0:MTLOC, KBM1))
         ALKDEN = 0.0 !
         Allocate (ALKREM(0:MTLOC, KBM1))
         ALKREM = 0.0 !
         Allocate (ALKNH4SED(0:MTLOC, KBM1))
         ALKNH4SED = 0.0 !only bottom value is different to zero
         Allocate (ALKNO3SED(0:MTLOC, KBM1))
         ALKNO3SED = 0.0 !only bottom value is different to zero
!
         Allocate (DICUPT_GL(MGL, KBM1))
         DICUPT_GL = 0.0 !DICUPT=-CP1-CP2  TDIC uptake
         Allocate (DICBMP_GL(MGL, KBM1))
         DICBMP_GL = 0.0 !DICBMP=DOR1+DOR2 TDIC source by basal met and photoresp
         Allocate (DICPRD_GL(MGL, KBM1))
         DICPRD_GL = 0.0 !DICPRD=DOP1+DOP2  "   "     by predation
         Allocate (DICMNL_GL(MGL, KBM1))
         DICMNL_GL = 0.0 !DICMNL=MNLLDOC+MNLRDOC
         Allocate (DICDEN_GL(MGL, KBM1))
         DICDEN_GL = 0.0 !DICDEN=DENIT
    !ALLOCATE(DICNIT_GL(MGL,KBM1));   DICNIT_GL  = 0.0  !DICNIT=-NT
         Allocate (DICGAS_GL(MGL, KBM1))
         DICGAS_GL = 0.0 !DICGAS=FLUXCO2  !only surface value is different to zero
         Allocate (DICSED_GL(MGL, KBM1))
         DICSED_GL = 0.0 !DICSED=BENDIC   !only bottom value is different to zero
!
         Allocate (ALKNH4_GL(MGL, KBM1))
         ALKNH4_GL = 0.0 !
         Allocate (ALKNO3_GL(MGL, KBM1))
         ALKNO3_GL = 0.0 !
         Allocate (ALKNIT_GL(MGL, KBM1))
         ALKNIT_GL = 0.0 !
         Allocate (ALKDEN_GL(MGL, KBM1))
         ALKDEN_GL = 0.0 !
         Allocate (ALKREM_GL(MGL, KBM1))
         ALKREM_GL = 0.0 !
         Allocate (ALKNH4SED_GL(MGL, KBM1))
         ALKNH4SED_GL = 0.0 !only bottom value is different to zero
         Allocate (ALKNO3SED_GL(MGL, KBM1))
         ALKNO3SED_GL = 0.0 !only bottom value is different to zero
!
!
!
!
      End Subroutine CO2SYSCONST_ALLOC
!
!
  !************************************************************************
  !**         S U B R O U T I N E   CO2SYSCONST_DEALLOC                  **
  !***********************************************************************0
      Subroutine CO2SYSCONST_DEALLOC
!
         If (ALLOCATED(KW)) DEALLOCATE (KW)
         If (ALLOCATED(K0)) DEALLOCATE (K0)
         If (ALLOCATED(K1)) DEALLOCATE (K1)
         If (ALLOCATED(K2)) DEALLOCATE (K2)
         If (ALLOCATED(TB)) DEALLOCATE (TB)
         If (ALLOCATED(KBb)) DEALLOCATE (KBb)
         If (ALLOCATED(TF)) DEALLOCATE (TF)
         If (ALLOCATED(KF)) DEALLOCATE (KF)
         If (ALLOCATED(TS)) DEALLOCATE (TS)
         If (ALLOCATED(KS)) DEALLOCATE (KS)
         If (ALLOCATED(IonS)) DEALLOCATE (IonS)
    !IF(ALLOCATED(TP))DEALLOCATE(TP)
    !IF(ALLOCATED(KP1))DEALLOCATE(KP1)
    !IF(ALLOCATED(KP2))DEALLOCATE(KP2)
    !IF(ALLOCATED(KP3))DEALLOCATE(KP3)
    !IF(ALLOCATED(TSi))DEALLOCATE(TSi)
    !IF(ALLOCATED(KSi))DEALLOCATE(KSi)
         If (ALLOCATED(FugFac)) DEALLOCATE (FugFac)
         If (ALLOCATED(VPFac)) DEALLOCATE (VPFac)
!
         If (ALLOCATED(pCO2atm)) DEALLOCATE (pCO2atm)
         If (ALLOCATED(CO2star_sat)) DEALLOCATE (CO2star_sat)
         If (ALLOCATED(CO2star_surf)) DEALLOCATE (CO2star_surf)
!
         If (ALLOCATED(DICUPT)) DEALLOCATE (DICUPT)
         If (ALLOCATED(DICBMP)) DEALLOCATE (DICBMP)
         If (ALLOCATED(DICPRD)) DEALLOCATE (DICPRD)
         If (ALLOCATED(DICMNL)) DEALLOCATE (DICMNL)
         If (ALLOCATED(DICDEN)) DEALLOCATE (DICDEN)
    !IF(ALLOCATED(DICNIT))DEALLOCATE(DICNIT)
         If (ALLOCATED(DICGAS)) DEALLOCATE (DICGAS)
         If (ALLOCATED(DICSED)) DEALLOCATE (DICSED)
         If (ALLOCATED(ALKNH4)) DEALLOCATE (ALKNH4)
         If (ALLOCATED(ALKNO3)) DEALLOCATE (ALKNO3)
         If (ALLOCATED(ALKNIT)) DEALLOCATE (ALKNIT)
         If (ALLOCATED(ALKDEN)) DEALLOCATE (ALKDEN)
         If (ALLOCATED(ALKREM)) DEALLOCATE (ALKREM)
         If (ALLOCATED(ALKNH4SED)) DEALLOCATE (ALKNH4SED)
         If (ALLOCATED(ALKNO3SED)) DEALLOCATE (ALKNO3SED)
!
         If (ALLOCATED(DICUPT_GL)) DEALLOCATE (DICUPT_GL)
         If (ALLOCATED(DICBMP_GL)) DEALLOCATE (DICBMP_GL)
         If (ALLOCATED(DICPRD_GL)) DEALLOCATE (DICPRD_GL)
         If (ALLOCATED(DICMNL_GL)) DEALLOCATE (DICMNL_GL)
         If (ALLOCATED(DICDEN_GL)) DEALLOCATE (DICDEN_GL)
    !IF(ALLOCATED(DICNIT_GL))DEALLOCATE(DICNIT_GL)
         If (ALLOCATED(DICGAS_GL)) DEALLOCATE (DICGAS_GL)
         If (ALLOCATED(DICSED_GL)) DEALLOCATE (DICSED_GL)
         If (ALLOCATED(ALKNH4_GL)) DEALLOCATE (ALKNH4_GL)
         If (ALLOCATED(ALKNO3_GL)) DEALLOCATE (ALKNO3_GL)
         If (ALLOCATED(ALKNIT_GL)) DEALLOCATE (ALKNIT_GL)
         If (ALLOCATED(ALKDEN_GL)) DEALLOCATE (ALKDEN_GL)
         If (ALLOCATED(ALKREM_GL)) DEALLOCATE (ALKREM_GL)
         If (ALLOCATED(ALKNH4SED_GL)) DEALLOCATE (ALKNH4SED_GL)
         If (ALLOCATED(ALKNO3SED_GL)) DEALLOCATE (ALKNO3SED_GL)
!
      End Subroutine CO2SYSCONST_DEALLOC
!
!
  !************************************************************************
  !**         S U B R O U T I N E   C O 2 S Y S C O N S T A N T S        **
  !******************************************************************************!
  !   Constants for CO2SYS calculations   (CO2SYS.m v1.1 sept 2011)              !
  !                                                                               !
  !   In this subroutine we calculate all the parameters from CO2SYS.m           !
  !   that depend on temperature, salinity and user inputs                         !
  !                                                                               !
  !   Laura Bianucci, Feb 2015 - laura.bianucci@pnnl.gov                         !
  !                                                                               !
  !                                                                               !
  !    Below, follows the 'help' from the matlab subroutine Constants:               !
  !                                                                               !
  !% SUB Constants, version 04.01, 10-13-97, written by Ernie Lewis.               !
  !% Inputs: pHScale%, WhichKs%, WhoseKSO4%, Sali, TempCi, Pdbar                   !
  !% Outputs: K0, K(), T(), fH, FugFac, VPFac                                       !
  !% This finds the Constants of the CO2 system in seawater or freshwater,       !
  !% corrects them for pressure, and reports them on the chosen pH scale.        !
  !% The process is as follows: the Constants (except KS, KF which stay on the   !
  !% free scale - these are only corrected for pressure) are                       !
  !%       1) evaluated as they are given in the literature                       !
  !%       2) converted to the SWS scale in mol/kg-SW or to the NBS scale        !
  !%       3) corrected for pressure                                                   !
  !%       4) converted to the SWS pH scale in mol/kg-SW                              !
  !%       5) converted to the chosen pH scale                                      !
  !%                                                                                  !
  !%       PROGRAMMER'S NOTE: all logs are log base e                                  !
  !%       PROGRAMMER'S NOTE: all Constants are converted to the pH scale           !
  !%               pHScale% (the chosen one) in units of mol/kg-SW               !
  !%               except KS and KF are on the free scale                           !
  !%               and KW is in units of (mol/kg-SW)^2                           !
  !******************************************************************************!
      Subroutine CO2SYSCONSTANTS
!
         Use MOD_HYDROVARS, Only: ZZ, D,ZZ2D
		 USE MOD_FILEINFO, ONLY: INK1K2
         Use MOD_WQM, Only: T, SALT, TALK, PO4 ! (Peng et al, PengCorrection of TALK)
         Implicit None
         Save
!
    !! need to update this one to define ph,pco2,tdic,talk and probably all variables !WQM defines T, SALT, D
!
!
         Real (SP) :: TempK, RT, logTempK, Pbar, sqrSal, S, SWStoTOT, &
        & FREEtoTOT, fH
         Real (SP) :: aux, aux2, A_1, B_1, C_1, deltaV, Kappa, lnK1fac, &
        & lnK2fac, lnKBfac, lnKWfac, lnKFfac, lnKSfac, lnKP1fac, &
        & lnKP2fac, lnKP3fac, lnKSifac
         Real (SP) :: Delta, pK2
!
         Real (SP) :: RGasConstant = 83.1451_SP ! ml bar-1 K-1 mol-1, DOEv2
    !     REAL(SP) :: RGasConstant = 83.14472 ! ml bar-1 K-1 mol-1, DOEv3 !(Dickson et al., 2007). !Commented in CO2SYS.m v1.1
         Real (SP) :: P1atm = 1.01325_SP ! in bar
!
         ! even though save is set above, reset it here, since these really need saved and someone might remove above save because save is "bad"
         Logical,save :: NEEDREAD = .TRUE.
         Integer,save :: WHICH_K1K2, WHOSEKSO4
         Integer :: I, K
!
!
    ! Read WHICH_K1K2 and WHOSEKSO4 from input file K1_K2_KSO4.dat
         IF (NEEDREAD) THEN
         Open (Unit=INK1K2, File='inputs/input_K1_K2_KSO4.dat', Status='ol&
        &d')
!
         Read (INK1K2, 1099) WHICH_K1K2, WHOSEKSO4
!
         Close (INK1K2)
         NEEDREAD = .FALSE.
         ENDIF ! NEEDREAD
!
1099     Format (/ (2 I2))
!
    !if(MSR)write(*,*)'WHICH_K1K2, WHOSEKSO4='
    !if(MSR)write(*,*)WHICH_K1K2, WHOSEKSO4
!
    ! Calculate constants inside a double DO loop
         Do K = 1, KBM1
            Do I = 1, MLOC
               TempK = T (I, K) + 273.15_SP
               RT = RGasConstant * TempK
               logTempK = Log (TempK)!LOG to provide result in double precision
               Pbar = - D (I) * ZZ2D (I,K) / 10.0_SP !depth of this layer (positive,from surface) divided 10 to convert from m (dbar) to bar.
               S = SALT (I, K)
               sqrSal = Sqrt (S)
!
          !if(MSR)write(*,*)'D(I),ZZ2D(I,K),Pbar='
          !if(MSR)write(*,*)D(I),ZZ2D(I,K),Pbar
!
          !------------------------------------------------------------------------------
          !Calculate TB - Total Borate
               If (WHICH_K1K2 .Eq. 8) Then
                  TB (I, K) = 0.0_SP !Pure water (Millero, 1979)
                  S = 0.0_SP
                  sqrSal = 0.0_SP
!
               Else If (WHICH_K1K2 .Eq. 6 .Or. WHICH_K1K2 .Eq. 7) Then
                  TB (I, K) = 0.0004106_SP * S / 35.0_SP ! in mol/kg-SW
             !    % this is .00001173 *Sali
             !    % this is about 1% lower than Uppstrom's value
             !    % Culkin, F., in Chemical Oceanography,
             !    % ed. Riley and Skirrow, 1965:
             !    % GEOSECS references this, but this value is not explicitly
             !    % given here
!
               Else
                  If (WHOSEKSO4 .Eq. 1 .Or. WHOSEKSO4 .Eq. 2) Then
                ! If user opted for Uppstrom's values:
                !        % Uppstrom, L., Deep-Sea Research 21:161-162, 1974:
                !        % this is .000416 *Sali /35. = .0000119 *Sali
                !        % TB = (0.000232 /10.811) *(S /1.80655); % in mol/kg-SW
                !!LB not sure why above is commented and version below used instead,
                !!although the difference is very small (1e-9)
                     TB (I, K) = 0.0004157_SP * S / 35.0_SP ! in mol/kg-SW
!
                  Else
                ! If user opted for the new Lee values:
                !        % Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.
                !         % Geochimica Et Cosmochimica Acta 74 (6): 1801â??1811.
                     TB (I, K) = 0.0004326_SP * S / 35.0_SP ! in mol/kg-SW
!
                  End If
               End If
!
          !------------------------------------------------------------------------------
          ! Calculate TF (fluoride), TS (sulfate), IonS (SO4), KO (solubility of CO2)
          ! TF:
          !% Riley, J. P., Deep-Sea Research 12:219-220, 1965:
          !% this is .000068 *Sali /35. = .00000195 *Sali
          ! TS:
          !% Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
          !% this is .02824 *Sali /35. = .0008067 *Sali
          ! K0:
          !% Weiss, R. F., Marine Chemistry 2:203-215, 1974.
          ! IonS:
          !% This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4
!
               TF (I, K) = (0.000067_SP/18.998_SP) * (S/1.80655_SP)! in mol/kg-SW
               TS (I, K) = (0.14_SP/96.062_SP) * (S/1.80655_SP)! in mol/kg-SW
               IonS (I, K) = 19.924_SP * S / (1000.0_SP-1.005_SP*S)
!
               aux = TempK / 100.0_SP
               aux2 = - 60.2409_SP + 93.4517_SP / aux + 23.3585_SP * &
              & Log (aux) + S * &
              & (0.023517_SP-0.023656_SP*aux+0.0047036_SP*aux**2)
               K0 (I, K) = Exp (aux2)! this is in mol/kg-SW/atm
!
!
          !------------------------------------------------------------------------------
          !Calculate KS -
               If (WHOSEKSO4 .Eq. 1 .Or. WHOSEKSO4 .Eq. 3) Then
             !    % Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
             !    % The goodness of fit is .021.
             !    % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
             !    % TYPO on p. 121: the constant e9 should be e8.
             !    % This is from eqs 22 and 23 on p. 123, and Table 4 on p 121:
!
                  aux = - 4276.1_SP / TempK + 141.328_SP - 23.093_SP * &
                 & logTempK + &
                 & (-13856.0_SP/TempK+324.57_SP-47.986_SP*logTempK) * &
                 & Sqrt (IonS(I, K)) + IonS (I, K) * &
                 & (35474.0_SP/TempK-771.54_SP+114.723_SP*logTempK) + &
                 & (-2698.0_SP/TempK) * Sqrt (IonS(I, K)) * IonS (I, K) &
                 & + (1776.0_SP/TempK) * IonS (I, K) ** 2
                  KS (I, K) = Exp (aux) * (1.0_SP-0.001005_SP*S)!convert to mol/kg-SW
!
               Else
             !    % Khoo, et al, Analytical Chemistry, 49(1):29-34, 1977
             !    % KS was found by titrations with a hydrogen electrode
             !    % of artificial seawater containing sulfate (but without F)
             !    % at 3 salinities from 20 to 45 and artificial seawater NOT
             !    % containing sulfate (nor F) at 16 salinities from 15 to 45,
             !    % both at temperatures from 5 to 40 deg C.
             !    % KS is on the Free pH scale (inherently so).
             !    % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
             !    % He finds log(beta) which = my pKS;
             !    % his beta is an association constant.
             !    % The rms error is .0021 in pKS, or about .5% in KS.
             !    % This is equation 20 on p. 33:
!
                  aux = 647.59_SP / TempK - 6.3451_SP + 0.019085_SP * &
                 & TempK - 0.5208_SP * Sqrt (IonS(I, K))
                  KS (I, K) = 10.0_SP ** (-aux) * &
                 & (1.0_SP-0.001005_SP*S)! % convert to mol/kg-SW
!
               End If
!
!
          !------------------------------------------------------------------------------
          ! Calculate KF:
          !% Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979
               aux = 1590.2_SP / TempK - 12.641_SP + 1.525_SP * Sqrt &
              & (IonS(I, K))
               KF (I, K) = Exp (aux) * (1.0_SP-0.001005_SP*S)! % convert to mol/kg-SW
!
          !% Another expression exists for KF: Perez and Fraga 1987. Not used here since ill defined for low salinity. (to be used f
          !% Nonetheless, P&F87 might actually be better than the fit of D&R79 above, which is based on only three salinities: [0 26
          !% lnKF = 874 /TempK - 9.68 + 0.111 *S**0.5;
          !% KF   = exp(lnKF);                   % this is on the free pH scale in mol/kg-SW
!
          !------------------------------------------------------------------------------
          ! Calculate pH Scale Conversion Factors:
          !        These are NOT pressure-corrected
               SWStoTOT = (1.0_SP+TS(I, K)/KS(I, K)) / (1.0_SP+TS(I, &
              & K)/KS(I, K)+TF(I, K)/KF(I, K))
!
               FREEtoTOT = 1.0_SP + TS (I, K) / KS (I, K)
!
!
          !------------------------------------------------------------------------------
          ! CalculatefH
          !% Use GEOSECS's value for cases 1,2,3,4,5 (and 6) to convert pH scales.
!
               If (WHICH_K1K2 .Eq. 8) Then
                  fH = 1.0_SP ! this shouldn't occur in the program for this case
!
               Else If (WHICH_K1K2 .Eq. 7) Then
             !    % Peng et al, Tellus 39B:439-458, 1987:
             !    % They reference the GEOSECS report, but round the value
             !    % given there off so that it is about .008 (1%) lower. It
             !    % doesn't agree with the check value they give on p. 456.
                  fH = 1.29_SP - 0.00204_SP * TempK + &
                 & (0.00046_SP-0.00000148_SP*TempK) * S ** 2
!
             !!add a note for the case when WHICH_K1K2=7: (we can't calculate PengCorrection)
                  Write (*,*) '!*****.......*****.......*****.......***&
                 &**.......*****.......'
                  Write (*,*) '          K1K2=7 --> TALK corrected with&
                 & PO4'
                  Write (*,*) '!*****.......*****.......*****.......***&
                 &**.......*****.......'
                  TALK (I, K) = TALK (I, K) - PO4 (I, K)
                  Write (*,*) 'LBnote: done TALK corrrection'
               Else
             !    % Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition,
             !    % v. 3, 1982 (p. 80);
                  fH = 1.2948_SP - 0.002036_SP * TempK + &
                 & (0.0004607_SP-0.000001475_SP*TempK) * S ** 2
!
               End If
!
!
          !------------------------------------------------------------------------------
          ! Calculate KB:
               If (WHICH_K1K2 .Eq. 8) Then
                  KBb (I, K) = 0.0_SP !Pure water (Millero, 1979)
!
               Else If (WHICH_K1K2 .Eq. 6 .Or. WHICH_K1K2 .Eq. 7) Then
             !    % This is for GEOSECS and Peng et al.
             !    % Lyman, John, UCLA Thesis, 1957
             !    % fit by Li et al, JGR 74:5507-5525, 1969:
!
                  aux = - 9.26_SP + 0.00886_SP * S + 0.01_SP * T (I, K)!Temp in Celcius here
                  KBb (I, K) = 10.0_SP ** (aux) / fH ! convert to the SWS scale
!
               Else
             !    % Dickson, A. G., Deep-Sea Research 37:755-766, 1990:
!
                  aux = - 8966.9_SP - 2890.53_SP * sqrSal - 77.942_SP * &
                 & S + 1.728_SP * sqrSal * S - 0.0996_SP * S ** 2
                  aux2 = aux / TempK + 148.0248_SP + 137.1942_SP * &
                 & sqrSal + 1.62142_SP * S + logTempK * &
                 & (-24.4344_SP-25.085_SP*sqrSal-0.2474_SP*S) + &
                 & 0.053105_SP * sqrSal * TempK
                  KBb (I, K) = Exp (aux2) / SWStoTOT ! convert to SWS pH scale
!
               End If
!
          !------------------------------------------------------------------------------
          ! Calculate KW:
               If (WHICH_K1K2 .Eq. 8) Then
             !    % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
             !    % refit data of Harned and Owen, The Physical Chemistry of
             !    % Electrolyte Solutions, 1958
                  aux = 148.9802_SP - 13847.26_SP / TempK - 23.6521_SP &
                 & * logTempK
!
!
               Else If (WHICH_K1K2 .Eq. 7) Then
             !    % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
                  aux = 148.9802_SP - 13847.26_SP / TempK - 23.6521_SP &
                 & * logTempK + sqrSal * &
                 & (-79.2447_SP+3298.72_SP/TempK+12.0408_SP*logTempK) - &
                 & 0.019813_SP * S
               Else
             !    % Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995.
             !    % his check value of 1.6 umol/kg-SW should be 6.2
                  aux = 148.9802_SP - 13847.26_SP / TempK - 23.6521_SP &
                 & * logTempK + sqrSal * &
                 & (-5.977_SP+118.67_SP/TempK+1.0495_SP*logTempK) - &
                 & 0.01615_SP * S
               End If
!
               KW (I, K) = Exp (aux)! % this is on the SWS pH scale in (mol/kg-SW)^2
!
               If (WHICH_K1K2 .Eq. 6) Then
                  KW (I, K) = 0.0_SP ! GEOSECS doesn't include OH effects
               End If
!
!
          !------------------------------------------------------------------------------
          !!**LB note: we currently ignore phosphate and silicate protonation           **
          !!**     (we do not model PO4 nor Si and consider effects on TALK negligible) **
          !! Calculate KP1, KP2, KP3, KSi:
          !
          !        IF(WHICH_K1K2.eq.7) THEN
          !                KP1 = 0.02_SP
          !!    % Peng et al don't include the contribution from this term,
          !!    % but it is so small it doesn't contribute. It needs to be
          !!    % kept so that the routines work ok.
          !!    % KP2, KP3 from Kester, D. R., and Pytkowicz, R. M.,
          !!    % Limnology and Oceanography 12:243-252, 1967.
          !!    % (these are only for sals 33 to 36 and are on the NBS scale)
          !!    % KSi from Sillen, Martell, and Bjerrum,  Stability Constants of metal-ion complexes,
          !!    % The Chemical Society (London), Special Publ. 17:751, 1964
          !
          !                KP2(I,K) = exp(-9.039_SP - 1450.0_SP /TempK)        & ! this is on the NBS scale
          !                            /fH                                ! convert to SWS scale
          !                KP3(I,K) = exp(4.466_SP- 7276.0_SP /TempK)         & ! this is on the NBS scale
          !                            /fH                                ! convert to SWS scale
          !                KSi(I,K) = 0.0000000004_SP                      & ! this is on the NBS scale
          !                            /fH                                ! convert to SWS scale
          !
          !        ELSE IF (WHICH_K1K2.eq.6.or.WHICH_K1K2.eq.8) THEN
          !!    % Neither the GEOSECS choice nor the freshwater choice
          !!    % include contributions from phosphate or silicate.
          !                 KP1(I,K) = 0.0_SP
          !                 KP2(I,K) = 0.0_SP
          !                 KP3(I,K) = 0.0_SP
          !                 KSi(I,K) = 0.0_SP
          !
          !        ELSE
          !!    % Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
          !!    % KP1, KP2, KP3 are on the SWS pH scale in mol/kg-SW.
          !!    % KSi was given on the SWS pH scale in molal units.
          !
          !              aux = -4576.752_SP /TempK + 115.54_SP - 18.453_SP * logTempK    &
          !                    + (-106.736_SP /TempK + 0.69171_SP) *sqrSal            &
          !                    + (-0.65643_SP/TempK - 0.01844_SP) *S
          !              KP1(I,K) = exp(aux)
          !
          !              aux = -8814.715_SP /TempK + 172.1033_SP- 27.927_SP*logTempK     &
          !                   + (-160.34_SP /TempK + 1.3566_SP) * sqrSal              &
          !                   + (0.37335_SP /TempK - 0.05778_SP) * S
          !              KP2(I,K) = exp(aux)
          !
          !              aux = -3070.75_SP /TempK - 18.126_SP +                       &
          !                    + (17.27039_SP /TempK + 2.81197_SP) *sqrSal            &
          !                    + (-44.99486_SP /TempK - 0.09984) *S
          !              KP3(I,K) = exp(aux)
          !
          !              aux = -8904.2_SP /TempK + 117.4_SP - 19.334_SP *logTempK        &
          !                    + (-458.79_SP /TempK + 3.5913_SP) *sqrt(IonS(I,K))     &
          !                    + (188.74_SP /TempK - 1.5998_SP) *IonS(I,K)            &
          !                    +(-12.1652_SP /TempK + 0.07871_SP) *IonS(I,K)**2
          !              KSi(I,K) = exp(aux)                        & ! this is on the SWS pH scale in mol/kg-H2O
          !                         *(1.0_SP - 0.001005_SP *S)               ! convert to mol/kg-SW
          !
          !        END IF
!
          !------------------------------------------------------------------------------
          ! Calculate K1, K2:
               If (WHICH_K1K2 .Eq. 1) Then
             !    % ROY et al, Marine Chemistry, 44:249-267, 1993
             !    % (see also: Erratum, Marine Chemistry 45:337, 1994
             !    % and Erratum, Marine Chemistry 52:183, 1996)
             !    % Typo: in the abstract on p. 249: in the eq. for lnK1* the
             !    % last term should have S raised to the power 1.5.
             !    % They claim standard deviations (p. 254) of the fits as
             !    % .0048 for lnK1 (.5% in K1) and .007 in lnK2 (.7% in K2).
             !    % They also claim (p. 258) 2s precisions of .004 in pK1 and
             !    % .006 in pK2. These are consistent, but Andrew Dickson
             !    % (personal communication) obtained an rms deviation of about
             !    % .004 in pK1 and .003 in pK2. This would be a 2s precision
             !    % of about 2% in K1 and 1.5% in K2.
             !    % T:  0-45  S:  5-45. Total Scale. Artificial sewater.
             !    % This is eq. 29 on p. 254 and what they use in their abstract:
!
                  aux = 2.83655_SP - 2307.1266_SP / TempK - &
                 & 1.5529413_SP * logTempK + &
                 & (-0.20760841_SP-4.0484_SP/TempK) * sqrSal + &
                 & 0.08468345_SP * S - 0.00654208_SP * sqrSal * S
                  K1 (I, K) = Exp (aux) * (1.0_SP-0.001005_SP*S) / &
                 & SWStoTOT ! convert to SWS pH scale
!
             !    % This is eq. 30 on p. 254 and what they use in their abstract:
                  aux = - 9.226508_SP - 3351.6106_SP / TempK - &
                 & 0.2005743_SP * logTempK + &
                 & (-0.106901773_SP-23.9722_SP/TempK) * sqrSal + &
                 & 0.1130822_SP * S - 0.00846934_SP * sqrSal * S
                  K2 (I, K) = Exp (aux) * (1.0_SP-0.001005_SP*S) / &
                 & SWStoTOT ! convert to SWS pH scale
!
               Else If (WHICH_K1K2 .Eq. 2) Then
             !    % GOYET AND POISSON, Deep-Sea Research, 36(11):1635-1654, 1989
             !    % The 2s precision in pK1 is .011, or 2.5% in K1.
             !    % The 2s precision in pK2 is .02, or 4.5% in K2.
             !    % This is in Table 5 on p. 1652 and what they use in the abstract:
!
                  aux = 812.27_SP / TempK + 3.356_SP - 0.00171_SP * S * &
                 & logTempK + 0.000091_SP * S ** 2
                  K1 (I, K) = 10.0_SP ** (-aux)! this is on the SWS pH scale in mol/kg-SW
!
             !    % This is in Table 5 on p. 1652 and what they use in the abstract:
                  aux = 1450.87_SP / TempK + 4.604_SP - 0.00385_SP * S &
                 & * logTempK + 0.000182_SP * S ** 2
                  K2 (I, K) = 10.0_SP ** (-aux)! this is on the SWS pH scale in mol/kg-SW
!
!
               Else If (WHICH_K1K2 .Eq. 3) Then
             !    % HANSSON refit BY DICKSON AND MILLERO
             !    % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
             !    % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
             !    % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973
             !    % and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973.
             !    % on the SWS pH scale in mol/kg-SW.
             !    % Hansson gave his results on the Total scale (he called it
             !    % the seawater scale) and in mol/kg-SW.
             !    % Typo in DM on p. 1739 in Table 4: the equation for pK2*
             !    % for Hansson should have a .000132 *S^2
             !    % instead of a .000116 *S^2.
             !    % The 2s precision in pK1 is .013, or 3% in K1.
             !    % The 2s precision in pK2 is .017, or 4.1% in K2.
             !    % This is from Table 4 on p. 1739.
!
                  aux = 851.4_SP / TempK + 3.237_SP - 0.0106_SP * S + &
                 & 0.000105_SP * S ** 2
                  K1 (I, K) = 10.0_SP ** (-aux)! this is on the SWS pH scale in mol/kg-SW
!
             !    % This is from Table 4 on p. 1739.
                  aux = - 3885.4_SP / TempK + 125.844_SP - 18.141_SP * &
                 & logTempK - 0.0192_SP * S + 0.000132_SP * S ** 2
                  K2 (I, K) = 10.0_SP ** (-aux)! this is on the SWS pH scale in mol/kg-SW
!
!
               Else If (WHICH_K1K2 .Eq. 4) Then
             !    % MEHRBACH refit BY DICKSON AND MILLERO
             !    % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
             !    % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
             !    % refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
             !    % on the SWS pH scale in mol/kg-SW.
             !    % Mehrbach et al gave results on the NBS scale.
             !    % The 2s precision in pK1 is .011, or 2.6% in K1.
             !    % The 2s precision in pK2 is .020, or 4.6% in K2.
             !     % Valid for salinity 20-40.
             !    % This is in Table 4 on p. 1739.
!
                  aux = 3670.7_SP / TempK - 62.008_SP + 9.7944_SP * &
                 & logTempK - 0.0118_SP * S + 0.000116_SP * S ** 2
                  K1 (I, K) = 10.0_SP ** (-aux)! this is on the SWS pH scale in mol/kg-SW
!
             !    % This is in Table 4 on p. 1739.
                  aux = 1394.7_SP / TempK + 4.777_SP - 0.0184_SP * S + &
                 & 0.000118_SP * S ** 2
                  K2 (I, K) = 10.0_SP ** (-aux)! his is on the SWS pH scale in mol/kg-SW
!
!
               Else If (WHICH_K1K2 .Eq. 5) Then
             !    % HANSSON and MEHRBACH refit BY DICKSON AND MILLERO
             !    % Dickson and Millero, Deep-Sea Research,34(10):1733-1743, 1987
             !    % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
             !    % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973,
             !    % Hansson, Acta Chemica Scandanavia, 27:931-944, 1973,
             !    % and Mehrbach et al, Limnol. Oceanogr.,18(6):897-907, 1973
             !    % on the SWS pH scale in mol/kg-SW.
             !    % Typo in DM on p. 1740 in Table 5: the second equation
             !    % should be pK2* =, not pK1* =.
             !    % The 2s precision in pK1 is .017, or 4% in K1.
             !    % The 2s precision in pK2 is .026, or 6% in K2.
             !    % Valid for salinity 20-40.
             !    % This is in Table 5 on p. 1740.
                  aux = 845.0_SP / TempK + 3.248_SP - 0.0098_SP * S + &
                 & 0.000087_SP * S ** 2
!
                  K1 (I, K) = 10.0_SP ** (-aux)! this is on the SWS pH scale in mol/kg-SW
!
             !    % This is in Table 5 on p. 1740.
                  aux = 1377.3_SP / TempK + 4.824_SP - 0.0185_SP * S + &
                 & 0.000122_SP * S ** 2
                  K2 (I, K) = 10.0_SP ** (-aux)! this is on the SWS pH scale in mol/kg-SW
!
!
               Else If (WHICH_K1K2 .Eq. 6 .Or. WHICH_K1K2 .Eq. 7) Then
             !    % GEOSECS and Peng et al use K1, K2 from Mehrbach et al,
             !    % Limnology and Oceanography, 18(6):897-907, 1973.
             !    % I.e., these are the original Mehrbach dissociation constants.
             !    % The 2s precision in pK1 is .005, or 1.2% in K1.
             !    % The 2s precision in pK2 is .008, or 2% in K2.
                  aux = - 13.7201_SP + 0.031334_SP * TempK + 3235.76_SP &
                 & / TempK + 1.3e-5_SP * S * TempK - 0.1032_SP * sqrSal
                  K1 (I, K) = 10.0_SP ** (-aux) / fH ! convert to SWS scale
!
                  aux = 5371.9645_SP + 1.671221_SP * TempK + 0.22913_SP &
                 & * S + 18.3802_SP * Log10 (S) - 128375.28_SP / TempK &
                 & - 2194.3055_SP * Log10 (TempK) - 8.0944e-4_SP * S * &
                 & TempK - 5617.11_SP * Log10 (S) / TempK + 2.136_SP * &
                 & S / TempK ! pK2 is not defined for S=0, since log10(0)=-inf
                  K2 (I, K) = 10.0_SP ** (-aux) / fH ! convert to SWS scale
!
               Else If (WHICH_K1K2 .Eq. 8) Then
             !    % PURE WATER CASE
             !    % Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
             !    % K1 from refit data from Harned and Davis,
             !    % J American Chemical Society, 65:2030-2037, 1943.
             !    % K2 from refit data from Harned and Scholes,
             !    % J American Chemical Society, 43:1706-1709, 1941.
             !    % This is only to be used for S=0 water (note the absence of S in the below formulations)
             !    % These are the thermodynamic Constants:
                  aux = 290.9097_SP - 14554.21_SP / TempK - 45.0575_SP &
                 & * logTempK
                  K1 (I, K) = Exp (aux)
!
                  aux = 207.6548_SP - 11843.79_SP / TempK - 33.6485_SP &
                 & * logTempK
                  K2 (I, K) = Exp (aux)
!
               Else If (WHICH_K1K2 .Eq. 9) Then
             !   % From Cai and Wang 1998, for estuarine use.
             !    % Data used in this work is from:
             !    % K1: Merhback (1973) for S>15, for S<15: Mook and Keone (1975)
             !    % K2: Merhback (1973) for S>20, for S<20: Edmond and Gieskes (1970)
             !    % Sigma of residuals between fits and above data: Â±0.015, +0.040 for K1 and K2, respectively.
             !    % S 0-40, Temp 0.2-30
             !   % Limnol. Oceanogr. 43(4) (1998) 657-668
             !    % On the NBS scale
             !    % Their check values for F1 don't work out, not sure if this was correctly published...
                  aux = 200.1_SP / TempK + 0.3220_SP
                  aux2 = 3404.71_SP / TempK + 0.032786_SP * TempK - &
                 & 14.8435_SP - 0.071692_SP * aux * sqrSal + &
                 & 0.0021487_SP * S
                  K1 (I, K) = 10.0_SP ** (-aux2) / fH ! convert to SWS scale (uncertain at low S due to junction potential);
!
                  aux = - 129.24_SP / TempK + 1.4381_SP
                  aux2 = 2902.39_SP / TempK + 0.02379_SP * TempK - &
                 & 6.4980_SP - 0.3191_SP * aux * sqrSal + 0.0198_SP * S
                  K2 (I, K) = 10.0_SP ** (-aux2) / fH ! convert to SWS scale (uncertain at low S due to junction potential);
!
               Else If (WHICH_K1K2 .Eq. 10) Then
             !   % From Lueker, Dickson, Keeling, 2000
             !    % This is Mehrbach's data refit after conversion to the total scale, for comparison with their equilibrator work.
             !   % Mar. Chem. 70 (2000) 105-119
             !   % Total scale and kg-sw
!
!
                  aux = 3633.86_SP / TempK - 61.2172_SP + 9.6777_SP * &
                 & logTempK - 0.011555_SP * S + 0.0001152_SP * S ** 2
                  K1 (I, K) = 10.0_SP ** (-aux) / SWStoTOT ! convert to SWS pH scale
!
                  aux = 471.78_SP / TempK + 25.929_SP - 3.16967_SP * &
                 & logTempK - 0.01781_SP * S + 0.0001122_SP * S ** 2
                  K2 (I, K) = 10.0_SP ** (-aux) / SWStoTOT ! convert to SWS pH scale
!
!
               Else If (WHICH_K1K2 .Eq. 11) Then
             !    % Mojica Prieto and Millero 2002. Geochim. et Cosmochim. Acta. 66(14) 2529-2540.
             !    % sigma for pK1 is reported to be 0.0056
             !    % sigma for pK2 is reported to be 0.010
             !    % This is from the abstract and pages 2536-2537
                  aux = - 43.6977_SP - 0.0129037_SP * S + 1.364e-4_SP * &
                 & S ** 2 + 2885.378_SP / TempK + 7.045159_SP * Log &
                 & (TempK)
                  K1 (I, K) = 10.0_SP ** (-aux)! this is on the SWS pH scale in mol/kg-SW
!
                  aux = - 452.0940_SP + 13.142162_SP * S - 8.101e-4_SP &
                 & * S ** 2 + 21263.61_SP / TempK + 68.483143_SP * Log &
                 & (TempK) + (-581.4428_SP*S+0.259601_SP*S**2) / TempK &
                 & - 1.967035_SP * S * Log (TempK)
                  K2 (I, K) = 10.0_SP ** (-aux)! this is on the SWS pH scale in mol/kg-SW
!
!
               Else If (WHICH_K1K2 .Eq. 12) Then
             !    % Millero et al., 2002. Deep-Sea Res. I (49) 1705-1723.
             !    % Calculated from overdetermined WOCE-era field measurements
             !    % sigma for pK1 is reported to be 0.005
             !    % sigma for pK2 is reported to be 0.008
             !    % This is from page 1715
                  aux = 6.359_SP - 0.00664_SP * S - 0.01322_SP * T (I, &
                 & K) + 4.989e-5_SP * T (I, K) ** 2 !T here is in Celcius
                  K1 (I, K) = 10.0_SP ** (-aux)! this is on the SWS pH scale in mol/kg-SW
!
                  aux = 9.867_SP - 0.01314_SP * S - 0.01904_SP * T (I, &
                 & K) + 2.448e-5_SP * T (I, K) ** 2 !T here is in Celcius
                  K2 (I, K) = 10.0_SP ** (-aux)! this is on the SWS pH scale in mol/kg-SW
!
!
               Else If (WHICH_K1K2 .Eq. 13) Then
             !   % From Millero 2006 work on pK1 and pK2 from titrations
             !    % Millero, Graham, Huang, Bustos-Serrano, Pierrot. Mar.Chem. 100 (2006) 80-94.
             !   % S=1 to 50, T=0 to 50. On seawater scale (SWS). From titrations in Gulf Stream seawater.
                  aux = - 126.34048_SP + 6320.813_SP / TempK + &
                 & 19.568224_SP * logTempK
                  A_1 = 13.4191_SP * sqrSal + 0.0331_SP * S - &
                 & 5.33e-5_SP * S ** 2
                  B_1 = - 530.123_SP * sqrSal - 6.103_SP * S
                  C_1 = - 2.06950_SP * sqrSal
                  aux2 = A_1 + B_1 / TempK + C_1 * logTempK + aux ! aux2=pK1; pK1 sigma = 0.0054
                  K1 (I, K) = 10.0_SP ** (-aux2)
!
                  aux = - 90.18333_SP + 5143.692_SP / TempK + &
                 & 14.613358_SP * logTempK
                  A_1 = 21.0894_SP * sqrSal + 0.1248_SP * S - &
                 & 3.687e-4_SP * S ** 2
                  B_1 = - 772.483_SP * sqrSal - 20.051_SP * S
                  C_1 = - 3.3336_SP * sqrSal
                  aux2 = A_1 + B_1 / TempK + C_1 * logTempK + aux !aux2=pK2;  pK2 sigma = 0.011
                  K2 (I, K) = 10.0_SP ** (-aux2)
!
               Else If (WHICH_K1K2 .Eq. 14) Then
             !    % From Millero, 2010, also for estuarine use.
             !    % Marine and Freshwater Research, v. 61, p. 139â??142.
             !    % Fits through compilation of real seawater titration results:
             !    % Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), Millero et al. (2006)
             !    % Constants for K's on the SWS;
             !    % This is from page 141 (pK1_0)
                  aux = - 126.34048_SP + 6320.813_SP / TempK + &
                 & 19.568224_SP * logTempK
             !    % This is from their table 2, page 140.
                  A_1 = 13.4038_SP * sqrSal + 0.03206_SP * S - &
                 & 5.242e-5_SP * S ** 2
                  B_1 = - 530.659_SP * sqrSal - 5.8210_SP * S
                  C_1 = - 2.0664_SP * sqrSal
                  aux2 = aux + A_1 + B_1 / TempK + C_1 * logTempK
                  K1 (I, K) = 10.0_SP ** (-aux2)
!
             !    % This is from page 141 (pK2_0)
                  aux = - 90.18333_SP + 5143.692_SP / TempK + &
                 & 14.613358_SP * logTempK
             !    % This is from their table 3, page 140.
                  A_1 = 21.3728_SP * sqrSal + 0.1218_SP * S - &
                 & 3.688e-4_SP * S ** 2
                  B_1 = - 788.289_SP * sqrSal - 19.189_SP * S
             !            C_1 = -3.374_SP   * sqrSal
                  C_1 = - 3.3738_SP * sqrSal !!LB: using coef from Millero's spreadsheet rather than the publication Millero et al (2010)
             !!(as suggested by Orr et al, Biogeosciences Disc, 2014 - footnote table 10)
                  pK2 = aux + A_1 + B_1 / TempK + C_1 * logTempK
                  K2 (I, K) = 10.0_SP ** (-pK2)
!
               End If
!
!
!
          !------------------------------------------------------------------------------
          !------------------------------------------------------------------------------
          !%***************************************************************************
          !%CorrectKsForPressureNow:
          !% Currently: For WhichKs% = 1 to 7, all Ks (except KF and KS, which are on
          !%       the free scale) are on the SWS scale.
          !%       For WhichKs% = 6, KW set to 0, KP1, KP2, KP3, KSi don't matter.
          !%       For WhichKs% = 8, K1, K2, and KW are on the "pH" pH scale
          !%       (the pH scales are the same in this case); the other Ks don't matter.
          !%
          !%
          !% No salinity dependence is given for the pressure coefficients here.
          !% It is assumed that the salinity is at or very near Sali = 35.
          !% These are valid for the SWS pH scale, but the difference between this and
          !% the total only yields a difference of .004 pH units at 1000 bars, much
          !% less than the uncertainties in the values.
          !%****************************************************************************
          !% The sources used are:
          !% Millero, 1995:
          !%       Millero, F. J., Thermodynamics of the carbon dioxide system in the
          !%       oceans, Geochemica et Cosmochemica Acta 59:661-677, 1995.
          !%       See table 9 and eqs. 90-92, p. 675.
          !%       TYPO: a factor of 10^3 was left out of the definition of Kappa
          !%       TYPO: the value of R given is incorrect with the wrong units
          !%       TYPO: the values of the a's for H2S and H2O are from the 1983
          !%                values for fresh water
          !%       TYPO: the value of a1 for B(OH)3 should be +.1622
          !%        Table 9 on p. 675 has no values for Si.
          !%       There are a variety of other typos in Table 9 on p. 675.
          !%       There are other typos in the paper, and most of the check values
          !%       given don't check.
          !% Millero, 1992:
          !%       Millero, Frank J., and Sohn, Mary L., Chemical Oceanography,
          !%       CRC Press, 1992. See chapter 6.
          !%       TYPO: this chapter has numerous typos (eqs. 36, 52, 56, 65, 72,
          !%               79, and 96 have typos).
          !% Millero, 1983:
          !%       Millero, Frank J., Influence of pressure on chemical processes in
          !%       the sea. Chapter 43 in Chemical Oceanography, eds. Riley, J. P. and
          !%       Chester, R., Academic Press, 1983.
          !%       TYPO: p. 51, eq. 94: the value -26.69 should be -25.59
          !%       TYPO: p. 51, eq. 95: the term .1700t should be .0800t
          !%       these two are necessary to match the values given in Table 43.24
          !% Millero, 1979:
          !%       Millero, F. J., The thermodynamics of the carbon dioxide system
          !%       in seawater, Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
          !%       See table 5 and eqs. 7, 7a, 7b on pp. 1656-1657.
          !% Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
          !%       TYPO: the pressure dependence of K2 should have a 16.4, not 26.4
          !%       This matches the GEOSECS results and is in Edmond and Gieskes.
          !% Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic acid,
          !%       boric acid, and the pH of seawater, Limnology and Oceanography
          !%       13:403-417, 1968.
          !% Edmond, John M. and Gieskes, J. M. T. M., The calculation of the degree of
          !%       seawater with respect to calcium carbonate under in situ conditions,
          !%       Geochemica et Cosmochemica Acta, 34:1261-1291, 1970.
          !%****************************************************************************
          !% These references often disagree and give different fits for the same thing.
          !% They are not always just an update either; that is, Millero, 1995 may agree
          !%       with Millero, 1979, but differ from Millero, 1983.
          !% For WhichKs% = 7 (Peng choice) I used the same factors for KW, KP1, KP2,
          !%       KP3, and KSi as for the other cases. Peng et al didn't consider the
          !%       case of P different from 0. GEOSECS did consider pressure, but didn't
          !%       include Phos, Si, or OH, so including the factors here won't matter.
          !% For WhichKs% = 8 (freshwater) the values are from Millero, 1983 (for K1, K2,
          !%       and KW). The other aren't used (TB = TS = TF = TP = TSi = 0.), so
          !%       including the factors won't matter.
          !%****************************************************************************
          !%       deltaVs are in cm3/mole
          !%       Kappas are in cm3/mole/bar
          !%****************************************************************************
!
          !------------------------------------------------------------------------------
          ! Correct K1, K2, KB For Pressure:
!
               If (WHICH_K1K2 .Eq. 8) Then
             !    %***PressureEffectsOnK1inFreshWater:
             !    %               This is from Millero, 1983.
                  deltaV = - 30.54_SP + 0.1849_SP * T (I, K) - &
                 & 0.0023366_SP * T (I, K) ** 2
                  Kappa = (-6.22_SP+0.1368_SP*T(I, K)-0.001233_SP*T(I, &
                 & K)**2) / 1000.0_SP
                  lnK1fac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
!
             !    %***PressureEffectsOnK2inFreshWater:
             !    %               This is from Millero, 1983.
                  deltaV = - 29.81_SP + 0.115_SP * T (I, K) - &
                 & 0.001816_SP * T (I, K) ** 2
                  Kappa = (-5.74_SP+0.093_SP*T(I, K)-0.001896_SP*T(I, &
                 & K)**2) / 1000.0_SP
                  lnK2fac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
                  lnKBfac = 0.0_SP !this doesn't matter since TB = 0 for this case
!
               Else If (WHICH_K1K2 .Eq. 6 .Or. WHICH_K1K2 .Eq. 7) Then
             !    %               GEOSECS Pressure Effects On K1, K2, KB (on the NBS scale)
             !    %               Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982 quotes
             !    %               Culberson and Pytkowicz, L and O 13:403-417, 1968:
             !    %               but the fits are the same as those in
             !    %               Edmond and Gieskes, GCA, 34:1261-1291, 1970
             !    %               who in turn quote Li, personal communication
                  lnK1fac = (24.2_SP-0.085_SP*T(I, K)) * Pbar / RT
                  lnK2fac = (16.4_SP-0.04_SP*T(I, K)) * Pbar / RT
             !    %               Takahashi et al had 26.4, but 16.4 is from Edmond and Gieskes
             !    %               and matches the GEOSECS results
                  lnKBfac = (27.5_SP-0.095_SP*T(I, K)) * Pbar / RT
!
               Else
             !    %***PressureEffectsOnK1:
             !    %               These are from Millero, 1995.
             !    %               They are the same as Millero, 1979 and Millero, 1992.
             !    %               They are from data of Culberson and Pytkowicz, 1968.
                  deltaV = - 25.5_SP + 0.1271_SP * T (I, K)
             !    %                 'deltaV = deltaV - .151 *(Sali - 34.8); % Millero, 1979
                  Kappa = (-3.08_SP+0.0877_SP*T(I, K)) / 1000.0_SP
             !    %                 'Kappa = Kappa  - .578 *(Sali - 34.8)/1000.; % Millero, 1979
                  lnK1fac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
             !    %               The fits given in Millero, 1983 are somewhat different.
!
             !    %***PressureEffectsOnK2:
             !    %               These are from Millero, 1995.
             !    %               They are the same as Millero, 1979 and Millero, 1992.
             !    %               They are from data of Culberson and Pytkowicz, 1968.
                  deltaV = - 15.82_SP - 0.0219_SP * T (I, K)
             !    %                  'deltaV = deltaV + .321 *(Sali - 34.8); % Millero, 1979
                  Kappa = (1.13_SP-0.1475_SP*T(I, K)) / 1000.0_SP
             !    %                 'Kappa = Kappa - .314 *(Sali - 34.8) /1000: % Millero, 1979
                  lnK2fac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
             !    %               The fit given in Millero, 1983 is different.
             !    %               Not by a lot for deltaV, but by much for Kappa. %
!
             !    %***PressureEffectsOnKB:
             !    %               This is from Millero, 1979.
             !    %               It is from data of Culberson and Pytkowicz, 1968.
                  deltaV = - 29.48_SP + 0.1622_SP * T (I, K) - &
                 & 0.002608_SP * T (I, K) ** 2
             !    %               Millero, 1983 has:
             !    %                 'deltaV = -28.56 + .1211 *TempCi - .000321 *TempCi *TempCi
             !    %               Millero, 1992 has:
             !    %                 'deltaV = -29.48 + .1622 *TempCi + .295 *(Sali - 34.8)
             !    %               Millero, 1995 has:
             !    %                 'deltaV = -29.48 - .1622 *TempCi - .002608 *TempCi *TempCi
             !    %                 'deltaV = deltaV + .295 *(Sali - 34.8); % Millero, 1979
                  Kappa = - 2.84_SP / 1000.0_SP ! Millero, 1979
             !    %               Millero, 1992 and Millero, 1995 also have this.
             !    %                 'Kappa = Kappa + .354 *(Sali - 34.8) /1000: % Millero,1979
             !    %               Millero, 1983 has:
             !    %                 'Kappa = (-3 + .0427 *TempCi) /1000
                  lnKBfac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
!
!
               End If
!
!
!
          !------------------------------------------------------------------------------
          ! Correct KW For Pressure:
               If (WHICH_K1K2 .Eq. 8) Then
             !    % PressureEffectsOnKWinFreshWater:
             !    %               This is from Millero, 1983.
                  deltaV = - 25.6_SP + 0.2324_SP * T (I, K) - &
                 & 0.0036246_SP * T (I, K) ** 2
                  Kappa = (-7.33_SP+0.1368_SP*T(I, K)-0.001233_SP*T(I, &
                 & K)**2) / 1000.0_SP
                  lnKWfac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
!
             !    %               NOTE the temperature dependence of KappaK1 and KappaKW
             !    %               for fresh water in Millero, 1983 are the same.
!
               Else
             !    % GEOSECS doesn't include OH term, so this won't matter.
             !    % Peng et al didn't include pressure, but here I assume that the KW correction
             !    %       is the same as for the other seawater cases.
             !    % PressureEffectsOnKW:
             !    %               This is from Millero, 1983 and his programs CO2ROY(T).BAS.
                  deltaV = - 20.02_SP + 0.1119_SP * T (I, K) - &
                 & 0.001409_SP * T (I, K) ** 2
             !    %               Millero, 1992 and Millero, 1995 have:
                  Kappa = (-5.13_SP+0.0794_SP*T(I, K)) / 1000.0_SP ! Millero, 1983
             !    %               Millero, 1995 has this too, but Millero, 1992 is different.
                  lnKWfac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
             !    %               Millero, 1979 does not list values for these.
!
               End If
!
          !------------------------------------------------------------------------------
          ! Pressure ffects On KF:
          !%       This is from Millero, 1995, which is the same as Millero, 1983.
          !%       It is assumed that KF is on the free pH scale.
               deltaV = - 9.78_SP - 0.009_SP * T (I, K) - 0.000942_SP * &
              & T (I, K) ** 2
               Kappa = (-3.91_SP+0.054_SP*T(I, K)) / 1000.0_SP
               lnKFfac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
!
          !------------------------------------------------------------------------------
          !% Pressure ffects On KS:
          !%       This is from Millero, 1995, which is the same as Millero, 1983.
          !%       It is assumed that KS is on the free pH scale.
               deltaV = - 18.03_SP + 0.0466_SP * T (I, K) + 0.000316_SP &
              & * T (I, K) ** 2
               Kappa = (-4.53_SP+0.09_SP*T(I, K)) / 1000.0_SP
               lnKSfac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
!
          !------------------------------------------------------------------------------
          ! Correct KP1, KP2, KP3, KSi For Pressure:
          !% These corrections don't matter for the GEOSECS choice (WhichKs% = 6) and
          !%       the freshwater choice (WhichKs% = 8). For the Peng choice I assume
          !%       that they are the same as for the other choices (WhichKs% = 1 to 5).
          !% The corrections for KP1, KP2, and KP3 are from Millero, 1995, which are the
          !%       same as Millero, 1983.
!
          !% PressureEffectsOnKP1:
               deltaV = - 14.51_SP + 0.1211_SP * T (I, K) - 0.000321_SP &
              & * T (I, K) ** 2
               Kappa = (-2.67_SP+0.0427_SP*T(I, K)) / 1000.0_SP
               lnKP1fac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
!
          !% PressureEffectsOnKP2:
               deltaV = - 23.12_SP + 0.1758_SP * T (I, K) - 0.002647_SP &
              & * T (I, K) ** 2
               Kappa = (-5.15_SP+0.09_SP*T(I, K)) / 1000.0_SP
               lnKP2fac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
!
          !% PressureEffectsOnKP3:
               deltaV = - 26.57_SP + 0.202_SP * T (I, K) - 0.003042_SP &
              & * T (I, K) ** 2
               Kappa = (-4.08_SP+0.0714_SP*T(I, K)) / 1000.0_SP
               lnKP3fac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
!
          !% PressureEffectsOnKSi:
          !%  The only mention of this is Millero, 1995 where it is stated that the
          !%    values have been estimated from the values of boric acid. HOWEVER,
          !%    there is no listing of the values in the table.
          !%    I used the values for boric acid from above.
               deltaV = - 29.48_SP + 0.1622_SP * T (I, K) - 0.002608_SP &
              & * T (I, K) ** 2
               Kappa = - 2.84_SP / 1000.0_SP
               lnKSifac = (-deltaV+0.5_SP*Kappa*Pbar) * Pbar / RT
!
!
          !------------------------------------------------------------------------------
          ! Correct Ks For Pressure Here:
!
               aux = Exp (lnK1fac)
               K1 (I, K) = K1 (I, K) * aux
!
               aux = Exp (lnK2fac)
               K2 (I, K) = K2 (I, K) * aux
!
               aux = Exp (lnKWfac)
               KW (I, K) = KW (I, K) * aux
!
               aux = Exp (lnKBfac)
               KBb (I, K) = KBb (I, K) * aux
!
               aux = Exp (lnKFfac)
               KF (I, K) = KF (I, K) * aux
!
               aux = Exp (lnKSfac)
               KS (I, K) = KS (I, K) * aux
!
          !aux = exp(lnKP1fac)           !not considering phosphate in ph !LB
          !KP1(I,K) = KP1(I,K) *aux
!
          !aux = exp(lnKP2fac)
          !KP2(I,K) = KP2(I,K) *aux
!
          !aux = exp(lnKP3fac)
          !KP3(I,K) = KP3(I,K) *aux
!
          !aux = exp(lnKSifac)
          !KSi(I,K) = KSi(I,K) *aux    !not modeling Silicate !LB
!
          !------------------------------------------------------------------------------
          !******************************************************************************
          !     At this point, all K constants are in SWS scale.
          !     If we wanted to calculate them in another scale, we would have to
          !     re-calculate SWStoTOT and FREEtoTOT with KS and KF pressure-corrected
          !    and then transform the Ks to ather scale
          !
          !    Since we choose to work with SWS pH scale in FVCOM-ICM, these transformations
          !     (available in CO2SYS.m) are not included here
          !******************************************************************************
          !------------------------------------------------------------------------------
!
          !------------------------------------------------------------------------------
          ! Calculate Fugacity Constants:
          !% This assumes that the pressure is at one atmosphere, or close to it.
          !% Otherwise, the Pres term in the exponent affects the results.
          !%       Weiss, R. F., Marine Chemistry 2:203-215, 1974.
          !%       Delta and aux in cm3/mol
          !**LB note: checking the publication above, the approximation is OK at
          !           pressures typical of natural environments (and up to 10 atm)
!
               Delta = 57.7_SP - 0.118_SP * TempK
               aux = - 1636.75_SP + 12.0408_SP * TempK - 0.0327957_SP * &
              & TempK ** 2 + 3.16528_SP * 0.00001_SP * TempK ** 3
!
          !% For a mixture of CO2 and air at 1 atm (at low CO2 concentrations);
               FugFac (I, K) = Exp ((aux+2.0_SP*Delta)*P1atm/RT)
!
               If (WHICH_K1K2 .Eq. 6 .Or. WHICH_K1K2 .Eq. 7) Then
             !% GEOSECS and Peng assume pCO2 = fCO2, or FugFac = 1
                  FugFac (I, K) = 1.0_SP
               End If
!
          !------------------------------------------------------------------------------
          ! Calculate VPFac:
          !% Weiss, R. F., and Price, B. A., Nitrous oxide solubility in water and
          !%       seawater, Marine Chemistry 8:347-359, 1980.
          !% They fit the data of Goff and Gratch (1946) with the vapor pressure
          !%       lowering by sea salt as given by Robinson (1954).
          !% This fits the more complicated Goff and Gratch, and Robinson equations
          !%       from 273 to 313 deg K and 0 to 40 Sali with a standard error
          !%       of .015%, about 5 uatm over this range.
          !% This may be on IPTS-29 since they didn't mention the temperature scale,
          !%       and the data of Goff and Gratch came before IPTS-48.
          !% The references are:
          !% Goff, J. A. and Gratch, S., Low pressure properties of water from -160 deg
          !%       to 212 deg F, Transactions of the American Society of Heating and
          !%       Ventilating Engineers 52:95-122, 1946.
          !% Robinson, Journal of the Marine Biological Association of the U. K.
          !%       33:449-455, 1954.
          !%       This is eq. 10 on p. 350.
          !%       This is in atmospheres.
               aux = Exp (24.4543_SP-67.4509_SP*(100.0_SP/TempK)-&
              & 4.8489_SP*Log(TempK/100.0_SP))!VPWP in CO2SYS.m
               aux2 = Exp (-0.000544_SP*S)!VPCorrWP in CO2SYS.m
               aux = aux * aux2 !VPSWWP   in CO2SYS.m
               VPFac (I, K) = 1.0_SP - aux ! this assumes 1 atmosphere
!
          !------------------------------------------------------------------------------
            End Do
         End Do
!
      End Subroutine CO2SYSCONSTANTS
!
End Module MOD_CO2SYS
