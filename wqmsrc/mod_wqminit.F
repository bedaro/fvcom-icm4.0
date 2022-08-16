!mod_wqminit.F
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
Module MOD_WQMINIT
  !
      Use MOD_PREC, Only: SP
      Use MOD_SIZES, Only: NCP, NHQP, NBCP, NDP, NOIP, NSSFP
  !
      Use MOD_LIMS, Only: KBM1, MTLOC
  !
      Implicit None
  !
      Save
  !***** Variable declarations
  !
      Integer :: RSODP, DLTDP, SNPDP, TFLDP, PLTDP, APLDP, OPLDP, DIADP
      Integer :: COURB, COURFS, COURBS, F, SB, S1LNMAX, S2LNMAX, &
     & S3LNMAX, DIFFFS, DIFFBS
  !
      Real (SP) :: NXSNP, NXPLT, NXAPL, NXTFL, NXKFL, NXTVD, NXOPL, &
     & NXMBL, NXDIA
  !
  !            REAL(SP)::  NO3T2I,        &!    Moved this to mod_sed.F
  !                        NH4T2I
  !
      Real (SP) :: MAXDLT, MXDLT
      Real (SP) :: TM1, TM2, TM3, MINSTEP
  !
      Character (Len=3) :: SNPC, RSIC, RSOC, BCC, S1C, S2C, MDC * 3, &
     & PLTC, FLC, MBLC, BFC, VBC, QPLTC, XYDFC, XYDFU, ZDFC, ICOC, &
     & ATMC, SAVLC, SEDC,SEDTC, AUTOC, SPLTC, TFLC, DIAC, STLC, APLTC, KFLC, &
     & OPLC, BFOC, BAOC, DFOC, S3C, SAVMC, SAVPLTC, DFLC, SFLC, SFLOXC, &
     & DFLOXC
  !
      Character :: EXT1 * 1, EXT2 * 2, EXT3 * 3, EXT4 * 4, EXT_R * 1, &
     & EXT_R2 * 2
  !
      Character (Len=8) :: SLC, HYDC, BNDTC, CONSC, ICIC
  !
      Character (Len=8) :: TIMEA, TIMEB, SPVARM, PRINTM
  !
      Character (Len=72) :: TITLE (6), OLDTITLE (6), FILENAME
  !
      Character (Len=100) :: FILENAME2
  !
      Character (Len=72) :: MAPFN, GEOFN, ICIFN, RSIFN, AGRFN, STLFN, &
     & MRLFN, KFLFN, ICOFN
  !
      Character (Len=72) :: SNPFN, PLTFN, APLFN, DIAFN, TFLFN, RSOFN, &
     & OPLFN, MBLFN, SFIFN, SFOFN
  !
      Logical :: RESTART_OUT, SNAPSHOTS, END_RUN, MODIFY_ICONC, &
     & VOLUME_BALANCE, QUICKEST, UPWIND, ICOND_OUT, UNI_ICON_IN, &
     & UNI_ICON_IN_SED_VAR, BIN_ICON_IN, RES_ICON_IN, AUTO_STEPPING, &
     & STOP_RUN, PLOTS, OXYGEN_PLOTS, ZOO_CALC, NEW_VOLUMES, RESTART_IN
  !
      Logical :: TEMPERATURE_CALC, ALGAE_CALC, CARBON_CALC, &
     & NITROGEN_CALC, PHOSPHORUS_CALC, COD_CALC, OXYGEN_CALC, &
     & SILICA_CALC
  !
      Logical :: LEFT_FLOWB (NHQP), RIGHT_FLOWB (NHQP), LEFTM1_BOUNDARY &
     & (NHQP), RIGHTP1_BOUNDARY (NHQP)
  !
      Logical :: IS_ATM_FCO2, CARBONATE_CALC !
  !      Logical :: SAVE_PH_OUTPUT, SAVE_PCO2_OUTPUT  !
  !
  !***** Dimension declarations
  !
      Real (SP), Dimension (NHQP, 2) :: DEN1, DEN2, DEN3, TP1, TP2, &
     & TP3, T2, SF2
  !
      Real (SP), Dimension (NHQP) :: T1, T3, SF1
  !
      Real (SP), Dimension (NHQP) :: TERM1, TERM2, TERM3, GRAD1, GRAD2, &
     & GRAD3, COUR
  !
      Real (SP) :: GRAD (NHQP, 3), TERM (NHQP, 3), IFLOWP (NBCP), IBT &
     & (NHQP), COURBSV (NHQP), COURVSV (NHQP)
  !
      Equivalence (GRAD1, GRAD(1, 1))
      Equivalence (GRAD2, GRAD(1, 2))
      Equivalence (GRAD3, GRAD(1, 3))
      Equivalence (TERM1, TERM(1, 1))
      Equivalence (TERM2, TERM(1, 2))
      Equivalence (TERM3, TERM(1, 3))
!
      Real (SP), Allocatable, Dimension (:, :, :) :: C1MIN, C1MAX
!
      Real (SP), Allocatable :: DOVDAYS (:, :, :)
      Real (SP) :: OINT (NOIP)
      Character (Len=3) :: ACC (NCP)
  !
      Real (SP) :: CIC (NCP), COUT (NCP)
  !
      Real (SP) :: SFEEDI (NSSFP)!WLong changed dimension size from 10 to NSSFP
  !
  !
  !
      Real (SP) :: DLTD (NDP), SNPD (NDP), RSOD (NDP), SNPF (NDP), &
     & DLTVAL (NDP), DLTMAX (NDP), DLTFTN (NDP), PLTD (NDP), PLTF &
     & (NDP), APLTD (NDP), APLF (NDP), TFLD (NDP), TFLF (NDP), KFLD &
     & (NDP), KFLF (NDP), OPLD (NDP), OPLF (NDP), DIAD (NDP), DIAF &
     & (NDP), MBLD (NDP), MBLF (NDP)
  !
  !***** Data declarations
      Real (SP), Parameter :: DOCLIT = 0.0, LPOCLIT = 0.0, POCLIT = &
     & 0.0, PBSLIT = 0.0, DSILLIT = 0.0
  !
Contains
  !
  !subroutines:
  !subroutine     WQMINIT_ALLOC()
  !
      Subroutine WQMINIT_ALLOC
         Implicit None
    !
         Allocate (C1MIN(MTLOC, KBM1, NCP))
         C1MIN = 0.0
         Allocate (C1MAX(MTLOC, KBM1, NCP))
         C1MAX = 0.0
         Allocate (DOVDAYS(MTLOC, KBM1, NOIP))
         DOVDAYS = 0.0
         PHOSPHORUS_CALC = .False.
         IS_ATM_FCO2 = .False. !
         Return
      End Subroutine WQMINIT_ALLOC
  !
      Subroutine WQMINIT_DEALLOC
         Implicit None
         If (ALLOCATED(C1MIN)) DEALLOCATE (C1MIN)
         If (ALLOCATED(C1MAX)) DEALLOCATE (C1MAX)
         If (ALLOCATED(DOVDAYS)) DEALLOCATE (DOVDAYS)
         Return
      End Subroutine WQMINIT_DEALLOC
  !
End Module MOD_WQMINIT
!
