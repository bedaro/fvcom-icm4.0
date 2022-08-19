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
Module MOD_KIN
      Use MOD_PREC, Only: SP, SMOD
      Use MOD_ZOOP, Only: CTSZ, CTLZ, KHCSZ, KHCLZ, MZEROSZ, MZEROLZ, &
     & UB1SZ, UB2SZ, UB3SZ, UB1LZ, UB2LZ, UB3LZ, ULDSZ, URDSZ, ULPSZ, &
     & URPSZ, ULLZ, URLZ, USZLZ, DOCRITSZ, DOCRITLZ, ANCSZ, ANCLZ, &
     & APCSZ, APCLZ, AOCRSZ, AOCRLZ, FRSASZ, FRSALZ, FLDOCSZ, FRDOCSZ, &
     & FLPOCSZ, FRPOCSZ, FLDONSZ, FRDONSZ, FLPONSZ, FRPONSZ, FLDOPSZ, &
     & FRDOPSZ, FLPOPSZ, FRPOPSZ, FNH4SZ, FPO4SZ, FLDOCLZ, FRDOCLZ, &
     & FLPOCLZ, FRPOCLZ, FLDONLZ, FRDONLZ, FLPONLZ, FRPONLZ, FLDOPLZ, &
     & FRDOPLZ, FLPOPLZ, FRPOPLZ, FNH4LZ, FPO4LZ, B1ASZ, B2ASZ, B3ASZ, &
     & LPOCASZ, RPOCASZ, PRASZ, B1ALZ, B2ALZ, B3ALZ, SZALZ, LPOCALZ, &
     & RPOCALZ, PRALZ, CLSZ, CLLZ, RSZ, RLZ, RMAXSZ, RMAXLZ, BMSZ, &
     & BMLZ, BMRSZ, BMRLZ, MSZ, MLZ, PRSZLZ, GSZ, GLZ, ESZ, ELZ, RFSZ, &
     & RFLZ, PRSZ, PRLZ, LDOCASZ, BPRSZ, BPRLZ, RDOCASZ, LDOCSZ, &
     & LPOCSZ, RPOCSZ, LDOCLZ, LPOCLZ, RPOCLZ, NH4SZ, LDONSZ, LPONSZ, &
     & RPONSZ, NH4LZ, LDONLZ, LPONLZ, RPONLZ, PO4SZ, LDOPSZ, LPOPSZ, &
     & RPOPSZ, PO4LZ, LDOPLZ, LPOPLZ, RPOPLZ, RDOCSZ, RDONSZ, RDOPSZ, &
     & RDOCLZ, RDONLZ, RDOPLZ, PIB1SZ, PIB2SZ, PIB3SZ, PIB1LZ, PIB2LZ, &
     & PIB3LZ, B1SZ, B2SZ, B3SZ, B1LZ, B2LZ, B3LZ, DOSZ, DOLZ, SASZ, &
     & SUSZ, SALZ, SULZ, FTLZ, FTSZ, FTBMSZ, FTBMLZ, FTPRSZ, FTPRLZ, &
     & SZ, LZ, TVPRSZ, TVPRLZ, GET_FT,GET_FTPR,GET_FTBM,TMSZ,KTGSZ1, &
     & KTGSZ2,KTBSZ,KTBLZ,TRSZ,TRLZ,TRPRSZ,KTPRSZ,TRPRLZ,KTPRLZ,KTGLZ1,&
     & KTGLZ2,TMLZ

!
  !************************************************************************
  !** additional paramters were made spatially variable these include:   **
  !**     KDC,    KLC,    KRC,    KDN,    KLN,    KRN,    KDP,           **
  !**     KLP,    KRP,   KSUA,   KCOD, KDCALG, KLCALG, KRCALG,           **
  !**  KDNALG, KLNALG, KRNALG, KDPALG, KLPALG, KRPALG,    NTM,           **
  !**      DL,      R        ***** MNOEL 2-20-93 *****                   **
  !**  **** took out kfl_%%% subroutines that performed writing          **
  !**       operations and included them in main program so they could   **
  !**       be averaged over time periods instead of instantaneous       **
  !**       output.  This process included vectorizing some variables.   **
  !**                        ***** MNOEL 7-28-93 *****                   **
  !**                                                                    **
  !**                Kinetics Subroutines for CE-QUAL-IC                 **
  !**                                                                    **
  !**                            Version 1.0                             **
  !**                         February 24, 1989                          **
  !**                                                                    **
  !**                    Water Quality Modeling Group                    **
  !**                    U.S. Army Corps of Engineers                    **
  !**                    Waterways Experiment Station                    **
  !**                    Vicksburg, Mississippi 39180                    **
  !**                                                                    **
  !************************************************************************
  !
Contains
  !Subroutines:
  !	Subroutine TEMPER()
  !	Subroutine SOLIDS()
  !	Subroutine ZOOPL()
  !	Subroutine CARBON()
  !	Subroutine NITROG()
  !	Subroutine PHOSPH()
  !	Subroutine CODMND()
  !	Subroutine OXYGEN()
  !	Subroutine SILICA()
  !	Subroutine BEN_FLUX()
  !
  !
  !************************************************************************
  !**                 S U B R O U T I N E   T E M P E R                  **
  !************************************************************************
  !
      Subroutine TEMPER (DTT)
    !
         Use MOD_WQM, Only: KT, T, TE !Equilibrium temperature (degC)
    !
         Use MOD_LIMS, Only: MLOC, KBM1, MTLOC
    !
         Use MOD_HYDROVARS, Only: DZ,DZ2D, D
    !
# if defined (WET_DRY)
	USE  MOD_WD, ONLY: ISWETN
# endif
         Implicit None
         Real (SP) :: RHO_, CP
         Real (SP) :: DTT (0:MTLOC, KBM1)
         Integer :: I
         Data RHO_ / 1.0E6 /, CP / 4.1796 /
    !
         DTT = 0.0
    !
         Do I = 1, MLOC
#  if !defined (WET_DRY)
			IF(D(I) > 0.0_SP)THEN
#  else
			IF(ISWETN(I) == 1)THEN
#  endif
				DTT (I, 1) = KT / (RHO_*CP*D(I)*DZ2D(I,1)) * (TE-T(I, 1))
			ELSE
				DTT(I,1)=0.0_SP
			ENDIF
         End Do
    !
         Return
      End Subroutine TEMPER
  !
  !************************************************************************
  !**                 S U B R O U T I N E   S O L I D S                  **
  !************************************************************************
  !
      Subroutine SOLIDS (DTSSI, FLXSSSI)
    !
         Use MOD_WQM, Only: SEDIMENT_CALC, SSI, V2, WSS, WSSHI, WSSNET 
    !
         Use MOD_LIMS, Only: MLOC, KBM1, MTLOC
# if defined (WET_DRY)
	USE  MOD_WD, ONLY: ISWETN
# endif
    !
         Use MOD_HYDROVARS, Only: DZ, DZ2D, D
    !
         Implicit None
         Integer :: I, K
         Real (SP), Dimension (0:MTLOC, KBM1) :: DTSSI, FLXSSSI
    !******* Set special settling rate for high solids concs
         DTSSI = 0.0
         FLXSSSI = 0.0
    !
         WSSHI = 0.0
         Do K = 1, KBM1
            Do I = 1, MLOC
               If (SSI(I, K) < 100.) Then !Adi 2/16/18
                  WSSHI (I, K) = WSS (I, K)
               Else
                 WSSHI (I, K) = 5.0
               End If
            End Do
         End Do
    !
    !******* Settle particulate fraction
    !
         Do I = 1, MLOC

#  if !defined (WET_DRY)
			IF(D(I) > 0.0_SP)THEN
#  else
			IF(ISWETN(I) == 1)THEN
#  endif
            FLXSSSI (I, 1) = WSSHI (I, 1) * SSI (I, 1) * V2 (I, 1) / &
           & (D(I)*DZ2D(I,1)*86400.)
            DTSSI (I, 1) = - WSSHI (I, 1) * SSI (I, 1) / &
           & (D(I)*DZ2D(I,1)*86400.)
		    ELSE
				FLXSSSI(I,1)=0.0_SP
				DTSSI(I,1)=0.0_SP
			ENDIF
         End Do
    !
         Do K = 2, KBM1
            Do I = 1, MLOC
#  if !defined (WET_DRY)
			IF(D(I) > 0.0_SP)THEN
#  else
			IF(ISWETN(I) == 1)THEN
#  endif
               FLXSSSI (I, K) = WSSHI (I, K) * SSI (I, K) * V2 (I, K) / &
              & (D(I)*DZ2D(I,K)*86400.)
               DTSSI (I, K) = (WSSHI(I, K-1)*SSI(I, K-1)-WSSHI(I, &
              & K)*SSI(I, K)) / (D(I)*DZ2D(I,K)*86400.)

			ELSE
				FLXSSSI(I,K)=0.0_SP
				DTSSI(I,K)=0.0_SP
			ENDIF
            End Do
         End Do
    !
    !******  Resuspension
    !
         If (SEDIMENT_CALC) Then
       !
            Do I = 1, MLOC
#  if !defined (WET_DRY)
				IF(D(I) > 0.0_SP)THEN
#  else
				IF(ISWETN(I) == 1)THEN
#  endif
					DTSSI (I, KBM1) = DTSSI (I, KBM1) + (WSSHI(I, &
              & KBM1)-WSSNET(I)) * SSI (I, KBM1) / (D(I)*DZ2D(I,KBM1)) / &
              & 86400.

				ELSE
					DTSSI (I, KBM1)=0.0_SP
				ENDIF
            End Do
         End If
    !
         Return
      End Subroutine SOLIDS
  !
  !
  !************************************************************************
  !**                  S U B R O U T I N E   Z O O P L                   **
  !************************************************************************
  !
      Subroutine ZOOPL (DTSZ, DTLZ)
    !
    !
         Use MOD_LIMS, Only: MLOC, KBM1, MTLOC
    !
         Use MOD_WQM, Only: ANC1, ANC2, ANC3, ASC1, ASC2, ASC3, B1, B2, &
        & B3, DOXG, JDAY, LDOC, LPOC, LPON, LPOP, Q1, Q2, Q3, RDOC, &
        & RPOC, RPON, RPOP, T 
    !
    !
         Implicit None
         Real (SP) :: NRRATE, DOREF, CRRATE, FRLD, FRRD, FRLP, FRRP, &
        & ANCPSZ, FRL, FRR, ANCPLZ, APCPSZ, PRRATE, APCPLZ, FRB1SZ, &
        & FRB2SZ, FRB3SZ, FRB1LZ, FRB2LZ, FRB3LZ
         Real (SP), Dimension (0:MTLOC, KBM1) :: DTSZ, DTLZ
         Integer :: I, J, K, ITEMP
    !
    ! DETERMINE JULIAN DAY
         DTSZ = 0.0
         DTLZ = 0.0
    !
         J = AINT(1.0 + SMOD (JDAY, 365.25_SP))
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
          ! COMPUTE PREY AVAILABLE TO ZOOPLANKTON
          !
               B1ASZ (I, K) = Max (0., B1(I, K)-CTSZ)
               B2ASZ (I, K) = Max (0., B2(I, K)-CTSZ)

               LDOCASZ (I, K) = Max (0., LDOC(I, K)-CTSZ)
               RDOCASZ (I, K) = Max (0., RDOC(I, K)-CTSZ)
               LPOCASZ (I, K) = Max (0., LPOC(I, K)-CTSZ)
               RPOCASZ (I, K) = Max (0., RPOC(I, K)-CTSZ)

			   PRASZ (I, K) = UB1SZ * B1ASZ (I, K) + UB2SZ * B2ASZ (I, &
              & K) + ULDSZ * LDOCASZ (I, K) + URDSZ * RDOCASZ (I, K) + &
			  & ULPSZ * LPOCASZ (I, K) + URPSZ * RPOCASZ (I, K)
          !
               B1ALZ (I, K) = Max (0., B1(I, K)-CTLZ)
               B2ALZ (I, K) = Max (0., B2(I, K)-CTLZ)

               SZALZ (I, K) = Max (0., SZ(I, K)-CTLZ)
               LPOCALZ (I, K) = Max (0., LPOC(I, K)-CTLZ)
               RPOCALZ (I, K) = Max (0., RPOC(I, K)-CTLZ)

			   PRALZ (I, K) = UB1LZ * B1ALZ (I, K) + UB2LZ * B2ALZ (I, &
              & K) + USZLZ * SZALZ (I, K) + ULLZ * LPOCALZ (I, K) + &
			  & URLZ * RPOCALZ (I, K)
          !
            End Do
         End Do
    !
    ! COMPUTE ZOOPLANKTON RATION AND METABOLISM
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
          ! NUTRIENT LIMITATION
          IF ((KHCSZ /= 0.0) .AND. (KHCLZ /= 0.0)) THEN
               CLSZ (I, K) = PRASZ (I, K) / (KHCSZ+PRASZ(I, K)+1.0E-30)
               CLLZ (I, K) = PRALZ (I, K) / (KHCLZ+PRALZ(I, K)+1.0E-30)
          ELSE
              WRITE(*,*)'Make sure the values of KHCSZ and KHCLZ are non-zero'
              STOP
          END IF
          !
          ! RATION (GM PREY CARBON/GM ZOOPL CARBON/DAY)
          !
               ITEMP = 10. * (T(I, K)+0.05)
               RSZ (I, K) = GET_FT (T(I, K),TMSZ,KTGSZ1,KTGSZ2) * CLSZ (I, K) * RMAXSZ (I, K)
               ! The above line uses the bi-model temperature function
               ! However, we replaced the temperature function according to the theory
                RLZ (I, K) = GET_FT (T(I, K),TMLZ,KTGLZ1,KTGLZ2) * CLLZ (I, K) * RMAXLZ (I, K)

              IF(RSZ(I, K) > 1.0E+10) THEN
                 WRITE(*,*)'RSZ goes to infinity....', CLSZ (I, K),RMAXSZ (I, K)
                 STOP
               END IF
          !
          ! BASAL METABOLISM (1/DAY)
          !
               BMSZ (I, K) = BMRSZ (I, K) * GET_FTBM (T(I, K), TRSZ, KTBSZ)

               BMLZ (I, K) = BMRLZ (I, K) * GET_FTBM (T(I, K), TRLZ, KTBLZ)
          ! LOSS TO PREDATION (GM C/M**3/DAY)
          !
            IF (PRALZ(I,K) /= 0.0) THEN
                PRSZLZ (I, K) = USZLZ * SZALZ (I, K) / (PRALZ(I, &
               & K)+1.0E-30) * RLZ (I, K) * LZ (I, K)
             ELSE
               PRSZLZ = 0.0
             END IF

               PRSZ (I, K) = BPRSZ (I, K) * SZ (I, K) * SZ (I, K) * &
              & GET_FTPR (T(I, K), TRPRSZ, KTPRSZ) * TVPRSZ (J)

              PRLZ (I, K) = BPRLZ (I, K) * LZ (I, K) * LZ (I, K) * &
              & GET_FTPR (T(I, K), TRPRLZ, KTPRLZ) * TVPRLZ (J)
          !
            End Do
         End Do
    !
    ! COMPUTE MORTALITY FACTORS
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
          ! DISSOLVED OXYGEN EFFECTS (1/DAY)
          !
               DOREF = Min (DOCRITSZ, DOXG(I, K))
               MSZ (I, K) = MZEROSZ * (1.-DOREF/DOCRITSZ)
               DOREF = Min (DOCRITLZ, DOXG(I, K))
               MLZ (I, K) = MZEROLZ * (1.-DOREF/DOCRITLZ)
          !
            End Do
         End Do
    !
    ! RATE OF CHANGE OF ZOOPLANKTON BIOMASS (GM C/M**3/DAY)
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
               GSZ (I, K) = ESZ (I, K) * (1-RFSZ(I, K)) * RSZ (I, K)
               GLZ (I, K) = ELZ (I, K) * (1-RFLZ(I, K)) * RLZ (I, K)
               DTSZ (I, K) = ((GSZ(I, K)-BMSZ(I, K)-MSZ(I, K))*SZ(I, &
              & K)-PRSZLZ(I, K)-PRSZ(I, K)) / 86400.
               DTLZ (I, K) = ((GLZ(I, K)-BMLZ(I, K)-MLZ(I, K))*LZ(I, &
              & K)-PRLZ(I, K)) / 86400.
          !
            End Do
         End Do
    !--------------------------------------------------------------------------------
    ! EFFECT OF ZOOPLANKTON ON ALGAE (GM C/M**3/DAY)
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
             IF(PRASZ(I,K) /= 0.0) THEN
               B1SZ (I, K) = UB1SZ * B1ASZ (I, K) / (PRASZ(I, &
              & K)+1.0E-30) * RSZ (I, K) * SZ (I, K)
               B2SZ (I, K) = UB2SZ * B2ASZ (I, K) / (PRASZ(I, &
              & K)+1.0E-30) * RSZ (I, K) * SZ (I, K)
               B3SZ (I, K) = UB3SZ * B3ASZ (I, K) / (PRASZ(I, &
              & K)+1.0E-30) * RSZ (I, K) * SZ (I, K)
            ELSE
              B1SZ (I, K) = 0.0
              B2SZ (I, K) = 0.0
              B3SZ (I, K) = 0.0
            END IF

            IF(PRALZ(I,K) /= 0.0) THEN
               B1LZ (I, K) = UB1LZ * B1ALZ (I, K) / (PRALZ(I, &
              & K)+1.0E-30) * RLZ (I, K) * LZ (I, K)
               B2LZ (I, K) = UB2LZ * B2ALZ (I, K) / (PRALZ(I, &
              & K)+1.0E-30) * RLZ (I, K) * LZ (I, K)
               B3LZ (I, K) = UB3LZ * B3ALZ (I, K) / (PRALZ(I, &
              & K)+1.0E-30) * RLZ (I, K) * LZ (I, K)
            ELSE
              B1LZ (I, K) = 0.0
              B2LZ (I, K) = 0.0
              B3LZ (I, K) = 0.0
            END IF

          !
            End Do
         End Do
    !---------------------------------------------------------------------------------
    ! EFECT OF ZOOPLANKTON ON CARBON (GM C/M**3/DAY)
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
               CRRATE = ((1.-ESZ(I, K))*RSZ(I, K)+MSZ(I, K)) * SZ (I, &
              & K) + PRSZ (I, K)
            IF(PRASZ(I, K) /= 0.0) THEN
               FRLD = ULDSZ * LDOCASZ (I, K) / (PRASZ(I, K)+1.0E-30)
               FRRD = URDSZ * RDOCASZ (I, K) / (PRASZ(I, K)+1.0E-30)
               FRLP = ULPSZ * LPOCASZ (I, K) / (PRASZ(I, K)+1.0E-30)
               FRRP = URPSZ * RPOCASZ (I, K) / (PRASZ(I, K)+1.0E-30)
             ELSE
               FRLD = 0.0
               FRRD = 0.0
               FRLP = 0.0
               FRRP = 0.0
             END IF

               LDOCSZ (I, K) = - FRLD * RSZ (I, K) * SZ (I, K) + CRRATE &
              & * FLDOCSZ
               RDOCSZ (I, K) = - FRRD * RSZ (I, K) * SZ (I, K) + CRRATE &
              & * FRDOCSZ
               LPOCSZ (I, K) = - FRLP * RSZ (I, K) * SZ (I, K) + CRRATE &
              & * FLPOCSZ
               RPOCSZ (I, K) = - FRRP * RSZ (I, K) * SZ (I, K) + CRRATE &
              & * FRPOCSZ
          !
               CRRATE = ((1.-ELZ(I, K))*RLZ(I, K)+MLZ(I, K)) * LZ (I, &
              & K) + PRLZ (I, K)
            IF(PRALZ(I, K) /= 0.0) THEN
               FRLP = ULLZ * LPOCALZ (I, K) / (PRALZ(I, K)+1.0E-30)
               FRRP = URLZ * RPOCALZ (I, K) / (PRALZ(I, K)+1.0E-30)
             ELSE
               FRLP = 0.0
               FRRP = 0.0
             END IF
               LDOCLZ (I, K) = CRRATE * FLDOCLZ
               RDOCLZ (I, K) = CRRATE * FRDOCLZ
               LPOCLZ (I, K) = - FRLP * RLZ (I, K) * LZ (I, K) + CRRATE &
              & * FLPOCLZ
               RPOCLZ (I, K) = - FRRP * RLZ (I, K) * LZ (I, K) + CRRATE &
              & * FRPOCLZ
          !
            End Do
         End Do
    !---------------------------------------------------------------------------------
    ! EFFECT OF ZOOPLANKTON ON NITROGEN (GM N/M**3/DAY)
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
            IF(PRASZ(I, K) /= 0.0) THEN
               ANCPSZ = (ANC1*UB1SZ*B1ASZ(I, K)+ANC2*UB2SZ*B2ASZ(I, &
              & K)+ANC3*UB3SZ*B3ASZ(I, K)+ULPSZ*LPON(I, K)*LPOCASZ(I, &
              & K)/(LPOC(I, K)+1.0E-30)+URPSZ*RPON(I, K)*RPOCASZ(I, &
              & K)/(RPOC(I, K)+1.0E-30)) / (PRASZ(I, K)+1.0E-30)
            ELSE
              ANCPSZ = 0.0
            END IF

               NRRATE = ((ANCPSZ-ANCSZ*ESZ(I, K)*(1.0-RFSZ(I, &
              & K)))*RSZ(I, K)+(BMSZ(I, K)+MSZ(I, K))*ANCSZ) * SZ (I, &
              & K) + PRSZ (I, K) * ANCSZ

            IF(PRASZ(I, K) /= 0.0) THEN
               FRL = ULPSZ * LPOCASZ (I, K) / (PRASZ(I, K)+1.0E-30)
               FRR = URPSZ * RPOCASZ (I, K) / (PRASZ(I, K)+1.0E-30)
            ELSE
               FRL = 0.0
               FRR = 0.0
            END IF

               NH4SZ (I, K) = NRRATE * FNH4SZ
               IF(NH4SZ (I, K) > 1.0E+10) THEN
                 WRITE(*,*)'NH4SZ goes to infinity',NRRATE,CLSZ(I, K),PRASZ(I, K),RSZ(I, K),FTSZ (ITEMP),RMAXSZ (I, K)
                 STOP
               END IF
               LDONSZ (I, K) = NRRATE * FLDONSZ
               RDONSZ (I, K) = NRRATE * FRDONSZ
               LPONSZ (I, K) = - FRL * LPON (I, K) / (LPOC(I, &
              & K)+1.0E-30) * RSZ (I, K) * SZ (I, K) + NRRATE * FLPONSZ
               RPONSZ (I, K) = - FRR * RPON (I, K) / (RPOC(I, &
              & K)+1.0E-30) * RSZ (I, K) * SZ (I, K) + NRRATE * FRPONSZ
          !
            IF(PRALZ(I, K) /= 0.0) THEN
               ANCPLZ = (ANC1*UB1LZ*B1ALZ(I, K)+ANC2*UB2LZ*B2ALZ(I, &
              & K)+ANC3*UB3LZ*B3ALZ(I, K)+ANCSZ*USZLZ*SZALZ(I, &
              & K)+ULLZ*LPON(I, K)*LPOCALZ(I, K)/(LPOC(I, &
              & K)+1.0E-30)+URLZ*RPON(I, K)*RPOCALZ(I, K)/(RPOC(I, &
              & K)+1.0E-30)) / (PRALZ(I, K)+1.0E-30)
            ELSE
              ANCPLZ = 0.0
            END IF

               NRRATE = ((ANCPLZ-ANCLZ*ELZ(I, K)*(1.0-RFLZ(I, &
              & K)))*RLZ(I, K)+(BMLZ(I, K)+MLZ(I, K))*ANCLZ) * LZ (I, &
              & K) + PRLZ (I, K) * ANCLZ
            IF(PRALZ(I, K) /= 0.0) THEN
               FRL = ULLZ * LPOCALZ (I, K) / (PRALZ(I, K)+1.0E-30)
               FRR = URLZ * RPOCALZ (I, K) / (PRALZ(I, K)+1.0E-30)
            ELSE
               FRL = 0.0
               FRR = 0.0
            END IF

               NH4LZ (I, K) = NRRATE * FNH4LZ
               LDONLZ (I, K) = NRRATE * FLDONLZ
               RDONLZ (I, K) = NRRATE * FRDONLZ
               LPONLZ (I, K) = - FRL * LPON (I, K) / (LPOC(I, &
              & K)+1.0E-30) * RLZ (I, K) * LZ (I, K) + NRRATE * FLPONLZ
               RPONLZ (I, K) = - FRR * RPON (I, K) / (RPOC(I, &
              & K)+1.0E-30) * RLZ (I, K) * LZ (I, K) + NRRATE * FRPONLZ
          !
            End Do
         End Do
    !--------------------------------------------------------------------------------------
    ! EFFECT OF ZOOPLANKTON ON PHOSPHORUS (GM P/M**3/DAY)
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
               APCPSZ = (Q1(I, K)*UB1SZ*B1ASZ(I, K)+Q2(I, &
              & K)*UB2SZ*B2ASZ(I, K)+Q3(I, K)*UB3SZ*B3ASZ(I, &
              & K)+ULPSZ*LPOP(I, K)*LPOCASZ(I, K)/(LPOC(I, &
              & K)+1.0E-30)+URPSZ*RPOP(I, K)*RPOCASZ(I, K)/(RPOC(I, &
              & K)+1.0E-30)) / (PRASZ(I, K)+1.0E-30)
               PRRATE = ((APCPSZ-APCSZ*ESZ(I, K)*(1.-RFSZ(I, &
              & K)))*RSZ(I, K)+(BMSZ(I, K)+MSZ(I, K))*APCSZ) * SZ (I, &
              & K) + PRSZ (I, K) * APCSZ
               FRL = ULPSZ * LPOCASZ (I, K) / (PRASZ(I, K)+1.0E-30)
               FRR = URPSZ * RPOCASZ (I, K) / (PRASZ(I, K)+1.0E-30)
               PO4SZ (I, K) = PRRATE * FPO4SZ
               LDOPSZ (I, K) = PRRATE * FLDOPSZ
               RDOPSZ (I, K) = PRRATE * FRDOPSZ
               LPOPSZ (I, K) = - FRL * LPOP (I, K) / (LPOC(I, &
              & K)+1.0E-30) * RSZ (I, K) * SZ (I, K) + PRRATE * FLPOPSZ
               RPOPSZ (I, K) = - FRR * RPOP (I, K) / (RPOC(I, &
              & K)+1.0E-30) * RSZ (I, K) * SZ (I, K) + PRRATE * FRPOPSZ
               PIB1SZ (I, K) = Q1 (I, K) * B1SZ (I, K)
               PIB2SZ (I, K) = Q2 (I, K) * B2SZ (I, K)
               PIB3SZ (I, K) = Q3 (I, K) * B3SZ (I, K)
          !
               APCPLZ = (Q1(I, K)*UB1LZ*B1ALZ(I, K)+Q2(I, &
              & K)*UB2LZ*B2ALZ(I, K)+Q3(I, K)*UB3LZ*B3ALZ(I, &
              & K)+USZLZ*APCSZ*SZALZ(I, K)+ULLZ*LPOP(I, K)*LPOCALZ(I, &
              & K)/(LPOC(I, K)+1.0E-30)+URLZ*RPOP(I, K)*RPOCALZ(I, &
              & K)/(RPOC(I, K)+1.0E-30)) / (PRALZ(I, K)+1.0E-30)
               PRRATE = ((APCPLZ-APCLZ*ELZ(I, K)*(1.-RFLZ(I, &
              & K)))*RLZ(I, K)+(BMLZ(I, K)+MLZ(I, K))*APCLZ) * LZ (I, &
              & K) + PRLZ (I, K) * APCLZ
               FRL = ULLZ * LPOCALZ (I, K) / (PRALZ(I, K)+1.0E-30)
               FRR = URLZ * RPOCALZ (I, K) / (PRALZ(I, K)+1.0E-30)
               PO4LZ (I, K) = PRRATE * FPO4LZ
               LDOPLZ (I, K) = PRRATE * FLDOPLZ
               RDOPLZ (I, K) = PRRATE * FRDOPLZ
               LPOPLZ (I, K) = - FRL * LPOP (I, K) / (LPOC(I, &
              & K)+1.0E-30) * RLZ (I, K) * LZ (I, K) + PRRATE * FLPOPLZ
               RPOPLZ (I, K) = - FRR * RPOP (I, K) / (RPOC(I, &
              & K)+1.0E-30) * RLZ (I, K) * LZ (I, K) + PRRATE * FRPOPLZ
               PIB1LZ (I, K) = Q1 (I, K) * B1LZ (I, K)
               PIB2LZ (I, K) = Q2 (I, K) * B2LZ (I, K)
               PIB3LZ (I, K) = Q3 (I, K) * B3LZ (I, K)
          !
            End Do
         End Do
    !-----------------------------------------------------------------------------------------
    ! EFFECT OF ZOOPLANKTON ON DISSOLVED OXYGEN (GM DO/M**3/DAY)
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
               DOSZ (I, K) = (ESZ(I, K)*RFSZ(I, K)*RSZ(I, K)+BMSZ(I, &
              & K)) * AOCRSZ * SZ (I, K)
               DOLZ (I, K) = (ELZ(I, K)*RFLZ(I, K)*RLZ(I, K)+BMLZ(I, &
              & K)) * AOCRLZ * LZ (I, K)
          !
            End Do
         End Do

    !----------------------------------------------------------------------------------------
    ! EFFECT OF ZOOPLANKTON ON SILICA (GM SI/M**3/DAY)
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
               FRB1SZ = UB1SZ * B1ASZ (I, K) / (PRASZ(I, K)+1.0E-30)
               FRB2SZ = UB2SZ * B2ASZ (I, K) / (PRASZ(I, K)+1.0E-30)
               FRB3SZ = UB3SZ * B3ASZ (I, K) / (PRASZ(I, K)+1.0E-30)
               SASZ (I, K) = (FRB1SZ*ASC1+FRB2SZ*ASC2+FRB3SZ*ASC3) * &
              & RSZ (I, K) * SZ (I, K) * FRSASZ
               SUSZ (I, K) = (FRB1SZ*ASC1+FRB2SZ*ASC2+FRB3SZ*ASC3) * &
              & RSZ (I, K) * SZ (I, K) * (1.-FRSASZ)
          !
               FRB1LZ = UB1LZ * B1ALZ (I, K) / (PRALZ(I, K)+1.0E-30)
               FRB2LZ = UB2LZ * B2ALZ (I, K) / (PRALZ(I, K)+1.0E-30)
               FRB3LZ = UB3LZ * B3ALZ (I, K) / (PRALZ(I, K)+1.0E-30)
               SALZ (I, K) = (FRB1LZ*ASC1+FRB2LZ*ASC2+FRB3LZ*ASC3) * &
              & RLZ (I, K) * LZ (I, K) * FRSALZ
               SULZ (I, K) = (FRB1LZ*ASC1+FRB2LZ*ASC2+FRB3LZ*ASC3) * &
              & RLZ (I, K) * LZ (I, K) * (1.-FRSALZ)
          !
            End Do
         End Do

!--------------------------------------------------------------------------------------
         Return
  End Subroutine ZOOPL
  !
  !
  !************************************************************************
  !**                 S U B R O U T I N E   C A R B O N                  **
  !************************************************************************
  !
      Subroutine CARBON (DTLDOC, DTRDOC, DTLPOC, DTRPOC, FLXSPOC)
    !
         Use MOD_LIMS, Only: MLOC, KBM1, MTLOC
    !
         Use MOD_WQM, Only: ALGDOC, ALGPOC, B1, B2, BM1, BM2, DOXG, &
        & FCLD1, FCLD2, FCLDP, FCLP1, FCLP2, FCLPP, FCRD1, FCRD2, &
        & FCRDP, FCRP1, FCRP2, FCRPP, FDOP, KDCALG, KHR1, KHR2, KLDC, &
        & P1, P2, PR1, PR2, PRSP1, PRSP2, AANOX, DENIT, FTHDR, FTMNL, &
        & HDRLPOC, HDRRPOC, KHCOAG, KHNDN, KHODOC, KLCALG, KLPC, &
        & KRCOAG, KRDC, KRPC, KTHDR, KTMNL, LDOC, LPOC, MNLLDOC, &
        & MNLRDOC, NO3, RATOX, RDOC, RPOC, SALT, T, TRHDR, TRMNL, B, &
        & BENDOC, DLALGC, SAV_CALC, SAV_LOADS, SEDIMENT_CALC, V2, WSL, &
        & WSLNET, WSR, WSRNET
    !
         Use MOD_HYDROVARS, Only: DZ,DZ2D,D
    !
         Use MOD_SAV, Only: NSAVCELL, SAVCELL, LDOCSAVW, RDOCSAVW, &
        & LPOCSAVW, RPOCSAVW, LDOCEPIW, RDOCEPIW, LPOCEPIW, RPOCEPIW
    !
         Implicit None
         Real (SP) :: LPOC1, LPOC2, LPOC3, LPSETL, KLDOC, KLPOC, LDOC1, &
        & LDOC2, LDOC3
         Integer :: I, K
         Real (SP) :: DOCBM1, DOCBM2, DOCBM3, CP1, CP2, CP3, RDOC1, &
        & RDOC2, RDOC3, RPOC1, RPOC2, RPOC3, ALGCAR, SALTC, COAG, FDO1, &
        & FDO2, FDO3
         Real (SP) :: DICBM, DICPR, RPSETL, FTN
         Real (SP), Dimension (0:MTLOC, KBM1) :: DTLDOC, DTRDOC, &
        & DTLPOC, DTRPOC, FLXSPOC
    !
         DTLDOC = 0.0
         DTRDOC = 0.0
         DTLPOC = 0.0
         DTRPOC = 0.0
         FLXSPOC = 0.0
    !Note: *** modified to take out all algae 3 calculations ***
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
          !********* Algal contribution
          !
               DOCBM1 = KHR1 / (KHR1+DOXG(I, K))
               DOCBM2 = KHR2 / (KHR2+DOXG(I, K))
          !       DOCBM3 = KHR3/(KHR3+DOXG(I,K))
               CP1 = (P1(I, K)*PRSP1+BM1(I, K)) * B1 (I, K)
               CP2 = (P2(I, K)*PRSP2+BM2(I, K)) * B2 (I, K)
          !       CP3 = (P3(I,K)*PRSP3+BM3(I,K))*B3(I,K)
               LDOC1 = FCLD1 * CP1 + (1.-FCLD1-FCRD1-FCLP1-FCRP1) * &
              & DOCBM1 * CP1 + (FCLDP+FDOP*DOCBM1) * PR1 (I, K)
               LDOC2 = FCLD2 * CP2 + (1.-FCLD2-FCRD2-FCLP2-FCRP2) * &
              & DOCBM2 * CP2 + (FCLDP+FDOP*DOCBM2) * PR2 (I, K)
          !       LDOC3 = FCLD3*CP3+(1.-FCLD3-FCRD3-FCLP3-FCRP3)*DOCBM3*CP3+ (FCLDP+FDOP*DOCBM3)*PR3(I,K)
               RDOC1 = FCRD1 * CP1 + FCRDP * PR1 (I, K)
               RDOC2 = FCRD2 * CP2 + FCRDP * PR2 (I, K)
          !       RDOC3 = FCRD3*CP3+FCRDP*PR3(I,K)
               LPOC1 = FCLP1 * CP1 + FCLPP * PR1 (I, K)
               LPOC2 = FCLP2 * CP2 + FCLPP * PR2 (I, K)
          !       LPOC3      = FCLP3*CP3+FCLPP*PR3(I,K)
               RPOC1 = FCRP1 * CP1 + FCRPP * PR1 (I, K)
               RPOC2 = FCRP2 * CP2 + FCRPP * PR2 (I, K)
          !       RPOC3      = FCRP3*CP3+FCRPP*PR3(I,K)
          !       ALGDOC(I,K)  = LDOC1+LDOC2+LDOC3+RDOC1+RDOC2+RDOC3
          !       ALGPOC(I,K)  = LPOC1+LPOC2+LPOC3+RPOC1+RPOC2+RPOC3
               ALGDOC (I, K) = LDOC1 + LDOC2 + RDOC1 + RDOC2
               ALGPOC (I, K) = LPOC1 + LPOC2 + RPOC1 + RPOC2
          !********* Mineralization and hydrolysis
          !
               ALGCAR = B1 (I, K) + B2 (I, K)!+B3(I,K)
               KLDOC = KLDC (I, K) + KDCALG (I, K) * ALGCAR
               KLPOC = KLPC (I, K) + KLCALG (I, K) * ALGCAR
               SALTC = Max (0., SALT(I, K))
               COAG = 0.5 * (1.+Tanh(SALTC-KHCOAG)) * KRCOAG (I, K) * &
              & RDOC (I, K)
               FTMNL (I, K) = Exp (KTMNL*(T(I, K)-TRMNL))
               FTHDR (I, K) = Exp (KTHDR*(T(I, K)-TRHDR))
               RATOX (I, K) = (DOXG(I, K)+AANOX*KHODOC) / &
              & (KHODOC+DOXG(I, K))
               DENIT (I, K) = KLDOC * FTMNL (I, K) * AANOX * KHODOC / &
              & (KHODOC+DOXG(I, K)) * NO3 (I, K) / (KHNDN+NO3(I, K)) * &
              & LDOC (I, K)
               MNLLDOC (I, K) = KLDOC * FTMNL (I, K) * DOXG (I, K) / &
              & (KHODOC+DOXG(I, K)) * LDOC (I, K)
               MNLRDOC (I, K) = KRDC (I, K) * FTMNL (I, K) * DOXG (I, &
              & K) / (KHODOC+DOXG(I, K)) * RDOC (I, K)
               HDRLPOC (I, K) = RATOX (I, K) * KLPOC * FTHDR (I, K) * &
              & LPOC (I, K)
               HDRRPOC (I, K) = RATOX (I, K) * KRPC (I, K) * FTHDR (I, &
              & K) * RPOC (I, K)
          !
          !********* Change in carbon species

               DTLDOC (I, K) = (LDOC1+LDOC2-MNLLDOC(I, K)-DENIT(I, &
              & K)+HDRLPOC(I, K)+HDRRPOC(I, K)+LDOCSZ(I, K)+LDOCLZ(I, &
              & K)) / 86400.
               DTRDOC (I, K) = (RDOC1+RDOC2-MNLRDOC(I, K)+RDOCSZ(I, &
              & K)+RDOCLZ(I, K)-COAG) / 86400.
               DTLPOC (I, K) = (LPOC1+LPOC2-HDRLPOC(I, K)+LPOCSZ(I, &
              & K)+LPOCLZ(I, K)) / 86400.
               DTRPOC (I, K) = (RPOC1+RPOC2-HDRRPOC(I, K)+RPOCSZ(I, &
              & K)+RPOCLZ(I, K)+COAG) / 86400.
          !
            End Do
         End Do
    !
    ! NET EFFECT OF ALGAE ON CARBON FOR MASS BALANCE
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
               FDO1 = DOXG (I, K) / (KHR1+DOXG(I, K))
               FDO2 = DOXG (I, K) / (KHR2+DOXG(I, K))

               DICBM = (P1(I, K)*PRSP1+BM1(I, K)) * &
              & (1.-FCLD1-FCRD1-FCLP1-FCRP1) * FDO1 * B1 (I, K) + &
              & (P2(I, K)*PRSP2+BM2(I, K)) * &
              & (1.-FCLD2-FCRD2-FCLP2-FCRP2) * FDO2 * B2 (I, K)!&

               DICPR = (PR1(I, K)+PR2(I, K)) * FDOP
               DLALGC (I, K) = P1 (I, K) * B1 (I, K) + P2 (I, K) * B2 &
              & (I, K) - DICBM - DICPR
            End Do
         End Do
    !
         Do I = 1, MLOC
       !
       !********* Settling
       !
            LPSETL = - WSL (I, 1) * LPOC (I, 1) / (D(I)*DZ2D(I,1))
            RPSETL = - WSR (I, 1) * RPOC (I, 1) / (D(I)*DZ2D(I,1))
            DTLPOC (I, 1) = DTLPOC (I, 1) + LPSETL / 86400.
            DTRPOC (I, 1) = DTRPOC (I, 1) + RPSETL / 86400.
       !
       !********* Settling flux
       !
            FLXSPOC (I, 1) = (WSL(I, 1)*LPOC(I, 1)+WSR(I, 1)*RPOC(I, &
           & 1)) * V2 (I, 1) / (D(I)*DZ2D(I,1)*86400.)
       !
         End Do
    !
         Do K = 2, KBM1
            Do I = 1, MLOC
          !
          !********* Settling
          !
               LPSETL = (WSL(I, K-1)*LPOC(I, K-1)-WSL(I, K)*LPOC(I, K)) &
              & / (D(I)*DZ2D(I,K))
               RPSETL = (WSR(I, K-1)*RPOC(I, K-1)-WSR(I, K)*RPOC(I, K)) &
              & / (D(I)*DZ2D(I,K))
               DTLPOC (I, K) = DTLPOC (I, K) + LPSETL / 86400.
               DTRPOC (I, K) = DTRPOC (I, K) + RPSETL / 86400.
          !
          !********* Settling flux
          !
               FLXSPOC (I, K) = (WSL(I, K)*LPOC(I, K)+WSR(I, K)*RPOC(I, &
              & K)) * V2 (I, K) / (D(I)*DZ2D(I,K)*86400.)
          !
            End Do
         End Do
    !
    !******* Benthic fluxes
    !
         Do I = 1, MLOC
            DTLDOC (I, KBM1) = DTLDOC (I, KBM1) + BENDOC (I) / &
           & (D(I)*DZ2D(I,KBM1)) / 86400.
         End Do
    !
    !******* Resuspension
    !
         If (SEDIMENT_CALC) Then
       !
            Do I = 1, MLOC
               DTLPOC (I, KBM1) = DTLPOC (I, KBM1) + (WSL(I, &
              & KBM1)-WSLNET(I)) * LPOC (I, KBM1) / (D(I)*DZ2D(I,KBM1)) / &
              & 86400.
               DTRPOC (I, KBM1) = DTRPOC (I, KBM1) + (WSR(I, &
              & KBM1)-WSRNET(I)) * RPOC (I, KBM1) / (D(I)*DZ2D(I,KBM1)) / &
              & 86400.
            End Do
         End If
    !
    !******* SAV and epiphytes
    !
         If (SAV_CALC) Then
            Do I = 1, NSAVCELL
               B = SAVCELL (I)

!------------------ temporary modification to accomadate SAV module --------------
              DO K=1,6
                    DTLDOC (B, K) = DTLDOC (B, K) + &
                   & ((LDOCSAVW(B)+LDOCEPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
                   & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

                    DTRDOC (B, K) = DTRDOC (B, K) + &
                   & ((RDOCSAVW(B)+RDOCEPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
                   & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

                    DTLPOC (B, K) = DTLPOC (B, K) + &
                   & ((LPOCSAVW(B)+LPOCEPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
                   & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

                    DTRPOC (B, K) = DTRPOC (B, K) + &
                   & ((RPOCSAVW(B)+RPOCEPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
                   & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))
!-----------------------------------------------------------------------------------------
              END DO
            End Do
         End If

         If (SAV_LOADS) Then
            Do I = 1, NSAVCELL
               B = SAVCELL (I)
               DTLDOC (B, KBM1) = DTLDOC (B, KBM1) + LDOCSAVW (B) / &
              & (D(B)*DZ2D(B,KBM1)) / 86400.
               DTRDOC (B, KBM1) = DTRDOC (B, KBM1) + RDOCSAVW (B) / &
              & (D(B)*DZ2D(B,KBM1)) / 86400.
               DTLPOC (B, KBM1) = DTLPOC (B, KBM1) + LPOCSAVW (B) / &
              & (D(B)*DZ2D(B,KBM1)) / 86400.
               DTRPOC (B, KBM1) = DTRPOC (B, KBM1) + RPOCSAVW (B) / &
              & (D(B)*DZ2D(B,KBM1)) / 86400.
            End Do
         End If
         Return
      End Subroutine CARBON
  !
  !************************************************************************
  !**                 S U B R O U T I N E   N I T R O G                  **
  !************************************************************************
  !
      Subroutine NITROG (DTNH4, DTNO3, DTLDON, DTRDON, DTLPON, DTRPON, &
     & FLXSPON)
    !
         Use MOD_LIMS, Only: MLOC, KBM1, MTLOC
    !
         Use MOD_WQM, Only: ANC1, ANC2, B1, B2, BM1, BM2, DOXG, FNI1, &
        & FNIP, KHN1, KHN2, KHNNT, KHONT, KTNT1, KTNT2, NFIX, NH4, NO3, &
        & NT, NTM, P1, P1NNF, P2, PN1, PN2, PRSP1, PRSP2, T, TMNT, &
        & ATMOS_LOADS, ALGDON, ALGNH4, ALGNO3, ALGPON, B3, FNI2, FNLD1, &
        & FNLD2, FNLDP, FNLP1, FNLP2, FNLPP, FNRD1, FNRD2, FNRDP, &
        & FNRP1, FNRP2, FNRPP, KDNALG, KHCOAG, KHNAVG, KLDN, KLNALG, &
        & KLPN, KRCOAG, PR1, PR2, RDON, SALT, ANDC, ATMNH4, B, BENDON, &
        & BENNH4, BENNO3, DENIT, DENNO3, FTHDR, FTMNL, HDRLPON, &
        & HDRRPON, KRDN, KRPN, LDON, LPON, MNLLDON, MNLRDON, PRECIP, &
        & RATOX, RPON, SAV_CALC, SEDIMENT_CALC, V2, WSL, WSLNET, WSR, &
        & WSRNET, ATMLDON, ATMNO3, ATMRDON
    !
         Use MOD_HYDROVARS, Only: DZ, DZ2D, D
         Use MOD_SAV, Only: NSAVCELL, SAVCELL, NH4SAVW, NO3SAVW, &
        & LDONSAVW, RDONSAVW, LPONSAVW, RPONSAVW, NH4EPIW, NO3EPIW, &
        & LDONEPIW, RDONEPIW, LPONEPIW, RPONEPIW
    !
         Implicit None
         Real (SP) :: KLDON, KLPON, LPON1, LPON2, LPON3, NH4A1, NH4A2, &
        & NH4A3, NO3A1, NO3A2, NO3A3, LPSETL, NP1, NP2, NP3, LDON1, &
        & LDON2, LDON3, NL1MOD
         Real (SP) :: FTN, RDON1, RDON2, RDON3, RPON1, RPON2, RPON3
         Real (SP) :: ALGCAR, SALTC, COAG, RPSETL
         Integer :: I, K
         Real (SP), Dimension (0:MTLOC, KBM1) :: DTNH4, DTNO3, DTLDON, &
        & DTRDON, DTLPON, DTRPON, FLXSPON !all in gN/m^3/sec
    !
         DTLDON = 0.0
         DTRDON = 0.0
         DTLPON = 0.0
         DTRPON = 0.0
         FLXSPON = 0.0
    !
    ! NITRIFICATION
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
          !********* Temperature effect
          !
               If (T(I, K) < TMNT) Then
                  FTN = Exp (-KTNT1*(T(I, K)-TMNT)**2)
               Else
                  FTN = Exp (-KTNT2*(TMNT-T(I, K))**2)
               End If
          !
               NT (I, K) = DOXG (I, K) / (KHONT+DOXG(I, K)) * NH4 (I, &
              & K) / (KHNNT+NH4(I, K)) * FTN * NTM (I, K)
          !
            End Do
         End Do
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
          !********* Algal nitrogen preference
               PN1 (I, K) = NH4 (I, K) * (NO3(I, K)/((KHN1+NH4(I, &
              & K))*(KHN1+NO3(I, K)))+KHN1/((1.E-30+NH4(I, K)+NO3(I, &
              & K))*(KHN1+NO3(I, K))))
               PN2 (I, K) = NH4 (I, K) * (NO3(I, K)/((KHN2+NH4(I, &
              & K))*(KHN2+NO3(I, K)))+KHN2/((1.E-30+NH4(I, K)+NO3(I, &
              & K))*(KHN2+NO3(I, K))))

          !********* Nitrogen Fixation
          !
               NFIX (I, K) = (P1(I, K)-P1NNF(I, K)) * ANC1 * B1 (I, K)
          !
          !********* Algal sources/sinks
          !
               NP1 = ANC1 * (P1(I, K)*PRSP1+BM1(I, K)) * B1 (I, K)
               NP2 = ANC2 * (P2(I, K)*PRSP2+BM2(I, K)) * B2 (I, K)
          !       NP3 = ANC3*(P3(I,K)*PRSP3+BM3(I,K))*B3(I,K)
               NH4A1 = FNI1 * NP1 - PN1 (I, K) * P1NNF (I, K) * ANC1 * &
              & B1 (I, K) + FNIP * PR1 (I, K) * ANC1
               NH4A2 = FNI2 * NP2 - PN2 (I, K) * P2 (I, K) * ANC2 * B2 &
              & (I, K) + FNIP * PR2 (I, K) * ANC2
          !       NH4A3 = FNI3*NP3-PN3(I,K)*P3(I,K)*ANC3*B3(I,K)+FNIP*PR3(I,K)*ANC3
               NO3A1 = (PN1(I, K)-1.) * P1NNF (I, K) * ANC1 * B1 (I, K)
               NO3A2 = (PN2(I, K)-1.) * P2 (I, K) * ANC2 * B2 (I, K)
          !       NO3A3      = (PN3(I,K)-1.)*P3(I,K)*ANC3*B3(I,K)
               LDON1 = FNLD1 * NP1 + FNLDP * PR1 (I, K) * ANC1
               LDON2 = FNLD2 * NP2 + FNLDP * PR2 (I, K) * ANC2
          !       LDON3      = FNLD3*NP3+FNLDP*PR3(I,K)*ANC3
               RDON1 = FNRD1 * NP1 + FNRDP * PR1 (I, K) * ANC1
               RDON2 = FNRD2 * NP2 + FNRDP * PR2 (I, K) * ANC2
          !       RDON3      = FNRD3*NP3+FNRDP*PR3(I,K)*ANC3
               LPON1 = FNLP1 * NP1 + FNLPP * PR1 (I, K) * ANC1
               LPON2 = FNLP2 * NP2 + FNLPP * PR2 (I, K) * ANC2
          !       LPON3     = FNLP3*NP3+FNLPP*PR3(I,K)*ANC3
               RPON1 = FNRP1 * NP1 + FNRPP * PR1 (I, K) * ANC1
               RPON2 = FNRP2 * NP2 + FNRPP * PR2 (I, K) * ANC2
          !       RPON3     = FNRP3*NP3+FNRPP*PR3(I,K)*ANC3
               ALGNH4 (I, K) = NH4A1 + NH4A2 !+NH4A3
               ALGNO3 (I, K) = NO3A1 + NO3A2 !+NO3A3
          !       ALGDON(I,K) = LDON1+LDON2+LDON3+RDON1+RDON2+RDON3
          !       ALGPON(I,K) = LPON1+LPON2+LPON3+RPON1+RPON2+RPON3
               ALGDON (I, K) = LDON1 + LDON2 + RDON1 + RDON2
               ALGPON (I, K) = LPON1 + LPON2 + RPON1 + RPON2
          !********* Mineralization and hydrolysis
          !
          !     ALGCAR = B1 (I, K) + B2 (I, K) + B3 (I, K)
			   ALGCAR = B1 (I, K) + B2 (I, K)! + B3 (I, K) !AN
               KLDON = KLDN (I, K) + KDNALG (I, K) * ALGCAR * KHNAVG / &
              & (KHNAVG+NH4(I, K)+NO3(I, K))
               KLPON = KLPN (I, K) + KLNALG (I, K) * ALGCAR * KHNAVG / &
              & (KHNAVG+NH4(I, K)+NO3(I, K))
               SALTC = Max (0., SALT(I, K))
               COAG = 0.5 * (1.+Tanh(SALTC-KHCOAG)) * KRCOAG (I, K) * &
              & RDON (I, K)
               MNLLDON (I, K) = RATOX (I, K) * KLDON * FTMNL (I, K) * &
              & LDON (I, K)
               MNLRDON (I, K) = RATOX (I, K) * KRDN (I, K) * FTMNL (I, &
              & K) * RDON (I, K)
               HDRLPON (I, K) = RATOX (I, K) * KLPON * FTHDR (I, K) * &
              & LPON (I, K)
               HDRRPON (I, K) = RATOX (I, K) * KRPN (I, K) * FTHDR (I, &
              & K) * RPON (I, K)
               DENNO3 (I, K) = - ANDC * DENIT (I, K)
          !
          !********* Change in nitrogen species

               DTNH4 (I, K) = (NH4A1+NH4A2+MNLLDON(I, K)+MNLRDON(I, &
              & K)-NT(I, K)+NH4SZ(I, K)+NH4LZ(I, K)) / 86400.
               DTNO3 (I, K) = (NT(I, K)-ANDC*DENIT(I, K)+NO3A1+NO3A2) / &
              & 86400.
               DTLDON (I, K) = (LDON1+LDON2-MNLLDON(I, K)+HDRLPON(I, &
              & K)+HDRRPON(I, K)+LDONSZ(I, K)+LDONLZ(I, K)) / 86400.
               DTRDON (I, K) = (RDON1+RDON2-MNLRDON(I, K)+RDONSZ(I, &
              & K)+RDONLZ(I, K)-COAG) / 86400.
               DTLPON (I, K) = (LPON1+LPON2-HDRLPON(I, K)+LPONSZ(I, &
              & K)+LPONLZ(I, K)) / 86400.
               DTRPON (I, K) = (RPON1+RPON2-HDRRPON(I, K)+RPONSZ(I, &
              & K)+RPONLZ(I, K)+COAG) / 86400.
          !
            End Do
         End Do
    !
         Do I = 1, MLOC
       !
       !********* Settling flux
       !
            FLXSPON (I, 1) = (WSL(I, 1)*LPON(I, 1)+WSR(I, 1)*RPON(I, &
           & 1)) * V2 (I, 1) / (D(I)*DZ2D(I,1)*86400.)
       !
       !********* Settling
       !
            LPSETL = - WSL (I, 1) * LPON (I, 1) / (D(I)*DZ2D(I,1))
            RPSETL = - WSR (I, 1) * RPON (I, 1) / (D(I)*DZ2D(I,1))
            DTLPON (I, 1) = DTLPON (I, 1) + LPSETL / 86400.
            DTRPON (I, 1) = DTRPON (I, 1) + RPSETL / 86400.
       !
         End Do
    !
         Do K = 2, KBM1
            Do I = 1, MLOC
          !
          !********* Settling flux
          !
               FLXSPON (I, K) = (WSL(I, K)*LPON(I, K)+WSR(I, K)*RPON(I, &
              & K)) * V2 (I, K) / (D(I)*DZ2D(I,K)*86400.)
          !
          !********* Settling
          !
               LPSETL = (WSL(I, K-1)*LPON(I, K-1)-WSL(I, K)*LPON(I, K)) &
              & / (D(I)*DZ2D(I,K))
               RPSETL = (WSR(I, K-1)*RPON(I, K-1)-WSR(I, K)*RPON(I, K)) &
              & / (D(I)*DZ2D(I,K))
               DTLPON (I, K) = DTLPON (I, K) + LPSETL / 86400.
               DTRPON (I, K) = DTRPON (I, K) + RPSETL / 86400.
          !
            End Do
         End Do
    !
    !******* Benthic fluxes
    !
         Do I = 1, MLOC
            DTNH4 (I, KBM1) = DTNH4 (I, KBM1) + BENNH4 (I) / &
           & (D(I)*DZ2D(I,KBM1)) / 86400. 								!gN/m^3/sec
            DTNO3 (I, KBM1) = DTNO3 (I, KBM1) + BENNO3 (I) / &
           & (D(I)*DZ2D(I,KBM1)) / 86400.
            DTLDON (I, KBM1) = DTLDON (I, KBM1) + BENDON (I) / &
           & (D(I)*DZ2D(I,KBM1)) / 86400.
         End Do
    !
    !******* Resuspension
    !
         If (SEDIMENT_CALC) Then
       !
            Do I = 1, MLOC
               DTLPON (I, KBM1) = DTLPON (I, KBM1) + (WSL(I, &
              & KBM1)-WSLNET(I)) * LPON (I, KBM1) / (D(I)*DZ2D(I,KBM1)) / &
              & 86400.
               DTRPON (I, KBM1) = DTRPON (I, KBM1) + (WSR(I, &
              & KBM1)-WSRNET(I)) * RPON (I, KBM1) / (D(I)*DZ2D(I,KBM1)) / &
              & 86400.
            End Do
         End If
    !
    !******* SAV and epiphytes
    !
         If (SAV_CALC) Then
            DO I = 1, NSAVCELL
               B = SAVCELL (I)

            DO K=1,6

               DTNH4 (B, K) = DTNH4 (B, K) + &
               & ((NH4SAVW(B)+NH4EPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
               & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

                DTNO3 (B, K) = DTNO3 (B, K) + &
               & ((NO3SAVW(B)+NO3EPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
               & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

                DTLDON (B, K) = DTLDON (B, K) + &
               & ((LDONSAVW(B)+LDONEPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
               & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

                DTRDON (B, K) = DTRDON (B, K) + &
               & ((RDONSAVW(B)+RDONEPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
               & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

                DTLPON (B, K) = DTLPON (B, K) + &
               & ((LPONSAVW(B)+LPONEPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
               & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

                DTRPON (B, K) = DTRPON (B, K) + &
               & ((RPONSAVW(B)+RPONEPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
               & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

            END DO
!---------------------------------------------------------------------------------
           END DO
         End If
    !
    !******* Atmospheric Loads
    !
         If (ATMOS_LOADS) Then ! 
            Do I = 1, MLOC
               DTNH4 (I, 1) = DTNH4 (I, 1) + PRECIP * ATMNH4 / &
              & (D(I)*DZ2D(I,1))
               DTNO3 (I, 1) = DTNO3 (I, 1) + PRECIP * ATMNO3 / &
              & (D(I)*DZ2D(I,1))
               DTLDON (I, 1) = DTLDON (I, 1) + PRECIP * ATMLDON / &
              & (D(I)*DZ2D(I,1))
               DTRDON (I, 1) = DTRDON (I, 1) + PRECIP * ATMRDON / &
              & (D(I)*DZ2D(I,1))
            End Do
         End If
    !
         Return
      End Subroutine NITROG
  !
  !************************************************************************
  !**                 S U B R O U T I N E   P H O S P H                  **
  !************************************************************************
  !
      Subroutine PHOSPH (DTPO4, DTLDOP, DTRDOP, DTLPOP, DTRPOP, DTPIP, &
     & DTPIB1, DTPIB2, DTPIB3, FLXSPO4, FLXSPOP, FLXSPIP)
    !
         Use MOD_LIMS, Only: MLOC, KBM1, MTLOC
    !
         Use MOD_WQM, Only: B1, B2, BM1, BM2, FPI1, FPI2, FPIP, FPLD1, &
        & FPLD2, FPLDP, FPLP1, FPLPP, FPRD1, FPRD2, FPRDP, KHP1, KHP2, &
        & P1, P2, PO4, PR1, PR2, PRSP1, PRSP2, Q01, Q02, Q1, Q2, VMAX1, &
        & VMAX2, ALGDOP, ALGPO4, ALGPOP, B3, ATMOS_LOADS, FPLP2, FPRP1, &
        & FPRP2, FPRPP, FTHDR, FTMNL, HDRLPOP, HDRRPOP, KADPO4, KDPALG, &
        & KHPAVG, KLDP, KLPALG, KLPP, KRDP, KRPP, LDOP, LPOP, MNLLDOP, &
        & MNLRDOP, RATOX, RDOP, RPOP, SSI, WSL, WSS, ATMLDOP, ATMPO4, &
        & ATMRDOP, B, BENDOP, BENPO4, PIB1, PIB2, PIB3, PIP, PIP_CALC, &
        & PRECIP, SAV_CALC, SEDIMENT_CALC, SOLIDS_CALC, V2, WS1, &
        & WS1NET, WS2, WS2NET, WS3, WSLNET, WSR, WSRNET, WSSHI, WSSNET
    !
         Use MOD_HYDROVARS, Only: DZ, DZ2D, D
         Use MOD_SAV, Only: NSAVCELL, SAVCELL, PO4SAVW, LDOPSAVW, &
        & RDOPSAVW, LPOPSAVW, RPOPSAVW, PO4EPIW, LDOPEPIW, RDOPEPIW, &
        & LPOPEPIW, RPOPEPIW
    !
         Implicit None
         Real (SP) PO41, PO42
         Real (SP) :: KLDOP, KLPOP, LPOP1, LPOP2, LPOP3, LPSETL, LDOP1, &
        & LDOP2, LDOP3
         Real (SP) :: DF, PO4AVL, PUP1, PUP2, PUP3, PP1, PP2, PP3, &
        & PO43, RDOP1, RDOP2, RDOP3, RPOP1, RPOP2, RPOP3
         Real (SP) :: ALGCAR, PFU, PPU, PFD, PPD, PO4SET, RPSETL, &
        & PIB1SETL, PIB2SETL, PIB3SETL, PF, PP
         Integer :: I, J, K
         Real (SP), Dimension (0:MTLOC, KBM1) :: DTPO4, DTLDOP, DTRDOP, &
        & DTLPOP, DTRPOP, DTPIP, DTPIB1, DTPIB2, DTPIB3, FLXSPO4, &
        & FLXSPOP, FLXSPIP
    !
         DTPO4 = 0.0
         DTLDOP = 0.0
         DTRDOP = 0.0
         DTLPOP = 0.0
         DTRPOP = 0.0
         DTPIP = 0.0
         DTPIB1 = 0.0
         DTPIB2 = 0.0
         DTPIB3 = 0.0
         FLXSPO4 = 0.0
         FLXSPOP = 0.0
         FLXSPIP = 0.0
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
          !********* Algae sources/sinks
          !
          !********* Uptake via Monod formulation
          !
               DF = 1. !1./(1.+KADPO4*SSI(I,K))
               PO4AVL = Max (0., DF*PO4(I, K))
               PUP1 = VMAX1 * B1 (I, K) * PO4AVL / (KHP1+PO4AVL+1.0E-6)
               PUP2 = VMAX2 * B2 (I, K) * PO4AVL / (KHP2+PO4AVL+1.0E-6)
          !       PUP3   = VMAX3*B3(I,K)*PO4AVL/(KHP3+PO4AVL+1.0E-6)
          !
               PP1 = Q01 * (P1(I, K)*PRSP1+BM1(I, K)) * B1 (I, K)
               PP2 = Q02 * (P2(I, K)*PRSP2+BM2(I, K)) * B2 (I, K)
          !       PP3  = Q03*(P3(I,K)*PRSP3+BM3(I,K))*B3(I,K)
               PO41 = FPI1 * PP1 - PUP1 + FPIP * PR1 (I, K) * Q1 (I, K)
               PO42 = FPI2 * PP2 - PUP2 + FPIP * PR2 (I, K) * Q2 (I, K)
          !       PO43 = FPI3*PP3-PUP3+FPIP*PR3(I,K)*Q3(I,K)
               LDOP1 = FPLD1 * PP1 + FPLDP * PR1 (I, K) * Q1 (I, K)
               LDOP2 = FPLD2 * PP2 + FPLDP * PR2 (I, K) * Q2 (I, K)
          !       LDOP3     = FPLD3*PP3+FPLDP*PR3(I,K)*Q3(I,K)
               RDOP1 = FPRD1 * PP1 + FPRDP * PR1 (I, K) * Q1 (I, K)
               RDOP2 = FPRD2 * PP2 + FPRDP * PR2 (I, K) * Q2 (I, K)
          !       RDOP3     = FPRD3*PP3+FPRDP*PR3(I,K)*Q3(I,K)
               LPOP1 = FPLP1 * PP1 + FPLPP * PR1 (I, K) * Q1 (I, K)
               LPOP2 = FPLP2 * PP2 + FPLPP * PR2 (I, K) * Q2 (I, K)
          !       LPOP3     = FPLP3*PP3+FPLPP*PR3(I,K)*Q3(I,K)
               RPOP1 = FPRP1 * PP1 + FPRPP * PR1 (I, K) * Q1 (I, K)
               RPOP2 = FPRP2 * PP2 + FPRPP * PR2 (I, K) * Q2 (I, K)
          !       RPOP3     = FPRP3*PP3+FPRPP*PR3(I,K)*Q3(I,K)
               ALGPO4 (I, K) = PO41 + PO42 !+PO43
               ALGDOP (I, K) = LDOP1 + LDOP2 + RDOP1 + RDOP2
               ALGPOP (I, K) = LPOP1 + LPOP2 + RPOP1 + RPOP2
          !       ALGDOP(I,K) = LDOP1+LDOP2+LDOP3+RDOP1+RDOP2+RDOP3
          !       ALGPOP(I,K) = LPOP1+LPOP2+LPOP3+RPOP1+RPOP2+RPOP3
          !
          !
          !********* Mineralization and hydrolysis
          !
          !     ALGCAR = B1 (I, K) + B2 (I, K) + B3 (I, K)
			   ALGCAR = B1 (I, K) + B2 (I, K) ! + B3 (I, K) !AN Commented out
               DF = 1. !1./(1.+KADPO4*SSI(I,K))
               PO4AVL = Max (1.E-6, DF*PO4(I, K))
               KLDOP = KLDP (I, K) + KDPALG (I, K) * ALGCAR * KHPAVG / &
              & (KHPAVG+PO4AVL)
               KLPOP = KLPP (I, K) + KLPALG (I, K) * ALGCAR * KHPAVG / &
              & (KHPAVG+PO4AVL)
               MNLLDOP (I, K) = RATOX (I, K) * KLDOP * FTMNL (I, K) * &
              & LDOP (I, K)
               MNLRDOP (I, K) = RATOX (I, K) * KRDP (I, K) * FTMNL (I, &
              & K) * RDOP (I, K)
               HDRLPOP (I, K) = RATOX (I, K) * KLPOP * FTHDR (I, K) * &
              & LPOP (I, K)
               HDRRPOP (I, K) = RATOX (I, K) * KRPP (I, K) * FTHDR (I, &
              & K) * RPOP (I, K)
          !
          !********* Change in phosphorus species

               DTPO4 (I, K) = (PO41+PO42+MNLLDOP(I, K)+MNLRDOP(I, &
              & K)+PO4SZ(I, K)+PO4LZ(I, K)) / 86400.
               DTLDOP (I, K) = (LDOP1+LDOP2-MNLLDOP(I, K)+HDRLPOP(I, &
              & K)+LDOPSZ(I, K)+LDOPLZ(I, K)) / 86400.
               DTRDOP (I, K) = (RDOP1+RDOP2-MNLRDOP(I, K)+HDRRPOP(I, &
              & K)+RDOPSZ(I, K)+RDOPLZ(I, K)) / 86400.
               DTLPOP (I, K) = (LPOP1+LPOP2-HDRLPOP(I, K)+LPOPSZ(I, &
              & K)+LPOPLZ(I, K)) / 86400.
               DTRPOP (I, K) = (RPOP1+RPOP2-HDRRPOP(I, K)+RPOPSZ(I, &
              & K)+RPOPLZ(I, K)) / 86400.
               DTPIB1 (I, K) = (PUP1-PP1-PR1(I, K)*Q1(I, K)-PIB1SZ(I, &
              & K)-PIB1LZ(I, K)) / 86400.
               DTPIB2 (I, K) = (PUP2-PP2-PR2(I, K)*Q2(I, K)-PIB2SZ(I, &
              & K)-PIB2LZ(I, K)) / 86400.
          !      DTPIB3(I,K) = (PUP3-PP3-PR3(I,K)*Q3(I,K)-PIB3SZ(I,K)-PIB3LZ(I,K))/86400.
          !
            End Do
         End Do
    !
         Do I = 1, MLOC
       !
       !********* Settling
       !
            PFD = KADPO4 * SSI (I, 1) / (1.+KADPO4*SSI(I, 1))
            PPD = PFD * PO4 (I, 1)
            PO4SET = - WSS (I, 1) * PPD / (D(I)*DZ2D(I,1))
            LPSETL = - WSL (I, 1) * LPOP (I, 1) / (D(I)*DZ2D(I,1))
            RPSETL = - WSR (I, 1) * RPOP (I, 1) / (D(I)*DZ2D(I,1))
            PIB1SETL = - WS1 (I, 1) * PIB1 (I, 1) / (D(I)*DZ2D(I,1))
            PIB2SETL = - WS2 (I, 1) * PIB2 (I, 1) / (D(I)*DZ2D(I,1))
       !     PIB3SETL = -WS3(I,1)*PIB3(I,1)/(D(I)*DZ2D(I,1))
            DTPO4 (I, 1) = DTPO4 (I, 1) + PO4SET / 86400.
            DTLPOP (I, 1) = DTLPOP (I, 1) + LPSETL / 86400.
            DTRPOP (I, 1) = DTRPOP (I, 1) + RPSETL / 86400.
            DTPIB1 (I, 1) = DTPIB1 (I, 1) + PIB1SETL / 86400.
            DTPIB2 (I, 1) = DTPIB2 (I, 1) + PIB2SETL / 86400.
       !     DTPIB3(I,1) = DTPIB3(I,1)+PIB3SETL/86400.
       !
       !********* Settling flux
       !
            FLXSPO4 (I, 1) = WSS (I, 1) * PPD * V2 (I, 1) / &
           & (D(I)*DZ2D(I,1)*86400.)

            FLXSPOP (I, 1) = (WSL(I, 1)*LPOP(I, 1)+WSR(I, 1)*RPOP(I, &
           & 1)+WS1(I, 1)*PIB1(I, 1)+WS2(I, 1)*PIB2(I, 1)) * V2 (I, 1) &
           & / (D(I)*DZ2D(I,1)*86400.)
         End Do
    !
         Do K = 2, KBM1
            Do I = 1, MLOC
          !
          !********* Settling below top depth
          !
               PFU = KADPO4 * SSI (I, K-1) / (1.+KADPO4*SSI(I, K-1))!Particle fraction of P in upperward grid (K-1)
               PPU = PFU * PO4 (I, K-1)!Particulate Inorganic P (PO4 in particles) (gP/m^3) at K-1
               PFD = KADPO4 * SSI (I, K) / (1.+KADPO4*SSI(I, K))!Particle fraction of P in downward grid(K)
               PPD = PFD * PO4 (I, K)!Particulate Inorganic P (PO4 in particles) (gP/m^3) at K
          !
               PO4SET = (WSS(I, K-1)*PPU-WSS(I, K)*PPD) / (D(I)*DZ2D(I,K))
          !
               LPSETL = (WSL(I, K-1)*LPOP(I, K-1)-WSL(I, K)*LPOP(I, K)) &
              & / (D(I)*DZ2D(I,K))
               RPSETL = (WSR(I, K-1)*RPOP(I, K-1)-WSR(I, K)*RPOP(I, K)) &
              & / (D(I)*DZ2D(I,K))
               PIB1SETL = (WS1(I, K-1)*PIB1(I, K-1)-WS1(I, K)*PIB1(I, &
              & K)) / (D(I)*DZ2D(I,K))
               PIB2SETL = (WS2(I, K-1)*PIB2(I, K-1)-WS2(I, K)*PIB2(I, &
              & K)) / (D(I)*DZ2D(I,K))
               PIB3SETL = (WS3(I, K-1)*PIB3(I, K-1)-WS3(I, K)*PIB3(I, &
              & K)) / (D(I)*DZ2D(I,K))
               DTPO4 (I, K) = DTPO4 (I, K) + PO4SET / 86400.
               DTLPOP (I, K) = DTLPOP (I, K) + LPSETL / 86400.
               DTRPOP (I, K) = DTRPOP (I, K) + RPSETL / 86400.
               DTPIB1 (I, K) = DTPIB1 (I, K) + PIB1SETL / 86400.
               DTPIB2 (I, K) = DTPIB2 (I, K) + PIB2SETL / 86400.
               DTPIB3 (I, K) = DTPIB3 (I, K) + PIB3SETL / 86400.
          !
          !********* Settling flux!
          !
               FLXSPO4 (I, K) = WSS (I, K) * PPD * V2 (I, K) / &
              & (D(I)*DZ2D(I,K)*86400.)
               FLXSPOP (I, K) = (WSL(I, K)*LPOP(I, K)+WSR(I, K)*RPOP(I, &
              & K)+WS1(I, K)*PIB1(I, K)+WS2(I, K)*PIB2(I, K)+WS3(I, &
              & K)*PIB3(I, K)) * V2 (I, K) / (D(I)*DZ2D(I,K)*86400.)
          !
            End Do
         End Do
    !
    !******* Benthic fluxes
    !
         Do I = 1, MLOC
            DTPO4 (I, KBM1) = DTPO4 (I, KBM1) + BENPO4 (I) / &
           & (D(I)*DZ2D(I,KBM1)) / 86400.
            DTLDOP (I, KBM1) = DTLDOP (I, KBM1) + BENDOP (I) / &
           & (D(I)*DZ2D(I,KBM1)) / 86400.
         End Do
    !
    !******* Resuspension
    !
         If (SEDIMENT_CALC) Then
       !
            Do I = 1, MLOC
               PF = KADPO4 * SSI (I, KBM1) / (1.+KADPO4*SSI(I, KBM1))
               PP = PF * PO4 (I, KBM1)
               DTPO4 (I, KBM1) = DTPO4 (I, KBM1) + (WSS(I, &
              & KBM1)-WSSNET(I)) * PP / (D(I)*DZ2D(I,KBM1)) / 86400.
               DTLPOP (I, KBM1) = DTLPOP (I, KBM1) + (WSL(I, &
              & KBM1)-WSLNET(I)) * LPOP (I, KBM1) / (D(I)*DZ2D(I,KBM1)) / &
              & 86400.
               DTRPOP (I, KBM1) = DTRPOP (I, KBM1) + (WSR(I, &
              & KBM1)-WSRNET(I)) * RPOP (I, KBM1) / (D(I)*DZ2D(I,KBM1)) / &
              & 86400.
               DTPIB1 (I, KBM1) = DTPIB1 (I, KBM1) + (WS1(I, &
              & KBM1)-WS1NET(I)) * PIB1 (I, KBM1) / (D(I)*DZ2D(I,KBM1)) / &
              & 86400.
               DTPIB2 (I, KBM1) = DTPIB2 (I, KBM1) + (WS2(I, &
              & KBM1)-WS2NET(I)) * PIB2 (I, KBM1) / (D(I)*DZ2D(I,KBM1)) / &
              & 86400.
          !       DTPIB3(I,KBM1) = DTPIB3(I,KBM1)+(WS3(I,KBM1)-WS3NET(I))*      &
          !                        PIB3(I,KBM1)/(D(I)*DZ2D(I,KBM1))/86400.
            End Do
         End If
    !
    !******* SAV and epiphytes
    !
         If (SAV_CALC) Then
            Do I = 1, NSAVCELL
               B = SAVCELL (I)
!----------------------------------------------------------------------------------
            DO K=1,6
                  DTPO4 (B, K) = DTPO4 (B, K) + &
                 & ((PO4SAVW(B)+PO4EPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
                 & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

                  DTLDOP (B, K) = DTLDOP (B, K) + &
                 & ((LDOPSAVW(B)+LDOPEPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
                 & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

                  DTRDOP (B, K) = DTRDOP (B, K) + &
                 & ((RDOPSAVW(B)+RDOPEPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
                 & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

                  DTLPOP (B, K) = DTLPOP (B, K) + &
                 & ((LPOPSAVW(B)+LPOPEPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
                 & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))

                  DTRPOP (B, K) = DTRPOP (B, K) + &
                 & ((RPOPSAVW(B)+RPOPEPIW(B)) / (D(B)*DZ2D(B,K)) / 86400.)* &
                 & (DZ2D(B,K)/SUM(DZ2D(B,1:6)))
            END DO

!----------------------------------------------------------------------------------
            End Do
         End If
    !
    !******* Atmospheric loads
    !
         If (ATMOS_LOADS) Then ! 
            Do I = 1, MLOC
               DTPO4 (I, 1) = DTPO4 (I, 1) + PRECIP * ATMPO4 / &
              & (D(I)*DZ2D(I,1))
               DTLDOP (I, 1) = DTLDOP (I, 1) + PRECIP * ATMLDOP / &
              & (D(I)*DZ2D(I,1))
               DTRDOP (I, 1) = DTRDOP (I, 1) + PRECIP * ATMRDOP / &
              & (D(I)*DZ2D(I,1))
            End Do
         End If
    !
    !******* Particulate Inorganic Phosphorus
    !
         If (SOLIDS_CALC .And. PIP_CALC) Then
       !
       !******  Settling
            Do I = 1, MLOC
               FLXSPIP (I, 1) = WSSHI (I, 1) * PIP (I, 1) * V2 (I, 1) / &
              & (D(I)*DZ2D(I,1)*86400.)!
               DTPIP (I, 1) = - WSSHI (I, 1) * PIP (I, 1) / &
              & (D(I)*DZ2D(I,1)*86400.)
            End Do
       !
            Do K = 2, KBM1
               Do I = 1, MLOC
                  FLXSPIP (I, K) = WSSHI (I, K) * PIP (I, K) * V2 (I, &
                 & K) / (D(I)*DZ2D(I,K)*86400.)
                  DTPIP (I, K) = (WSSHI(I, K-1)*PIP(I, K-1)-WSSHI(I, &
                 & K)*PIP(I, K)) / (D(I)*DZ2D(I,K)*86400.)
               End Do
            End Do
       !
       !******  Resuspension
       !
            If (SEDIMENT_CALC) Then
               Do I = 1, MLOC
                  DTPIP (I, KBM1) = DTPIP (I, KBM1) + (WSSHI(I, &
                 & KBM1)-WSSNET(I)) * PIP (I, KBM1) / (D(I)*DZ2D(I,KBM1)) / &
                 & 86400.
               End Do
            End If
       !
         End If
    !
         Return
      End Subroutine PHOSPH
  !
  !
  !************************************************************************
  !**                S U B R O U T I N E   C O D M N D                   **
  !************************************************************************
  !
      Subroutine CODMND (DTCOD)
    !
         Use MOD_LIMS, Only: MLOC, MTLOC, KBM1
    !
         Use MOD_WQM, Only: COD, DOXG, FTCOD, KCOD, KHOCOD, KTCOD, T, &
        & TRCOD
         Use MOD_HYDROVARS, Only: DZ,DZ2D,D
         Implicit None
         Integer :: I, K
         Real (SP) :: DTCOD (0:MTLOC, KBM1)
    !
         DTCOD = 0.0
    !
    !******* Change in chemical oxygen demand
    !
         Do K = 1, KBM1
            Do I = 1, MLOC

               FTCOD (I, K) = KCOD (I, K) * Exp (KTCOD*(T(I, K)-TRCOD))
               DTCOD (I, K) = (-DOXG(I, K)/(KHOCOD+DOXG(I, K))*FTCOD(I, &
              & K)*COD(I, K)) / 86400.
               DTCOD (I, K) = - FTCOD (I, K) * COD (I, K) / 86400.
            End Do
         End Do
    !
    !******* Sediment demand
         Do I = 1, MLOC
            DTCOD (I, KBM1) = DTCOD (I, KBM1)
         End Do
    !
         Return
      End Subroutine CODMND
  !
  !************************************************************************
  !**                 S U B R O U T I N E   O X Y G E N                  **
  !************************************************************************
  !
      Subroutine OXYGEN (DTDO)
    !
         Use MOD_LIMS, Only: MLOC, MTLOC, KBM1
    !
         Use MOD_WQM, Only: AOCR, AONT, B1, B2, BM1, BM2, DCOD, DOPR, &
        & DOXG, FCLD1, FCLD2, FCLP1, FCLP2, FCRD1, FCRD2, FCRP1, FCRP2, &
        & FDOP, KHR1, KHR2, NITRIF, NT, P1, P2, PN1, PN2, PR1, PR2, &
        & PRSP1, PRSP2, AREAR, B, BENDO, BREAR, COD, CREAR, DDOC, &
        & DORALG, DOS, FTCOD, KHOCOD, MNLLDOC, MNLRDOC, REAERDO, RESP, &
        & SALT, SAV_CALC, SAV_LOADS, T, WMS, JDAY
    !
         Use MOD_HYDROVARS, Only: DZ,DZ2D, D
         Use MOD_SAV, Only: NSAVCELL, SAVCELL, DOSAVW, DOEPIW, SEDDOSAV
    !
         Use MOD_WQMINIT, Only: COD_CALC
!
         Use MOD_CO2SYS, Only: AIRSEA_OPTION, & 
        & JDstart1
!
         Implicit None
         Real (SP) :: KRDO
         Integer :: I, K
         Real (SP) :: FRDO1, FRDO2, FRDO3, CP1, CP2, CP3, ALGDO, DOR1, &
        & DOR2, DOR3, DOP1, DOP2, DOP3
         Real (SP) :: FACTOR, CHLORN, TDOS, RNU, Sc, AA, TS
         Real (SP) :: DTDO (0:MTLOC, KBM1)
    !
         DTDO = 0.0
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
          !********* Nitrification
          !
               NITRIF (I, K) = AONT * NT (I, K)
          !
          !********* Uptake/respiration
          !
               FRDO1 = (1.-FCLD1-FCRD1-FCLP1-FCRP1) * DOXG (I, K) / &
              & (KHR1+DOXG(I, K))
               FRDO2 = (1.-FCLD2-FCRD2-FCLP2-FCRP2) * DOXG (I, K) / &
              & (KHR2+DOXG(I, K))
          !       FRDO3  = (1.-FCLD3-FCRD3-FCLP3-FCRP3)*DOXG(I,K)/(KHR3+DOXG(I,K))
               CP1 = P1 (I, K) * PRSP1 + BM1 (I, K)
               CP2 = P2 (I, K) * PRSP2 + BM2 (I, K)
          !       CP3    = P3(I,K)*PRSP3+BM3(I,K)
          !       ALGDO  = (FRDO1*CP1*B1(I,K)+FRDO2*CP2*B2(I,K)                       &
          !                +FRDO3*CP3*B3(I,K))*AOCR
               ALGDO = (FRDO1*CP1*B1(I, K)+FRDO2*CP2*B2(I, K)) * AOCR
               DOR1 = ((1.3-0.3*PN1(I, K))*P1(I, K)-FRDO1*CP1) * AOCR * &
              & B1 (I, K)
               DOR2 = ((1.3-0.3*PN2(I, K))*P2(I, K)-FRDO2*CP2) * AOCR * &
              & B2 (I, K)
          !       DOR3   = ((1.3-0.3*PN3(I,K))*P3(I,K)-FRDO3*CP3)*AOCR*B3(I,K)
               DOP1 = FDOP * PR1 (I, K) * AOCR * DOXG (I, K) / &
              & (KHR1+DOXG(I, K))
               DOP2 = FDOP * PR2 (I, K) * AOCR * DOXG (I, K) / &
              & (KHR2+DOXG(I, K))
          !       DOP3   = FDOP*PR3(I,K)*AOCR*DOXG(I,K)/(KHR3+DOXG(I,K))
               DOPR (I, K) = DOP1 + DOP2 !+DOP3
          !
          !LB: debug DOXG
               If (DOXG(I, K) == 0.0) Then
                  Write (*,*) 'I,K,DOXG=', I, K, DOXG (I, K)
                  Write (*,*) 'DCOD(I-1,K)=', DCOD (I-1, K)
                  Write (*,*) 'FTCOD(I,K),COD(I,K)=', FTCOD (I, K), COD &
                 & (I, K)
               End If
          !
               If (COD_CALC) Then
                  DCOD (I, K) = DOXG (I, K) / (KHOCOD+DOXG(I, K)) * &
                 & FTCOD (I, K) * COD (I, K)
               End If
          !
               DDOC (I, K) = AOCR * (MNLLDOC(I, K)+MNLRDOC(I, K))
               DORALG (I, K) = DOR1 + DOR2 !+DOR3
          !
          !********* Change in dissolved oxygen
          !
          !     DTDO(I,K)  = (DOR1+DOR2+DOR3-DOP1-DOP2-DOP3-DDOC(I,K)-DCOD(I,K) &
          !        &-NITRIF(I,K)-DOSZ(I,K)-DOLZ(I,K))/86400.
		       DTDO (I, K)  = (DOR1 + DOR2 - DOP1 - DOP2 - DDOC(I, K) &
              &- NITRIF(I, K) - DOSZ(I, K) - DOLZ(I, K))/86400. 
			   RESP (I, K) = ALGDO + DOPR (I, K) + DDOC (I, K) + NITRIF &
              & (I, K) + DOSZ (I, K) + DOLZ (I, K)
               If (COD_CALC) Then
                  DTDO (I, K) = DTDO (I, K) - DCOD (I, K) / 86400.
                  RESP (I, K) = RESP (I, K) + DCOD (I, K)
               End If
          !
            End Do
         End Do
    !
    !
    !******* Reaeration
!
         If (AIRSEA_OPTION == 1) Then
            FACTOR = AREAR * (BREAR*WMS) ** CREAR
            Do I = 1, MLOC
               CHLORN = SALT (I, 1) / 1.80655
               TDOS = T (I, 1)
               RNU = 0.54 + 0.7 * TDOS / 30 - 0.07 * SALT (I, 1) / 35.
          !     KRDO      = 0.157*RNU*(1.5*WMS)**1.5
               KRDO = FACTOR * RNU
               DOS = 14.5532 + TDOS * (0.0054258*TDOS-0.38217) - CHLORN &
              & * (0.1665+TDOS*(9.796E-5*TDOS-5.866E-3))
               DTDO (I, 1) = DTDO (I, 1) + KRDO / (D(I)*DZ2D(I,1)) * &
              & (DOS-DOXG(I, 1)) / 86400.

               REAERDO (I, 1) = KRDO / (D(I)*DZ2D(I,1)) * (DOS-DOXG(I, 1))!gO2/m^3/day rate of areation in surface layer
          !Stumm and Morgan (third edition) pp242
            End Do
         !Else
		 Else If(AIRSEA_OPTION>1 .AND. AIRSEA_OPTION<=5) Then
!
            If (AIRSEA_OPTION == 2) Then 		!Wanninkhof (1992)
               FACTOR = 0.31 * WMS ** 2.0 		!0.31 for winds in steady state;  0.39 for long-term averaged winds
            Else If (AIRSEA_OPTION == 3) Then 	!Nightingale (2000)
               FACTOR = (0.33*WMS+0.222*WMS**2.0)
            Else If (AIRSEA_OPTION == 4) Then 	!Ho (2006)
               FACTOR = 0.266 * WMS ** 2.0
            Else If (AIRSEA_OPTION == 5) Then 	!Wanninkhof et al (2013): wind in equation is <U^2>,
												!meaning the time mean of U^2
               FACTOR = 0.251 * WMS ** 2.0      !over the time interval of flux determination
												!(they compute monthly fluxes with 6hourly wind data)
			   !FACTOR = AREAR * (BREAR*WMS) ** CREAR
            End If
!
            Do I = 1, MLOC
!
               Sc = 1953.4 - 128.0 * T (I, 1) + 3.9918 * T (I, 1) ** &
              & 2.0 - 0.050091 * T (I, 1) ** 3.0 !Schmidt number for O2 (W92)

               KRDO = FACTOR * Sqrt (660.0/Sc)!gas transfer velocity (W12) in cm/hr
               KRDO = KRDO * 24.0 / 100.0 !m/day

               TS = Log ((298.15-T(I, 1))/(273.15+T(I, 1)))
               AA = 2.00907 + TS * (3.22014+TS*(4.05010+TS*(4.94457+&
              & TS*(-0.256847+TS*3.88767)))) + SALT (I, 1) * (-&
              & 0.00624523+TS*(-0.00737614+TS*(-0.0103410+TS*-&
              & 0.00817083))) + (-0.000000488682) * SALT (I, 1) * SALT &
              & (I, 1)

		  !  Converting from ml/l to mg/l (1mg/l=22.391/31.998 ml/l ~ 0.7ml/l)
               DOS = Exp (AA) * 1.429
!
               DTDO (I, 1) = DTDO (I, 1) + KRDO / (D(I)*DZ(1)) * &
              & (DOS-DOXG(I, 1)) / 86400.
!
          !debugging areation
               REAERDO (I, 1) = KRDO / (D(I)*DZ(1)) * (DOS-DOXG(I, 1))!gO2/m^3/day rate of areation in surface layer
          !Stumm and Morgan (third edition) pp242
            End Do
		 Else If(AIRSEA_OPTION == 6) Then ! Combination of Option 1 method and Wanninkhof et al (2013)
		    FACTOR = AREAR * (BREAR*WMS) ** CREAR
!
			Do I = 1, MLOC
!
		    If (JDAY <= JDstart1(1) .OR. JDAY >= JDstart1(2)) Then
			   CHLORN = SALT (I, 1) / 1.80655
               TDOS = T (I, 1)
               RNU = 0.54 + 0.7 * TDOS / 30 - 0.07 * SALT (I, 1) / 35.

               KRDO = FACTOR * RNU
               DOS = 14.5532 + TDOS * (0.0054258*TDOS-0.38217) - CHLORN &
              & * (0.1665+TDOS*(9.796E-5*TDOS-5.866E-3))
!
			Else If (JDAY > JDstart1(1) .AND. JDAY < JDstart1(2)) Then

          !Schmidt number for O2 (W92)
               Sc = 1953.4 - 128.0 * T (I, 1) + 3.9918 * T (I, 1) ** &
              & 2.0 - 0.050091 * T (I, 1) ** 3.0

               KRDO = FACTOR * Sqrt (660.0/Sc)!gas transfer velocity (W12) in cm/hr
               KRDO = KRDO * 24.0 / 100.0 !m/day

               TS = Log ((298.15-T(I, 1))/(273.15+T(I, 1)))
               AA = 2.00907 + TS * (3.22014+TS*(4.05010+TS*(4.94457+&
              & TS*(-0.256847+TS*3.88767)))) + SALT (I, 1) * (-&
              & 0.00624523+TS*(-0.00737614+TS*(-0.0103410+TS*-&
              & 0.00817083))) + (-0.000000488682) * SALT (I, 1) * SALT &
              & (I, 1)
!
               DOS = Exp (AA) * 1.429 !  Converting from ml/l to mg/l (1mg/l=22.391/31.998 ml/l ~ 0.7ml/l)

			 End If
!
               DTDO (I, 1) = DTDO (I, 1) + KRDO / (D(I)*DZ2D(I,1)) * &
              & (DOS-DOXG(I, 1)) / 86400.
!
          !Debugging aeration

               REAERDO (I, 1) = KRDO / (D(I)*DZ2D(I,1)) * (DOS-DOXG(I, 1))!gO2/m^3/day rate of areation
																		  !in surface layer

          !Stumm and Morgan (third edition) pp242

            End Do
         End If
!
    !******* Sediment oxygen demand
    !
         Do I = 1, MLOC
            DTDO (I, KBM1) = DTDO (I, KBM1) + BENDO (I) / &
           & (D(I)*DZ2D(I,KBM1)) / 86400. !gO2/m^3/sec
         End Do
    !BENDO needs to be gO2/m2/day
    !
    !
    !******* SAV and epiphytes
    !
         If (SAV_CALC) Then
            Do I = 1, NSAVCELL
               B = SAVCELL (I)
               !DTDO (B, KBM1) = DTDO (B, KBM1) + &
              !& (DOSAVW(B)+DOEPIW(B)+SEDDOSAV(B)) / (D(B)*DZ2D(B,KBM1)) / &
            !  & 86400.

              DTDO (B, KBM1) = DTDO (B, KBM1) + &
             & (SEDDOSAV(B)) / (D(B)*DZ2D(B,KBM1)) / &
             & 86400.

              DO K=1,6
                DTDO (B, K) = DTDO (B, K) + ((DOSAVW(B)+DOEPIW(B))/(D(B)*DZ2D(B,K)) / &
                & 86400.)*(DZ2D(B,K)/SUM(DZ2D(B,1:6)))
              END DO
!------------------------------------------------------------------------------------

            End Do
         End If

         If (SAV_LOADS) Then
            Do I = 1, NSAVCELL
               B = SAVCELL (I)
               DTDO (B, KBM1) = DTDO (B, KBM1) + DOSAVW (B) / &
              & (D(B)*DZ2D(B,KBM1)) / 86400.
            End Do
         End If
    !
         Return
      End Subroutine OXYGEN
  !
  !************************************************************************
  !**                  S U B R O U T I N E   S I L I C A                 **
  !************************************************************************
  !
      Subroutine SILICA (DTSIUPB, DTSIAT, FLXSSI)
    !
         Use MOD_LIMS, Only: MLOC, MTLOC, KBM1
    !
         Use MOD_WQM, Only: ALGRES, ALGUP, ASC1, ASC2, ASC3, B1, B2, &
        & B3, BM1, BM2, BM3, FSAP, KADSA, KSUA, KTSUA, P1, P2, P3, PR1, &
        & PR2, PR3, PRSP1, PRSP2, PRSP3, PSD, SIUPB, SAP, T, TRSUA, &
        & BENSA, SIAT, SEDIMENT_CALC, SSI, V2, WSS, WSSNET, WSU, WSUNET
    !
         Use MOD_HYDROVARS, Only: DZ, D, DZ2D
    !
         Implicit None
         Integer :: I, K
         Real (SP) :: SUP, SP1, SP2, SP3, PFU, PF, SASETL, SUSETL
         Real (SP), Dimension (0:MTLOC, KBM1) :: DTSIUPB, DTSIAT, &
        & FLXSSI
    !
         DTSIUPB = 0.0
         DTSIAT = 0.0
         FLXSSI = 0.0
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
          !********* Particulate silica dissolution
          !
               PSD (I, K) = KSUA (I, K) * SIUPB (I, K) * Exp &
              & (KTSUA*(T(I, K)-TRSUA))
          !
          !********* Predation
          !
               SAP (I, K) = FSAP * (PR1(I, K)*ASC1+PR2(I, &
              & K)*ASC2+PR3(I, K)*ASC3)
               SUP = (1.-FSAP) * (PR1(I, K)*ASC1+PR2(I, K)*ASC2+PR3(I, &
              & K)*ASC3)
          !
          !********* Algal uptake and respiration
          !
          !     ALGUP (I, K) = P1 (I, K) * ASC1 * B1 (I, K) + P2 (I, K) &
          !    & * ASC2 * B2 (I, K) + P3 (I, K) * ASC3 * B3 (I, K)
			   ALGUP (I, K) = P1 (I, K) * ASC1 * B1 (I, K) + P2 (I, K) &
              & * ASC2 * B2 (I, K)!  + P3 (I, K) * ASC3 * B3 (I, K) !AN
          !
               SP1 = ASC1 * (P1(I, K)*PRSP1+BM1(I, K)) * B1 (I, K)
               SP2 = ASC2 * (P2(I, K)*PRSP2+BM2(I, K)) * B2 (I, K)
          !     SP3 = ASC3 * (P3(I, K)*PRSP3+BM3(I, K)) * B3 (I, K)
          !     ALGRES (I, K) = SP1 + SP2 + SP3
			   ALGRES (I, K) = SP1 + SP2 ! + SP3 !AN
          !
          !********* Change in silica
          !
               DTSIUPB (I, K) = (SUP-PSD(I, K)+ALGRES(I, K)+SUSZ(I, &
              & K)+SULZ(I, K)) / 86400.
               DTSIAT (I, K) = (SAP(I, K)+PSD(I, K)-ALGUP(I, K)+SASZ(I, &
              & K)+SALZ(I, K)) / 86400.
          !
            End Do
         End Do
    !
         Do I = 1, MLOC
       !
       !********* Settling
       !
            PF = KADSA * SSI (I, 1) / (1.+KADSA*SSI(I, 1))
            SASETL = - WSS (I, 1) * PF * SIAT (I, 1) / (D(I)*DZ2D(I,1))
            SUSETL = - WSU (I, 1) * SIUPB (I, 1) / (D(I)*DZ2D(I,1))
       !
       !********* Change in silica
       !
            DTSIUPB (I, 1) = DTSIUPB (I, 1) + SUSETL / 86400.
            DTSIAT (I, 1) = DTSIAT (I, 1) + SASETL / 86400.
       !
       !********* Settling fluxes
       !
       !
       !PF*SIAT is the concentration of particulate silica bound on SSI (avaiable for algae and zooplankton to graze)
       !but not part of algae and zooplankton yet
       !
            FLXSSI (I, 1) = (WSS(I, 1)*PF*SIAT(I, 1)+WSU(I, 1)*SIUPB(I, &
           & 1)) * V2 (I, 1) / (D(I)*DZ2D(I,1)*86400.)

         End Do
    !
         Do K = 2, KBM1
            Do I = 1, MLOC
          !
          !********* Settling (source term of layer k due to settling)
          !
               PFU = KADSA * SSI (I, K-1) / (1.+KADSA*SSI(I, K-1))
               PF = KADSA * SSI (I, K) / (1.+KADSA*SSI(I, K))
               SASETL = (WSS(I, K-1)*PFU*SIAT(I, K-1)-WSS(I, &
              & K)*PF*SIAT(I, K)) / (D(I)*DZ2D(I,K))
               SUSETL = (WSU(I, K-1)*SIUPB(I, K-1)-WSU(I, K)*SIUPB(I, &
              & K)) / (D(I)*DZ2D(I,K))
          !
          !********* Change in silica
          !
               DTSIUPB (I, K) = DTSIUPB (I, K) + SUSETL / 86400.
               DTSIAT (I, K) = DTSIAT (I, K) + SASETL / 86400.
          !
          !********* Settling fluxes
          !
               FLXSSI (I, K) = (WSS(I, K)*PF*SIAT(I, K)+WSU(I, &
              & K)*SIUPB(I, K)) * V2 (I, K) / (D(I)*DZ2D(I,K)*86400.)
          !
            End Do
         End Do
    !
    !******* Benthic fluxes
    !
         Do I = 1, MLOC
            DTSIAT (I, KBM1) = DTSIAT (I, KBM1) + BENSA (I) / &
           & (D(I)*DZ2D(I,KBM1)) / 86400.
         End Do
    !
    !******* Resuspension
    !
         If (SEDIMENT_CALC) Then
       !
            Do I = 1, MLOC
               PF = KADSA * SSI (I, KBM1) / (1.+KADSA*SSI(I, KBM1))
               DTSIAT (I, KBM1) = DTSIAT (I, KBM1) + (WSS(I, &
              & KBM1)-WSSNET(I)) * PF * SIAT (I, KBM1) / &
              & (D(I)*DZ2D(I,KBM1)) / 86400.
               DTSIUPB (I, KBM1) = DTSIUPB (I, KBM1) + (WSU(I, &
              & KBM1)-WSUNET(I)) * SIUPB (I, KBM1) / (D(I)*DZ2D(I,KBM1)) / &
              & 86400.
            End Do
         End If
    !
         Return
      End Subroutine SILICA
  !
  !************************************************************************
  !**               S U B R O U T I N E   B E N C O M P                  **
  !************************************************************************
  !
      Subroutine BEN_FLUX !simple bottom boundary flux calcualtin (parallel to sediment flux model)
    !modulate the input benthic flux to give some temperature control
    !
         Use MOD_LIMS, Only: MLOC, MTLOC, KBM1
    !
         Use MOD_WQM, Only: BENCOD, BENCODB, BENDO, BENDOB, BENNH4, &
        & BENNH4B, BENNO3, BENNO3B, BENPO4, BENPO4B, BENSA, DOXG, KHSO, &
        & KSDOC, KSNH4, KSNO3, KSO, KSPO4, KSSA, MTCNO3, NO3, SEDNO3, &
        & T, TRSDOC, TRSNH4, TRSNO3, TRSO, TRSPO4, TRSSA, B1, B2, &
        & BENDOC, BENDOCB, BENSAB, LPOC, LPON, LPOP, PCFWS, PNFWS, &
        & PPFWS, PSFWS, SSFWS, SIUPB, Q1, Q2, RPOC, RPON, RPOP, WS1, &
        & WS2, WSL, WSR, WSU, ANC1, ANC2, ASC1, ASC2, SIAT, KADPO4, &
        & KADSA, PO4, SSI, WSS
    !
         Use MOD_HYDROVARS, Only: DZ, D

# if defined (WET_DRY)
	USE  MOD_WD, ONLY: ISWETN
# endif
    !
         Implicit None
         Integer :: I, K
         Real (SP) :: FTSOD, FTNH4, FTNO3, FTPO4, FTSA, FTDOC
         Real (SP) :: PF, PPO4
    !
         Do I = 1, MLOC
       !
       !******* Calculate temperature effects
       !
            FTSOD = Exp (KSO*(T(I, KBM1)-TRSO))
            FTNH4 = Exp (KSNH4*(T(I, KBM1)-TRSNH4))
            FTNO3 = Exp (KSNO3*(T(I, KBM1)-TRSNO3))
            FTPO4 = Exp (KSPO4*(T(I, KBM1)-TRSPO4))
            FTSA = Exp (KSSA*(T(I, KBM1)-TRSSA))
            FTDOC = Exp (KSDOC*(T(I, KBM1)-TRSDOC))
       !
       !********* Dissolved oxygen effect on sediment oxygen demand
       !
            BENDO (I) = BENDOB (I) * FTSOD * DOXG (I, KBM1) / &
           & (KHSO+DOXG(I, KBM1))
            BENCOD (I) = BENCODB (I) - BENDOB (I) * FTSOD * KHSO / &
           & (KHSO+DOXG(I, KBM1))!
       !
       !********* Balance of nutrients (g/m2/day)
       !
            BENNH4 (I) = BENNH4B (I) * FTNH4
            BENNO3 (I) = BENNO3B (I) + FTNO3 * MTCNO3 * (SEDNO3-NO3(I, &
           & KBM1))
            BENPO4 (I) = BENPO4B (I) * FTPO4
            BENSA (I) = BENSAB (I) * FTSA
            BENDOC (I) = BENDOCB (I) * FTDOC
       !
         End Do

		 DO I=1,MLOC
#  if !defined (WET_DRY)
			IF(D(I) <= 0.0_SP) THEN
#  else
			IF(ISWETN(I) == 0)THEN
#  endif
				BENDO (I)  = 0.0_SP
				BENCOD (I) = 0.0_SP
				BENNH4 (I) = 0.0_SP
				BENNO3 (I) = 0.0_SP
				BENPO4 (I) = 0.0_SP
				BENSA (I)  = 0.0_SP
				BENDOC (I) = 0.0_SP
            ENDIF
		 ENDDO
    !
    !******* Compute particle flux to sediments for mass-balance purposes
    !
         Do I = 1, MLOC

       !
       !******* First zero flux accumlators
       !
            PCFWS (I) = 0.
            PNFWS (I) = 0.
            PPFWS (I) = 0.
            PSFWS (I) = 0.
            SSFWS (I) = 0.
       !
       !******* Now accumulate fluxes of labile and refractory particles as well as suspended solids
       !
!!!!!
            PPFWS (I) = PPFWS (I) - WSL (I, KBM1) * LPOP (I, KBM1) - &
           & WSR (I, KBM1) * RPOP (I, KBM1)
            PNFWS (I) = PNFWS (I) - WSL (I, KBM1) * LPON (I, KBM1) - &
           & WSR (I, KBM1) * RPON (I, KBM1)
            PCFWS (I) = PCFWS (I) - WSL (I, KBM1) * LPOC (I, KBM1) - &
           & WSR (I, KBM1) * RPOC (I, KBM1)
            PSFWS (I) = PSFWS (I) - WSU (I, KBM1) * SIUPB (I, KBM1)
            SSFWS (I) = SSFWS (I) - WSS (I, KBM1) * SSI (I, KBM1)
       !
       !******* Now accumulate fluxes of algal biomass

       !Fixed bug  in PCFWS - twp WS2's instead of WS2 and WS3
            PCFWS (I) = PCFWS (I) - WS1 (I, KBM1) * B1 (I, KBM1) - WS2 &
           & (I, KBM1) * B2 (I, KBM1)
       !
       ! below per tykim replaced ANC3 with ANC2 in second equation
       !
            PNFWS (I) = PNFWS (I) - WS1 (I, KBM1) * ANC1 * B1 (I, KBM1) &
           & - WS2 (I, KBM1) * ANC2 * B2 (I, KBM1)

            PPFWS (I) = PPFWS (I) - Q1 (I, KBM1) * WS1 (I, KBM1) * B1 &
           & (I, KBM1) - Q2 (I, KBM1) * WS2 (I, KBM1) * B2 (I, KBM1)

            PSFWS (I) = PSFWS (I) - WS1 (I, KBM1) * ASC1 * B1 (I, KBM1) &
           & - WS2 (I, KBM1) * ASC2 * B2 (I, KBM1)

       !******* Accumulate adsorbed phosphate and silica flux (things that move with suspended solids****
       !
            PF = KADPO4 * SSI (I, KBM1) / (1.+KADPO4*SSI(I, KBM1))
            PPO4 = PF * PO4 (I, KBM1)
            PPFWS (I) = PPFWS (I) - WSS (I, KBM1) * PPO4
            PF = KADSA * SSI (I, KBM1) / (1.+KADSA*SSI(I, KBM1))
            PSFWS (I) = PSFWS (I) - WSS (I, KBM1) * PF * SIAT (I, KBM1)
       !
         End Do

		DO I=1,MLOC
#  if !defined (WET_DRY)
			IF(D(I) <= 0.0_SP)THEN
#  else
			IF(ISWETN(I) == 0)THEN
#  endif
				BENDO (I) = 0.0_SP
				BENCOD (I) = 0.0_SP
				BENNH4 (I) = 0.0_SP
				BENNO3 (I) = 0.0_SP
				BENPO4 (I) = 0.0_SP
				BENSA (I) = 0.0_SP
				BENDOC (I) = 0.0_SP
            ENDIF
		 ENDDO

         Return
      End Subroutine BEN_FLUX
  !
End Module MOD_KIN
