!mod_algal.F
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
!** This is a new version of algal input and kinetics.  Created for    **
!** The Tributary Refinements portion of the Chesapeake Bay study.     **
!**                                                                    **
!** Refinements include creation of an algal input routine, a new      **
!** P vs. I curve, a silica fraction for all algal groups, option to   **
!** spatially vary inputs, lookup table for temperature functions,     **
!** switch to uE from langleys, option to specify minimum or product   **
!** function for growth limitations.                                   **
!**                                                                    **
!**                  Algal  Subroutines for CE-QUAL-ICM                **
!**                                                                    **
!**                  Chesapeake Bay Tributary Refinements              **
!**                         February 7, 1996                           **
!**                                                                    **
!**                    Water Quality Modeling Group                    **
!**                    U.S. Army Corps of Engineers                    **
!**                    Waterways Experiment Station                    **
!**                    Vicksburg, Mississippi 39180                    **
!**                                                                    **
!************************************************************************
!
Module MOD_ALGAL
  !
      Use MOD_PREC, Only: SP
      Use MOD_SIZES, Only: MGL
      Use MOD_LIMS, Only: KBM1, MTLOC, MLOC
  !
      Use MOD_WQM, Only: ACHLC1, ACHLC2, ACHLC3, ALOFN, ALPHMIN1, &
     & ALPHMIN2, ALPHMIN3, ALPHRAT1, ALPHRAT2, ALPHRAT3, ANC1, ANC2, &
     & ANC3, ASC1, ASC2, ASC3, ASRAT, B1, B2, B2GR, B3, BCHLC1, BCHLC2, &
     & BCHLC3, BM1, BM2, BMR1, BMR2, BMR3, BPR1, BPR2, BPR3, CCHL1, &
     & CCHL2, CCHL3, CCHLC1, CCHLC2, CCHLC3, CFIX, CHLCMN1, CHLCMN2, &
     & CHLCMN3, FCLD1, FCLD2, FCLD3, FCLDP, FCLP1, FCLP2, FCLP3, FCLPP, &
     & FCRD1, FCRD2, FCRD3, FCRDP, FCRP1, FCRP2, FCRP3, FCRPP, FDOP, &
     & FI1, FI2, FNFIX, FNI1, FNI2, FNI3, FNIP, FNLD1, FNLD2, FNLD3, &
     & FNLDP, FNLP1, FNLP2, FNLP3, FNLPP, FNRD1, FNRD2, FNRD3, FNRDP, &
     & FNRP1, FNRP2, FNRP3, FNRPP, FNUP, FPI1, FPI2, FPI3, FPIP, FPLD1, &
     & FPLD2, FPLD3, FPLDP, FPLP1, FPLP2, FPLP3, FPLPP, FPRD1, FPRD2, &
     & FPRD3, FPRDP, FPRP1, FPRP2, FPRP3, FPRPP, FSAP, GPP, IK1, IK2, &
     & IT, JDAY, KHN1, KHN2, KHN3, KHNAVG, KHNFIX, KHP1, KHP2, KHP3, &
     & KHPAVG, KHR1, KHR2, KHR3, KHS1, KHS2, KHS3, KHST1, KHST2, KHST3, &
     & MINPROD, NH4, NL1, NL2, NO3, NPP, P1, P1NNF, P2, PL1, PL2, PM1, &
     & PM2, PM3, PN1, PN2, PO4, PR1, PR2, PRSP1, PRSP2, PRSP3, Q01, &
     & Q02, Q03, SEDIMENT_CALC, SNFIX, STF1, STF2, STF3, T, &
     & TOTAL_NETPP, TVPR, V2, VMAX1, VMAX2, VMAX3, WS1, WS1NET, WS2, &
     & WS2NET, WS3, WS3NET
      Use MOD_FILEINFO, Only: AGR, ALO
  !
      Use MOD_HYDROVARS, Only: DZ, DZ2D, D
  !Took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: SERIAL, MSR, PAR
  !
      Use MOD_OWQ, Only: IAVG, KESS
      Use MOD_ZOOP, Only: B1SZ, B2SZ, B1LZ, B2LZ
      Implicit None
  !
  !LB: moved the definitions below to this location to have acces to TRPR, TMP1, TR1, TMP2, TR2, TMP3, TR3 withing get_f*** function
      Real (SP) :: KTG11, KTG12, KTB1, KTG21, KTG22, KTB2, KTG31
      Real (SP) :: KTG32, KTB3, KTPR
      Real (SP) :: TRPR, TMP1, TR1, TMP2, TR2, TMP3, TR3
!
  !   MOVE ALL VARIABLES OUT OF GLOBAL (HERE)
  !   AND DOWN TO WHERE THEY BELONG
  !
  !
Contains
  !
  !Subroutines:
  ! subroutine	ALG_READ()
  ! subroutine	ALGAE()
  !
  !
  !************************************************************************
  !**             S U B R O U T I N E   A L G _ R E A D                  **
  !************************************************************************
  !
      Subroutine ALG_READ
#  if defined (MULTIPROCESSOR)
         Use MOD_PAR, Only: NGID, NHN, HN_LST
    !
#  endif
         Implicit None
    !
         Integer :: I, J, K
         Real (SP), Allocatable :: RTMP21 (:, :), RTMP22 (:, :), RTMP23 &
        & (:, :)
         Character (Len=72) :: ALG_TITLE (6)
         Character (Len=8) :: SPVAR1, SPVAR2, SPVAR3
         Character (Len=8) :: PRINT1, PRINT2, PRINT3
         Character (Len=8) :: TPVAR, TPRINT, TB2GR, PRINTB2

         Read (AGR, 1010) (ALG_TITLE(J), J=1, 6)
    !
    ! READ PARAMETERS WHICH ARE IDENTICAL FOR ALL GROUPS
    !
         Read (AGR, 1040) MINPROD
         Read (AGR, 1030) TRPR, KTPR
         Read (AGR, 1030) FNIP, FNUP, FNLDP, FNRDP, FNLPP, FNRPP 
         Read (AGR, 1030) FPIP, FPLDP, FPRDP, FPLPP, FPRPP
         Read (AGR, 1030) FDOP, FCLDP, FCRDP, FCLPP, FCRPP
         Read (AGR, 1030) FSAP
    ! READ SPATIALLY-INVARIANT PARAMETERS FOR ALGAL GROUP 1
    !
         Read (AGR, 1030) ANC1, Q01, ASC1, STF1 !
         Read (AGR, 1030) CHLCMN1, ACHLC1, BCHLC1, CCHLC1 !CHL to C ratios etc
         Read (AGR, 1030) KHN1, KHP1, KHS1, KHR1, KHST1, KHNFIX, FNFIX !KH half constants
         Read (AGR, 1030) ALPHMIN1, ALPHRAT1, PRSP1, VMAX1
         Read (AGR, 1030) TMP1, TR1 !TMP1 is optimal temperature for grwoth
         Read (AGR, 1030) KTG11, KTG12, KTB1
         Read (AGR, 1030) FNI1, FNLD1, FNRD1, FNLP1, FNRP1 !

         Read (AGR, 1030) FPI1, FPLD1, FPRD1, FPLP1, FPRP1
         Read (AGR, 1030) FCLD1, FCRD1, FCLP1, FCRP1
    !
    ! READ SPATIALLY-INVARIANT PARAMETERS FOR ALGAL GROUP 2
    !
         Read (AGR, 1030) ANC2, Q02, ASC2, STF2
         Read (AGR, 1030) CHLCMN2, ACHLC2, BCHLC2, CCHLC2
         Read (AGR, 1030) KHN2, KHP2, KHS2, KHR2, KHST2
         Read (AGR, 1030) ALPHMIN2, ALPHRAT2, PRSP2, VMAX2
         Read (AGR, 1030) TMP2, TR2 !TMP2: optimal temperature for growth
         Read (AGR, 1030) KTG21, KTG22, KTB2
         Read (AGR, 1030) FNI2, FNLD2, FNRD2, FNLP2, FNRP2
         Read (AGR, 1030) FPI2, FPLD2, FPRD2, FPLP2, FPRP2
         Read (AGR, 1030) FCLD2, FCRD2, FCLP2, FCRP2
    !
    ! READ SPATIALLY-INVARIANT PARAMETERS FOR ALGAL GROUP 3
    !
         Read (AGR, 1030) ANC3, Q03, ASC3, STF3
         Read (AGR, 1030) CHLCMN3, ACHLC3, BCHLC3, CCHLC3
         Read (AGR, 1030) KHN3, KHP3, KHS3, KHR3, KHST3
         Read (AGR, 1030) ALPHMIN3, ALPHRAT3, PRSP3, VMAX3
         Read (AGR, 1030) TMP3, TR3
         Read (AGR, 1030) KTG31, KTG32, KTB3
         Read (AGR, 1030) FNI3, FNLD3, FNRD3, FNLP3, FNRP3
         Read (AGR, 1030) FPI3, FPLD3, FPRD3, FPLP3, FPRP3
         Read (AGR, 1030) FCLD3, FCRD3, FCLP3, FCRP3
    !
         KHNAVG = (KHN1+KHN2+KHN3) / 3.
         KHPAVG = (KHP1+KHP2+KHP3) / 3.
    !
    ! ARE REMAINING GROUP 1 PARAMETERS SPATIALLY VARYING?
    !
         Read (AGR, 1040) SPVAR1, PRINT1
         If (SPVAR1 == 'CONSTANT') Then
            Read (AGR, 1030) PM1 (1, 1), BMR1 (1, 1), BPR1 (1, 1)
            Do K = 1, KBM1
               Do I = 1, MTLOC
                  PM1 (I, K) = PM1 (1, 1)
                  BMR1 (I, K) = BMR1 (1, 1)
                  BPR1 (I, K) = BPR1 (1, 1)
               End Do
            End Do
         Else
            Allocate (RTMP21(MGL, KBM1))
            RTMP21 = 0.0
            Allocate (RTMP22(MGL, KBM1))
            RTMP22 = 0.0
            Allocate (RTMP23(MGL, KBM1))
            RTMP23 = 0.0
            Do K = 1, KBM1
               Do I = 1, MGL
                  Read (AGR, 1032) RTMP21 (I, K), RTMP22 (I, K), RTMP23 &
                 & (I, K)
               End Do
            End Do
            If (SERIAL) Then
               PM1 = RTMP21
               BMR1 = RTMP22
               BPR1 = RTMP23
            End If
#      if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do K = 1, KBM1
                  Do I = 1, MLOC
                     PM1 (I, K) = RTMP21 (NGID(I), K)
                     BMR1 (I, K) = RTMP22 (NGID(I), K)
                     BPR1 (I, K) = RTMP23 (NGID(I), K)
                  End Do
                  Do I = 1, NHN
                     PM1 (I+MLOC, K) = RTMP21 (HN_LST(I), K)
                     BMR1 (I+MLOC, K) = RTMP22 (HN_LST(I), K)
                     BPR1 (I+MLOC, K) = RTMP23 (HN_LST(I), K)
                  End Do
               End Do
            End If
#      endif
            Deallocate (RTMP21, RTMP22, RTMP23)
         End If
    !
    ! ARE REMAINING GROUP 2 PARAMETERS SPATIALLY VARYING?
    !
         Read (AGR, 1040) SPVAR2, PRINT2
         If (SPVAR2 == 'CONSTANT') Then
            Read (AGR, 1030) PM2 (1, 1), BMR2 (1, 1), BPR2 (1, 1)
            Do K = 1, KBM1
               Do I = 1, MTLOC
                  PM2 (I, K) = PM2 (1, 1)
                  BMR2 (I, K) = BMR2 (1, 1)
                  BPR2 (I, K) = BPR2 (1, 1)
               End Do
            End Do
         Else
            Allocate (RTMP21(MGL, KBM1))
            RTMP21 = 0.0
            Allocate (RTMP22(MGL, KBM1))
            RTMP22 = 0.0
            Allocate (RTMP23(MGL, KBM1))
            RTMP23 = 0.0
            Do K = 1, KBM1
               Do I = 1, MGL
                  Read (AGR, 1032) RTMP21 (I, K), RTMP22 (I, K), RTMP23 &
                 & (I, K)
               End Do
            End Do
            If (SERIAL) Then
               PM2 = RTMP21
               BMR2 = RTMP22
               BPR2 = RTMP23
            End If
#      if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do K = 1, KBM1
                  Do I = 1, MLOC
                     PM2 (I, K) = RTMP21 (NGID(I), K)
                     BMR2 (I, K) = RTMP22 (NGID(I), K)
                     BPR2 (I, K) = RTMP23 (NGID(I), K)
                  End Do
                  Do I = 1, NHN
                     PM2 (I+MLOC, K) = RTMP21 (HN_LST(I), K)
                     BMR2 (I+MLOC, K) = RTMP22 (HN_LST(I), K)
                     BPR2 (I+MLOC, K) = RTMP23 (HN_LST(I), K)
                  End Do
               End Do
            End If
#      endif
            Deallocate (RTMP21, RTMP22, RTMP23)
         End If
    !
    ! ARE REMAINING GROUP 3 PARAMETERS SPATIALLY VARYING?
    !
         Read (AGR, 1040) SPVAR3, PRINT3
         If (SPVAR3 == 'CONSTANT') Then
            Read (AGR, 1030) PM3 (1, 1), BMR3 (1, 1), BPR3 (1, 1)
            Do K = 1, KBM1
               Do I = 1, MTLOC
                  PM3 (I, K) = PM3 (1, 1)
                  BMR3 (I, K) = BMR3 (1, 1)
                  BPR3 (I, K) = BPR3 (1, 1)
               End Do
            End Do
         Else
            Allocate (RTMP21(MGL, KBM1))
            RTMP21 = 0.0
            Allocate (RTMP22(MGL, KBM1))
            RTMP22 = 0.0
            Allocate (RTMP23(MGL, KBM1))
            RTMP23 = 0.0
            Do K = 1, KBM1
               Do I = 1, MGL
                  Read (AGR, 1032) RTMP21 (I, K), RTMP22 (I, K), RTMP23 &
                 & (I, K)
               End Do
            End Do
            If (SERIAL) Then
               PM3 = RTMP21
               BMR3 = RTMP22
               BPR3 = RTMP23
            End If
#      if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do K = 1, KBM1
                  Do I = 1, MLOC
                     PM3 (I, K) = RTMP21 (NGID(I), K)
                     BMR3 (I, K) = RTMP22 (NGID(I), K)
                     BPR3 (I, K) = RTMP23 (NGID(I), K)
                  End Do
                  Do I = 1, NHN
                     PM3 (I+MLOC, K) = RTMP21 (HN_LST(I), K)
                     BMR3 (I+MLOC, K) = RTMP22 (HN_LST(I), K)
                     BPR3 (I+MLOC, K) = RTMP23 (HN_LST(I), K)
                  End Do
               End Do
            End If
#      endif
            Deallocate (RTMP21, RTMP22, RTMP23)
         End If
    !
    ! TIME DEPENDENCE OF PREDATION BY HIGHER TROPHIC LEVELS
    !
         Read (AGR, 1040) TPVAR, TPRINT
         If (TPVAR == 'CONSTANT') Then
            Read (AGR, 1060) TVPR (1)
            Do J = 2, 366
               TVPR (J) = TVPR (1)
            End Do
         Else
            Read (AGR, 1060) (TVPR(J), J=1, 366)
         End If
    !
    ! TIME DEPENDENCE OF GROWTH BY SPRING ALGAL GROUP
    !
         Read (AGR, 1040) TB2GR, PRINTB2
         If (TB2GR == 'CONSTANT') Then
            Read (AGR, 1060) B2GR (1)
            Do J = 2, 366
               B2GR (J) = B2GR (1)
            End Do
         Else
            Read (AGR, 1060) (B2GR(J), J=1, 366)!daily series of B2GR
         End If
    !
    !***** Input FORMAT statements
    !
1010     Format (A72)
1030     Format (/ / (8 X, 9 F8.0))
1032     Format (8 X, 9 F8.0)
1040     Format (/ / 8 X, 8 A8)
1060     Format (/ / (16 X, F8.0))
    !
    ! OUTPUT WHAT WAS INPUT
         Open (ALO, File=ALOFN)
         If (MSR) WRITE (ALO, 2002) (ALG_TITLE(J), J=1, 6)
    !
    ! WRITE ALGAL PROPERTIES WHICH ARE NOT SPATIALLY VARYING
    !
         If (MSR) WRITE (ALO, 3000)
         If (MSR) WRITE (ALO, 3010) CHLCMN1, CHLCMN2, CHLCMN3
         If (MSR) WRITE (ALO, 3012) ACHLC1, ACHLC2, ACHLC3
         If (MSR) WRITE (ALO, 3014) BCHLC1, BCHLC2, BCHLC3
         If (MSR) WRITE (ALO, 3016) CCHLC1, CCHLC2, CCHLC3
         If (MSR) WRITE (ALO, 3020) ANC1, ANC2, ANC3
         If (MSR) WRITE (ALO, 3022) Q01, Q02, Q03
         If (MSR) WRITE (ALO, 3030) ASC1, ASC2, ASC3
         If (MSR) WRITE (ALO, 3072) STF1, STF2, STF3
         If (MSR) WRITE (ALO, 3040) KHN1, KHN2, KHN3
         If (MSR) WRITE (ALO, 3050) KHP1, KHP2, KHP3
         If (MSR) WRITE (ALO, 3042) KHNFIX
         If (MSR) WRITE (ALO, 3043) FNFIX
         If (MSR) WRITE (ALO, 3060) KHS1, KHS2, KHS3
         If (MSR) WRITE (ALO, 3070) KHR1, KHR2, KHR3
         If (MSR) WRITE (ALO, 3074) KHST1, KHST2, KHST3
         If (MSR) WRITE (ALO, 3076) VMAX1, VMAX2, VMAX3
         If (MSR) WRITE (ALO, 3080) ALPHMIN1, ALPHMIN2, ALPHMIN3
         If (MSR) WRITE (ALO, 3090) ALPHRAT1, ALPHRAT2, ALPHRAT3
         If (MSR) WRITE (ALO, 4000) PRSP1, PRSP2, PRSP3
         If (MSR) WRITE (ALO, 4010) TMP1, TMP2, TMP3
         If (MSR) WRITE (ALO, 5020) TR1, TR2, TR3
         If (MSR) WRITE (ALO, 5030) KTG11, KTG21, KTG31
         If (MSR) WRITE (ALO, 5040) KTG12, KTG22, KTG32
         If (MSR) WRITE (ALO, 5050) KTB1, KTB2, KTB3
    !
    ! WRITE DISTRIBUTION OF BIOMASS UPON MORTALITY
    !
         If (MSR) WRITE (ALO, 4020)
         If (MSR) WRITE (ALO, 4030) FNI1, FNLD1, FNRD1, FNLP1, FNRP1
         If (MSR) WRITE (ALO, 4040) FPI1, FPLD1, FPRD1, FPLP1, FPRP1
         If (MSR) WRITE (ALO, 4050) FCLD1, FCRD1, FCLP1, FCRP1
         If (MSR) WRITE (ALO, 4060)
         If (MSR) WRITE (ALO, 4030) FNI2, FNLD2, FNRD2, FNLP2, FNRP2
         If (MSR) WRITE (ALO, 4040) FPI2, FPLD2, FPRD2, FPLP2, FPRP2
         If (MSR) WRITE (ALO, 4050) FCLD2, FCRD2, FCLP2, FCRP2
         If (MSR) WRITE (ALO, 4070)
         If (MSR) WRITE (ALO, 4030) FNI3, FNLD3, FNRD3, FNLP3, FNRP3
         If (MSR) WRITE (ALO, 4040) FPI3, FPLD3, FPRD3, FPLP3, FPRP3
         If (MSR) WRITE (ALO, 4050) FCLD3, FCRD3, FCLP3, FCRP3
         If (MSR) WRITE (ALO, 4090)
         If (MSR) WRITE (ALO, 4030) FNIP, FNLDP, FNRDP, FNLPP, FNRPP
         If (MSR) WRITE (ALO, 4040) FPIP, FPLDP, FPRDP, FPLPP, FPRPP
         If (MSR) WRITE (ALO, 4080) FDOP, FCLDP, FCRDP, FCLPP, FCRPP
         If (MSR) WRITE (ALO, 4082) FSAP, (1.-FSAP)
    !
         If (MINPROD == ' MINIMUM' .Or. MINPROD == ' PRODUCT') Then
            If (MSR) WRITE (ALO, 2010) MINPROD
         Else
            If (MSR) WRITE (ALO, 2020) MINPROD
            Stop
         End If
    !
    ! WRITE SPATIALLY-VARYING PARAMETERS FOR GROUP 1
    !
         If (MSR) WRITE (ALO, 2000)
         If (SPVAR1 == 'CONSTANT') Then
            If (MSR) WRITE (ALO, 2085)
         Else
            If (MSR) WRITE (ALO, 2086)
         End If
         If (MSR) WRITE (ALO, 5060)
         If (PRINT1 /= '     ALL') Then
            If (MSR) WRITE (ALO, 5070) PM1 (1, 1), BMR1 (1, 1), BPR1 &
           & (1, 1)
         Else
            Do K = 1, KBM1
               Do I = 1, MLOC
                  If (MSR) WRITE (ALO, 5080) I, K, PM1 (I, K), BMR1 (I, &
                 & K), BPR1 (I, K)
               End Do
            End Do
         End If
    !
    ! WRITE SPATIALLY-VARYING PARAMETERS FOR GROUP 2
    !
         If (MSR) WRITE (ALO, 2005)
         If (SPVAR2 == 'CONSTANT') Then
            If (MSR) WRITE (ALO, 2085)
         Else
            If (MSR) WRITE (ALO, 2086)
         End If
         If (MSR) WRITE (ALO, 5060)
         If (PRINT2 /= '     ALL') Then
            If (MSR) WRITE (ALO, 5070) PM2 (1, 1), BMR2 (1, 1), BPR2 &
           & (1, 1)
         Else
            Do K = 1, KBM1
               Do I = 1, MLOC
                  If (MSR) WRITE (ALO, 5080) I, K, PM2 (I, K), BMR2 (I, &
                 & K), BPR2 (I, K)
               End Do
            End Do
         End If
    !
    ! WRITE SPATIALLY-VARYING PARAMETERS FOR GROUP 3
    !
         If (MSR) WRITE (ALO, 2006)
         If (SPVAR3 == 'CONSTANT') Then
            If (MSR) WRITE (ALO, 2085)
         Else
            If (MSR) WRITE (ALO, 2086)
         End If
         If (MSR) WRITE (ALO, 5060)
         If (PRINT3 /= '     ALL') Then
            If (MSR) WRITE (ALO, 5070) PM3 (1, 1), BMR3 (1, 1), BPR3 &
           & (1, 1)
         Else
            Do K = 1, KBM1
               Do I = 1, MLOC
                  If (MSR) WRITE (ALO, 5080) I, K, PM3 (I, K), BMR3 (I, &
                 & K), BPR3 (I, K)
               End Do
            End Do
         End If
    !
    ! WRITE TERMS FOR PREDATION BY HIGHER TROPHIC LEVELS
    !
         If (TPVAR == 'CONSTANT') Then
            If (MSR) WRITE (ALO, 2087)
         Else
            If (MSR) WRITE (ALO, 2088)
         End If
         If (MSR) WRITE (ALO, 2040) TRPR, KTPR
         If (MSR) WRITE (ALO, 2092)
         If (TPRINT /= '     ALL') Then
            If (MSR) WRITE (ALO, 2090) (J, TVPR(J), J=1, 1)
         Else
            If (MSR) WRITE (ALO, 2090) (J, TVPR(J), J=1, 366)
         End If
    !
    ! WRITE TIME DEPENDENCE OF GROWTH OF SPRING ALGAL GROUP
    !
         If (TB2GR == 'CONSTANT') Then
            If (MSR) WRITE (ALO, 2094)
         Else
            If (MSR) WRITE (ALO, 2095)
         End If
         If (MSR) WRITE (ALO, 2093)
         If (PRINTB2 /= '     ALL') Then
            If (MSR) WRITE (ALO, 2090) (J, B2GR(J), J=1, 1)
         Else
            If (MSR) WRITE (ALO, 2090) (J, B2GR(J), J=1, 366)
         End If
    !
    !***** Output FORMAT statements
    !
2000     Format (/ ' ALGAL GROUP 1')
2002     Format (1 X, A72)
2005     Format (/ ' ALGAL GROUP 2')
2006     Format (/ ' ALGAL GROUP 3')
2010     Format (/ A8, ' FORMULATION SELECTED FOR LIGHT AND NUTRIENT', &
        & ' LIMITATIONS')
2020     Format (/ ' MINPROD INCORRECTLY SPECIFIED AS ', A8)
2040     Format (' PREDATION SPECIFIED AT ', F8.2, ' C.'/' TEMPERATURE &
        &EFFECT = ', F8.3, ' PER DEGREE')
2085     Format (/ ' REMAINING PARAMETERS ARE SPATIALLY-INVARIANT')
2086     Format (/ ' REMAINING PARAMETERS ARE SPATIALLY-VARYING')
2087     Format (/ ' HIGHER-LEVEL PREDATION IS TEMPORALLY-INVARIANT')
2088     Format (/ ' HIGHER-LEVEL PREDATION VARIES TEMPORALLY')
2090     Format (I8, F8.3)
2092     Format (/ '    DAY     TVPR')
2093     Format (/ '    DAY     B2GR')
2094     Format (/ ' BASE GROUP 2 GROWTH IS TEMPORALLY-INVARIANT')
2095     Format (/ ' BASE GROUP 2 VARIES TEMPORALLY')
3000     Format (/ ' ALGAL PROPERTIES', 30 X, 'GROUP 1   GROUP 2   GROU&
        &P 3')
3010     Format (/ ' MINIMUM GM CHL/GM C    ', 20 X, 3 F10.4)
3012     Format (' ACHCL (GM CHL/GM C)    ', 20 X, 3 F10.4)
3014     Format (' BCHCL (1/DEGREE C)     ', 20 X, 3 F10.4)
3016     Format (' CCHCL (M**2 DAY/E)     ', 20 X, 3 F10.4)
3020     Format (' GM N/GM C              ', 20 X, 3 F10.4)
3022     Format (' MINIMUM GM P/GM C      ', 20 X, 3 F10.4)
3030     Format (' GM SI/GM C             ', 20 X, 3 F10.4)
3040     Format (' KHN (GM N/M**3)        ', 20 X, 3 F10.4)
3042     Format (' KHNFIX (GM N/M**3)     ', 20 X, F10.4)
3043     Format (' FRACTION N FIXERS      ', 20 X, F10.4)
3050     Format (' KHP (GM P/M**3)        ', 20 X, 3 F10.4)
3060     Format (' KHS (GM SI/M**3)       ', 20 X, 3 F10.4)
3070     Format (' KHR (GM DO/M**3)       ', 20 X, 3 F10.4)
3072     Format (' SALT TOXICITY (1/DAY)  ', 20 X, 3 F10.4)
3074     Format (' KHST (ppt)             ', 20 X, 3 F10.4)
3076     Format (' VMAX (GM P/GM C/DAY)   ', 20 X, 3 F10.4)
3080     Format (' ALPHA (GM C/GM CHL/DAY)/(uE/M**2/SEC)      ', 3 &
        & F10.3)
3090     Format (' ALPHA TO P RATIO (1/(uE/M**2/SEC))', 9 X, 3 F10.3)
4000     Format (' PHOTORESPIRATION FRACTION     ', 13 X, 3 F10.3)
4010     Format (' OPTIMAL TEMPERATURE FOR PRODUCTION (C)', 5 X, 3 &
        & F10.3)
4020     Format (/ ' DISTRIBUTION OF ALGAE UPON MORTALITY' // ' GROUP 1&
        & RESPIRATION  DIS INORG  LAB DISS  REF DISS  LAB PART  REF PAR&
        &T')
4030     Format (' NITROGEN            ', 5 F10.3)
4040     Format (' PHOSPHORUS          ', 5 F10.3)
4050     Format (' CARBON              ', 10 X, 4 F10.3)
4060     Format (/ ' GROUP 2 RESPIRATION')
4070     Format (/ ' GROUP 3 RESPIRATION')
4080     Format (' CARBON              ', 5 F10.3)
4082     Format (' SILICA              ', F10.3, 10 X, F10.3)
4090     Format (/ ' PREDATION')
5000     Format (' PHOSPHORUS     ', 4 F10.3)
5010     Format (' CARBON         ', 4 F10.3)
5020     Format (' REFERENCE TEMPERATURE FOR RESPIRATION (C)  ', 3 &
        & F10.3)
5030     Format (' EFFECT OF SUBOPTIMAL TEMP ON PRODUCTION    ', 3 &
        & F10.4)
5040     Format (' EFFECT OF SUPEROPTIMAL TEMP ON PRODUCTION  ', 3 &
        & F10.4)
5050     Format (' EFFECT OF TEMPERATURE ON RESPIRATION       ', 3 &
        & F10.4)
5060     Format ('       CELL     PMAX     METAB     PRDTN',/'         &
        &    C/CHL/DAY   1/DAY     1/DAY' /)
5070     Format (10 X, 3 F10.3)
5080     Format (2 I10, 3 F10.3)
    !
         Return
      End Subroutine ALG_READ
  !
  !
  !
  !************************************************************************
  !**                  S U B R O U T I N E   A L G A E                   **
  !************************************************************************
  !
      Subroutine ALGAE (DTB1, DTB2, DTB3, FLXS1, FLXS2, FLXS3)
         Implicit None
         Integer :: I, J, K, ITEMP !, F !LBcleanup
         Real (SP), Dimension (0:MTLOC, KBM1) :: DTB1, DTB2, DTB3
         Real (SP), Dimension (0:MTLOC, KBM1) :: FLXS1, FLXS2, FLXS3
         Real (SP), Allocatable :: XL (:), STOX1 (:, :), STOX2 (:, :), &
        & STOX3 (:, :)
         Real (SP), Allocatable :: ITBOT (:, :), NL1NNF (:, :)
         Real (SP) :: GPP1, GPP2, GPP3 !FI01, FI02, FI03, , TREC!LBcleanup
         Real (SP) :: NETP1, NETP2, NETP3, IK, ITTOP, ITAVG, NUTLIM
         Real (SP) :: DIN, ALPHA, OPTDEPTH !SALTOX, DSIL,   !LBcleanup
    !
    !
         Allocate (XL(MLOC))
         XL = 0.0
         Allocate (STOX1(MLOC, KBM1))
         STOX1 = 0.0
         Allocate (STOX2(MLOC, KBM1))
         STOX2 = 0.0
         Allocate (STOX3(MLOC, KBM1))
         STOX3 = 0.0
         Allocate (NL1NNF(MLOC, KBM1))
         NL1NNF = 0.0
         Allocate (ITBOT(MLOC, KBM1))
         ITBOT = 0.0
    !
    ! DETERMINE JULIAN DAY
    !
         J = AINT(1.0 + AMOD (JDAY, 365.25))

!!! SALINITY TOXICITY
    ! TEMPERATURE EFFECTS ON PRODUCTION AND RESPIRATION
    !
         Do K = 1, KBM1
            Do I = 1, MLOC

               P1 (I, K) = PM1 (I, K) * GET_FT (T(I, K), TMP1, KTG11, &
              & KTG12)
               P2 (I, K) = PM2 (I, K) * GET_FT (T(I, K), TMP2, KTG21, &
              & KTG22) * B2GR (J)

               BM1 (I, K) = BMR1 (I, K) * GET_FTBM (T(I, K), TR1, KTB1) &
              & + STOX1 (I, K)!ALG1 metabolism
               BM2 (I, K) = BMR2 (I, K) * GET_FTBM (T(I, K), TR2, KTB2) &
              & + STOX2 (I, K)!ALG2 metabolism

               PR1 (I, K) = BPR1 (I, K) * B1 (I, K) * B1 (I, K) * &
              & GET_FTPR (T(I, K), TRPR, KTPR) * TVPR (J)

               PR2 (I, K) = BPR2 (I, K) * B2 (I, K) * B2 (I, K) * &
              & GET_FTPR (T(I, K), TRPR, KTPR) * TVPR (J)

            End Do
         End Do
    !
    ! NUTRIENT LIMITATIONS - NITROGEN
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
               DIN = NH4 (I, K) + NO3 (I, K)

               PN1 (I, K) = NH4 (I, K) * (NO3(I, K)/((KHN1+NH4(I, &
              & K))*(KHN1+NO3(I, K)))+KHN1/((1.E-30+NH4(I, K)+NO3(I, &
              & K))*(KHN1+NO3(I, K))))
               PN2 (I, K) = NH4 (I, K) * (NO3(I, K)/((KHN2+NH4(I, &
              & K))*(KHN2+NO3(I, K)))+KHN2/((1.E-30+NH4(I, K)+NO3(I, &
              & K))*(KHN2+NO3(I, K))))

               NL1 (I, K) = (NH4(I, K)+DIN-PN1(I, K)*DIN) / &
              & (KHN1+NH4(I, K)+DIN-PN1(I, K)*DIN+1E-30)
               NL1NNF (I, K) = (NH4(I, K)+DIN-PN1(I, K)*DIN) / &
              & (KHN1+NH4(I, K)+DIN-PN1(I, K)*DIN+1E-30)
               NL2 (I, K) = (NH4(I, K)+DIN-PN2(I, K)*DIN) / &
              & (KHN2+NH4(I, K)+DIN-PN2(I, K)*DIN+1E-30)

            End Do
         End Do
    !
    !
    ! PHOSPHORUS
    !! DROOP FORMULA FOR P LIMIT
    !! 
         Do K = 1, KBM1
            Do I = 1, MLOC

               PL1 (I, K) = PO4 (I, K) / (KHP1+PO4(I, K))
               PL2 (I, K) = PO4 (I, K) / (KHP2+PO4(I, K))

            End Do
         End Do
    !
    !! SILICA

    ! PHOTOSYNTHESIS VS IRRADIANCE

         Do K = 1, KBM1
            Do I = 1, MLOC
               ALPHA = ALPHMIN1 + ALPHRAT1 * P1 (I, K)
          !Changed IK to reflect Pmax/alpha - not modified P
               IK = P1 (I, K) / (ALPHA+1.0E-10)!

               IK1 (I, K) = IK

               FI1 (I, K) = IAVG (I, K) / (Sqrt(IK*IK+IAVG(I, &
              & K)*IAVG(I, K)+1.0E-10))

               ALPHA = ALPHMIN2 + ALPHRAT2 * P2 (I, K)

               IK = P2 (I, K) / (ALPHA+1.0E-10)!
          !
               IK2 (I, K) = IK
          !
               FI2 (I, K) = IAVG (I, K) / (Sqrt(IK*IK+IAVG(I, &
              & K)*IAVG(I, K)+1.0E-10))

            End Do
         End Do
    !
    !     *** Compute carbon to chlorophyll ratio in each cell
    ! THESE ARE SOMETIME ALL ZERO - MUCH SIMPLER CODE
         If (ACHLC1 .Eq. 0 .And. ACHLC2 .Eq. 0 .And. ACHLC3 .Eq. 0) &
        & Then
            Do K = 1, KBM1
               Do I = 1, MLOC
                  CCHL1 (I, K) = 1. / CHLCMN1
                  CCHL2 (I, K) = 1. / CHLCMN2
                  CCHL3 (I, K) = 1. / CHLCMN3
               End Do
            End Do
       !
         Else ! 14 SEPT 2009
            Do I = 1, MLOC
               ITTOP = IT
               OPTDEPTH = KESS (I, 1) * D (I) * DZ2D (I,1)
               ITBOT (I, 1) = ITTOP * Exp (-OPTDEPTH)
               ITAVG = (ITTOP-ITBOT(I, 1)) / OPTDEPTH
               NUTLIM = AMIN1 (NL1(I, 1), PL1(I, 1))

               CCHL1 (I, 1) = 1. / (CHLCMN1+ACHLC1*Exp(BCHLC1*T(I, &
              & 1))*NUTLIM/Exp(CCHLC1*ITAVG))
               NUTLIM = AMIN1 (NL2(I, 1), PL2(I, 1))

               CCHL2 (I, 1) = 1. / (CHLCMN2+ACHLC2*Exp(BCHLC2*T(I, &
              & 1))*NUTLIM/Exp(CCHLC2*ITAVG))

            End Do
       !
            Do K = 2, KBM1
               Do I = 1, MLOC
                  ITTOP = ITBOT (I, K-1)
                  OPTDEPTH = KESS (I, K) * D (I) * DZ2D (I,K)
                  ITBOT (I, K) = ITTOP * Exp (-OPTDEPTH)
                  ITAVG = (ITTOP-ITBOT(I, K)) / OPTDEPTH
                  NUTLIM = AMIN1 (NL1(I, K), PL1(I, K))

                  CCHL1 (I, K) = 1. / (CHLCMN1+ACHLC1*Exp(BCHLC1*T(I, &
                 & K))*NUTLIM/Exp(CCHLC1*ITAVG))
                  NUTLIM = AMIN1 (NL2(I, K), PL2(I, K))

                  CCHL2 (I, K) = 1. / (CHLCMN2+ACHLC2*Exp(BCHLC2*T(I, &
                 & K))*NUTLIM/Exp(CCHLC2*ITAVG))

               End Do
            End Do
         End If ! ENDIF 14 SEPT 2009
    !
    !END output section
    ! COMPUTE ASSIMILATION RATIO AT WATER SURFACE (GM C/GM CHL/DAY),
    ! CARBON FIXATION (GM C/M**3/DAY)  COPIED FROM CHES BAY CODE MARCH 14, 2007
    !
         Do I = 1, MLOC

            ASRAT (I) = (P1(I, 1)*AMIN1(NL1(I, 1), PL1(I, 1))*B1(I, &
           & 1)+P2(I, 1)*AMIN1(NL2(I, 1), PL2(I, 1))*B2(I, 1)) / (B1(I, &
           & 1)+B2(I, 1)+1.0e-6)

       ! Below without P or Si limitation
            NETP1 = (P1(I, 1)*AMIN1(NL1(I, 1), PL1(I, &
           & 1))*(1.-PRSP1)/CCHL1(I, 1)-BM1(I, 1)) * B1 (I, 1)
            NETP2 = (P2(I, 1)*AMIN1(NL2(I, 1), PL2(I, &
           & 1))*(1.-PRSP2)/CCHL2(I, 1)-BM2(I, 1)) * B2 (I, 1)

            CFIX (I) = NETP1 + NETP2
         End Do
    !
    ! EFFECTS OF LIGHT AND NUTRIENTS ON PRODUCTION
    ! CONVERT FROM GM C/GM CHL/DAY TO SPECIFIC GROWTH RATE
    !
         If (MINPROD == ' MINIMUM') Then
       !
            Do K = 1, KBM1
               Do I = 1, MLOC
             !Below without P or Si limitation
                  P1NNF (I, K) = P1 (I, K) * AMIN1 (FI1(I, K), &
                 & NL1NNF(I, K), PL1(I, K)) / CCHL1 (I, K)
                  P1 (I, K) = P1 (I, K) * AMIN1 (FI1(I, K), NL1(I, K), &
                 & PL1(I, K)) / CCHL1 (I, K)
                  P2 (I, K) = P2 (I, K) * AMIN1 (FI2(I, K), NL2(I, K), &
                 & PL2(I, K)) / CCHL2 (I, K)

               End Do
            End Do
       !
         Else
       !
            Do K = 1, KBM1
               Do I = 1, MLOC

             !Below without P or Si limitation
                  P1NNF (I, K) = P1 (I, K) * FI1 (I, K) * AMIN1 &
                 & (NL1NNF(I, K), PL1(I, K)) / CCHL1 (I, K)
                  P1 (I, K) = P1 (I, K) * FI1 (I, K) * AMIN1 (NL1(I, &
                 & K), PL1(I, K)) / CCHL1 (I, K)
                  P2 (I, K) = P2 (I, K) * FI2 (I, K) * AMIN1 (NL2(I, &
                 & K), PL2(I, K)) / CCHL2 (I, K)

               End Do
            End Do
       !
         End If
    ! RATE OF CHANGE DUE TO PRODUCTION, RESPIRATION, PREDATION
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
          !
               NETP1 = (P1(I, K)*(1.-PRSP1)-BM1(I, K)) * B1 (I, K)
               NETP2 = (P2(I, K)*(1.-PRSP2)-BM2(I, K)) * B2 (I, K)

               GPP1 = P1 (I, K) * B1 (I, K)   !gC/m^3/day
               GPP2 = P2 (I, K) * B2 (I, K)

               GPP (I, K) = (GPP1+GPP2) * D (I) * DZ2D (I,K)   ! gC/m^3/day * m ==> gC/m^2/day
               NPP (I, K) = (NETP1+NETP2) * D (I) * DZ2D (I,K)
          ! is this a rate or a total? looks like a rate. when added back to B1?
               DTB1 (I, K) = (NETP1-PR1(I, K)-B1SZ(I, K)-B1LZ(I, K)) / &
              & 86400.

               DTB2 (I, K) = (NETP2-PR2(I, K)-B2SZ(I, K)-B2LZ(I, K)) / &
              & 86400.

            End Do
         End Do
         TOTAL_NETPP (:) = 0.0
         Do K = 1, KBM1
            Do I = 1, MLOC
               TOTAL_NETPP (I) = TOTAL_NETPP (I) + NPP (I, K) !gC/m^2/day
            End Do
         End Do
    ! NITROGEN FIXATION IN SURFACE BOXES (G N/G CHL/DAY)
    !
         Do I = 1, MLOC
            SNFIX (I) = (P1(I, 1)-P1NNF(I, 1)) * ANC1 * CCHL1 (I, 1)
         End Do
    !
    ! RATE OF CHANGE DUE TO SETTLING
         Do I = 1, MLOC
       ! CHANGED X/Y/86400 TO X/(Y*86400)
            DTB1 (I, 1) = DTB1 (I, 1) - WS1 (I, 1) * B1 (I, 1) / &
           & (D(I)*DZ2D(I,1)*86400.)
            DTB2 (I, 1) = DTB2 (I, 1) - WS2 (I, 1) * B2 (I, 1) / &
           & (D(I)*DZ2D(I,1)*86400.)

         End Do

         Do K = 2, KBM1
            Do I = 1, MLOC
          ! CHANGED X/Y/86400 TO X/(Y*86400)
               DTB1 (I, K) = DTB1 (I, K) + (WS1(I, K-1)*B1(I, &
              & K-1)-WS1(I, K)*B1(I, K)) / (D(I)*DZ2D(I,K)*86400.)
               DTB2 (I, K) = DTB2 (I, K) + (WS2(I, K-1)*B2(I, &
              & K-1)-WS2(I, K)*B2(I, K)) / (D(I)*DZ2D(I,K)*86400.)

            End Do
         End Do

    ! SETTLING FLUX FOR MASS BALANCE
    !
         Do K = 1, KBM1
            Do I = 1, MLOC
               FLXS1 (I, K) = WS1 (I, K) * B1 (I, K) * V2 (I, K) / &
              & (D(I)*DZ2D(I,K)*86400.)
               FLXS2 (I, K) = WS2 (I, K) * B2 (I, K) * V2 (I, K) / &
              & (D(I)*DZ2D(I,K)*86400.)

            End Do
         End Do
    !
    !  RESUSPENSION
    !
         If (SEDIMENT_CALC) Then
       ! MAKE THIS SIMPLER LOOP
       !
       !
       !  resuspension rate should be (WS1-WS1NET ) in m/d
       !  and resuspension is counter-acting on the sinking of WS1
       !  so the net sinking is WS1NET = WS1 - (WS1-WS1NET) , i.e. sinking - resuspension
       !
       !
            Do I = 1, MLOC
               XL (I) = 1.0 / (D(I)*DZ2D(I,KBM1)*86400)
               DTB1 (I, KBM1) = DTB1 (I, KBM1) + (WS1(I, &
              & KBM1)-WS1NET(I)) * B1 (I, KBM1) * XL (I)
               DTB2 (I, KBM1) = DTB2 (I, KBM1) + (WS2(I, &
              & KBM1)-WS2NET(I)) * B2 (I, KBM1) * XL (I)

            End Do

         End If
    !
         Deallocate (XL, STOX1, STOX2, STOX3, NL1NNF, ITBOT)
    !
    !
         Return
      End Subroutine ALGAE
  !
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
End Module MOD_ALGAL
