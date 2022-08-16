!mod_wd.F
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
!**                     Adopted from CE-QUAL-ICM  Model                 **
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
Module MOD_ZOOP
  !
      Use MOD_PREC, Only: SP
      Use MOD_LIMS, Only: MTLOC, KBM1
  !
      Use MOD_SIZES, Only: MGL
  !
      Use MOD_WQM, Only: ZFOFN, ZOOFN !,         &!
  !
  !Zooplankton parameters
  !
  !    Integer :: KTBSZ, KTBLZ
      Real (SP) ::  KTBSZ, KTBLZ
  !
      Real (SP), Dimension (366) :: TVPRSZ, TVPRLZ
  !
  !
      Real (SP) :: CTSZ, CTLZ, KHCSZ, KHCLZ, MZEROSZ, MZEROLZ !mortality at zero dissolved oxygen (1/day) for mesozooplankton
  !
      Real (SP) :: UB1SZ, UB2SZ, UB3SZ, UB1LZ, UB2LZ, UB3LZ, UDSZ, &
     & ULDSZ, URDSZ, ULPSZ, URPSZ, ULLZ, URLZ, USZLZ, TRSZ, TRLZ, &
     & DOCRITSZ, DOCRITLZ, ANCSZ, ANCLZ, APCSZ, APCLZ, AOCRSZ, AOCRLZ, &
     & FRSASZ, FRSALZ !fraction of mesozooplankton silica recycled to dissolved silica pool, range [0,1]
  !
      Real (SP) :: FLDOCSZ, FRDOCSZ, FLPOCSZ, FRPOCSZ, FLDONSZ, &
     & FRDONSZ, FLPONSZ, FRPONSZ, FLDOPSZ, FRDOPSZ, FLPOPSZ, FRPOPSZ, &
     & FNH4SZ, FPO4SZ, FLDOCLZ, FRDOCLZ, FLPOCLZ, FRPOCLZ, FLDONLZ, &
     & FRDONLZ, FLPONLZ, FRPONLZ, FLDOPLZ, FRDOPLZ, FLPOPLZ, FRPOPLZ, &
     & FNH4LZ, FPO4LZ, FUREASZ, FUREALZ
  ! FLDOCSZ,&    !fraction of microzooplankton carbon released to LDOC, range [0,1]
  ! FRDOCSZ,&    !fraction of microzooplankton carbon released to RDOC, range [0,1]
  ! FLPOCSZ,&    !fraction of microzooplankton carbon released to LPOC, range [0,1]
  ! FRPOCSZ,&    !fraction of microzooplankton carbon released to RPOC, range [0,1]
  ! FLDONSZ,&    !fraction of microzooplankton nitrogen released to LDON, range [0,1]
  ! FRDONSZ,&    !fraction of microzooplankton nitrogen released to RDON, range [0,1]
  ! FLPONSZ,&    !fraction of microzooplankton nitrogen released to LPON, range [0,1]
  ! FRPONSZ,&    !fraction of microzooplankton nitrogen released to RPON, range [0,1]
  ! FLDOPSZ,&    !fraction of microzooplankton phosphorus released to LDOP range [0,1]
  ! FRDOPSZ,&    !fraction of microzooplankton phosphorus released to RDOP range [0,1]
  ! FLPOPSZ,&    !fraction of microzooplankton phosphorus released to LPOP range [0,1]
  ! FRPOPSZ,&    !fraction of microzooplankton phosphorus released to RPOP range [0,1]
  ! FNH4SZ, &    !fraction of microzooplankton nitrogen recycled to DIN as NH4, range [0,1]
  ! FPO4SZ, &    !fraction of microzooplankton phosphorus recycled to DIP as PO4, range [0,1]
  !
  ! FLDOCLZ,&    !fraction of mesozooplankton carbon released to LDOC, range [0,1]
  ! FRDOCLZ,&    !fraction of mesozooplankton carbon released to RDOC, range [0,1]
  ! FLPOCLZ,&    !fraction of mesozooplankton carbon released to LPOC, range [0,1]
  ! FRPOCLZ,&    !fraction of mesozooplankton carbon released to RPOC, range [0,1]
  ! FLDONLZ,&    !fraction of mesozooplankton nitrogen released to LDON, range [0,1]
  ! FRDONLZ,&    !fraction of mesozooplankton nitrogen released to RDON, range [0,1]
  ! FLPONLZ,&    !fraction of mesozooplankton nitrogen released to LPON, range [0,1]
  ! FRPONLZ,&    !fraction of mesozooplankton nitrogen released to RPON, range [0,1]
  ! FLDOPLZ,&    !fraction of mesozooplankton phosphorus released to LDOP range [0,1]
  ! FRDOPLZ,&    !fraction of mesozooplankton phosphorus released to RDOP range [0,1]
  ! FLPOPLZ,&    !fraction of mesozooplankton phosphorus released to LPOP range [0,1]
  ! FRPOPLZ,&    !fraction of mesozooplankton phosphorus released to RPOP range [0,1]
  ! FNH4LZ, &    !fraction of mesozooplankton nitrogen recycled to DIN as NH4, range [0,1]
  ! FPO4LZ, &    !fraction of mesozooplankton phosphorus recycled to DIP as PO4, range [0,1]
  ! FUREASZ,&    !Never used !!!
  ! FUREALZ      !Never used!!!
!
      Real (SP), Allocatable :: B1ASZ (:, :), B2ASZ (:, :), B3ASZ (:, &
     & :), LPOCASZ (:, :), RPOCASZ (:, :), PRASZ (:, :), B1ALZ (:, :), &
     & B2ALZ (:, :), B3ALZ (:, :), SZALZ (:, :), LPOCALZ (:, :), &
     & RPOCALZ (:, :), PRALZ (:, :), CLSZ (:, :), CLLZ (:, :), RSZ (:, &
     & :), RLZ (:, :), RMAXSZ (:, :), RMAXLZ (:, :), BMSZ (:, :), BMLZ &
     & (:, :), BMRSZ (:, :), BMRLZ (:, :), MSZ (:, :), MLZ (:, :), &
     & PRSZLZ (:, :), GSZ (:, :), GLZ (:, :), ESZ (:, :), ELZ (:, :), &
     & RFSZ (:, :), RFLZ (:, :), PRSZ (:, :), PRLZ (:, :), LDOCASZ (:, &
     & :), BPRSZ (:, :), BPRLZ (:, :), RDOCASZ (:, :), LDOCSZ (:, :), &
     & LPOCSZ (:, :), RPOCSZ (:, :), LDOCLZ (:, :), LPOCLZ (:, :), &
     & RPOCLZ (:, :), NH4SZ (:, :), LDONSZ (:, :), LPONSZ (:, :), &
     & RPONSZ (:, :), NH4LZ (:, :), LDONLZ (:, :), LPONLZ (:, :), &
     & RPONLZ (:, :), PO4SZ (:, :), LDOPSZ (:, :), LPOPSZ (:, :), &
     & RPOPSZ (:, :), PO4LZ (:, :), LDOPLZ (:, :), LPOPLZ (:, :), &
     & RPOPLZ (:, :), RDOCSZ (:, :), RDONSZ (:, :), RDOPSZ (:, :), &
     & RDOCLZ (:, :), RDONLZ (:, :), RDOPLZ (:, :), PIB1SZ (:, :), &
     & PIB2SZ (:, :), PIB3SZ (:, :), PIB1LZ (:, :), PIB2LZ (:, :), &
     & PIB3LZ (:, :), B1SZ (:, :), B2SZ (:, :), B3SZ (:, :), B1LZ (:, &
     & :), B2LZ (:, :), B3LZ (:, :), DOSZ (:, :), DOLZ (:, :), SASZ (:, &
     & :), SUSZ (:, :), SALZ (:, :), SULZ (:, :)
  !
  !
      Real (SP), Allocatable :: ACLSZ (:, :), ACLLZ (:, :), ARSZ (:, &
     & :), ARLZ (:, :), ABMSZ (:, :), ABMLZ (:, :), AMSZ (:, :), AMLZ &
     & (:, :), APRSZLZ (:, :), AGSZ (:, :), AGLZ (:, :), ADOCSZ (:, :), &
     & APOCSZ (:, :), ADOCLZ (:, :), APOCLZ (:, :), ANH4SZ (:, :), &
     & ADONSZ (:, :), APONSZ (:, :), ANH4LZ (:, :), ADONLZ (:, :), &
     & APONLZ (:, :), APO4SZ (:, :), ADOPSZ (:, :), APOPSZ (:, :), &
     & APO4LZ (:, :), ADOPLZ (:, :), APOPLZ (:, :), APRSZ (:, :), APRLZ &
     & (:, :), APISZ (:, :), APILZ (:, :)
  !
      Real (SP), Allocatable :: AB1SZ (:, :), AB2SZ (:, :), AB3SZ (:, &
     & :), AB1LZ (:, :), AB2LZ (:, :), AB3LZ (:, :), ADOSZ (:, :), &
     & ADOLZ (:, :), ASASZ (:, :), ASUSZ (:, :), ASALZ (:, :), ASULZ &
     & (:, :)
  !
      Real (SP), Dimension (-50:400) :: FTLZ, FTSZ, FTBMSZ, FTBMLZ, &
     & FTPRSZ, FTPRLZ
  !
      Real (SP), Allocatable :: SZ (:, :), LZ (:, :)
      Real (SP) :: TMSZ, KTGSZ1, KTGSZ2, TRPRSZ, KTPRSZ
      Real (SP) :: TMLZ, KTGLZ1, KTGLZ2, TRPRLZ, KTPRLZ
!
Contains
!
  !Subroutines
  !	subroutine 	ZOOP_READ()
  !	subroutine	ZOOP_ALLOC()
  !	subroutine	ZOOP_DEALLOC()
  !
  !
  !************************************************************************
  !**             S U B R O U T I N E   Z O OP _ R E A D                 **
  !************************************************************************
  !
      Subroutine ZOOP_READ
    !
         Use MOD_FILEINFO, Only: ZOO, ZFO
         Use MOD_LIMS, Only: MLOC
    !
    !Took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
         Use MOD_CONTROL, Only: SERIAL, MSR, PAR
#if defined (MULTIPROCESSOR)
         Use MOD_PAR, Only: NGID, NHN, HN_LST
#endif
         Implicit None
         Character (Len=72) :: TITLE (6)
    !
         Character (Len=8) :: SPVARSZ, SPVARLZ, PRINTSZ, PRINTLZ, &
        & TVARSZ, TPRINTSZ, TVARLZ, TPRINTLZ

         Real (SP) :: FTB_ZOOP (0:35)
    !
         Integer :: I, J, K
    !
         Real (SP) :: TLOOK !look up able temperature for temperature control
    !
    !temporary variables for reading data
         Real (SP), Allocatable :: RTMP21 (:, :), RTMP22 (:, :), RTMP23 &
        & (:, :), RTMP24 (:, :), RTMP25 (:, :)
    !
    ! TITLE CARDS of zooplanton input file
    !
         Open (ZOO, File=ZOOFN, Status='OLD')
         Read (ZOO, 1010) (TITLE(J), J=1, 6)
    !
    ! READ SPATIALLY-INVARIANT PARAMETERS FOR GROUP 1
    !
         Read (ZOO, 1030) CTSZ, KHCSZ, DOCRITSZ, MZEROSZ
         Read (ZOO, 1030) ANCSZ, APCSZ, AOCRSZ
         Read (ZOO, 1030) UB1SZ, UB2SZ, UB3SZ, ULDSZ, URDSZ, ULPSZ, &
        & URPSZ
         Read (ZOO, 1030) TMSZ, KTGSZ1, KTGSZ2, TRSZ, KTBSZ, TRPRSZ, &
        & KTPRSZ
         Read (ZOO, 1030) FLDOCSZ, FRDOCSZ, FLPOCSZ, FRPOCSZ
         Read (ZOO, 1030) FNH4SZ, FUREASZ, FLDONSZ, FRDONSZ, FLPONSZ, &
        & FRPONSZ
         Read (ZOO, 1030) FPO4SZ, FLDOPSZ, FRDOPSZ, FLPOPSZ, FRPOPSZ
         Read (ZOO, 1030) FRSASZ
    !
    ! READ SPATIALLY-INVARIANT PARAMETERS FOR GROUP 2
    !
         Read (ZOO, 1030) CTLZ, KHCLZ, DOCRITLZ, MZEROLZ
         Read (ZOO, 1030) ANCLZ, APCLZ, AOCRLZ
         Read (ZOO, 1030) UB1LZ, UB2LZ, UB3LZ, USZLZ, ULLZ, URLZ
         Read (ZOO, 1030) TMLZ, KTGLZ1, KTGLZ2, TRLZ, KTBLZ, TRPRLZ, &
        & KTPRLZ
         Read (ZOO, 1030) FLDOCLZ, FRDOCLZ, FLPOCLZ, FRPOCLZ
         Read (ZOO, 1030) FNH4LZ, FUREALZ, FLDONLZ, FRDONLZ, FLPONLZ, &
        & FRPONLZ
         Read (ZOO, 1030) FPO4LZ, FLDOPLZ, FRDOPLZ, FLPOPLZ, FRPOPLZ
         Read (ZOO, 1030) FRSALZ
    !
    ! ARE REMAINING GROUP 1 PARAMETERS SPATIALLY VARYING?
    !
         Read (ZOO, 1040) SPVARSZ, PRINTSZ
         If (SPVARSZ == 'CONSTANT') Then
            Read (ZOO, 1030) RMAXSZ (1, 1), ESZ (1, 1), RFSZ (1, 1), &
           & BMRSZ (1, 1), BPRSZ (1, 1)
            Do K = 1, KBM1
               Do I = 1, MTLOC
                  RMAXSZ (I, K) = RMAXSZ (1, 1)
                  ESZ (I, K) = ESZ (1, 1)
                  RFSZ (I, K) = RFSZ (1, 1)
                  BMRSZ (I, K) = BMRSZ (1, 1)
                  BPRSZ (I, K) = BPRSZ (1, 1)
               End Do
            End Do
         Else
            Allocate (RTMP21(MGL, KBM1))
            RTMP21 = 0.0
            Allocate (RTMP22(MGL, KBM1))
            RTMP22 = 0.0
            Allocate (RTMP23(MGL, KBM1))
            RTMP23 = 0.0
            Allocate (RTMP24(MGL, KBM1))
            RTMP24 = 0.0
            Allocate (RTMP25(MGL, KBM1))
            RTMP25 = 0.0
            Do K = 1, KBM1
               Do I = 1, MGL
                  Read (ZOO, 1030) RTMP21 (I, K), RTMP22 (I, K), RTMP23 &
                 & (I, K), RTMP24 (I, K), RTMP25 (I, K)
               End Do
            End Do
            If (SERIAL) Then
               RMAXSZ = RTMP21
               ESZ = RTMP22
               RFSZ = RTMP23
               BMRSZ = RTMP24
               BPRSZ = RTMP25
            End If
#if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do K = 1, KBM1
                  Do I = 1, MLOC
                     RMAXSZ (I, K) = RTMP21 (NGID(I), K)
                     ESZ (I, K) = RTMP22 (NGID(I), K)
                     RFSZ (I, K) = RTMP23 (NGID(I), K)
                     BMRSZ (I, K) = RTMP24 (NGID(I), K)
                     BPRSZ (I, K) = RTMP25 (NGID(I), K)
                  End Do
                  Do I = 1, NHN
                     RMAXSZ (I+MLOC, K) = RTMP21 (HN_LST(I), K)
                     ESZ (I+MLOC, K) = RTMP22 (HN_LST(I), K)
                     RFSZ (I+MLOC, K) = RTMP23 (HN_LST(I), K)
                     BMRSZ (I+MLOC, K) = RTMP24 (HN_LST(I), K)
                     BPRSZ (I+MLOC, K) = RTMP25 (HN_LST(I), K)
                  End Do
               End Do
            End If
#endif
            Deallocate (RTMP21, RTMP22, RTMP23, RTMP24, RTMP25)
         End If
    !
    ! ARE REMAINING GROUP 2 PARAMETERS SPATIALLY VARYING?
    !
         Read (ZOO, 1040) SPVARLZ, PRINTLZ
         If (SPVARLZ == 'CONSTANT') Then
            Read (ZOO, 1030) RMAXLZ (1, 1), ELZ (1, 1), RFLZ (1, 1), &
           & BMRLZ (1, 1), BPRLZ (1, 1)
            Do K = 1, KBM1
               Do I = 1, MLOC
                  RMAXLZ (I, K) = RMAXLZ (1, 1)
                  ELZ (I, K) = ELZ (1, 1)
                  RFLZ (I, K) = RFLZ (1, 1)
                  BMRLZ (I, K) = BMRLZ (1, 1)
                  BPRLZ (I, K) = BPRLZ (1, 1)
               End Do
            End Do
         Else
            Allocate (RTMP21(MGL, KBM1))
            RTMP21 = 0.0
            Allocate (RTMP22(MGL, KBM1))
            RTMP22 = 0.0
            Allocate (RTMP23(MGL, KBM1))
            RTMP23 = 0.0
            Allocate (RTMP24(MGL, KBM1))
            RTMP24 = 0.0
            Allocate (RTMP25(MGL, KBM1))
            RTMP25 = 0.0
            Do K = 1, KBM1
               Do I = 1, MGL
                  Read (ZOO, 1030) RTMP21 (I, K), RTMP22 (I, K), RTMP23 &
                 & (I, K), RTMP24 (I, K), RTMP25 (I, K)
               End Do
            End Do
            If (SERIAL) Then
               RMAXLZ = RTMP21
               ELZ = RTMP22
               RFLZ = RTMP23
               BMRLZ = RTMP24
               BPRLZ = RTMP25
            End If
#if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do K = 1, KBM1
                  Do I = 1, MLOC
                     RMAXLZ (I, K) = RTMP21 (NGID(I), K)
                     ELZ (I, K) = RTMP22 (NGID(I), K)
                     RFLZ (I, K) = RTMP23 (NGID(I), K)
                     BMRLZ (I, K) = RTMP24 (NGID(I), K)
                     BPRLZ (I, K) = RTMP25 (NGID(I), K)
                  End Do
                  Do I = 1, NHN
                     RMAXLZ (I+MLOC, K) = RTMP21 (HN_LST(I), K)
                     ELZ (I+MLOC, K) = RTMP22 (HN_LST(I), K)
                     RFLZ (I+MLOC, K) = RTMP23 (HN_LST(I), K)
                     BMRLZ (I+MLOC, K) = RTMP24 (HN_LST(I), K)
                     BPRLZ (I+MLOC, K) = RTMP25 (HN_LST(I), K)
                  End Do
               End Do
            End If
#endif
            Deallocate (RTMP21, RTMP22, RTMP23, RTMP24, RTMP25)
         End If
    !
    ! TIME DEPENDENCE OF PREDATION ON GROUP 1
    !
         Read (ZOO, 1040) TVARSZ, TPRINTSZ
         If (TVARSZ == 'CONSTANT') Then
            Read (ZOO, 1060) TVPRSZ (1)
            Do J = 2, 366
               TVPRSZ (J) = TVPRSZ (1)
            End Do
         Else
            Read (ZOO, 1060) (TVPRSZ(J), J=1, 366)
         End If
    !
    ! TIME DEPENDENCE OF PREDATION ON GROUP 2
    !
         Read (ZOO, 1040) TVARLZ, TPRINTLZ
         If (TVARLZ == 'CONSTANT') Then
            Read (ZOO, 1060) TVPRLZ (1)
            Do J = 2, 366
               TVPRLZ (J) = TVPRLZ (1)
            End Do
         Else
            Read (ZOO, 1060) (TVPRLZ(J), J=1, 366)
         End If
    !
    ! CREATE LOOK-UP TABLE FOR TEMPERATURE EFFECTS
    !
         Do I = - 50, 350
            TLOOK = FLOAT (I) / 10.
            If (TLOOK < TMSZ) Then
               FTSZ (I) = Exp (-KTGSZ1*(TLOOK-TMSZ)**2)
            Else
               FTSZ (I) = Exp (-KTGSZ2*(TMSZ-TLOOK)**2)
            End If
       !       IF (T(B) < TMLZ) THEN
       !         FTLZ = EXP(-KTGLZ1*(T(B)-TMLZ)**2)
       !       ELSE
       !         FTLZ = EXP(-KTGLZ2*(TMLZ-T(B))**2)
       !       ENDIF
            FTBMSZ (I) = Exp (KTBSZ*(TLOOK-TRSZ))
            FTBMLZ (I) = Exp (KTBLZ*(TLOOK-TRLZ))
            FTPRSZ (I) = Exp (KTPRSZ*(TLOOK-TRPRSZ))
            FTPRLZ (I) = Exp (KTPRLZ*(TLOOK-TRPRLZ))
         End Do
    !
    ! READ BI-MODAL TEMPERATURE FUNCTION FOR GROUP TWO
    !
         Read (ZOO, 1050) (FTB_ZOOP(I), I=0, 35)
         Do I = - 50, 0
            FTLZ (I) = 0.
         End Do
         Do I = 0, 34
            Do J = 0, 10
               K = 10 * I + J
               FTLZ (K) = FTB_ZOOP (I) + FLOAT (J) * &
              & (FTB_ZOOP(I+1)-FTB_ZOOP(I)) / 10.
            End Do
         End Do
    !
         Close (38)
    !
    !
    !***** Input FORMAT statements
    !
1010     Format (A72)
1030     Format (/ / (8 X, 9 F8.0))
1040     Format (/ / 8 X, 8 A8)
1050     Format (/ / (8 X, F8.0))
1060     Format (/ / (16 X, F8.0))
    !
    !
    ! OUTPUT WHAT WAS INPUT
    !
         If (MSR) OPEN (ZFO, FILE=ZFOFN, STATUS='UNKNOWN')
         If (MSR) WRITE (ZFO, 2002) (TITLE(J), J=1, 6)
    !
    ! WRITE SPATIALLY-INVARIANT PARAMETERS FOR GROUP 1
    !
         If (MSR) WRITE (ZFO, 2000)
         If (MSR) WRITE (ZFO, 2010) CTSZ, KHCSZ, DOCRITSZ, MZEROSZ
         If (MSR) WRITE (ZFO, 2020) ANCSZ, APCSZ, AOCRSZ
         If (MSR) WRITE (ZFO, 2030) UB1SZ, UB2SZ, UB3SZ, ULDSZ, URDSZ, &
        & ULPSZ, URPSZ
         If (MSR) WRITE (ZFO, 2040) TMSZ, KTGSZ1, KTGSZ2, TRSZ, KTBSZ, &
        & TRPRSZ, KTPRSZ
         If (MSR) WRITE (ZFO, 2045)
         If (MSR) WRITE (ZFO, 2050) FLDOCSZ, FRDOCSZ, FLPOCSZ, FRPOCSZ
         If (MSR) WRITE (ZFO, 2060) FNH4SZ, FUREASZ, FLDONSZ, FRDONSZ, &
        & FLPONSZ, FRPONSZ
         If (MSR) WRITE (ZFO, 2070) FPO4SZ, FLDOPSZ, FRDOPSZ, FLPOPSZ, &
        & FRPOPSZ
         If (MSR) WRITE (ZFO, 2080) FRSASZ
    !
    ! WRITE SPATIALLY-VARYING PARAMETERS FOR GROUP 1
    !
         If (SPVARSZ == 'CONSTANT') Then
            If (MSR) WRITE (ZFO, 2085)
         Else
            If (MSR) WRITE (ZFO, 2086)
         End If
         If (MSR) WRITE (ZFO, 2090)
         If (PRINTSZ /= '     ALL') Then
            If (MSR) WRITE (ZFO, 3000) RMAXSZ (1, 1), ESZ (1, 1), RFSZ &
           & (1, 1), BMRSZ (1, 1), BPRSZ (1, 1)
         Else
            Do K = 1, KBM1
               Do I = 1, MLOC
                  If (MSR) WRITE (ZFO, 3010) I, K, RMAXSZ (I, K), ESZ &
                 & (I, K), RFSZ (I, K), BMRSZ (I, K), BPRSZ (I, K)
               End Do
            End Do
         End If
    !
    ! WRITE TEMPORALLY-VARYING PARAMETERS FOR GROUP 1
    !
         If (TVARSZ == 'CONSTANT') Then
            If (MSR) WRITE (ZFO, 2087)
         Else
            If (MSR) WRITE (ZFO, 2088)
         End If
         If (MSR) WRITE (ZFO, 2092)
         If (TPRINTSZ /= '     ALL') Then
            If (MSR) WRITE (ZFO, 3040) (J, TVPRSZ(J), J=1, 1)
         Else
            If (MSR) WRITE (ZFO, 3040) (J, TVPRSZ(J), J=1, 366)
         End If
    !
    ! WRITE SPATIALLY-INVARIANT PARAMETERS FOR GROUP 2
    !
         If (MSR) WRITE (ZFO, 2005)
         If (MSR) WRITE (ZFO, 2010) CTLZ, KHCLZ, DOCRITLZ, MZEROLZ
         If (MSR) WRITE (ZFO, 2020) ANCLZ, APCLZ, AOCRLZ
         If (MSR) WRITE (ZFO, 2035) UB1LZ, UB2LZ, UB3LZ, USZLZ, ULLZ, &
        & URLZ
         If (MSR) WRITE (ZFO, 2040) TMLZ, KTGLZ1, KTGLZ2, TRLZ, KTBLZ, &
        & TRPRLZ, KTPRLZ
         If (MSR) WRITE (ZFO, 2045)
         If (MSR) WRITE (ZFO, 2050) FLDOCLZ, FRDOCLZ, FLPOCLZ, FRPOCLZ
         If (MSR) WRITE (ZFO, 2060) FNH4LZ, FUREALZ, FLDONLZ, FRDONLZ, &
        & FLPONLZ, FRPONLZ
         If (MSR) WRITE (ZFO, 2070) FPO4LZ, FLDOPLZ, FRDOPLZ, FLPOPLZ, &
        & FRPOPLZ
         If (MSR) WRITE (ZFO, 2080) FRSALZ
    !
    ! WRITE SPATIALLY-VARYING PARAMETERS FOR GROUP 2
    !
         If (SPVARLZ == 'CONSTANT') Then
            If (MSR) WRITE (ZFO, 2085)
         Else
            If (MSR) WRITE (ZFO, 2086)
         End If
         If (MSR) WRITE (ZFO, 2090)
         If (PRINTLZ /= '     ALL') Then
            If (MSR) WRITE (ZFO, 3000) RMAXLZ (1, 1), ELZ (1, 1), RFLZ &
           & (1, 1), BMRLZ (1, 1), BPRLZ (1, 1)
         Else
            Do K = 1, KBM1
               Do I = 1, MLOC
                  If (MSR) WRITE (ZFO, 3010) I, K, RMAXLZ (I, K), ELZ &
                 & (I, K), RFLZ (I, K), BMRLZ (I, K), BPRLZ (I, K)
               End Do
            End Do
         End If
    !
    ! WRITE TEMPORALLY-VARYING PARAMETERS FOR GROUP 2
    !
         If (TVARLZ == 'CONSTANT') Then
            If (MSR) WRITE (ZFO, 2087)
         Else
            If (MSR) WRITE (ZFO, 2088)
         End If
         If (MSR) WRITE (ZFO, 2092)
         If (TPRINTLZ /= '     ALL') Then
            If (MSR) WRITE (ZFO, 3040) (J, TVPRLZ(J), J=1, 1)
         Else
            If (MSR) WRITE (ZFO, 3040) (J, TVPRLZ(J), J=1, 366)
         End If
    !
    ! WRITE GROUP 2 TEMPERATURE FUNCTION
    !
         If (MSR) WRITE (ZFO, 3020)
         Do I = 0, 350, 10
            If (MSR) WRITE (ZFO, 3030) I / 10, FTLZ (I)
         End Do
    !
    ! KURT GLAESEMANN - only close file on MSR
         If (MSR) CLOSE (39)
    !
    !
    !***** Output FORMAT statements
    !
2000     Format (/ ' ZOOPLANKTON GROUP 1')
2002     Format (1 X, A72)
2005     Format (/ ' ZOOPLANKTON GROUP 2')
2010     Format (/ ' CARBON THRESHOLD = ', F8.2, ' GM C/M**3' / ' HALF-&
        &SATURATION CONCENTRATION FOR CARBON UPTAKE = ', F8.2, ' GM C/M&
        &**3' / ' DO AT WHICH MORTALITY COMMENCES = ', F8.2, ' GM DO/M*&
        &*3' / ' MORTALITY UNDER ANOXIC CONDITIONS = ', F8.2, ' /DAY')
2020     Format (' NITROGEN TO CARBON RATIO = ', F8.3, ' GM N/GM C'/' P&
        &HOSPHORUS TO CARBON RATIO = ', F8.3, ' GM P/GM C'/' RESPIRATIO&
        &N RATIO = ', F8.3, ' GM DO/GM C')
2030     Format (' UTILIZATION OF ALGAL GROUP 1   ', F8.2/'            &
        &    ALGAL GROUP 2   ', F8.2/'                ALGAL GROUP 3   '&
        & , F8.2/'                LABILE DOC      ', F8.2/'            &
        &    REFRACTORY DOC  ', F8.2/'                LABILE PARTICLES'&
        & , F8.2/'                REFRACTORY PTCLS', F8.2)
2035     Format (' UTILIZATION OF ALGAL GROUP 1   ', F8.2/'            &
        &    ALGAL GROUP 2   ', F8.2/'                ALGAL GROUP 3   '&
        & , F8.2/'                ZOOPL GROUP 1   ', F8.2/'            &
        &    LABILE PARTICLES', F8.2/'                REFRACTORY PTCLS'&
        & , F8.2)
2040     Format (' OPTIMUM PRODUCTION SPECIFIED AT ', F8.2, ' C.'/' KT1&
        & = ', F8.3, ' KT2 = ', F8.3, ' PER DEGREE**2'/' METABOLISM SPE&
        &CIFIED AT ', F8.2, ' C.'/' TEMPERATURE EFFECT = ', F8.3, ' PER&
        & DEGREE'/' PREDATION SPECIFIED AT ', F8.2, ' C.'/' TEMPERATURE&
        & EFFECT = ', F8.3, ' PER DEGREE')
2045     Format (/ ' RECYCLING BY ZOOPLANKTON' / '            DIS INORG&
        &  LAB DISS  REF DISS  LAB PART  REF PART')
2050     Format (' CARBON    ', 10 X, 4 F10.3)
2060     Format (' NITROGEN  ', 6 F10.3)
2070     Format (' PHOSPHORUS', 5 F10.3)
2080     Format (' SILICA    ', F10.3)
2085     Format (/ ' REMAINING PARAMETERS ARE SPATIALLY-INVARIANT')
2086     Format (/ ' REMAINING PARAMETERS VARY SPATIALLY')
2087     Format (/ ' HIGHER-LEVEL PREDATION IS TEMPORALLY-INVARIANT')
2088     Format (/ ' HIGHER-LEVEL PREDATION VARIES TEMPORALLY')
2090     Format (/ '    BOX   RATION   EFNCY  RESPFR   METAB  PRLOSS' / &
        & '           1/DAY                   1/DAY M**3/GM C/' / '    &
        &                                        DAY')
2092     Format (/ '    DAY     TVPR')
3000     Format (8 X, 5 F8.3)
3010     Format (2 I8, 5 F8.3)
3020     Format (/ ' TEMPERATURE FUNCTION FOR GROUP 2' / '    TEMP     &
        & f(T)' /)
3030     Format (I8, F10.4)
3040     Format (I8, F8.3)
    !
         Return
      End Subroutine ZOOP_READ
    !
    !
    !
      Subroutine ZOOP_ALLOC !Added this for allocating zooplankton related variables
      !
         Allocate (B1ASZ(MTLOC, KBM1))
         B1ASZ = 0.0
         Allocate (B2ASZ(MTLOC, KBM1))
         B2ASZ = 0.0
         Allocate (B3ASZ(MTLOC, KBM1))
         B3ASZ = 0.0
         Allocate (LPOCASZ(MTLOC, KBM1))
         LPOCASZ = 0.0
         Allocate (RPOCASZ(MTLOC, KBM1))
         RPOCASZ = 0.0
         Allocate (PRASZ(MTLOC, KBM1))
         PRASZ = 0.0
         Allocate (B1ALZ(MTLOC, KBM1))
         B1ALZ = 0.0
         Allocate (B2ALZ(MTLOC, KBM1))
         B2ALZ = 0.0
         Allocate (B3ALZ(MTLOC, KBM1))
         B3ALZ = 0.0
         Allocate (SZALZ(MTLOC, KBM1))
         SZALZ = 0.0
         Allocate (LPOCALZ(MTLOC, KBM1))
         LPOCALZ = 0.0
         Allocate (RPOCALZ(MTLOC, KBM1))
         RPOCALZ = 0.0
         Allocate (PRALZ(MTLOC, KBM1))
         PRALZ = 0.0
         Allocate (CLSZ(MTLOC, KBM1))
         CLSZ = 0.0
         Allocate (CLLZ(MTLOC, KBM1))
         CLLZ = 0.0
         Allocate (RSZ(MTLOC, KBM1))
         RSZ = 0.0
         Allocate (RLZ(MTLOC, KBM1))
         RLZ = 0.0
         Allocate (RMAXSZ(MTLOC, KBM1))
         RMAXSZ = 0.0
         Allocate (RMAXLZ(MTLOC, KBM1))
         RMAXLZ = 0.0
         Allocate (BMSZ(MTLOC, KBM1))
         BMSZ = 0.0
         Allocate (BMLZ(MTLOC, KBM1))
         BMLZ = 0.0
         Allocate (BMRSZ(MTLOC, KBM1))
         BMRSZ = 0.0
         Allocate (BMRLZ(MTLOC, KBM1))
         BMRLZ = 0.0
         Allocate (MSZ(MTLOC, KBM1))
         MSZ = 0.0
         Allocate (MLZ(MTLOC, KBM1))
         MLZ = 0.0
         Allocate (PRSZLZ(MTLOC, KBM1))
         PRSZLZ = 0.0
         Allocate (GSZ(MTLOC, KBM1))
         GSZ = 0.0
         Allocate (GLZ(MTLOC, KBM1))
         GLZ = 0.0
         Allocate (ESZ(MTLOC, KBM1))
         ESZ = 0.0
         Allocate (ELZ(MTLOC, KBM1))
         ELZ = 0.0
         Allocate (RFSZ(MTLOC, KBM1))
         RFSZ = 0.0
         Allocate (RFLZ(MTLOC, KBM1))
         RFLZ = 0.0
         Allocate (PRSZ(MTLOC, KBM1))
         PRSZ = 0.0
         Allocate (PRLZ(MTLOC, KBM1))
         PRLZ = 0.0
         Allocate (LDOCASZ(MTLOC, KBM1))
         LDOCASZ = 0.0
         Allocate (BPRSZ(MTLOC, KBM1))
         BPRSZ = 0.0
         Allocate (BPRLZ(MTLOC, KBM1))
         BPRLZ = 0.0
         Allocate (RDOCASZ(MTLOC, KBM1))
         RDOCASZ = 0.0
      !
         Allocate (LDOCSZ(MTLOC, KBM1))
         LDOCSZ = 0.0
         Allocate (LPOCSZ(MTLOC, KBM1))
         LPOCSZ = 0.0
         Allocate (RPOCSZ(MTLOC, KBM1))
         RPOCSZ = 0.0
         Allocate (LDOCLZ(MTLOC, KBM1))
         LDOCLZ = 0.0
         Allocate (LPOCLZ(MTLOC, KBM1))
         LPOCLZ = 0.0
         Allocate (RPOCLZ(MTLOC, KBM1))
         RPOCLZ = 0.0
         Allocate (NH4SZ(MTLOC, KBM1))
         NH4SZ = 0.0
         Allocate (LDONSZ(MTLOC, KBM1))
         LDONSZ = 0.0
         Allocate (LPONSZ(MTLOC, KBM1))
         LPONSZ = 0.0
         Allocate (RPONSZ(MTLOC, KBM1))
         RPONSZ = 0.0
         Allocate (NH4LZ(MTLOC, KBM1))
         NH4LZ = 0.0
         Allocate (LDONLZ(MTLOC, KBM1))
         LDONLZ = 0.0
         Allocate (LPONLZ(MTLOC, KBM1))
         LPONLZ = 0.0
         Allocate (RPONLZ(MTLOC, KBM1))
         RPONLZ = 0.0
         Allocate (PO4SZ(MTLOC, KBM1))
         PO4SZ = 0.0
         Allocate (LDOPSZ(MTLOC, KBM1))
         LDOPSZ = 0.0
         Allocate (LPOPSZ(MTLOC, KBM1))
         LPOPSZ = 0.0
         Allocate (RPOPSZ(MTLOC, KBM1))
         RPOPSZ = 0.0
         Allocate (PO4LZ(MTLOC, KBM1))
         PO4LZ = 0.0
         Allocate (LDOPLZ(MTLOC, KBM1))
         LDOPLZ = 0.0
         Allocate (LPOPLZ(MTLOC, KBM1))
         LPOPLZ = 0.0
         Allocate (RPOPLZ(MTLOC, KBM1))
         RPOPLZ = 0.0
         Allocate (RDOCSZ(MTLOC, KBM1))
         RDOCSZ = 0.0
         Allocate (RDONSZ(MTLOC, KBM1))
         RDONSZ = 0.0
         Allocate (RDOPSZ(MTLOC, KBM1))
         RDOPSZ = 0.0
         Allocate (RDOCLZ(MTLOC, KBM1))
         RDOCLZ = 0.0
         Allocate (RDONLZ(MTLOC, KBM1))
         RDONLZ = 0.0
         Allocate (RDOPLZ(MTLOC, KBM1))
         RDOPLZ = 0.0
         Allocate (PIB1SZ(MTLOC, KBM1))
         PIB1SZ = 0.0
         Allocate (PIB2SZ(MTLOC, KBM1))
         PIB2SZ = 0.0
         Allocate (PIB3SZ(MTLOC, KBM1))
         PIB3SZ = 0.0
         Allocate (PIB1LZ(MTLOC, KBM1))
         PIB1LZ = 0.0
         Allocate (PIB2LZ(MTLOC, KBM1))
         PIB2LZ = 0.0
         Allocate (PIB3LZ(MTLOC, KBM1))
         PIB3LZ = 0.0
      !
         Allocate (B1SZ(MTLOC, KBM1))
         B1SZ = 0.0
         Allocate (B2SZ(MTLOC, KBM1))
         B2SZ = 0.0
         Allocate (B3SZ(MTLOC, KBM1))
         B3SZ = 0.0
         Allocate (B1LZ(MTLOC, KBM1))
         B1LZ = 0.0
         Allocate (B2LZ(MTLOC, KBM1))
         B2LZ = 0.0
         Allocate (B3LZ(MTLOC, KBM1))
         B3LZ = 0.0
         Allocate (DOSZ(MTLOC, KBM1))
         DOSZ = 0.0
         Allocate (DOLZ(MTLOC, KBM1))
         DOLZ = 0.0
         Allocate (SASZ(MTLOC, KBM1))
         SASZ = 0.0
         Allocate (SUSZ(MTLOC, KBM1))
         SUSZ = 0.0
         Allocate (SALZ(MTLOC, KBM1))
         SALZ = 0.0
         Allocate (SULZ(MTLOC, KBM1))
         SULZ = 0.0
      !
      !Moved the following from mod_wqm.F
         Allocate (ACLSZ(MTLOC, KBM1))
         ACLSZ = 0.0
         Allocate (ACLLZ(MTLOC, KBM1))
         ACLLZ = 0.0
         Allocate (ARSZ(MTLOC, KBM1))
         ARSZ = 0.0
         Allocate (ARLZ(MTLOC, KBM1))
         ARLZ = 0.0
         Allocate (ABMSZ(MTLOC, KBM1))
         ABMSZ = 0.0
         Allocate (ABMLZ(MTLOC, KBM1))
         ABMLZ = 0.0
         Allocate (AMSZ(MTLOC, KBM1))
         AMSZ = 0.0
         Allocate (AMLZ(MTLOC, KBM1))
         AMLZ = 0.0
         Allocate (APRSZLZ(MTLOC, KBM1))
         APRSZLZ = 0.0
         Allocate (AGSZ(MTLOC, KBM1))
         AGSZ = 0.0
         Allocate (AGLZ(MTLOC, KBM1))
         AGLZ = 0.0
         Allocate (ADOCSZ(MTLOC, KBM1))
         ADOCSZ = 0.0
         Allocate (APOCSZ(MTLOC, KBM1))
         APOCSZ = 0.0
         Allocate (ADOCLZ(MTLOC, KBM1))
         ADOCLZ = 0.0
         Allocate (APOCLZ(MTLOC, KBM1))
         APOCLZ = 0.0
         Allocate (ANH4SZ(MTLOC, KBM1))
         ANH4SZ = 0.0
         Allocate (ADONSZ(MTLOC, KBM1))
         ADONSZ = 0.0
         Allocate (APONSZ(MTLOC, KBM1))
         APONSZ = 0.0
         Allocate (ANH4LZ(MTLOC, KBM1))
         ANH4LZ = 0.0
         Allocate (ADONLZ(MTLOC, KBM1))
         ADONLZ = 0.0
         Allocate (APONLZ(MTLOC, KBM1))
         APONLZ = 0.0
         Allocate (APO4SZ(MTLOC, KBM1))
         APO4SZ = 0.0
         Allocate (ADOPSZ(MTLOC, KBM1))
         ADOPSZ = 0.0
         Allocate (APOPSZ(MTLOC, KBM1))
         APOPSZ = 0.0
         Allocate (APO4LZ(MTLOC, KBM1))
         APO4LZ = 0.0
         Allocate (ADOPLZ(MTLOC, KBM1))
         ADOPLZ = 0.0
         Allocate (APOPLZ(MTLOC, KBM1))
         APOPLZ = 0.0
         Allocate (APRSZ(MTLOC, KBM1))
         APRSZ = 0.0
         Allocate (APRLZ(MTLOC, KBM1))
         APRLZ = 0.0
         Allocate (APISZ(MTLOC, KBM1))
         APISZ = 0.0
         Allocate (APILZ(MTLOC, KBM1))
         APILZ = 0.0
         Allocate (AB1SZ(MTLOC, KBM1))
         AB1SZ = 0.0
         Allocate (AB2SZ(MTLOC, KBM1))
         AB2SZ = 0.0
         Allocate (AB3SZ(MTLOC, KBM1))
         AB3SZ = 0.0
         Allocate (AB1LZ(MTLOC, KBM1))
         AB1LZ = 0.0
         Allocate (AB2LZ(MTLOC, KBM1))
         AB2LZ = 0.0
         Allocate (AB3LZ(MTLOC, KBM1))
         AB3LZ = 0.0
         Allocate (ADOSZ(MTLOC, KBM1))
         ADOSZ = 0.0
         Allocate (ADOLZ(MTLOC, KBM1))
         ADOLZ = 0.0
         Allocate (ASASZ(MTLOC, KBM1))
         ASASZ = 0.0
         Allocate (ASUSZ(MTLOC, KBM1))
         ASUSZ = 0.0
         Allocate (ASALZ(MTLOC, KBM1))
         ASALZ = 0.0
         Allocate (ASULZ(MTLOC, KBM1))
         ASULZ = 0.0
      !
      !Moved here from mod_wqm.F
         Allocate (SZ(0:MTLOC, KBM1))
         SZ = 0.0
         Allocate (LZ(0:MTLOC, KBM1))
         LZ = 0.0
      !
      End Subroutine ZOOP_ALLOC
    !
    !
      Subroutine ZOOP_DEALLOC !Added this for de-allocating zooplankton related variables
      !
         If (ALLOCATED(B1ASZ)) DEALLOCATE (B1ASZ)
         If (ALLOCATED(B2ASZ)) DEALLOCATE (B2ASZ)
         If (ALLOCATED(B3ASZ)) DEALLOCATE (B3ASZ)
         If (ALLOCATED(LPOCASZ)) DEALLOCATE (LPOCASZ)
         If (ALLOCATED(RPOCASZ)) DEALLOCATE (RPOCASZ)
         If (ALLOCATED(PRASZ)) DEALLOCATE (PRASZ)
         If (ALLOCATED(B1ALZ)) DEALLOCATE (B1ALZ)
         If (ALLOCATED(B2ALZ)) DEALLOCATE (B2ALZ)
         If (ALLOCATED(B3ALZ)) DEALLOCATE (B3ALZ)
         If (ALLOCATED(SZALZ)) DEALLOCATE (SZALZ)
         If (ALLOCATED(LPOCALZ)) DEALLOCATE (LPOCALZ)
         If (ALLOCATED(RPOCALZ)) DEALLOCATE (RPOCALZ)
         If (ALLOCATED(PRALZ)) DEALLOCATE (PRALZ)
         If (ALLOCATED(CLSZ)) DEALLOCATE (CLSZ)
         If (ALLOCATED(CLLZ)) DEALLOCATE (CLLZ)
         If (ALLOCATED(RSZ)) DEALLOCATE (RSZ)
         If (ALLOCATED(RLZ)) DEALLOCATE (RLZ)
         If (ALLOCATED(RMAXSZ)) DEALLOCATE (RMAXSZ)
         If (ALLOCATED(RMAXLZ)) DEALLOCATE (RMAXLZ)
         If (ALLOCATED(BMSZ)) DEALLOCATE (BMSZ)
         If (ALLOCATED(BMLZ)) DEALLOCATE (BMLZ)
         If (ALLOCATED(BMRSZ)) DEALLOCATE (BMRSZ)
         If (ALLOCATED(BMRLZ)) DEALLOCATE (BMRLZ)
         If (ALLOCATED(MSZ)) DEALLOCATE (MSZ)
         If (ALLOCATED(MLZ)) DEALLOCATE (MLZ)
         If (ALLOCATED(PRSZLZ)) DEALLOCATE (PRSZLZ)
      !
         If (ALLOCATED(GSZ)) DEALLOCATE (GSZ)
         If (ALLOCATED(GLZ)) DEALLOCATE (GLZ)
         If (ALLOCATED(ESZ)) DEALLOCATE (ESZ)
         If (ALLOCATED(ELZ)) DEALLOCATE (ELZ)
         If (ALLOCATED(RFSZ)) DEALLOCATE (RFSZ)
         If (ALLOCATED(RFLZ)) DEALLOCATE (RFLZ)
         If (ALLOCATED(PRSZ)) DEALLOCATE (PRSZ)
         If (ALLOCATED(PRLZ)) DEALLOCATE (PRLZ)
         If (ALLOCATED(LDOCASZ)) DEALLOCATE (LDOCASZ)
         If (ALLOCATED(BPRSZ)) DEALLOCATE (BPRSZ)
         If (ALLOCATED(BPRLZ)) DEALLOCATE (BPRLZ)
         If (ALLOCATED(RDOCASZ)) DEALLOCATE (RDOCASZ)
      !
         If (ALLOCATED(LDOCSZ)) DEALLOCATE (LDOCSZ)
         If (ALLOCATED(LPOCSZ)) DEALLOCATE (LPOCSZ)
         If (ALLOCATED(RPOCSZ)) DEALLOCATE (RPOCSZ)
         If (ALLOCATED(LDOCLZ)) DEALLOCATE (LDOCLZ)
         If (ALLOCATED(LPOCLZ)) DEALLOCATE (LPOCLZ)
         If (ALLOCATED(RPOCLZ)) DEALLOCATE (RPOCLZ)
         If (ALLOCATED(NH4SZ)) DEALLOCATE (NH4SZ)
         If (ALLOCATED(LDONSZ)) DEALLOCATE (LDONSZ)
         If (ALLOCATED(LPONSZ)) DEALLOCATE (LPONSZ)
         If (ALLOCATED(RPONSZ)) DEALLOCATE (RPONSZ)
         If (ALLOCATED(NH4LZ)) DEALLOCATE (NH4LZ)
         If (ALLOCATED(LDONLZ)) DEALLOCATE (LDONLZ)
         If (ALLOCATED(LPONLZ)) DEALLOCATE (LPONLZ)
         If (ALLOCATED(RPONLZ)) DEALLOCATE (RPONLZ)
         If (ALLOCATED(PO4SZ)) DEALLOCATE (PO4SZ)
         If (ALLOCATED(LDOPSZ)) DEALLOCATE (LDOPSZ)
         If (ALLOCATED(LPOPSZ)) DEALLOCATE (LPOPSZ)
         If (ALLOCATED(RPOPSZ)) DEALLOCATE (RPOPSZ)
         If (ALLOCATED(PO4LZ)) DEALLOCATE (PO4LZ)
         If (ALLOCATED(LDOPLZ)) DEALLOCATE (LDOPLZ)
         If (ALLOCATED(LPOPLZ)) DEALLOCATE (LPOPLZ)
         If (ALLOCATED(RPOPLZ)) DEALLOCATE (RPOPLZ)
         If (ALLOCATED(RDOCSZ)) DEALLOCATE (RDOCSZ)
         If (ALLOCATED(RDONSZ)) DEALLOCATE (RDONSZ)
         If (ALLOCATED(RDOPSZ)) DEALLOCATE (RDOPSZ)
         If (ALLOCATED(RDOCLZ)) DEALLOCATE (RDOCLZ)
         If (ALLOCATED(RDONLZ)) DEALLOCATE (RDONLZ)
         If (ALLOCATED(RDOPLZ)) DEALLOCATE (RDOPLZ)
         If (ALLOCATED(PIB1SZ)) DEALLOCATE (PIB1SZ)
         If (ALLOCATED(PIB2SZ)) DEALLOCATE (PIB2SZ)
         If (ALLOCATED(PIB3SZ)) DEALLOCATE (PIB3SZ)
         If (ALLOCATED(PIB1LZ)) DEALLOCATE (PIB1LZ)
         If (ALLOCATED(PIB2LZ)) DEALLOCATE (PIB2LZ)
         If (ALLOCATED(PIB3LZ)) DEALLOCATE (PIB3LZ)
      !
      !Moved here from wqm_main.F
      !
         If (ALLOCATED(B1SZ)) DEALLOCATE (B1SZ)
         If (ALLOCATED(B2SZ)) DEALLOCATE (B2SZ)
         If (ALLOCATED(B3SZ)) DEALLOCATE (B3SZ)
         If (ALLOCATED(B1LZ)) DEALLOCATE (B1LZ)
         If (ALLOCATED(B2LZ)) DEALLOCATE (B2LZ)
         If (ALLOCATED(B3LZ)) DEALLOCATE (B3LZ)
         If (ALLOCATED(DOSZ)) DEALLOCATE (DOSZ)
         If (ALLOCATED(DOLZ)) DEALLOCATE (DOLZ)
         If (ALLOCATED(SASZ)) DEALLOCATE (SASZ)
         If (ALLOCATED(SUSZ)) DEALLOCATE (SUSZ)
         If (ALLOCATED(SALZ)) DEALLOCATE (SALZ)
         If (ALLOCATED(SULZ)) DEALLOCATE (SULZ)
      !
         If (ALLOCATED(ACLSZ)) DEALLOCATE (ACLSZ)
         If (ALLOCATED(ACLLZ)) DEALLOCATE (ACLLZ)
         If (ALLOCATED(ARSZ)) DEALLOCATE (ARSZ)
         If (ALLOCATED(ARLZ)) DEALLOCATE (ARLZ)
         If (ALLOCATED(ABMSZ)) DEALLOCATE (ABMSZ)
         If (ALLOCATED(ABMLZ)) DEALLOCATE (ABMLZ)
         If (ALLOCATED(AMSZ)) DEALLOCATE (AMSZ)
         If (ALLOCATED(AMLZ)) DEALLOCATE (AMLZ)
         If (ALLOCATED(APRSZLZ)) DEALLOCATE (APRSZLZ)
         If (ALLOCATED(AGSZ)) DEALLOCATE (AGSZ)
         If (ALLOCATED(AGLZ)) DEALLOCATE (AGLZ)
         If (ALLOCATED(ADOCSZ)) DEALLOCATE (ADOCSZ)
         If (ALLOCATED(APOCSZ)) DEALLOCATE (APOCSZ)
         If (ALLOCATED(ADOCLZ)) DEALLOCATE (ADOCLZ)
         If (ALLOCATED(APOCLZ)) DEALLOCATE (APOCLZ)
         If (ALLOCATED(ANH4SZ)) DEALLOCATE (ANH4SZ)
         If (ALLOCATED(ADONSZ)) DEALLOCATE (ADONSZ)
         If (ALLOCATED(APONSZ)) DEALLOCATE (APONSZ)
         If (ALLOCATED(ANH4LZ)) DEALLOCATE (ANH4LZ)
         If (ALLOCATED(ADONLZ)) DEALLOCATE (ADONLZ)
         If (ALLOCATED(APONLZ)) DEALLOCATE (APONLZ)
         If (ALLOCATED(APO4SZ)) DEALLOCATE (APO4SZ)
         If (ALLOCATED(ADOPSZ)) DEALLOCATE (ADOPSZ)
         If (ALLOCATED(APOPSZ)) DEALLOCATE (APOPSZ)
         If (ALLOCATED(APO4LZ)) DEALLOCATE (APO4LZ)
         If (ALLOCATED(ADOPLZ)) DEALLOCATE (ADOPLZ)
         If (ALLOCATED(APOPLZ)) DEALLOCATE (APOPLZ)
         If (ALLOCATED(APRSZ)) DEALLOCATE (APRSZ)
         If (ALLOCATED(APRLZ)) DEALLOCATE (APRLZ)
         If (ALLOCATED(APISZ)) DEALLOCATE (APISZ)
         If (ALLOCATED(APILZ)) DEALLOCATE (APILZ)
         If (ALLOCATED(AB1SZ)) DEALLOCATE (AB1SZ)
         If (ALLOCATED(AB2SZ)) DEALLOCATE (AB2SZ)
         If (ALLOCATED(AB3SZ)) DEALLOCATE (AB3SZ)
         If (ALLOCATED(AB1LZ)) DEALLOCATE (AB1LZ)
         If (ALLOCATED(AB2LZ)) DEALLOCATE (AB2LZ)
         If (ALLOCATED(AB3LZ)) DEALLOCATE (AB3LZ)
         If (ALLOCATED(ADOSZ)) DEALLOCATE (ADOSZ)
         If (ALLOCATED(ADOLZ)) DEALLOCATE (ADOLZ)
         If (ALLOCATED(ASASZ)) DEALLOCATE (ASASZ)
         If (ALLOCATED(ASUSZ)) DEALLOCATE (ASUSZ)
         If (ALLOCATED(ASALZ)) DEALLOCATE (ASALZ)
         If (ALLOCATED(ASULZ)) DEALLOCATE (ASULZ)
      !
      !Moved this from wqm_main.F
         If (ALLOCATED(SZ)) DEALLOCATE (SZ)
         If (ALLOCATED(LZ)) DEALLOCATE (LZ)
      !
      End Subroutine ZOOP_DEALLOC


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

      Function GET_FTPR (TEMPVAL, TRPR, KTPR)!m/d
         Real (SP) :: TEMPVAL, TRPR, KTPR !temperature (degC)
         Real (SP) :: GET_FTPR
         GET_FTPR = Exp (KTPR*(TEMPVAL-TRPR))
         Return
      End Function GET_FTPR

      Function GET_FTBM (TEMPVAL, TR, KTB)!m/d
         Real (SP) :: TEMPVAL, TR, KTB !temperature (degC)
         Real (SP) :: GET_FTBM
         GET_FTBM = Exp (KTB*(TEMPVAL-TR))
         Return
      End Function GET_FTBM

End Module MOD_ZOOP
