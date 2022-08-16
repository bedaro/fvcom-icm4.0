!tvds.F
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
!Subroutine TVDS()  !Time varying data series
!
!************************************************************************
!**                  S U B R O U T I N E   T V D S                     **
!************************************************************************
!
Subroutine TVDS (NXTVD)
  !
      Use MOD_PREC, Only: SP
      Use MOD_LIMS, Only: MLOC, MTLOC !Added for reading time-varying pCO2atm 12feb2016
  !
      Use MOD_SIZES, Only: MGL, NCP !,       &!
  !
      Use MOD_WQM, Only: BB, BENDOC, BENTHIC_FLUXES, BFLUX, BFLUXNX, &
     & DIAGNOSTICS, FD, FILGTH, IT, ITNX, JDAY, KT, METFN, METPTR, &
     & NXBFI, NXMET, TE, TMEND, TTSS, WMS, BENCOD, BENCODB, BENDO, &
     & BENDOB, BENDOCB, BENDON, BENDOP, BENNH4, BENNH4B, BENNO3, &
     & BENNO3B, BENPO4, BENPO4B, BENSA, BENSAB, BFIFN, BFIPTR, &
     & BFLUXNX_GL, KSDOC, ATMFN, ATMLDON, ATMLDOP, ATMNH4, ATMNO3, &
     & ATMOS_LOADS, ATMPO4, ATMPTR, ATMRDON, ATMRDOP, KHSO, KSNH4, &
     & KSNO3, KSO, KSPO4, KSSA, MTCNO3, NXATM, PRECIP, SAV_LOADS, &
     & SEDNO3, TRSDOC, TRSNH4, TRSNO3, TRSO, TRSPO4, TRSSA, B, NBB, &
     & NS1P, NS2P, NS3P, NXSAV, SAVPTR, SVIFN !,          &!
  !
      Use MOD_HYDROVARS, Only: ART1
      Use MOD_CONTROL, Only: SERIAL, PAR
  !
      Use MOD_FILEINFO, Only: DIA, BFI, BAI, MET, ATM, SVI,AIRC
  !
      Use MOD_SAV, Only: SPNAM, EPINAM, KHSTOX, KESAV, KEEPI, &
     & KHEP, NLIMEPI, KHNLEAF, KHNROOT, KHPLEAF, KHPROOT, KHNEPI, &
     & KHPEPI, KHNPSAV, KHNPEPI, PMSAV, PRSPSAVB, BMSAV, BMTBRREF, &
     & FDOSR, SALMAX, ALPHSAV, ALAC, SLSAV, ACDWSAV, ANDWSAV, APDWSAV, &
     & ANCSAV, APCSAV, HCAN, ACAN, BCAN, FNISAV, FNLDSAV, FNRDSAV, &
     & FNLPSAV, FNRPSAV, FPISAV, FPLDSAV, FPRDSAV, FPLPSAV, FPRPSAV, &
     & FCLDSAV, FCRDSAV, FCLPSAV, FCRPSAV, FDOSAV, WSSSAV, WSLSAV, &
     & WSRSAV, WS1SAV, WS2SAV, WS3SAV, WSUSAV, NDDI, NSAVCELL, &
     & NSAVSPCM, PMEPI, BMEPI, PREPI, PRSPEPI, ALPHEPI, CCHLEPI, &
     & ANCEPI, APCEPI, ADWCEPI, FCLDEPI, FCRDEPI, FNIEPI, FNLDEPI, &
     & FNRDEPI, FNLPEPI, FNRPEPI, FPIEPI, FPLDEPI, FPRDEPI, FPLPEPI, &
     & FPRPEPI, FCLDPEP, FCRDPEP, FCLPPEP, FCRPPEP, FNIPEP, FNLDPEP, &
     & FNRDPEP, FNLPPEP, FNRPPEP, FPIPEP, FPLDPEP, FPRDPEP, FPLPPEP, &
     & FPRPPEP, FPSR, FPLEAF, FPROOT, FPSTEM, FPTUBER, TRTBRLF, & ! AN FTPSAV,
     & FRPOCSAV, FRPONSAV, FRPOPSAV, & !AN FTRSAV, FTPEP, FTREP, FTPREP,
     & LEAF, ROOT, STEM, TUBER, EP, SAVAREA, SAVCELL, NSAVSPC, SAVDPH, &
     & SAVFRAC, PLEAF, BMLEAF, BMSTEM, BMROOT, BMTUBER, PEP, BMEP, &
     & PREP, SLSH, NLSAV, PLSAV, FNSEDSAV, FPSEDSAV, NLEPI, PLEPI, FHS, &
     & EPATN, SAVATN, WATATN, FISH, FIEP, NPPSAV, NPPEPI, DOSAVW, &
     & LDOCSAVW, RDOCSAVW, LPOCSAVW, RPOCSAVW, NH4SAVW, NO3SAVW, &
     & LDONSAVW, RDONSAVW, LPONSAVW, RPONSAVW, PO4SAVW, LDOPSAVW, &
     & RDOPSAVW, LPOPSAVW, RPOPSAVW, DOEPIW, LDOCEPIW, RDOCEPIW, &
     & LPOCEPIW, RPOCEPIW, NH4EPIW, NO3EPIW, LDONEPIW, RDONEPIW, &
     & LPONEPIW, RPONEPIW, PO4EPIW, LDOPEPIW, RDOPEPIW, LPOPEPIW, &
     & RPOPEPIW, SEDPOCSAV, SEDPONSAV, SEDPOPSAV, SEDDOSAV, SEDNH4SAV, &
     & SEDNO3SAV, SEDPO4SAV, DGRDAYS, ALEAF, AROOT, ASTEM, ATUBER, AEP, &
     & APLEAF, ABMLEAF, ABMTUBER, APEP, ABMEP, APREP, ASLSH, ANLSAV, &
     & APLSAV, ANLEPI, APLEPI, AFNSED, AFPSED, AFHS, AEPATN, AWATATN, &
     & AFISH, AFIEP, ANPPSAV, ANPPEPI, ADOSAVW, ADOCSAVW, APOCSAVW, &
     & ANH4SAVW, ANO3SAVW, ADONSAVW, APONSAVW, APO4SAVW, ADOPSAVW, &
     & APOPSAVW, ADOEPIW, ANH4EPIW, ANO3EPIW, APO4EPIW, ADOCEPIW, &
     & APOCEPIW, ADONEPIW, APONEPIW, ADOPEPIW, APOPEPIW, ASEDDOSAV, &
     & ASEDPOCSAV, ASEDPONSAV, ASEDPOPSAV, ASEDNH4SAV, ASEDPO4SAV
#   if defined (MULTIPROCESSOR)
      Use MOD_PAR, Only: NGID, NHN, HN_LST
#   endif
!
      Use MOD_CO2SYS, Only: MAPCO2ATM, pco2atmNX, NXPCO2, pCO2atm !Added to read time-varying pCO2atm 12feb2016
      Use MOD_CONTROL, Only: MSR 
!
      Implicit None
      Save
      Integer :: Max, I, JCON !Index for constituents
  !
  !    REAL(SP) :: FACTOR,FDNX,TENX,WMSNX,BFLUX,BFLUXNX,BFLUXNX_GL,PRECIPNX, &
  !     ANH4NX,ANO3NX,ALDONNX,ARDONNX,APO4NX,ALDOPNX,ARDOPNX,SDOCNX,  &
  !     SLPOCNX,SRPOCNX,SDONX
  !
      Real (SP) :: FACTOR, FDNX, TENX, WMSNX, PRECIPNX, ANH4NX, ANO3NX, &
     & ALDONNX, ARDONNX, APO4NX, ALDOPNX, ARDOPNX, SLDOCNX, SRDOCNX, &
     & SLPOCNX, SRPOCNX, SDONX !DO flux to water column due to SAV    (gC/m^2/day)
      Real (SP) :: S2LNX, KTNX, NXTVD, S1LNX, S3LNX
!
  !
  !Moved BFLUXNX, BFLUXNX_GL, BFLUX into wqm_module.f module name MOD_WQM!!!
  !
      Dimension S1LNX (NS1P, NCP), S2LNX (NS2P, NCP), S3LNX (NS3P, NCP)
  !
  !     DIMENSION S1LNX(NS1P,NCP), S2LNX(NS2P,NCP), S3LNX(NS3P,NCP),    &
  !               BFLUXNX(MTLOC,NCP),BFLUXNX_GL(MGL,NCP),BFLUX(MTLOC,NCP)
  !
  !******* Meteorologic data
  !
	  NXTVD = TMEND
      !FDNX = 1 !! commented out it bugs!
      !TENX = 0 !! commented out it bugs!
	  !WMSNX = 0 ! AN
10    Continue
      Do WHILE (JDAY >= NXMET)

	 !   Read (MET, 1010, END=11) NXMET, KTNX, TENX, ITNX, FDNX, WMSNX !  Wind speed (m/s)

         KT = KTNX				!coefficeint of heat exchange (surface)
         FD = FDNX				!fraction of day
     !
     ! CONVERT LANGLEYS TO EINSTEINS/M**2
     ! Commented conversion out since SJRWMD data already in einsteins/m**2
     !         IT = 0.093*ITNX
     ! convert w/m**2 to PAR w/m**2 - need to convert to Ein? - states PAR (E/m2/day)

     !     IT = ITNX*0.143  !If not in E/m2/day, then divide by 0.143 to convert from watt/m^2 to E/m^2/day
         IT = ITNX          !THIS is already in E/m2/day in input file

     ! need to change this
     ! TIME TO SUNRISE (SECONDS OF THE DAY COUNTING FROM 00:00 (MIDNIGHT) )
         TTSS = 86400. * (1-FD) / 2.
     !
	    !If(MSR) Write (*,*) 'FD, JDAY in TVDS= ', FD, JDAY !AN
        !If(MSR) Write (*,*) 'IO input in TVDS= ', IT !AN
	 !
         TE = Max (TENX, 0.0)
         WMS = WMSNX

     !   Move up AN
	     Read (MET, 1010, END=11) NXMET, KTNX, TENX, ITNX, FDNX, WMSNX !  Wind speed (m/s)
     !
         NXMET = (METPTR-1) * FILGTH + NXMET
		 NXMET = AINT(NXMET) !AN: This code fixes the daily estimate I0

		 !If(MSR) Write (*,*) 'JDAY, NXMET, METPTR  inside loop = ', JDAY, NXMET, METPTR !AN

      End Do
      Go To 12
  !
  !******* Open next data file
  !
11    Continue
      METPTR = METPTR + 1
#if defined (KILLTREES)  /*put diagnostics also in killtrees */
      If (DIAGNOSTICS) WRITE (DIA,*) 'Opening meteorologic file ', &
     & METPTR, ' at day ', JDAY
#endif
      Close (MET)
      Open (MET, File=METFN(METPTR), Status='OLD')
      Read (MET, 1000)
  !
  ! INITIALIZE DEGREE DAYS TO ZERO
  ! ong: this is now initialized by ALLOC_SAV() in wqm_sav.F
  !    DO I=1,MLOC
  !      DGRDAYS(I)=0  !this should be moved to wqm_sav.F
  !    ENDDO
  !
      Read (MET, 1010) NXMET, KTNX, TENX, ITNX, FDNX, WMSNX
      NXMET = (METPTR-1) * FILGTH + NXMET
      If (JDAY >= NXMET) Go To 10
12    Continue
      NXTVD = Min (NXTVD, NXMET)
	  !If(MSR) Write (*,*) 'NXTVD in TVDS= ', NXTVD, 'NXMET= ', NXMET !AN
  !
  !******* Boundary inflow concentrations
  !
  !JQI    IF (BOUNDARY_CONC) THEN
  !JQI      DO WHILE(JDAY >= NXCBC)THEN
  !JQI        OLDNXCBC = NXCBC
  !
  !*********** Reduce/increase concentrations
  !
  !JQI        DO JCON=4,12
  !JQI          DO JCB=1,NCB(JCON)
  !JQI            CBNX(JCB,JCON) = CBNX(JCB,JCON)*REDCBC
  !JQI          ENDDO
  !JQI        ENDDO
  !
  !JQI        DO JCON=JCON,19
  !JQI          DO JCB=1,NCB(JCON)
  !JQI            CBNX(JCB,JCON) = CBNX(JCB,JCON)*REDCBN
  !JQI          ENDDO
  !JQI        ENDDO
  !JQI
  !JQI        DO JCON=20,25
  !JQI          DO JCB=1,NCB(JCON)
  !JQI            CBNX(JCB,JCON) = CBNX(JCB,JCON)*REDCBP
  !JQI          ENDDO
  !JQI        ENDDO
  !
  !JQI        DO JCON=30,32
  !JQI          DO JCB=1,NCB(JCON)
  !JQI            CBNX(JCB,JCON) = CBNX(JCB,JCON)*REDCBP
  !JQI          ENDDO
  !JQI        ENDDO
  !JQI
  !JQI        DO JCON=1,NCP
  !JQI          DO JCB=1,NCB(JCON)
  !JQI            CBOLD(JCB,JCON) = CBNX(JCB,JCON)
  !JQI            CB(JCB,JCON)    = CBNX(JCB,JCON)
  !JQI          ENDDO
  !JQI        ENDDO
  !
  !******* Read in next set of boundary conditions
  !JQI
  !JQI        DO JCON=1,NCP
  !JQI          READ (CBC,1020,END=16) NXCBC,(CBNX(JCB,JCON),JCB=1,NCB(JCON))
  !JQI        ENDDO
  !JQI        NXCBC = (CBCPTR-1)*FILGTH+NXCBC
  !JQI      ENDDO
  !JQI      GO TO 18
  !
  !********* Open next data file
  !
  !JQI16    CONTINUE
  !JQI      CBCPTR = CBCPTR+1
  !JQI      IF (DIAGNOSTICS) WRITE (DIA,*) 'Opening boundary concentra',   &
  !JQI                                     'tion file ',CBCPTR,' at ',     &
  !JQI                                     'day ',JDAY
  !JQI      CLOSE (CBC)
  !JQI      OPEN  (CBC,FILE=CBCFN(CBCPTR),STATUS='OLD')
  !JQI      READ  (CBC,1080)
  !JQI      READ  (CBC,1030) (NCB(JCON),JCON=1,NCP)
  !JQI      READ  (CBC,1080)
  !JQI      IF(NXCBC >= OLDNXCBC) THEN
  !JQI        OLDNXCBC = NXCBC
  !JQI      ENDIF
  !JQI      DO JCON=1,NCP
  !JQI        READ  (CBC,1020) NXCBC,(CBNX(JCB,JCON),JCB=1,NCB(JCON))
  !JQI      ENDDO
  !JQI      NXCBC = (CBCPTR-1)*FILGTH+NXCBC
  !JQI18    CONTINUE
  !JQI      NXTVD = MIN(NXTVD,NXCBC)
  !
  !JQI    ENDIF
  !
  !******* Source One loads              !MNOEL   1-25-93
  !
  !JQI    IF (SOURCE_ONE) THEN
  !JQI      DO WHILE(JDAY >= NXS1)
  !JQI        DO JCON=1,NCP
  !JQI          DO JS1=1,S1LN(JCON)
  !JQI            S1L(JS1,JCON) = S1LNX(JS1,JCON)/86.4
  !JQI          ENDDO
  !JQI        ENDDO
  !
  !*********** Reduce/increase concentrations
  !
  !JQI        DO JCON=4,12
  !JQI          DO JS1=1,S1LN(JCON)
  !JQI            S1L(JS1,JCON) = S1L(JS1,JCON)*REDS1C
  !JQI          ENDDO
  !JQI        ENDDO
  !JQI        DO JCON=13,19
  !JQI          DO JS1=1,S1LN(JCON)
  !JQI            S1L(JS1,JCON) = S1L(JS1,JCON)*REDS1N
  !JQI          ENDDO
  !JQI        ENDDO
  !JQI        DO JCON=20,25
  !JQI          DO JS1=1,S1LN(JCON)
  !JQI            S1L(JS1,JCON) = S1L(JS1,JCON)*REDS1P
  !JQI          ENDDO
  !JQI        ENDDO
  !JQI        DO JCON=30,32
  !JQI          DO JS1=1,S1LN(JCON)
  !JQI            S1L(JS1,JCON) = S1L(JS1,JCON)*REDS1P
  !JQI          ENDDO
  !JQI        ENDDO
  !
  !*********** Read in next set of loads
  !
  !JQI        DO JCON=1,NCP
  !JQI          READ (S1,1020,END=112) NXS1,(S1LNX(JS1,JCON),JS1=1,S1LN(JCON))
  !JQI        ENDDO
  !JQI        NXS1 = (S1PTR-1)*FILGTH+NXS1
  !JQI      ENDDO
  !JQI      GO TO 115
  !
  !********* Open next data file
  !
  !JQI112   CONTINUE
  !JQI      S1PTR = S1PTR+1
  !JQI      IF (DIAGNOSTICS) WRITE (DIA,*) 'Opening Source One file ',    &
  !JQI                                      S1PTR,' at day ',JDAY
  !JQI      CLOSE (S1)
  !JQI      OPEN  (S1,FILE=S1FN(S1PTR),STATUS='OLD')
  !JQI      READ  (S1,1080)
  !JQI      READ  (S1,1030) (S1LN(JCON),JCON=1,NCP)
  !JQI      DO JCON=1,NCP
  !JQI        READ (S1,1030) (S1LB(JP,JCON),JP=1,S1LN(JCON))
  !JQI      ENDDO
  !JQI      READ (S1,1080)
  !JQI      DO JCON=1,NCP
  !JQI        READ (S1,1020) NXS1,(S1LNX(JS1,JCON),JS1=1,S1LN(JCON))
  !JQI      ENDDO
  !JQI      NXS1 = (S1PTR-1)*FILGTH+NXS1
  !JQI115   CONTINUE
  !JQI      NXTVD = MIN(NXTVD,NXS1)
  !JQI    ENDIF
  !
  !******* Source Two loads              !MNOEL   1-25-93
  !
  !JQI    IF (SOURCE_TWO) THEN
  !JQI      DO WHILE(JDAY >= NXS2)
  !JQI        DO JCON=1,NCP
  !JQI          DO JS2=1,S2LN(JCON)
  !JQI            S2L(JS2,JCON) = S2LNX(JS2,JCON)/86.4
  !JQI          ENDDO
  !JQI        ENDDO
  !
  !*********** Reduce/increase concentrations
  !
  !JQI        DO JCON=4,12
  !JQI          DO JS2=1,S2LN(JCON)
  !JQI            S2L(JS2,JCON) = S2L(JS2,JCON)*REDS2C
  !JQI          ENDDO
  !JQI        ENDDO
  !JQI        DO JCON=13,19
  !JQI          DO JS2=1,S2LN(JCON)
  !JQI            S2L(JS2,JCON) = S2L(JS2,JCON)*REDS2N
  !JQI          ENDDO
  !JQI        ENDDO
  !JQI        DO JCON=20,25
  !JQI          DO JS2=1,S2LN(JCON)
  !JQI            S2L(JS2,JCON) = S2L(JS2,JCON)*REDS2P
  !JQI          ENDDO
  !JQI        ENDDO
  !JQI        DO JCON=30,32
  !JQI          DO JS2=1,S2LN(JCON)
  !JQI            S2L(JS2,JCON) = S2L(JS2,JCON)*REDS2P
  !JQI          ENDDO
  !JQI        ENDDO
  !
  !*********** Read in next set of loads
  !
  !JQI        DO JCON=1,NCP
  !JQI          READ (S2,1020,END=119) NXS2,(S2LNX(JS2,JCON),JS2=1,S2LN(JCON))
  !JQI        ENDDO
  !JQI        NXS2 = (S2PTR-1)*FILGTH+NXS2
  !JQI      ENDDO
  !JQI      GO TO 122
  !
  !********* Open next data file
  !
  !JQI119   CONTINUE
  !JQI      S2PTR = S2PTR+1
  !JQI      IF (DIAGNOSTICS) WRITE (DIA,*) 'Opening Source Two ',        &
  !JQI                                     'file ',S2PTR,' at day ',     &
  !JQI                                      JDAY
  !JQI      CLOSE (S2)
  !JQI      OPEN  (S2,FILE=S2FN(S2PTR),STATUS='OLD')
  !JQI      READ  (S2,1080)
  !JQI      READ  (S2,1030) (S2LN(JCON),JCON=1,NCP)
  !JQI      DO JCON=1,NCP
  !JQI        READ (S2,1030) (S2LB(JNP,JCON),JNP=1,S2LN(JCON))
  !JQI      ENDDO
  !JQI      READ (S2,1080)
  !JQI      DO JCON=1,NCP
  !JQI        READ (S2,1020) NXS2,(S2LNX(JS2,JCON),JS2=1,S2LN(JCON))
  !JQI      ENDDO
  !JQI      NXS2 = (S2PTR-1)*FILGTH+NXS2
  !JQI122   CONTINUE
  !JQI      NXTVD = MIN(NXTVD,NXS2)
  !JQI    ENDIF
  !
  !******* Source Three loads
  !
  !JQI    IF (SOURCE_THR) THEN
  !JQI      DO WHILE(JDAY >= NXS3)
  !JQI        DO JCON=1,NCP
  !JQI          DO JS3=1,S3LN(JCON)
  !JQI            S3L(JS3,JCON) = S3LNX(JS3,JCON)/86.4
  !JQI          ENDDO
  !JQI        ENDDO
  !
  !*********** Reduce/increase concentrations
  !
  !JQI        DO JCON=4,12
  !JQI          DO JS3=1,S3LN(JCON)
  !JQI            S3L(JS3,JCON) = S3L(JS3,JCON)*REDS3C
  !JQI          ENDDO
  !JQI        ENDDO
  !JQI        DO JCON=13,19
  !JQI          DO JS3=1,S3LN(JCON)
  !JQI            S3L(JS3,JCON) = S3L(JS3,JCON)*REDS3N
  !JQI          ENDDO
  !JQI        ENDDO
  !JQI        DO JCON=20,25
  !JQI          DO JS3=1,S3LN(JCON)
  !JQI            S3L(JS3,JCON) = S3L(JS3,JCON)*REDS3P
  !JQI          ENDDO
  !JQI        ENDDO
  !JQI        DO JCON=30,32
  !JQI          DO JS3=1,S3LN(JCON)
  !JQI            S3L(JS3,JCON) = S3L(JS3,JCON)*REDS3P
  !JQI          ENDDO
  !JQI        ENDDO
  !
  !*********** Read in next set of loads
  !
  !JQI        DO JCON=1,NCP
  !JQI          READ (S3,1020,END=124) NXS3,(S3LNX(JS3,JCON),JS3=1,S3LN(JCON))
  !JQI        ENDDO
  !JQI        NXS3 = (S3PTR-1)*FILGTH+NXS3
  !JQI      ENDDO
  !JQI      GO TO 125
  !
  !********* Open next data file
  !
  !JQI124   CONTINUE
  !JQI      S3PTR = S3PTR+1
  !JQI      IF (DIAGNOSTICS) WRITE (DIA,*) 'Opening Source Three ',     &
  !JQI                                     'file ',S3PTR,' at day ',    &
  !JQI                                      JDAY
  !JQI      CLOSE (S3)
  !JQI      OPEN  (S3,FILE=S3FN(S3PTR),STATUS='OLD')
  !JQI      READ  (S3,1080)
  !JQI      READ  (S3,1030) (S3LN(JCON),JCON=1,NCP)
  !JQI      DO JCON=1,NCP
  !JQI        READ (S3,1030) (S3LB(JNP,JCON),JNP=1,S3LN(JCON))
  !JQI      ENDDO
  !JQI      READ (S3,1080)
  !JQI      DO JCON=1,NCP
  !JQI        READ (S3,1020) NXS3,(S3LNX(JS3,JCON),JS3=1,S3LN(JCON))
  !JQI      ENDDO
  !JQI      NXS3 = (S3PTR-1)*FILGTH+NXS3
  !JQI125   CONTINUE
  !JQI      NXTVD = MIN(NXTVD,NXS3)
  !JQI    ENDIF
  !
  !
  !******* Benthic fluxes
  !
      If (BENTHIC_FLUXES) Then
126      Continue
         Do WHILE (JDAY >= NXBFI)
            Do JCON = 1, 9
               Do BB = 1, MLOC
                  BFLUX (BB, JCON) = BFLUXNX (BB, JCON)
               End Do
            End Do
        !
        !Store previously read results (without temperature effects)
        !
            BENDOC = BFLUX (:, 1)
            BENNH4 = BFLUX (:, 2)
            BENNO3 = BFLUX (:, 3)
            BENDON = BFLUX (:, 4)
            BENPO4 = BFLUX (:, 5)
            BENDOP = BFLUX (:, 6)
            BENCOD = BFLUX (:, 7)
            BENDO = BFLUX (:, 8)
            BENSA = BFLUX (:, 9)
        !
            Do BB = 1, MLOC
           !JQI          BFLUXB(BB,1) = BFLUX(BB,1)                              ! DOC
           !JQI          BFLUXB(BB,2) = BFLUX(BB,2)                              ! NH4
           !JQI          BFLUXB(BB,3) = BFLUX(BB,3)                              ! NO3
           !JQI          BFLUXB(BB,4) = BFLUX(BB,5)                              ! PO4
           !JQI          BFLUXB(BB,5) = BFLUX(BB,7)                              ! COD
           !JQI          BFLUXB(BB,6) = BFLUX(BB,8)                              ! DO
           !JQI          BFLUXB(BB,7) = BFLUX(BB,9)                              ! SIAT
           !store data without the temperature effects, temperature effects
           !will be applied in wqm_kin.F using these values
               BENDOCB (BB) = BFLUX (BB, 1)! DOC
               BENNH4B (BB) = BFLUX (BB, 2)! NH4
               BENNO3B (BB) = BFLUX (BB, 3)! NO3
               BENPO4B (BB) = BFLUX (BB, 5)! PO4
               BENCODB (BB) = BFLUX (BB, 7)! COD
               BENDOB (BB) = BFLUX (BB, 8)! DO
               BENSAB (BB) = BFLUX (BB, 9)! SIAT
            End Do
        !
            Do BB = 1, MGL
               Read (BFI, 1103, END=127) NXBFI, (BFLUXNX_GL(BB, JCON), &
              & JCON=1, 9)
           !             WRITE(*,1104)NXBFI,(BFLUXNX_GL(1,JCON),JCON=1,9)
            End Do
        !
        ! ong        READ(BFI,1101) NXBFI
        ! ong        DO JCON=1,9
        ! ong          READ (BFI,1020,END=127) NXBFI,(BFLUXNX_GL(BB,JCON),BB=1,MGL)
        ! ong        ENDDO
        !
            If (SERIAL) BFLUXNX = BFLUXNX_GL
        !
#   if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do JCON = 1, 9
                  Do I = 1, MLOC
                     BFLUXNX (I, JCON) = BFLUXNX_GL (NGID(I), JCON)
                  End Do
                  Do I = 1, NHN
                     BFLUXNX (I+MLOC, JCON) = BFLUXNX_GL (HN_LST(I), &
                    & JCON)
                  End Do
               End Do
            End If
#   endif
        !
        !--moved to before the reading of BFI by WLong------
        !         DO BB=1,MLOC
        !  !JQI          BFLUXB(BB,1) = BFLUX(BB,1)                              ! DOC
        !  !JQI          BFLUXB(BB,2) = BFLUX(BB,2)                              ! NH4
        !  !JQI          BFLUXB(BB,3) = BFLUX(BB,3)                              ! NO3
        !  !JQI          BFLUXB(BB,4) = BFLUX(BB,5)                              ! PO4
        !  !JQI          BFLUXB(BB,5) = BFLUX(BB,7)                              ! COD
        !  !JQI          BFLUXB(BB,6) = BFLUX(BB,8)                              ! DO
        !  !JQI          BFLUXB(BB,7) = BFLUX(BB,9)                         ! SIAT
        !  !store data without the temperature effects, temperature effects
        !  !will be applied in wqm_kin.F using these values
        !            BENDOCB(BB) = BFLUX(BB,1)                              ! DOC
        !            BENNH4B(BB) = BFLUX(BB,2)                              ! NH4
        !            BENNO3B(BB) = BFLUX(BB,3)                              ! NO3
        !            BENPO4B(BB) = BFLUX(BB,5)                              ! PO4
        !            BENCODB(BB) = BFLUX(BB,7)                              ! COD
        !            BENDOB(BB)  = BFLUX(BB,8)                              ! DO
        !           BENSAB(BB)  = BFLUX(BB,9)                              ! SIAT
        !          ENDDO
        !---------------------------------------------------
        !
            NXBFI = (BFIPTR-1) * FILGTH + NXBFI
         End Do
         Go To 128
     !
     !********* Open next data file
     !
127      Continue
         BFIPTR = BFIPTR + 1
#if defined (KILLTREES)
         If (DIAGNOSTICS) WRITE (DIA,*) 'Opening benthic flux file ', &
        & BFIPTR, ' at day ', JDAY
#endif
         Close (BFI)
         Open (BFI, File=BFIFN(BFIPTR), Status='OLD')
     !
         Read (BFI, 1041)!skip 3 title lines
         Read (BFI, 1050) KSDOC, KSNH4, KSNO3, KSPO4, KSO, KSSA
     !       WRITE(*,1051) KSDOC, KSNH4, KSNO3, KSPO4, KSO, KSSA
         Read (BFI, 1050) TRSDOC, TRSNH4, TRSNO3, TRSPO4, TRSO, TRSSA
     !       WRITE(*,1051) TRSDOC,TRSNH4,TRSNO3,TRSPO4,TRSO,TRSSA
         Read (BFI, 1050) MTCNO3, SEDNO3, KHSO
     !       WRITE(*,1051) MTCNO3, SEDNO3, KHSO
         Read (BFI, 1100)
     !
         Do BB = 1, MGL
            Read (BFI, 1103, END=127) NXBFI, (BFLUXNX_GL(BB, JCON), &
           & JCON=1, 9)
        !          WRITE(*,1104)NXBFI,(BFLUXNX_GL(1,JCON),JCON=1,9)
         End Do
     !
         NXBFI = (BFIPTR-1) * FILGTH + NXBFI
     !
         Do JCON = 1, 9
            If (SERIAL) BFLUXNX (:, JCON) = BFLUXNX_GL (:, JCON)
        !
#   if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do I = 1, MLOC
                  BFLUXNX (I, JCON) = BFLUXNX_GL (NGID(I), JCON)
               End Do
               Do I = 1, NHN
                  BFLUXNX (I+MLOC, JCON) = BFLUXNX_GL (HN_LST(I), JCON)
               End Do
            End If
#   endif
         End Do
         If (JDAY >= NXBFI) Go To 126
128      Continue
         NXTVD = Min (NXTVD, NXBFI)
     ! ong debug
     !      WRITE(*,*)'Found record in file number:',BFIPTR
     !      WRITE(*,1104)NXBFI,(BFLUXNX_GL(1,JCON),JCON=1,9)
     !
      End If
  !
  !
  !******* Atmospheric Loads
  !
      If (ATMOS_LOADS) Then
133      Continue
         Do WHILE (JDAY >= NXATM)
            PRECIP = PRECIPNX / 8640000.
            ATMNH4 = ANH4NX
            ATMNO3 = ANO3NX
            ATMLDON = ALDONNX
            ATMRDON = ARDONNX
            ATMPO4 = APO4NX
            ATMLDOP = ALDOPNX
            ATMRDOP = ARDOPNX
            Read (ATM, 1010, END=134) NXATM, PRECIPNX, ANH4NX, ANO3NX, &
           & ALDONNX, ARDONNX, APO4NX, ALDOPNX, ARDOPNX
            NXATM = (ATMPTR-1) * FILGTH + NXATM
         End Do
         Go To 135
     !
     !********* Open next data file
     !
134      Continue
         ATMPTR = ATMPTR + 1
#if defined (KILLTREES)
         If (DIAGNOSTICS) WRITE (DIA,*) 'Opening atmospheric load', 'in&
        &g file ', ATMPTR, ' at ', 'day ', JDAY
#endif
         Close (ATM)
         Open (ATM, File=ATMFN(ATMPTR), Status='OLD')
         Read (ATM, 1000)
         Read (ATM, 1010) NXATM, PRECIPNX, ANH4NX, ANO3NX, ALDONNX, &
        & ARDONNX, APO4NX, ALDOPNX, ARDOPNX
         NXATM = (ATMPTR-1) * FILGTH + NXATM
         If (JDAY >= NXATM) Go To 133
135      Continue
         NXTVD = Min (NXTVD, NXATM)
      End If
  !
  !
  !******* Submerged Aquatic Vegetation
  !
      If (SAV_LOADS) Then
136      Continue
         Do WHILE (JDAY >= NXSAV)
            Do B = 1, MLOC
           !             FACTOR = SAVAREA(B)/SFA(B)            !Moved SFA is not calculated at all!
               FACTOR = SAVAREA (B) / ART1 (B)!replaced SFA by ART1 
               LDOCSAVW (B) = SLDOCNX * FACTOR
               RDOCSAVW (B) = SRDOCNX * FACTOR
               LPOCSAVW (B) = SLPOCNX * FACTOR
               RPOCSAVW (B) = SRPOCNX * FACTOR
               DOSAVW (B) = SDONX * FACTOR
            End Do
            Read (SVI, 1010, END=137) NXSAV, SLDOCNX, SRDOCNX, SLPOCNX, &
           & SRPOCNX, SDONX
            NXSAV = (SAVPTR-1) * FILGTH + NXSAV
         End Do
         Go To 138
     !
     !********* Open next data file
137      Continue
         SAVPTR = SAVPTR + 1
#if defined (KILLTREES)
         If (DIAGNOSTICS) WRITE (DIA,*) 'Opening aquatic vegetation', '&
        & file ', SAVPTR, ' at day ', JDAY
#endif
         Close (SVI)
         Open (SVI, File=SVIFN(SAVPTR), Status='OLD')
         Read (SVI, 1080)
         Read (SVI, 1026) (SAVAREA(B), B=1, NBB)
         Read (SVI, 1080)
         Read (SVI, 1010) NXSAV, SLDOCNX, SRDOCNX, SLPOCNX, SRPOCNX, &
        & SDONX
         NXSAV = (SAVPTR-1) * FILGTH + NXSAV
         If (JDAY >= NXSAV) Go To 136
138      Continue
         NXTVD = Min (NXTVD, NXSAV)
      End If
  !
  !
  !  ******* Atmospheric pCO2  ()
  !          NOTE: time-varying pCO2atm only valid for spatially uniform pCO2atm
      If (MAPCO2ATM == 'TIMEVAR') Then
!
         Do WHILE (JDAY >= NXPCO2)
            Do I = 0, MTLOC
               pCO2atm (I) = pco2atmNX
            End Do
            Read (AIRC, 1031) NXPCO2, pco2atmNX
         End Do
         NXTVD = Min (NXTVD, NXPCO2)
!
         If (MSR) WRITE (*,*) 'time-varying pCO2atm=', pCO2atm (1), 'ne&
        &xt: NXPCO2=', NXPCO2, 'd pco2atmNX=', pco2atmNX
      End If
  !
  !
  !******* Input FORMATs
  !
1000  Format (/ //)
1010  Format (10 F8.0, :/(:8 X, 9 F8.0))
1020  Format (8 X, 9 F8.0, :/(:16 X, 8 F8.0))
1025  Format (16 X, F8.0)
1026  Format (/ / (8 X, 9 F8.0))
1030  Format (/ / (8 X, 9 I8))
1040  Format (/)
1041  Format (/ /)!advance two lines
1050  Format (://8 X, 9 F8.0)!go down 2 lines and read
  !but do not advance to next line when all data
  !items are processed
1051  Format (://8 X, 9 F8.4)!go down 2 lines and write, but do not
  !continue going down when writing is finished
1060  Format (8 X, 9 F8.0)
1070  Format (/ / (:8 X, 6 F8.0))
1080  Format (/)
1082  Format (/ / 8 X, 2 A8 //)
1100  Format (/)
1101  Format (1 F8.0)
1102  Format (1 F8.0)
1103  Format (10 F8.0)
1104  Format (10 F8.4)
1031  Format (F8.1, F8.3)
  !
End Subroutine TVDS
!
