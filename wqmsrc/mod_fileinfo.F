!mod_fileinfo.F
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
Module MOD_FILEINFO
      Use MOD_SIZES, Only: NCP
  !
      Implicit None
      Save
      Integer :: DIA, CBC, S1, S2, S3, BFI, BAI, DFI, WCL, MET, BFO, &
     & BAO, DFO, KEI, ATM, STL, AGR, SVI, SVO, KFL, ZOO, ZFO, ALO, CON, &
     & RSO, SNP, PLT, APL, TFL, OPL, SFI, SFO, ICI, ICO, MRL, MBL, RSI, &
     & UNIT_LINKAGE, UNIT_STN, UNIT_HIS, AIRC,IREPORT,INRIV,INOBC,INPT,INK1K2, INCTR
  !DIA,
  !CBC,
  !S1,
  !S2,
  !S3,
  !BFI, 		  !Benthic Flux data file or Sediment Diagenesis module input file
  !BAI,                   !Benthic algae input file
  !DFI,			  !Deposition feeder input file
  !WCL,			  !Sediment overlying water column input file
  !MET,
  !BFO,
  !BAO,                   !Benthic algae output file
  !DFO,			  !Deposition feeder output file
  !KEI,
  !ATM,
  !STL,
  !AGR,
  !SVI,
  !SVO,
  !KFL,
  !ZOO,
  !ZFO,
  !ALO,
  !CON,
  !RSO,
  !SNP,
  !PLT,
  !APL,
  !TFL,
  !OPL,
  !SFI,
  !SFO,
  !ICI,
  !ICO,
  !MRL,
  !MBL,
  !RSI,
  !UNIT_LINKAGE,
  !UNIT_STN,
  !UNIT_HIS
!
      Character (Len=24) :: CNAME (NCP)!Wen Long this should be moved to wqm_alg.F and also wqm_alg and wqm_kin should be
  !combined, otherwise, we would need to separate out zooplankton etc
      Character (Len=18) :: CONFN = 'wqm_con.npt'
  !
  !--------------------end of data declarations---------------------------------c
  !
Contains
  !
  !subroutines:
  ! subroutine	INIT_FILE_INFO()
  !
      Subroutine INIT_FILE_INFO ()
    !
	     IREPORT=6            !file handle for screen output (default fortran screen output handle)
         UNIT_LINKAGE = 7 !file handle for LINKAGE input control
         UNIT_STN = 200 !file handle for station output control
         UNIT_HIS = 100 !file handle for history output control
		!actual handle will be UNIT_HIS+k, where k is layer number
		 INCTR=8        !control file
         CON = 10
		 INOBC=13  !file number of open boundary concentration
         BFI = 21
         BAI = 22  !Benthic alage input
         S1  = 23
         S2  = 24
         ATM = 25
		 STL = 26
         AGR = 27
         SVI = 28
         MET = 29
         KEI = 30
         SVO = 31
		 DFI = 32

         S3  = 34
         ICI = 35
         ICO = 36
         MRL = 37
         ZOO = 38
         ZFO = 39
         DIA = 40
         BFO = 41
		 RSO = 42
         SNP = 43
         PLT = 44
         APL = 45
         TFL = 46
         KFL = 47
         OPL = 48
         MBL = 49
         ALO = 50

		 SFI = 58
         SFO = 59
         CBC = 60
		 INRIV = 61  !file number for river flow
         INPT = 63	 !file number for river point source (concentration)
		 RSI = 75
         BAO = 76 	!Benthic algae output
         DFO = 77 	!Deposition feeder output

		 INK1K2=150 	!file handle for K1, K2, KSO4 constants for pH module
		 AIRC =151		!file handle for air sea flux of CO2 for pH  module

    !
    !The following SSNAMEs moved to wqm_sed.F
    !SSNAME(1)  = 'Sediment Temperature'
    !SSNAME(2)  = 'Sediment POP        '
    !SSNAME(3)  = 'Sediment PON        '
    !SSNAME(4)  = 'Sediment POC        '
    !SSNAME(5)  = 'Sediment PBS        '
    !SSNAME(6)  = 'Sediment PO4        '
    !SSNAME(7)  = 'Sediment NH4        '
    !SSNAME(8)  = 'Sediment NO3        '
    !SSNAME(9)  = 'Sediment HS         '
    !SSNAME(10) = 'Sediment CH4        '
    !SSNAME(11) = 'Sediment CH4        '
    !SSNAME(12) = 'Sediment SO4        '
    !SSNAME(13) = 'Sediment DSIL       '
    !SSNAME(14) = 'Benthic Stress      '
    !SSNAME(15) = 'Benthic Algae       '
    !SSNAME(16) = 'Deposit Feeders     '
    !SSNAME(17) = 'Suspension Feeders  '
    !
         CNAME (1) = 'Temperature'
         CNAME (2) = 'Salinity'
         CNAME (3) = 'Inorganic Solids'
         CNAME (4) = 'Algal Group 1'
         CNAME (5) = 'Algal Group 2'
         CNAME (6) = 'Algal Group 3'
         CNAME (7) = 'Zooplankton Group 1'
         CNAME (8) = 'Zooplankton Group 2'
         CNAME (9) = 'Labile DOC'
         CNAME (10) = 'Refractory DOC'
         CNAME (11) = 'Labile POC'
         CNAME (12) = 'Refractory POC'
         CNAME (13) = 'Ammonium'
         CNAME (14) = 'Nitrate-nitrite'
         CNAME (15) = 'Urea'
         CNAME (16) = 'Labile DON'
         CNAME (17) = 'Refractory DON'
         CNAME (18) = 'Labile PON'
         CNAME (19) = 'Refractory PON'
         CNAME (20) = 'Total phosphate'
         CNAME (21) = 'Labile DOP'
         CNAME (22) = 'Refractory DOP'
         CNAME (23) = 'Labile POP'
         CNAME (24) = 'Refractory POP'
         CNAME (25) = 'Particulate Inorganic P'
         CNAME (26) = 'COD'
         CNAME (27) = 'Dissolved oxygen'
         CNAME (28) = 'Particulate silica  '
         CNAME (29) = 'Dissolved silica'
         CNAME (30) = 'Internal P Group 1' !Phosphorus in algae 1
         CNAME (31) = 'Internal P Group 2' !Phosphorus in algae 2
         CNAME (32) = 'Internal P Group 3' !Phosphorus in algae 3
         CNAME (33) = 'Dissolved Inorganic Carbon' !
         CNAME (34) = 'Total Alkalinity' !
    !CNAME(35) =  'pH'    							 !
    !CNAME(36) =  'Seawater pCO2'    					         !
!
         Return
      End Subroutine INIT_FILE_INFO
  !
End Module MOD_FILEINFO
!
