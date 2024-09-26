!mod_DIC_TALK.f
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
Module MOD_DIC_TALK
!
      Use MOD_PREC, Only: SP
      Use MOD_CONTROL, Only: MSR
!
      Use MOD_LIMS, Only: MLOC, KBM1, MTLOC
#if defined (WET_DRY)
	  USE MOD_WD, ONLY: ISWETN
#endif
      Implicit None
      Integer :: AIRSEA_OPTION !choice of parameterization to calculate O2 and CO2 air-sea exchange
!
      Real (SP), Dimension (2) :: JDstart1 !AN

!
Contains
     !=====================================================================|
     !    This file has the following subroutines:
     !
     !   PCO2AIR_READ_INITIAL: Read atmospheric pCO2 data (initial record)
     !                and the choice of parameterization to compute air-sea
     !                exchange of O2 and CO2
     !   GET_AIRSEA_METHOD: Determine what AIRSEA_OPTION should be used
     !   TOTINC      : Calculation of changes in TDIC (DTTIC)
     !   ALKALINE    : Calculation of changes in TALK (DTALK)
     !
     !=====================================================================|
!
  !************************************************************************
  !**          S U B R O U T I N E   P C O 2 A I R _ R E A D             **
  !************************************************************************
!
      Subroutine PCO2AIR_READ_INITIAL
!
         Use MOD_SIZES, Only: MGL
         Use MOD_CONTROL, Only: PAR, SERIAL
         USE MOD_FILEINFO, ONLY: AIRC
!
         Use MOD_CO2SYS, Only: K0, FugFac, pCO2atm, &
            & MAPCO2ATM, UNITCO2ATM, NXPCO2, pco2atmNX

#    if defined (MULTIPROCESSOR)
         Use MOD_PAR
#    endif
!
         Implicit None

         Integer :: I
         Real (SP), Allocatable :: RTMP21 (:)
    ! Read file

         Open (AIRC, File='inputs/input_pCO2atm.dat', Status='OLD')

         Read (AIRC, 1009) AIRSEA_OPTION
         If (MSR) write (*,*) 'Parameterization air-sea exchange: (1=re&
        &aeration, 2=W92, 3=N00, 4=Ho06, 5=W12, 6=PNNL)'
         If (MSR) write (*,*) '  --> AIRSEA_OPTION=', AIRSEA_OPTION

         Read (AIRC, 1010) MAPCO2ATM

         Read (AIRC, 1011) UNITCO2ATM
!
         If (MAPCO2ATM == 'UNIFORM') Then
            Read (AIRC, 1030) pCO2atm (1)
            NXPCO2 = - 999.9
!
            Do I = 0, MTLOC
               pCO2atm (I) = pCO2atm (1)
            End Do
!
         Else If (MAPCO2ATM == 'TIMEVAR') Then
            Read (AIRC, 1031) NXPCO2, pco2atmNX
!
            Do I = 0, MTLOC
               pCO2atm (I) = pco2atmNX
            End Do
!
         Else 
            Allocate (RTMP21(MGL))
            RTMP21 = 0.0
!
            Do I = 1, MGL
               Read (AIRC, 1030) RTMP21 (I)
            End Do
!
            If (SERIAL) Then
               pCO2atm = RTMP21
            End If
#      if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do I = 1, MLOC
                  pCO2atm (I) = RTMP21 (NGID(I))
               End Do
               Do I = 1, NHN !
			     !
                  pCO2atm (I+MLOC) = RTMP21 (HN_LST(I))
               End Do
            End If
#      endif
            Deallocate (RTMP21)
         End If

         Read (AIRC, 1032) JDstart1(1), JDstart1(2)
         If (MSR) write (*,*) '  JDstart1=', JDstart1(1), JDstart1(2)

    !***** Input FORMAT statements
!
1009     Format (/ // (I1))
1010     Format (A7)
1011     Format (A4)
1030     Format (F8.0)
1031     Format (F8.1, F8.3)
1032     Format (F8.0, F8.0) !AN
!
         Return
      End Subroutine PCO2AIR_READ_INITIAL
!

      Integer Function GET_AIRSEA_METHOD()
        Use MOD_PREC, Only: SMOD
        Use MOD_WQM, Only: JDAY

        Real(SP) :: doy
        Integer :: temp

        If (AIRSEA_OPTION == 6) Then
          doy = AINT(1.0 + SMOD (JDAY, 365.25_SP))
			    If (doy <= JDstart1(1) .OR. doy >= JDstart1(2)) Then
            temp = 1
          Else
            temp = 5
          End If
        Else
          temp = AIRSEA_OPTION
        End If
        GET_AIRSEA_METHOD = temp
      End Function GET_AIRSEA_METHOD

  !************************************************************************
  !**                 S U B R O U T I N E   T O T I N C                  **
  !************************************************************************
!
      Subroutine TOTINC (DTTIC)
!
         Use MOD_HYDROVARS, Only: DZ,DZ2D, D, RHO !Added RHO
!
    USE MOD_WQM, ONLY :        &!
         B1,             &!
         B2,             &!
         BM1,            &!
         BM2,            &!
         DOXG,           &!
         FCLD1,          &!
         FCLD2,          &!
         FCLDP,          &!
         FCLP1,          &!
         FCLP2,          &!
         FCLPP,          &!
         FCRD1,          &!
         FCRD2,          &!
         FCRDP,          &!
         FCRP1,          &!
         FCRP2,          &!
         FDOP,           &!
         KHR1,           &!
         KHR2,           &!
         P1,             &!
         P2,             &!
         PR1,            &!
         PR2,            &!
         PRSP1,          &!
         PRSP2,          &!
         DENIT,          &!
         KHODOC,         &!
         MNLLDOC,        &!
         MNLRDOC,        &!
         SALT,           &!
         T,              &!
         AREAR,           &!
         BREAR,           &!
         CREAR,           &!
         WMS,             &
         NT,              &
         FTCOD,           &!
         KHOCOD,          &
         BENTHIC_FLUXES,  &
     & SAV_CALC,      & ! for SAV-DIC interaction
#if defined(CSOD_CALC)
		 CSODflux,        &
#elif defined(MASS_CALC)
		 PCFWS,           &
         BENCH4A,         &
         BENCH4G,         &
#endif
         COD
!
#if defined(MASS_CALC)
         Use MOD_SED, Only: BURIALC
#endif
!
         USE MOD_ZOOP, ONLY: AOCRSZ, AOCRLZ, DOSZ, DOLZ, B1SZ,B2SZ, &
                              & B1LZ,B2LZ ! AN activates it for zoopl, added zooplankton predation
!
         Use MOD_CO2SYS, Only: CO2star_sat, pCO2atm, CO2star_surf, &
        & UNITCO2ATM, K0, FugFac, DICUPT, DICBMP, &
        & DICPRD, DICMNL, DICDEN, DICGAS, DICSED !DICNIT

         Use MOD_WQMINIT, Only: COD_CALC, ZOO_CALC 
         USE MOD_SAV, ONLY : BM_LEAF, LEAF, ACDWSAV,PLEAF,SAVCELL,NSAVCELL     !  added for SAV

         Implicit None
!
         Integer :: I, J, K, SAVID,II
!
         Real (SP) :: CP1, CP2, FRDO1, FRDO2, DOR1, DOR2, DOP1, DOP2, &
        & CPCOD, CPSAV
         Real (SP) :: FLUXCO2, TDOS, RNU, KRDO, KRCO2, FACTOR, Sc !parameters in air-sea CO2 flux calculation
         Real (SP), Dimension (0:MTLOC, KBM1) :: DTTIC
         Real (SP) :: auxLB
!
#if defined(CSOD_CALC)
         Real (SP) :: RCCO
		 Real (SP) :: BENDIC !mmoleC/m^3/day
#elif defined(MASS_CALC)
         Real (SP) :: BENDIC !mmoleC/m^3/day
#endif
         Integer :: airsea_method
!
         DTTIC = 0.0
!
    !** note: coded for the sake of completeness
    !           all equations related to Algae 3
         CPSAV = 0.0
         Do K = 1, KBM1
            Do I = 1, MLOC
!
          !********* Uptake
               CP1 = P1 (I, K) * B1 (I, K)   !  1/d * gC/m3 =gC/m3/d
               CP2 = P2 (I, K) * B2 (I, K)
          !    CP3    = P3(I,K)*B3(I,K)
          !--------------------- uptake by the SAV ----------------------------------
          IF((SAV_CALC) .AND. (K<=6)) THEN
               CPSAV=0.0
               DO II=1,NSAVCELL
                  IF(I == SAVCELL(II)) THEN
                     CPSAV = 0.14*PLEAF(I,1)*LEAF(I,1)/(D(I)*DZ2D(I,K))   ! gC/m3   ! Added a factor 0.14 to compute net DIC uptake (hardcoded based on SSM calibration). 0.14 is arbitrary and needs to be updated by the user  
                     EXIT                                                           ! This accounts for DIC generation due to SAV respirartion that is not include in this version
                  END IF
               END DO
          ELSE
               CPSAV = 0.0
          END IF
          !--------------------------------------------------------------------------
          !********* Respiration
          !-- Fracion of losses by basal metabolism and photoresp that go to TDIC
               FRDO1 = (1.-FCLD1-FCRD1-FCLP1-FCRP1) * DOXG (I, K) / &
              & (KHR1+DOXG(I, K))
               FRDO2 = (1.-FCLD2-FCRD2-FCLP2-FCRP2) * DOXG (I, K) / &
              & (KHR2+DOXG(I, K))

               DOR1 = (P1(I, K)*PRSP1+BM1(I, K)) * FRDO1 * B1 (I, K)
               DOR2 = (P2(I, K)*PRSP2+BM2(I, K)) * FRDO2 * B2 (I, K)

!
          !-- Fracion of losses by predation that go to TDIC
               DOP1 = FDOP * (PR1 (I, K)+B1SZ(I, K)+B1LZ(I, K)) * DOXG (I, K) / (KHR1+DOXG(I, &
              & K))
               DOP2 = FDOP * (PR2 (I, K)+B2SZ(I, K)+B2LZ(I, K)) * DOXG (I, K) / (KHR2+DOXG(I, &
              & K))

               If (COD_CALC) Then
             !-- Production of TDIC by the oxydation of sulfide to sulfate (?) (COD)
                  CPCOD = DOXG (I, K) / (KHOCOD+DOXG(I, K)) * FTCOD (I, &
                 & K) * COD (I, K) * 12. / 64. ! CHANGE O2 EQUIVALENTS TO C
               End If !

#  if !defined (WET_DRY)
				IF(D(I) > 0.0_SP)THEN
#  else
				IF(ISWETN(I) == 1)THEN
#  endif

          !: accumulate fluxes per DThistoryFile(t_his_dlt) period (integration will be finalized in wqm_main.F, wqm_output subrou
          IF(SAV_CALC) THEN
               DICUPT (I, K) = DICUPT (I, K) + (-CP1-CP2-CPSAV) / 12. * & ! This term was added but it was not used in DTTIC in line 344 and 350 
              & 1000.0 !TDIC uptake mmolC/m3/d				  ! The term CPSAV is explicilty added in lines 344 and 350 **LK and TK 3/18/20 
          ELSE
               DICUPT (I, K) = DICUPT (I, K) + (-CP1-CP2) / 12. * &
              & 1000.0 !TDIC uptake mmolC/m3/d
          END IF
               DICBMP (I, K) = DICBMP (I, K) + (DOR1+DOR2) / 12. * &
              & 1000.0 !TDIC source by basal metabolism and photorespiration
               DICPRD (I, K) = DICPRD (I, K) + (DOP1+DOP2) / 12. * &
              & 1000.0 !TDIC source by predation
               DICMNL (I, K) = DICMNL (I, K) + (MNLLDOC(I, &
              & K)+MNLRDOC(I, K)) / 12. * 1000.0 !TDIC source by mineralization of LDOC and RDOC
               DICDEN (I, K) = DICDEN (I, K) + DENIT (I, K) / 12. * &
              & 1000.0 !TDIC source by water column denitrification
	      !
          !********* Change in TIC
          !
	           If (ZOO_CALC) Then

			       DTTIC(I,K) = (-CP1-CP2-CPSAV+DOR1+DOR2 +DOP1+DOP2    &   !-CP3+DOR3+DOP3 removed !CPSAV term added
                  & + MNLLDOC(I,K) +MNLRDOC(I,K)                &
                  & + DENIT(I,K)                                &   !: DENIT was missing in C.Cerco's code
                  & + DOSZ(I,K)/AOCRSZ+DOLZ(I,K)/AOCRLZ)        &   !: Added for Zooplankton
                  & /12./86400. * 1000.0
			   Else
                   DTTIC (I, K) = (-CP1-CP2-CPSAV+DOR1+DOR2+DOP1+DOP2+MNLLDOC(I, &
                  & K)+MNLRDOC(I, K)+DENIT(I, K)) / 12. / 86400. * 1000.0
			   End If
!
               If (COD_CALC) Then
                  DTTIC (I, K) = DTTIC (I, K) + CPCOD / 12. / 86400. * &
                 & 1000.0
               End If

				ENDIF
            End Do
         End Do
!
!
        !********  Atmospheric Exchange
        ! Transform pCO2atm OR fCO2atm to CO2star_sat:  CO2star_sat = K0 * fCO2atm = K0 * pCO2atm * FugFac    [K0 units: mol/kg-SW/atm]
        If (UNITCO2ATM == 'pCO2') Then !if we provided pCO2atm in microatm or ppmv
          Do I = 1, MLOC

            CO2star_sat (I) = K0 (I, 1) * pCO2atm (I) * FugFac (I, &
              & 1) * RHO (I, 1) * 1.E-3_SP !mmol/m3
          End Do
        Else
          Do I = 1, MLOC
            auxLB = K0 (I, 1) * pCO2atm (I) * RHO (I, 1) * 1.E-3_SP
               CO2star_sat (I) = auxLB

          End Do
        End If
        airsea_method = GET_AIRSEA_METHOD()

        If (airsea_method == 1) Then 
          FACTOR = AREAR * (BREAR*WMS) ** CREAR
          Do I = 1, MLOC
            TDOS = T (I, 1)
            RNU = 0.54 + 0.7 * TDOS / 30.0 - 0.07 * SALT (I, 1) / &
              & 35.0
			      ! KRDO      = 0.157*RNU*(1.5*WMS)**1.5
            KRDO = FACTOR * RNU !up to here, same as in calculation of air-sea O2 fluxes (m/day)
            KRCO2 = (32.0/44.0) ** 0.25 * KRDO !conversion of oxygen reareation into CO2 reareagion rate by Chapra (1997)

#  if !defined (WET_DRY)
				    IF(D(I) > 0.0_SP)THEN
#  else
				    IF(ISWETN(I) == 1)THEN
#  endif
					    FLUXCO2 = KRCO2 / (D(I)*DZ2D(I,1)) * &
                & (CO2star_sat(I)-CO2star_surf(I))

				    ELSE
					    FLUXCO2=0.0_SP
				    ENDIF
!
            DTTIC (I, 1) = DTTIC (I, 1) + FLUXCO2 / 86400.0

            DICGAS (I, 1) = DICGAS (I, 1) + FLUXCO2 !TDIC flux from air-sea gas exchange, accumulated
          End Do
!
        Else
          If (airsea_method == 2) Then !Wanninkhof (1992)
            FACTOR = 0.31 * WMS ** 2.0 !0.31 for winds in steady staate;  0.39 for long-term averaged winds
          Else If (airsea_method == 3) Then !Nightingale (2000)
            FACTOR = (0.33*WMS+0.222*WMS**2.0)
          Else If (airsea_method == 4) Then !Ho (2006)
            FACTOR = 0.266 * WMS ** 2.0
          Else If (airsea_method == 5) Then !Wanninkhof et al (2013): wind in equation is <U^2>,meaning the time mean of U^2
            FACTOR = 0.251 * WMS ** 2.0 !over the time interval of flux determination (they compute monthly fluxes with 6hourly wind data)
          End If
!
          Do I = 1, MLOC
            Sc = 2073.1 - 125.62 * T (I, 1) + 3.6276 * T (I, 1) ** &
              & 2.0 - 0.043219 * T (I, 1) ** 3.0 !Schmidt number for CO2 (W92)

            KRCO2 = FACTOR * Sqrt (660.0/Sc)!gas transfer velocity in cm/hr
            KRCO2 = KRCO2 * 24.0 / 100.0 !m/day

#  if !defined (WET_DRY)
            IF(D(I) > 0.0_SP)THEN
#  else
            IF(ISWETN(I) == 1)THEN
#  endif
!
              FLUXCO2 = KRCO2 / (D(I)*DZ2D(I,1)) * &
                  & (CO2star_sat(I)-CO2star_surf(I))!mmolC/m3/d
				    ELSE
					    FLUXCO2=0.0_SP
				    ENDIF
!
            DTTIC (I, 1) = DTTIC (I, 1) + FLUXCO2 / 86400.0 !mmolC/m3/sec

            DICGAS (I, 1) = DICGAS (I, 1) + FLUXCO2 !TDIC flux from air-sea gas exchange, accumulated
          End Do
        End If

	  !******** Sediment flux
         If ( .Not. BENTHIC_FLUXES) Then
#if defined(CSOD_CALC)
       !---Use QUAL2K approach:    DICflux    =    CSOD    *    RCCO    *  DZbottomLayerWaterColumn^(-1)
       !                         mmolC/m^3/d     gO2/m^2/d    mmolC/gO2        1/m

            RCCO = 1.0_SP / 2.667_SP !gC/gO2
            RCCO = RCCO / 12.0_SP    !molC/gO2  (1 molC = 12.0107 gC) : USE 12 instead of 12.0107 to match round up gN:molN to 14.0
            RCCO = RCCO * 1000.0_SP !mmolC/gO2
            Do I = 1, MLOC

#  if !defined (WET_DRY)
			IF(D(I) > 0.0_SP)THEN
#  else
			IF(ISWETN(I) == 1)THEN
#  endif
               BENDIC = CSODflux (I) * RCCO / (D(I)*DZ2D(I,KBM1))	!CSODflux = CSOD + XJCNO3 in mod_sed.F

			ELSE
			   BENDIC=0.0_SP
			ENDIF
               DTTIC (I, KBM1) = DTTIC (I, KBM1) + BENDIC / 86400.0_SP
			   !DICSED(I,KBM1) = BENDIC  !TDIC flux from sediments, snapshot
               DICSED (I, KBM1) = DICSED (I, KBM1) + BENDIC 		!TDIC flux from sediments, accumulated
            End Do

!
#elif defined(MASS_CALC)

       !---Use mass-balance approach:  DICflux  =      SettlingPOC_to_seafloor        -   CH4outfluxFromSediments  -  POCburial
       !                             mmolC/m^3/d    [JPCOC(I,1)+JPCOC(I,2)+JPCOC(I,3)    -BENCH4A(I) - BENCH4G()    -BURIALC(I) ]
       !                                            (3 labilities)   mgC/m^2/day                  gO2*/m^2/d            mgC/m^2/d
       !                                             or PCFWS(I)  in gC/m^2/d
       !                                    *1000.0/2.667  to convert from gO2/m^2/day to mgC/m^2/day and /12.0107 to change mgC to
       !! Note: DIAGENC(I) = depth integrated POC diagenesis rate, gC/m^2/d = JCX/1000
       !! JCX = basically the decomposition of POC in sediments aided by respiration of deposition feeder
!
       !BURIALC only considers burial in layer 2, right? Is this simplification a problem?
       !BURIALC(I) = (POC1+POC2+POC3)*W2
       !POC1,POC2,POC3 = !sediment POC G1,G2,G3 (mgC/m^3) at new time step
       !POC1 = (JPOC(I,1)*DLTS_H2+POC1TM1) /(1.+(XKPOC1+W2)*DLTS_H2)
       !W2 = burial velocity (m/d)

            Do I = 1, MLOC

#  if !defined (WET_DRY)
			IF(D(I) > 0.0_SP)THEN
#  else
			IF(ISWETN(I) == 1)THEN
#  endif
               BENDIC = (PCFWS(I)*1000.0_SP-(BENCH4A(I)+&
              & BENCH4G(I))*1000.0_SP/2.667_SP-BURIALC(I)) / 12.0_SP / &
              & (D(I)*DZ2D(I,KBM1))!: USE 12.0 instead of 12.0107 because rounding up gN:molN to 14.0
			ELSE
				BENDIC=0.0_SP
			ENDIF
               DTTIC (I, KBM1) = DTTIC (I, KBM1) + BENDIC / 86400.0_SP
			   !DICSED(I,KBM1) = BENDIC  !TDIC flux from sediments, snapshot
               DICSED (I, KBM1) = DICSED (I, KBM1) + BENDIC !TDIC flux from sediments, accumulated
            End Do
!
!
#endif
!
         End If
!
    !!NOTE: so far, if BENTHIC_FLUXES=true, we yet don't have BENDIC addition (no TDIC from sediments)
!

		 DO I=1, MLOC  !Set the DTTIC DICSED, DICGAS etc to zero for dry areas
#  if !defined (WET_DRY)
			IF(D(I) <= 0.0_SP)THEN
#  else
			IF(ISWETN(I) == 0)THEN
#  endif
				DO K=1,KBM1

					DTTIC(I,K)=0.0_SP 					!stopped change to DTTIC
					DICSED(I,K)=DICSED(I,K)+0.0_SP
					DICGAS(I,K)=DICGAS(I,K)+0.0_SP
					DICUPT (I, K) =DICUPT(I,K)+ 0.0_SP
					DICBMP (I, K) =DICBMP(I,K)+ 0.0_SP
					DICPRD (I, K) =DICPRD(I,K)+ 0.0_SP
					DICMNL (I, K) =DICMNL(I,K)+ 0.0_SP
					DICDEN (I, K) =DICDEN(I,K)+ 0.0_SP

				ENDDO
			ENDIF
		 ENDDO
         Return
      End Subroutine TOTINC
!
!
  !************************************************************************
  !**                 S U B R O U T I N E   A L K A L I N                **
  !************************************************************************
!
      Subroutine ALKALIN (DTTALK)
!
         Use MOD_HYDROVARS, Only: DZ,DZ2D, D !CURRENT DEPTH
!
         Use MOD_WQM, Only: B1, B2, ANC1, ANC2, P1NNF, P2, PN1, PN2, &
        & NT, DENNO3, MNLRDON, MNLLDON, BENNH4, BENNO3, PRECIP, &
        & ATMOS_LOADS, ATMNH4, ATMNO3
!
         Use MOD_CO2SYS, Only: ALKNH4, ALKNO3, ALKNIT, ALKDEN, ALKREM, &
        & ALKNH4SED, ALKNO3SED

         Implicit None

         Integer :: I, K
!
         Real (SP), Dimension (0:MTLOC, KBM1) :: DTTALK
         Real (SP) :: alk_nh4, alk_no3, alk_nit, alk_den, alk_rem
         Real (SP) :: alks_nh4sed, alks_no3sed
         Real (SP) :: NH4A1, NH4A2, NO3A1, NO3A2
!
         DTTALK = 0.0
!
         Do K = 1, KBM1
            Do I = 1, MLOC

#  if !defined (WET_DRY)
			IF(D(I) > 0.0_SP)THEN
#  else
			IF(ISWETN(I) == 1)THEN
#  endif
		  !
          !********* Algal nitrogen utilization
!
               NH4A1 = PN1 (I, K) * P1NNF (I, K) * ANC1 * B1 (I, K)
               NH4A2 = PN2 (I, K) * P2 (I, K) * ANC2 * B2 (I, K)
               NO3A1 = (1.-PN1(I, K)) * P1NNF (I, K) * ANC1 * B1 (I, K)
               NO3A2 = (1.-PN2(I, K)) * P2 (I, K) * ANC2 * B2 (I, K)

               alk_nh4 = - (NH4A1+NH4A2) / 14. * 1000.0 !ALGNH4 is the balance between NH4 uptake and sources by Predation and photorespiration
               alk_no3 = (NO3A1+NO3A2) / 14. * 1000.0
!
          !********* Nitrification  14 g N/mol N
!
               alk_nit = - 2. * NT (I, K) / 14. * 1000.0
!
          !********* Denitrification   !this one was missing in Cerco's model

               alk_den = - DENNO3 (I, K) / 14. * 1000.0 !OR ANDC*DENIT(I,K)/14.
!
          !********* Remineralization of organic DON   !this one was missing in Cerco's model
               alk_rem = (MNLRDON(I, K)+MNLLDON(I, K)) / 14. * 1000.0

          !********* Change in TALK

		ELSE
				alk_nh4 = 0.0_SP
				alk_no3 = 0.0_SP
				alk_nit = 0.0_SP
				alk_rem = 0.0_SP
				alk_den = 0.0_SP
		ENDIF

               DTTALK (I, K) = &
              & (alk_nh4+alk_no3+alk_nit+alk_den+alk_rem) / 86400.

               ALKNH4 (I, K) = ALKNH4 (I, K) + alk_nh4
               ALKNO3 (I, K) = ALKNO3 (I, K) + alk_no3
               ALKNIT (I, K) = ALKNIT (I, K) + alk_nit
               ALKDEN (I, K) = ALKDEN (I, K) + alk_den
               ALKREM (I, K) = ALKREM (I, K) + alk_rem

            End Do
         End Do
!
    !******** Sediment flux    of NO3- and NH4+
!
         Do I = 1, MLOC


#  if !defined (WET_DRY)
			IF(D(I) > 0.0_SP)THEN
#  else
			IF(ISWETN(I) == 1)THEN
#  endif

				alks_nh4sed = BENNH4 (I) / (D(I)*DZ2D(I,KBM1)) / 14. * 1000.0 !mmolN/m^3/day     !BENNH4 is in gN/m^2/day
				alks_no3sed = - BENNO3 (I) / (D(I)*DZ2D(I,KBM1)) / 14. * 1000.0 !mmolN/m^3/day     !BENNO3 is in gN/m^2/day

			ELSE
				alks_nh4sed=0.0_SP
				alks_no3sed=0.0_SP
			ENDIF

            DTTALK (I, KBM1) = DTTALK (I, KBM1) + alks_nh4sed / 86400. !mmolN/m^3/sec
            DTTALK (I, KBM1) = DTTALK (I, KBM1) + alks_no3sed / 86400. !mmolN/m^3/sec

            ALKNH4SED (I, KBM1) = ALKNH4SED (I, KBM1) + alks_nh4sed !accumulated
            ALKNO3SED (I, KBM1) = ALKNO3SED (I, KBM1) + alks_no3sed !accumulated

         End Do
!
    !******** Atm deposition of NO3- and NH4+ : not included yet
         If (ATMOS_LOADS) Then
            Do I = 1, MLOC
#  if !defined (WET_DRY)
				IF(D(I) > 0.0_SP)THEN
#  else
				IF(ISWETN(I) == 1)THEN
#  endif
					DTTALK (I, 1) = DTTALK (I, 1) + PRECIP * ATMNH4 / &
              & (D(I)*DZ2D(I,1)) / 14. * 1000.0
					DTTALK (I, 1) = DTTALK (I, 1) - PRECIP * ATMNO3 / &
              & (D(I)*DZ2D(I,1)) / 14. * 1000.0
				ELSE

				ENDIF
            End Do
         End If

         Return
      End Subroutine ALKALIN
!
End Module MOD_DIC_TALK
