!mod_owq.F
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
Module MOD_OWQ
  !
  !********************************************************************************
  !**                             OWQ Model Setup (Ocean Water Color)            **
  !********************************************************************************
  !
      Use MOD_LIMS, Only: MTLOC, KBM1, MLOC
  !
      Use MOD_SIZES, Only: MGL
      Use MOD_FILEINFO, Only: KEI
  !
      Use MOD_CONTROL, Only: PAR, SERIAL
# if defined (MULTIPROCESSOR)
  !
      Use MOD_PAR, Only: NHN, NGID, HN_LST
# endif
  !
      Use MOD_PREC, Only: SP, SMOD
  !
      Use MOD_ZOOP, Only: SZ, LZ
  !
      Real (SP) :: INITKE, INTKE, NXKEI
  !
      Real (SP) :: LAT, LONG
  !
      Real (SP) :: GSTAR440, PHISTAR676, BNSTARNTU, BNSTARCHL, &
     & PSTARINT, PSTAR440, PSTARCOLOR, SCDOM, EZEROINT, TCHL2CHL, &
     & DOFFSET, DEPTHKE
  !
  !
      Integer :: NWAVEL
  !
  !
      Real (SP), Allocatable, Dimension (:, :) :: COLOR, TURB
  !
      Real (SP), Save :: OPZ (9), G1 (9), G2 (9)
  !
      Real (SP), Dimension (100) :: ALAMB, GLAMB, PLAMB, PHILAMB, &
     & BLAMB, WAVEL, EZERO
  !
      Real (SP), Dimension (10) :: GMUZERO
  !
      Real (SP), Allocatable, Dimension (:, :) :: IAVG, IATBOT, KESS  !IAVG unit: Einstein/m^2/day
	  !Real (SP), Allocatable, Dimension (:, :) :: PARAD   
	  Real (SP), Allocatable, Dimension (:, :) :: TSS, PARAD 
  !
  !
Contains
  !
  !FUNCTION   OWQ(COLORS,CHL,TURBS,DEPTH,AVREFCOS)
  !SUBROUTINE LGHT_ATTN()			!calculate light attenuation
  !SUBROUTINE OWQ_ALLOC()
  !SUBROUTINE OWQ_DEALLOC()
  !SUBROUTINE OWQ_READ()
  !
  !
  !=========================================================================
  !
  !=========================================================================
      Function OWQ (COLORS, CHL, TURBS, DEPTH, AVREFCOS) Result &
     & (OWQ_OUT)!
    !
         Implicit None
         Real (SP) SUMI
         Real (SP), Dimension (100) :: ABWAT, ABCDOM, ABCHL, ABNAPM, &
        & ABTOT, SCAT, KDLAMB, EZ
         Real (SP) :: AG440, APHI676, AP440, B440
         Real (SP) :: COLORS, CHL, TURBS, AVREFCOS, PLR, DEPTH, &
        & REFSCAT, EZINT
         Integer :: I, J
         Real (SP) :: OWQ_OUT !
    !
         AG440 = GSTAR440 * COLORS
         APHI676 = PHISTAR676 * CHL
         AP440 = PSTARINT * (TURBS**PSTAR440) * (COLORS**PSTARCOLOR)
         B440 = TURBS * BNSTARNTU + CHL * BNSTARCHL
    !
    ! COMPUTE ABSORPTION AND SCATTTERING FOR EACH WAVELENGTH
    !
         Do I = 1, NWAVEL
            ABWAT (I) = ALAMB (I)
            ABCDOM (I) = GLAMB (I) * AG440
            ABCHL (I) = PHILAMB (I) * APHI676
            ABNAPM (I) = PLAMB (I) * AP440
            ABTOT (I) = ABWAT (I) + ABCDOM (I) + ABCHL (I) + ABNAPM (I)
            SCAT (I) = BLAMB (I) * B440
         End Do
    !
    ! COMPUTE ATTENUATION FOR EACH WAVELENGTH
    !
         Do I = 1, NWAVEL
       !
       ! FIRST ESTIMATE ATTENUATION BASED ON DEPTH TO 1% LIGHT PENETRATION
            KDLAMB (I) = Sqrt &
           & (ABTOT(I)*ABTOT(I)+GMUZERO(2)*ABTOT(I)*SCAT(I)) / AVREFCOS
       !
       ! PERCENT LIGHT REMAINING
            PLR = Exp (-KDLAMB(I)*DEPTH)
       !
       ! USE LOOKUP TABLE TO REFINE ESTIMATE OF SCATTERING FUNCTION BASED ON
       ! PERCENT LIGHT REMAINING
       !
            REFSCAT = GMUZERO (9)
            Do J = 1, 8
               If (PLR >= OPZ(J) .And. PLR < OPZ(J+1)) Then
                  REFSCAT = GMUZERO (J)
                  Go To 1
               End If
            End Do
       !
1           Continue
       !
       ! RE-EVALUATE ATTENUATION BASED ON REFINED SCATERING FUNCTION
            KDLAMB (I) = Sqrt &
           & (ABTOT(I)*ABTOT(I)+REFSCAT*ABTOT(I)*SCAT(I)) / AVREFCOS
       !
       ! COMPUTE REMAINING SPECTRAL IRRADIANCE
            EZ (I) = EZERO (I) * Exp (-KDLAMB(I)*DEPTH)
       !
         End Do
    !
    ! INTEGRATE PAR REMAINING AT THE REFERENCE DEPTH
         SUMI = 0
         Do I = 2, NWAVEL - 1
            SUMI = SUMI + EZ (I)
         End Do
         EZINT = 2.5 * (EZ(1)+EZ(NWAVEL)) + 5. * SUMI
    !
    ! COMPUTE DIFFUSE ATTENUATION COEFFICIENT
         OWQ_OUT = - LOG (EZINT/EZEROINT) / DEPTH
    !
         Return
      End Function OWQ
  !
  !
  !=========================================================================
  !
  !=========================================================================
  !
      Subroutine OWQ_READ
    !USE MOD_LIMS, ONLY: MLOC
         Use MOD_WQM, Only: SPVARKE, PRINTKE
         Implicit None
    !
    !
         Real (SP), Allocatable :: RTMP2 (:, :), RTMP21 (:, :)
         Real (SP) :: SUMI
    !
         Integer :: I, K
    !
         Read (KEI, 1032)
         Read (KEI, 1038) INTKE, INITKE
         Read (KEI, 1037) GSTAR440, PHISTAR676, BNSTARNTU, BNSTARCHL, &
        & PSTARINT, PSTAR440, PSTARCOLOR, SCDOM
         Read (KEI, 1036) LAT, LONG, NWAVEL, TCHL2CHL, DOFFSET, DEPTHKE
    !
    !
    !
    ! INITIALIZE KEI INTERVAL, KESS
    !
         NXKEI = INTKE / 24.
         Do K = 1, KBM1
            Do I = 1, MLOC
               KESS (I, K) = INITKE
            End Do
         End Do
    !
    ! CONVERT FROM DECIMAL DEGREES TO RADIANS
         LAT = 6.2832 * LAT / 360.
         LONG = 6.2832 * LONG / 360.
    !
         Read (KEI, 1032)
         Read (KEI,*) (WAVEL(I), ALAMB(I), GLAMB(I), PHILAMB(I), &
        & PLAMB(I), BLAMB(I), EZERO(I), I=1, NWAVEL)
    !
         Read (KEI, 1082) SPVARKE, PRINTKE
         If (SPVARKE == 'CONSTANT') Then
            Read (KEI, 1087) COLOR (1, 1), TURB (1, 1)
            Do K = 1, KBM1
               Do I = 1, MTLOC
                  COLOR (I, K) = COLOR (1, 1)
                  TURB (I, K) = TURB (1, 1)
               End Do
            End Do
         Else
            Allocate (RTMP2(MGL, KBM1))
            RTMP2 = 0.0
            Allocate (RTMP21(MGL, KBM1))
            RTMP21 = 0.0
            Do K = 1, KBM1
               Read (KEI, 1087) (COLOR(I, K), TURB(I, K), I=1, MGL)
            End Do
       !
            If (SERIAL) Then
               COLOR = RTMP2
               TURB = RTMP21
            End If
#        if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do K = 1, KBM1
                  Do I = 1, MLOC
                     COLOR (I, K) = RTMP2 (NGID(I), K)
                     TURB (I, K) = RTMP21 (NGID(I), K)
                  End Do
                  Do I = 1, NHN
                     COLOR (I+MLOC, K) = RTMP2 (HN_LST(I), K)
                     TURB (I+MLOC, K) = RTMP21 (HN_LST(I), K)
                  End Do
               End Do
            End If
#        endif
            Deallocate (RTMP2, RTMP21)
         End If
    !
    !
    ! INTEGRATE AVERAGE QUANTUM SPECTRUM OVER VISIBLE WAVELENGTHS
    !
         SUMI = 0
         Do I = 2, NWAVEL - 1
            SUMI = SUMI + EZERO (I)
         End Do
         EZEROINT = 2.5 * (EZERO(1)+EZERO(NWAVEL)) + 5. * SUMI
1032     Format (/)
1036     Format (/ / 8 X, 2 F8.1, I8, 3 F8.1)
1037     Format (/ / 8 X, 8 F8.1)
1038     Format (/ / 8 X, 8 F8.1)
1082     Format (/ / 8 X, 2 A8 //)
1087     Format (8 X, 4 F8.0)
    !
      End Subroutine OWQ_READ
  !
      Subroutine LGHT_ATTN
         Use MOD_WQM, Only: B1, B2, B3, CHLCMN1, CHLCMN2, CHLCMN3, I0, &
        & JDAY, LDOC, LPOC, RDOC, RPOC, SREFLECT, SSI ,JDAY!,
!
         Use MOD_HYDROVARS, Only: DZ, DZ2D, D
    !
         Implicit None
         Real (SP) :: ITHTAA
         Integer :: I, K
         Real (SP), Parameter :: TWOPI = 6.2832

         Real (SP) :: DOY, HOUR, SOLDATE, DECL, TAU, CTHTAA, THTAW, &
        & AVREFCOS
         !Real (SP) :: TOTCHL, TOTPOC, TSS, TOTDOC, IATTOP, OPTDEPTH
		 Real (SP) :: TOTCHL, TOTPOC, TOTDOC, IATTOP, OPTDEPTH 
    !

    !
         OPZ (1:9) = (/ 0.0, 0.01, 0.025, 0.033, 0.05, 0.1, 0.2, 0.286, &
        & 0.4 /)
         G1 (1:9) = (/ 0.425, 0.425, 0.364, 0.3383, 0.2996, 0.2175, &
        & 0.1232, 0.073167, 0.01546 /)
         G2 (1:9) = (/ 0.19, 0.19, 0.1396, 0.1187, 0.0875, 0.02288, - &
        & 0.048, - 0.08367, - 0.12386 /)
    !
    ! CONVERT MODEL JULIAN DAY TO DAY AND HOUR
    !
         DOY = SMOD (JDAY, 365._SP)
         HOUR = 24. * (DOY-AINT(DOY))
         DOY = AINT (DOY) + DOFFSET
    !
    ! COMPUTE COSIN OF IN-WATER SOLAR ANGLE
    !
    ! SOLAR DATE
         SOLDATE = TWOPI * DOY / 365.
    !
    ! DECLINATION
         DECL = TWOPI * (0.39637-(22.9133*Cos(SOLDATE))+&
        & (4.02543*Sin(SOLDATE))-(0.3872*Cos(2.*SOLDATE))+&
        & (0.052*Sin(2*SOLDATE))) / 360.
    !
    ! TIME OF DAY
         TAU = TWOPI * HOUR / 24.
    !
    ! COSIN OF THETA IN AIR
         CTHTAA = ((Sin(LAT))*(Sin(DECL))) - &
        & ((Cos(LAT))*(Cos(DECL))*(Cos(TAU)))
    !
    ! INCIDENCE OF THETA IN AIR
         ITHTAA = Acos (CTHTAA)
    !
    ! THETA IN WATER
    ! Here 1.33 is the water refraction index of light
         THTAW = Asin (Sin(ITHTAA)/1.33)
    !
    ! REFLECTANCE
         If (ITHTAA > TWOPI/4.) Then
            SREFLECT = 100.
         Else
            SREFLECT = Exp (0.6148+ITHTAA**3)
         End If
    !
    ! SCATTERING FUNCTION
         AVREFCOS = Cos (THTAW)
         GMUZERO (1) = 0.198
         Do I = 2, 9
            GMUZERO (I) = AVREFCOS * G1 (I) - G2 (I)
         End Do
    !
    ! EVALUATE ATTENUATION FOR EACH SURFACE CELL
    !
       !  Do I = 1, MLOC
       !     TOTCHL = (B1(I, 1)*CHLCMN1+B2(I, 1)*CHLCMN2+B3(I, &
       !    & 1)*CHLCMN3) * TCHL2CHL
       !!
       !! EMPIRICAL FUNCTIONS FOR COLOR AND TURBIDITY
       !!
       !     TOTPOC = B1 (I, 1) + B2 (I, 1) + B3 (I, 1) + SZ (I, 1) + LZ (I, 1) + LPOC (I, 1) + RPOC (I, 1)
       !     TSS(I, 1) = Max (SSI(I, 1)+2.5*TOTPOC, 1.0) !AN
       !     !TSS = Max (SSI(I, 1)+2.5*TOTPOC, 1.0)
       !     TURB (I, 1) = Min (0.4*TSS(I, 1), 12.) !AN
       !     !TURB (I, 1) = Min (0.4*TSS, 12.)
       !     TOTDOC = Min (35., LDOC(I, 1)+RDOC(I, 1))!  used to be 25
       !     COLOR (I, 1) = Max (12.2*Exp(0.14*TOTDOC), 10.)
       !!
       !     KESS (I, 1) = OWQ (COLOR(I, 1), TOTCHL, TURB(I, 1), DEPTHKE, AVREFCOS)
       !!
       !! ASSIGN SURFACE ATTENUATION TO SUB-SURFACE CELLS
       !!
       !     Do K = 2, KBM1
       !        KESS (I, K) = KESS (I, 1)
       !     End Do
       !  End Do

         Do I = 1, MLOC
       !     TOTCHL = (B1(I, 1)*CHLCMN1+B2(I, 1)*CHLCMN2+B3(I, &
       !    & 1)*CHLCMN3) * TCHL2CHL
			TOTCHL = (B1(I, 1)*CHLCMN1+B2(I, 1)*CHLCMN2) * TCHL2CHL 
       !
       ! EMPIRICAL FUNCTIONS FOR COLOR AND TURBIDITY
       !
            Do K = 1, KBM1
               !TOTPOC = B1 (I, K) + B2 (I, K) + B3 (I, K) + SZ (I, K) + LZ (I, K) + LPOC (I, K) + RPOC (I, K)
			   TOTPOC = B1 (I, K) + B2 (I, K) + SZ (I, K) + LZ (I, K) + LPOC (I, K) + RPOC (I, K) 
               TSS (I, K) = Max (SSI(I, K)+2.5*TOTPOC, 1.0) !original
               !TSS (I, K) = Max (SSI(I, K)+2.5*TOTPOC, 0.) 
               TURB (I, K) = Min (0.4*TSS(I, K), 12.)  !original
               !TURB (I, K) = Max (0.4*TSS(I, K), 0.)  !Avoid limit
               !TURB (I, K) = Max ((TSS(I, K)/1.6495)**(1/1.0055), 0.) 
               TOTDOC = Min (35., LDOC(I, K)+RDOC(I, K))!  used to be 25
               COLOR (I, K) = Max (12.2*Exp(0.14*TOTDOC), 10.) !original
               !COLOR (I, K) = Max (12.2*Exp(0.14*TOTDOC), 5.) 
       !
       ! ASSIGN SURFACE ATTENUATION TO SUB-SURFACE CELLS
       !
              if (K==1) then
               KESS (I, 1) = OWQ (COLOR(I, 1), TOTCHL, TURB(I, 1), DEPTHKE, AVREFCOS)
              end if
               KESS (I, K) = KESS (I, 1)
            End Do
         End Do
    !
    ! *** Compute irradiance in each cell
    !
         Do I = 1, MLOC
            IATTOP = I0 * (1.-SREFLECT/100.)
       !Added print below for irradiance at top of wc
       !    IF(MSR.and.I.eq.MLOC) WRITE(*,*)'IATTOP', IATTOP
            OPTDEPTH = KESS (I, 1) * D (I) * DZ2D (I,1)
            IATBOT (I, 1) = IATTOP * Exp (-OPTDEPTH)
            IAVG (I, 1) = (IATTOP-IATBOT(I, 1)) / OPTDEPTH
            PARAD(I,1)=IATTOP;
       !
            Do K = 2, KBM1
               IATTOP = IATBOT (I, K-1)
               OPTDEPTH = KESS (I, K) * D (I) * DZ2D (I,K)
               IATBOT (I, K) = IATTOP * Exp (-OPTDEPTH)
               IAVG (I, K) = (IATTOP-IATBOT(I, K)) / OPTDEPTH
			   PARAD(I,K)=IATTOP;
            End Do
         End Do
    !
         Return
      End Subroutine LGHT_ATTN
  !
  !=======================================================================
  !
  !=======================================================================
  !
      Subroutine OWQ_ALLOC
    !
         Allocate (COLOR(MTLOC, KBM1))
         COLOR = 0.0
         Allocate (TURB(MTLOC, KBM1))
         TURB = 0.0
		 Allocate (TSS(MTLOC, KBM1))
         TSS = 0.0 
		 Allocate (PARAD(MTLOC, KBM1)) 
         PARAD = 0.0 
         Allocate (IAVG(MTLOC, KBM1))
         IAVG = 0.0 !Moved this from mod_wqm.F
         Allocate (IATBOT(MTLOC, KBM1))
         IATBOT = 0.0 !Moved this from mod_wqm.F
         Allocate (KESS(MTLOC, KBM1))
         KESS = 0.0 !Moved this from mod_wqm.F
      End Subroutine OWQ_ALLOC
  !
  !=======================================================================
  !
  !=======================================================================
  !
      Subroutine OWQ_DEALLOC
         If (ALLOCATED(COLOR)) DEALLOCATE (COLOR)
         If (ALLOCATED(TURB)) DEALLOCATE (TURB)
		 If (ALLOCATED(TSS)) DEALLOCATE (TSS) 
		 If (ALLOCATED(PARAD)) DEALLOCATE (PARAD) 
         If (ALLOCATED(IAVG)) DEALLOCATE (IAVG)!Moved this to here from wqm_main.F
         If (ALLOCATED(IATBOT)) DEALLOCATE (IATBOT)!Moved this to here from wqm_main.F
         If (ALLOCATED(KESS)) DEALLOCATE (KESS)!Moved this to here from wqm_main.F
      End Subroutine OWQ_DEALLOC
  !
  !
End Module MOD_OWQ
