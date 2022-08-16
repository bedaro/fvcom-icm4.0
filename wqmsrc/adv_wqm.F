!adv_wqm.F
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!=============================================================================!
!
!Subroutine ADV_WQM()
!
!
Subroutine ADV_WQM

# if defined (WET_DRY)
	USE  MOD_WD, ONLY: ISWETN,ISWETNT
# endif
      Use MOD_LIMS, Only: KB, MYID, NPROCS, NCV, NUMPNT, MLOC, MTLOC, &
     & NTLOC, KBM1, NCV_I
      Use MOD_PREC, Only: SP
      Use MOD_TGE, Only: ISONB, NTSN, NBSN, NIEC, NTRG, XIJE, YIJE, &
     & DLTXE, DLTYE
!
      Use MOD_HYDROVARS, Only: ZERO, VX, VY, ART1, ART2, DZ,DZ2D, DT, DT1, &
     & DTFA, UU, VV, WTS, VISCOFH, THOUR, XFLUX_OBC_WQM !LB added XFLUX_OBC_WQM may17,2016
  !Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: PAR, MSR !	&           !!TRUE IF MULTIPROCESSOR RUN
!
      Use MOD_WQM, Only: NAC, AC, C2, C1, C2F, DLT, DTM, NCP, XYDF, &
     & ELTMS, ADVFLUX !LB added ADVFLUX jan26 2016
      Use MOD_WQMINIT, Only: XYDFU
      Use MOD_BCS, Only: INOPNT, PNT_TM, WVQDIST, PWQDIS, PDWQDIS, &
     & PQDIS, PDQDIS !!DISCHARGE at non-point source node
!
      Use MOD_BCMAP, Only: IOBCN, I_OBC_N !Added by LB 17may2016
!
# if defined (MULTIPROCESSOR)
      Use MOD_PAR, Only: NBN, BN_LOC, BN_MLT, BNC, NODE_MATCH_ALL !, & !
!
# endif
      Implicit None
      Real (SP), Allocatable, Dimension (:, :, :) :: XFLUX, XFLUX_ADV
      Real (SP), Dimension (0:MTLOC) :: PUPX, PUPY, PVPX, PVPY !LB 6may2016 turn MLOC into 0:MTLOC
      Real (SP), Dimension (0:MTLOC) :: PFPX, PFPY, PFPXD, PFPYD, &
     & VISCOFF !LB 6may2016 turn MLOC into 0:MTLOC
      Real (SP), Dimension (3*(NTLOC)) :: DTIJ !LB commented
  !REAL(SP), DIMENSION(3*(NTLOC),KBM1) :: DTIJ  !LB added to make as in HYD code -04/29/2016
      Real (SP), Dimension (3*(NTLOC), KBM1) :: UVN
      Real (SP) :: FFD, FF1, X11, Y11, X22, Y22, X33, Y33, XI, YI
      Real (SP) :: DXA, DYA, DXB_A (NCV_I), DYB_A (NCV_I), DXB_B &
     & (NCV_I), DYB_B (NCV_I), FIJ1, FIJ2, UN, DYB, DXB
      Real (SP) :: TXX, TYY, FXX, FYY, VISCOF, EXFLUX, TEMP
      Real (SP) :: FACT, FM1
      Real (SP) :: TT, TTIME, STPOINT, STPOINT1, STPOINT2
      Integer :: I, I1, I2, IA, IB, IC, J, J1, J2, JTMP, K, JJ, JCON, &
     & II
      Real (SP) :: C21MIN, C21MAX, C22MIN, C22MAX

      Logical XYDFUISON
      Real (SP) :: TMP1 (MLOC), TMP2 (MLOC)
  !
      Integer :: NJ1, NJ2, NJ3, EXTRA
  !TYKIM for nonpoint sources
      Integer :: L1, L2, IERR
      Real (SP) :: UFACT, DFACT
      Real (SP) :: auxLBadv
!
  !!ADVFLUX = 0.0  !LB initialize
  !------------------------------------------------------------------------------
  !
  !  Make XFLUX allocated
      Allocate (XFLUX(0:MTLOC, KB, NCP))
      Allocate (XFLUX_ADV(0:MTLOC, KB, NCP))
  !  Convert this to logic below
  !   FACT = 0.0
  !   FM1  = 1.0
  !   IF(XYDFU == ' ON')THEN
  !     FACT = 1.0
  !     FM1  = 0.0
  !   ENDIF
      XYDFUISON = .False.
      If (XYDFU == ' ON') XYDFUISON = .True.
  !
  ! Precalculate 
      Do I = 1, MLOC
         TMP1 (I) = DLT / DT (I) / ART1 (I)
         TMP2 (I) = DT (I) / DTFA (I)

      End Do
  !
  ! Precalculate more 
      Do I = 1, NCV_I
         IA = NIEC (I, 1)
         IB = NIEC (I, 2)
         XI = 0.5_SP * (XIJE(I, 1)+XIJE(I, 2))
         YI = 0.5_SP * (YIJE(I, 1)+YIJE(I, 2))
         DXB_A (I) = XI - VX (IA)
         DYB_A (I) = YI - VY (IA)
         DXB_B (I) = XI - VX (IB)
         DYB_B (I) = YI - VY (IB)
      End Do
  !
  !
  !--Initialize Fluxes-----------------------------------------------------------
  !
      XFLUX = 0.0
      XFLUX_ADV = 0.0
  !
  !--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------
  !
      Do I = 1, NCV
         I1 = NTRG (I)
         DTIJ (I) = DT1 (I1)!!!!commented by LB, as it is in FVCOM-HYD code -Arp29,2016
         Do K = 1, KBM1
            UVN (I, K) = VV (I1, K) * DLTXE (I) - UU (I1, K) * DLTYE &
           & (I)
         End Do
      End Do
  !
  !
  !--Calculate the Advection and Horizontal Diffusion Terms----------------------
  !
      Do JCON = 1, NAC
         II = AC (JCON)
         Do K = 1, KBM1
            PFPX = 0.0
            PFPY = 0.0

            PFPXD = 0.0 !LBtest
            PFPYD = 0.0
        !
            Do I = 1, MLOC
               Do J = 1, NTSN (I) - 1
                  I1 = NBSN (I, J)
                  I2 = NBSN (I, J+1)

#if defined (WET_DRY)

				IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN

					FF1 = 0.5_SP * (C2(I, K, II)+C2(I2, K, II))
				ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
					FF1 = 0.5_SP * (C2(I1, K, II)+C2(I, K, II))

				ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
					FF1 = 0.5_SP * (C2(I, K, II)+C2(I, K, II))
				ELSE
					FF1 = 0.5_SP * (C2(I1, K, II)+C2(I2, K, II))
				ENDIF
#else

				FF1 = 0.5_SP * (C2(I1, K, II)+C2(I2, K, II))
#endif

                  PFPX (I) = PFPX (I) + FF1 * (VY(I1)-VY(I2))
                  PFPY (I) = PFPY (I) + FF1 * (VX(I2)-VX(I1))

               End Do
               PFPX (I) = PFPX (I) / ART2 (I)
               PFPY (I) = PFPY (I) / ART2 (I)

            End Do
        !
            Do I = 1, MLOC
               VISCOFF (I) = VISCOFH (I, K)
            End Do
        !
            Do I = 1, NCV_I
               IA = NIEC (I, 1)
               IB = NIEC (I, 2)
           !
           ! 20% faster total run time -
			!  split UN > 0 and UN < 0 cases, and explicit loops for MAXVAL
           !DXA=XI-VX(IA)
           !DYA=YI-VY(IA)
           !DXB=XI-VX(IB)
           !DYB=YI-VY(IB)
           !FIJ1=C2(IA,K,II)+DXA*PFPX(IA)+DYA*PFPY(IA)
           !FIJ2=C2(IB,K,II)+DXB*PFPX(IB)+DYB*PFPY(IB)
!
           !C21MIN=MINVAL(C2(NBSN(IA,1:NTSN(IA)-1),K,II))
           !C21MIN=MIN(C21MIN, C2(IA,K,II))
           !C21MAX=MAXVAL(C2(NBSN(IA,1:NTSN(IA)-1),K,II))
           !C21MAX=MAX(C21MAX, C2(IA,K,II))
           !C22MIN=MINVAL(C2(NBSN(IB,1:NTSN(IB)-1),K,II))
           !C22MIN=MIN(C22MIN, C2(IB,K,II))
           !C22MAX=MAXVAL(C2(NBSN(IB,1:NTSN(IB)-1),K,II))
           !C22MAX=MAX(C22MAX, C2(IB,K,II))
           !IF(FIJ1 < C21MIN) FIJ1=C21MIN
           !IF(FIJ1 > C21MAX) FIJ1=C21MAX
           !IF(FIJ2 < C22MIN) FIJ2=C22MIN
           !IF(FIJ2 > C22MAX) FIJ2=C22MAX
!
           !UN=UVN(I,K)
!
           !VISCOF=XYDF*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5 + FM1)
!
           !TXX=0.5*(PFPXD(IA)+PFPXD(IB))*VISCOF
           !TYY=0.5*(PFPYD(IA)+PFPYD(IB))*VISCOF
!
           !!FXX=-DTIJ(I)*TXX*DLTYE(I)
           !!FYY= DTIJ(I)*TYY*DLTXE(I)
           !!LB added commented above and added below to have DTIJ as in HYD code -04/29/2016
           !FXX=-DTIJ(I,K)*TXX*DLTYE(I)
           !FYY= DTIJ(I,K)*TYY*DLTXE(I)
!
           !!EXFLUX=-UN*DTIJ(I)*((1.0+SIGN(1.0,UN))*FIJ2+                &
                !!        (1.0-SIGN(1.0,UN))*FIJ1)*0.5+FXX+FYY
           !!LB added commented above and added below to have DTIJ as in HYD code -04/29/2016
           !EXFLUX=-UN*DTIJ(I,K)*((1.0+SIGN(1.0,UN))*FIJ2+                &
                !        (1.0-SIGN(1.0,UN))*FIJ1)*0.5+FXX+FYY
           !
!
               UN = UVN (I, K)
           !
           !  Convert FACT and FM1 to logic
               If (XYDFUISON) Then
                  VISCOF = XYDF * ((VISCOFF(IA)+VISCOFF(IB))*0.5)
               Else
                  VISCOF = XYDF
               End If
           !   VISCOF=XYDF*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5 + FM1)
#if defined (WET_DRY)
			!Allow no diffusive flux between a wet node and a dry node
			IF(ISWETN(IA)==1.AND.ISWETN(IB)==1)THEN
# endif

               TXX = 0.5 * (PFPX(IA)+PFPX(IB)) * VISCOF
               TYY = 0.5 * (PFPY(IA)+PFPY(IB)) * VISCOF
# if defined (WET_DRY)
			ELSE
				TXX = 0.0_SP
				TYY = 0.0_SP
			ENDIF
# endif
           !
               FXX = - DTIJ (I) * TXX * DLTYE (I)
               FYY = DTIJ (I) * TYY * DLTXE (I)

               If (UN .Ge. 0) Then
                  IC = IB
                  DYB = DYB_B (I)
                  DXB = DXB_B (I)
               Else
                  IC = IA
                  DYB = DYB_A (I)
                  DXB = DXB_A (I)
               End If
               C22MIN = C2 (IC, K, II)
               C22MAX = C2 (IC, K, II)
               Do J = 1, NTSN (IC) - 1
                  TEMP = C2 (NBSN(IC, J), K, II)
                  C22MIN = Min (C22MIN, TEMP)
                  C22MAX = Max (C22MAX, TEMP)
               End Do
               FIJ2 = C2 (IC, K, II) + DXB * PFPX (IC) + DYB * PFPY &
              & (IC)
               If (FIJ2 < C22MIN) FIJ2 = C22MIN
               If (FIJ2 > C22MAX) FIJ2 = C22MAX
           !caps above for minmax criteria
               EXFLUX = - UN * DTIJ (I) * FIJ2 + FXX + FYY
           !!LB added commented above and added below to have DTIJ as in HYD code -04/29/2016
           !EXFLUX = - UN * DTIJ (I,K) * FIJ2 + FXX + FYY
!
           ! 
               XFLUX (IA, K, II) = XFLUX (IA, K, II) + EXFLUX
               XFLUX (IB, K, II) = XFLUX (IB, K, II) - EXFLUX
           !
               XFLUX_ADV (IA, K, II) = XFLUX_ADV (IA, K, II) + &
              & (EXFLUX-FXX-FYY)
               XFLUX_ADV (IB, K, II) = XFLUX_ADV (IB, K, II) - &
              & (EXFLUX-FXX-FYY)
            End Do
         End Do
      End Do
  !
  !
  !-Accumulate Fluxes at Boundary Nodes
  !
# if defined (MULTIPROCESSOR)
      If (PAR) Then
     ! DO MULTIPLE VARIABLES AT ONCE
     !    DO JCON=1,NAC
     !      II=AC(JCON)
     !      CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MTLOC,KB,MYID,NPROCS,XFLUX(:,:,II))
     !    ENDDO
#if(0)
         EXTRA = Mod (NAC, 3)
         If (EXTRA .Eq. 1) Then
            NJ1 = AC (1)
            Call NODE_MATCH (0, NBN, BN_MLT, BN_LOC, BNC, MTLOC, KB, &
           & MYID, NPROCS, XFLUX(:, :, NJ1))
         Else If (EXTRA .Eq. 2) Then
            NJ1 = AC (1)
            NJ2 = AC (2)
            Call NODE_MATCH (0, NBN, BN_MLT, BN_LOC, BNC, MTLOC, KB, &
           & MYID, NPROCS, XFLUX(:, :, NJ1), XFLUX(:, :, NJ2))
         End If
     !
         Do JCON = EXTRA + 1, NAC - 2, 3
            NJ1 = AC (JCON+0)
            NJ2 = AC (JCON+1)
            NJ3 = AC (JCON+2)
            Call NODE_MATCH (0, NBN, BN_MLT, BN_LOC, BNC, MTLOC, KB, &
           & MYID, NPROCS, XFLUX(:, :, NJ1), XFLUX(:, :, NJ2), XFLUX(:, &
           & :, NJ3))
         End Do
#else
     !  DO ALL AT ONCE
         Call NODE_MATCH_ALL (0, NBN, BN_MLT, BN_LOC, BNC, MTLOC, KB, &
        & NCP, MYID, NPROCS, XFLUX)
#endif
      End If
# endif

  !
  !!Added by LB in May 17, 2016, to use this xflux_obc_wqm in bcond_wqm.F (as done in HYD)
      Do JCON = 1, NAC
         II = AC (JCON)
         Do K = 1, KBM1
            If (IOBCN > 0) Then
               Do I = 1, IOBCN
                  I1 = I_OBC_N (I)
                  XFLUX_OBC_WQM (I, K, II) = XFLUX_ADV (I1, K, II)
               End Do
            End If
         End Do
      End Do
  !!end addition LB
!

  !NON-POINT SOURCE INTERPOLATION
      If (NUMPNT > 0) Then
         Allocate (PWQDIS(NUMPNT, NCP))
         PWQDIS = ZERO
         Allocate (PQDIS(NUMPNT))
         PQDIS = ZERO
         Call BRACKET (PNT_TM, THOUR, L1, L2, DFACT, UFACT, IERR)
         PWQDIS (:, :) = UFACT * PDWQDIS (:, :, L1) + DFACT * PDWQDIS &
        & (:, :, L2)
         PQDIS (:) = UFACT * PDQDIS (:, L1) + DFACT * PDQDIS (:, L2)

         PQDIS (:) = PQDIS (:) * Tanh (ELTMS/FLOAT(86400))

      End If
  !
      Do JCON = 1, NAC
         II = AC (JCON)
     !
     !
     !--Calculate the Vertical Terms------------------------------------------------
     !
         Do I=1,MLOC
#if defined (WET_DRY)
			IF(iswetn(i)*iswetnt(i) == 1) then	!Added wet_dry here
#endif
            Do K = 1, KBM1
               If (K == 1) Then
                  TEMP = - WTS (I, K+1) * (C2(I, K, II)*DZ2D(I,K+1)+C2(I, &
                 & K+1, II)*DZ2D(I,K)) / (DZ2D(I,K)+DZ2D(I,K+1))
               Else If (K == KBM1) Then
                  TEMP = WTS (I, K) * (C2(I, K, II)*DZ2D(I,K-1)+C2(I, K-1, &
                 & II)*DZ2D(I,K)) / (DZ2D(I,K)+DZ2D(I,K-1))
               Else
                  TEMP = WTS (I, K) * (C2(I, K, II)*DZ2D(I,K-1)+C2(I, K-1, &
                 & II)*DZ2D(I,K)) / (DZ2D(I,K)+DZ2D(I,K-1)) - WTS (I, K+1) * (C2(I, &
                 & K, II)*DZ2D(I,K+1)+C2(I, K+1, II)*DZ2D(I,K)) / &
                 & (DZ2D(I,K)+DZ2D(I,K+1))
               End If
           !
           !
           !--Total Fluxes ---------------------------------------------------------------
           !
               If (ISONB(I) == 2) Then
                  XFLUX (I, K, II) = TEMP * ART1 (I) / DZ2D (I,K)
               Else
                  XFLUX (I, K, II) = XFLUX (I, K, II) + TEMP * ART1 (I) &
                 & / DZ2D (I,K)
               End If
            End Do
#if defined (WET_DRY)
			END IF
#endif
         End Do
     !
     !TYKIM
     !===============================================================================
     !POINT SOURCE FLUX
     !===============================================================================
         If (NUMPNT > 0) Then
            Do J = 1, NUMPNT
               JJ = INOPNT (J)
               Do K = 1, KBM1
                  STPOINT = PWQDIS (J, II)

                  XFLUX (JJ, K, II) = XFLUX (JJ, K, II) - PQDIS (J) * &
                 & WVQDIST (J, K) * STPOINT / DZ2D (JJ,K)

               End Do
            End Do
         End If
     !
     !
     !--Update Water Quality Variables--------------------------------
     !
         Do I = 1, MLOC
#if defined (WET_DRY)
            IF(iswetn(i)*iswetnt(i) == 1 )THEN
#endif
				Do K = 1, KBM1

					C1 (I, K, II) = (C2(I, K, II)-XFLUX(I, K, II)*TMP1(I)) * &
              & TMP2 (I)

				C2F (I, K, II) = Max (C1(I, K, II), 0.0)
!
				If (II .Eq. 4 .Or. II .Eq. 5 .Or. II .Eq. 6) C2F (I, K, &
              & II) = Max (C1(I, K, II), 0.003)!make sure there are seeds for ALG1,2 and 3
!
				!LB: save flux of horiz,vert adv + horiz diffusion in term ADVFLUX in units mass/m3/d
				auxLBadv = - XFLUX (I, K, II) * TMP1 (I) * TMP2 (I)
				auxLBadv = auxLBadv / DLT * 86400.0
				!ADVFLUX(I, K, II) =  auxLBadv !LB: snapshot
				ADVFLUX (I, K, II) = ADVFLUX (I, K, II) + auxLBadv !LB:accumulated
				End Do
#if defined (WET_DRY)
			ELSE
				DO K=1,KBM1
						C2F(I,K,II)=C1(I,K,II)		!No change for dry areas
						C2F (I, K, II) = Max (C1(I, K, II), 0.0)
						If (II .Eq. 4 .Or. II .Eq. 5 .Or. II .Eq. 6) C2F (I, K, &
					& II) = Max (C1(I, K, II), 0.003)!make sure there are seeds for ALG1,2 and 3
				ENDDO
			ENDIF
#endif
         End Do !End of I loop
      End Do !End of JCON  loop (i.e. II)


      If (NUMPNT > 0) DEALLOCATE (PQDIS, PWQDIS)
      Deallocate (XFLUX)
      Deallocate (XFLUX_ADV)
      Return
End Subroutine ADV_WQM
!==============================================================================!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
