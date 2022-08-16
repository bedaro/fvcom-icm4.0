!bcond_wqm.F
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
!subroutine BCOND_WQM()
!subroutine BCOND_NUT_PERTURBATION()
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
Subroutine BCOND_WQM
  !==============================================================================|
  !   Set Boundary Conditions on Water Quality                                   |
  !==============================================================================|
  !
  !------------------------------------------------------------------------------|
      Use MOD_SIZES, Only: NCP, NOBTY
      Use MOD_PREC, Only: SP
      Use MOD_LIMS, Only: KBM1, NUMQBC, IINT
  !
      Use MOD_TGE, Only: NTSN, NBSN
  !
      Use MOD_BCMAP, Only: IOBCN, I_OBC_N, WQOBC, NUT_TM
  !
      Use MOD_HYDROVARS, Only: GRAV, ZERO, ART1, DZ,DZ2D, H, D, DT, &
     & UARD_OBCN, XFLUX_OBC, THOUR, RHO, XFLUX_OBC_WQM !
  !
  !Took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: DTI, INFLOW_TYPE, POINT_ST_TYPE, MSR
      Use MOD_BCS, Only: INODEQ, ICELLQ, QBC_TM, DWDIS, WDIS, N_ICELLQ
      Use MOD_WQM, Only: C2, C2F, NAC, AC, JDAY, SALT, DLT !
      Use MOD_PREC, Only: SP
      Use MOD_OBCS, Only: NEXT_OBC, DLTN_OBC
      Implicit None
      Integer :: I, J, K, J1, J11, J22, NCON2, JCON, II
      Real (SP) :: UFACT, FACT
      Real (SP) :: T2D, T2D_NEXT, T2D_OBC, XFLUX2D, TMP, TMAX, TMIN, &
     & RAMP_TS
      Real (SP) :: TAasFsalt !
      Integer :: stepsIn1day !
      Integer :: L1, L2, IERR
      Real (SP) :: CC, CP, PERT, PERT_NEXT
      Real (SP), Allocatable :: TTMP (:, :)
      Real (SP) :: TY2D_N, TYC2F
      Real (SP), Allocatable :: WQOBC_TMP (:, :, :)
  !------------------------------------------------------------------------------|
  !
  !
      stepsIn1day = 86400.0_SP / DLT
  !
  !--SET CONDITIONS FOR FRESH WATER INFLOW---------------------------------------|
  !
      If (NUMQBC > 0) Then
         Allocate (WDIS(NUMQBC, NCP))
         WDIS = ZERO
         Call BRACKET (QBC_TM, JDAY, L1, L2, FACT, UFACT, IERR)
         WDIS (:, :) = DWDIS (:, :, L1)!UFACT*DWDIS(:,L1) + FACT*DWDIS(:,L2)
      End If
  !
      If (POINT_ST_TYPE == 'specified') Then
         If (NUMQBC > 0) Then
            If (INFLOW_TYPE == 'node') Then
               Do JCON = 1, NAC
                  II = AC (JCON)
                  Do K = 1, KBM1
                     Do I = 1, NUMQBC
                        J11 = INODEQ (I)
                        C2F (J11, K, II) = WDIS (I, II)
                     End Do
                  End Do
               End Do
            Else If (INFLOW_TYPE == 'edge') Then
               Do JCON = 1, NAC
                  II = AC (JCON)
                  Do K = 1, KBM1
                     Do I = 1, NUMQBC
                        J11 = N_ICELLQ (I, 1)!first node of the element edge
                        J22 = N_ICELLQ (I, 2)!second node of the element edge
                        C2F (J11, K, II) = WDIS (I, II)!
                    !	it is distirbuted evenly by the two nodes connecting the edge?
                        C2F (J22, K, II) = WDIS (I, II)! !
                    !	it is distirbuted evenly by the two nodes connecting the edge?
                     End Do
                  End Do
               End Do
            End If
         End If
      End If
  !
      If (NUMQBC > 0) DEALLOCATE (WDIS)
  !
      If (IOBCN > 0) Then
     !
     !  SET WATER QUALITY CONDITIONS ON OUTER BOUNDARY
     !
         Allocate (WQOBC_TMP(1:IOBCN, 1:KBM1, 1:NCP))
         Call BRACKET (NUT_TM, THOUR, L1, L2, FACT, UFACT, IERR)
         If (IERR ==-1) Then
            WQOBC_TMP = 0.0_SP
         Else
        !
            WQOBC_TMP (1:IOBCN, 1:KBM1, 1:NCP) = UFACT * WQOBC &
           & (1:IOBCN, 1:KBM1, L1, 1:NCP) + FACT * WQOBC (1:IOBCN, &
           & 1:KBM1, L2, 1:NCP)
         End If
     !
#ifdef TA_OBC_INI_FROM_S
     !---------Make TALK (variable=34) a function of Salinity at open boundary - 
     !Write (*,*) '!  TALK at OBC      :    function of salinity'
         Do K = 1, KBM1
            Do I = 1, IOBCN
               J = I_OBC_N (I)
               TAasFsalt = (47.7_SP*SALT(J, K)+646.7_SP)!umol/kg-SW
               TAasFsalt = TAasFsalt * RHO (J, K) / 1000.0_SP !mmol/m3
           !
               WQOBC_TMP (I, K, 34) = TAasFsalt
            End Do
         End Do
     !---------End 
#endif
     !
     !RAMP_TS = Tanh (FLOAT(IINT)/FLOAT(8640+1)) !this is only the same as HYD if DLT=10sec
         RAMP_TS = Tanh (FLOAT(IINT)/FLOAT(stepsIn1day+1))!

         Do JCON = 1, NAC
            II = AC (JCON)
            Allocate (TTMP(NOBTY, KBM1))
            TTMP = 0.0_SP
            Do I = 1, IOBCN
               J = I_OBC_N (I)
               J1 = NEXT_OBC (I)
               T2D = 0.0_SP
               T2D_NEXT = 0.0_SP
               XFLUX2D = 0.0_SP
               Do K = 1, KBM1
                  T2D = T2D + C2 (J, K, II) * DZ2D (J,K)
                  T2D_NEXT = T2D_NEXT + C2F (J1, K, II) * DZ2D (J1,K)

                  XFLUX2D = XFLUX2D + XFLUX_OBC_WQM (I, K, II) * DZ2D (J,K)!Added  XFLUX_OBC_WQM is calculated in adv_wqm.F, analogous to FVCOM

               End Do
           !
               If (UARD_OBCN(I) > 0.0) Then
                  TMP = XFLUX2D + T2D * UARD_OBCN (I)

                  T2D_OBC = (T2D*DT(J)-TMP*DTI/ART1(J)) / D (J)

                  CC = Sqrt (GRAV*H(J)) * DTI / DLTN_OBC (I)
                  CP = CC + 1.0_SP
                  Do K = 1, KBM1
                     PERT_NEXT = C2F (J1, K, II) - T2D_NEXT
                     TYC2F = C2F (J1, K, II)
                     TY2D_N = T2D_NEXT
                     PERT = C2 (J, K, II) - T2D
                 !
                     C2F (J, K, II) = T2D_OBC + &
                    & (CC*PERT_NEXT+PERT*(1.0_SP-DTI/10800.0_SP)) / CP

                  End Do
              ! The end of Perturbation
              !
                  Do K = 1, KBM1
                     TMAX = MAXVAL (C2(NBSN(J, 1:NTSN(J)), K, II))
                     TMIN = MINVAL (C2(NBSN(J, 1:NTSN(J)), K, II))
                 !
                     If (K == 1) Then
                        TMAX = Max (TMAX, (C2(J, K, II)*DZ2D(J,K+1)+C2(J, &
                       & K+1, II)*DZ2D(J,K))/(DZ2D(J,K)+DZ2D(J,K+1)))
                        TMIN = Min (TMIN, (C2(J, K, II)*DZ2D(J,K+1)+C2(J, &
                       & K+1, II)*DZ2D(J,K))/(DZ2D(J,K)+DZ2D(J,K+1)))
                     Else If (K == KBM1) Then
                        TMAX = Max (TMAX, (C2(J, K, II)*DZ2D(J,K-1)+C2(J, &
                       & K-1, II)*DZ2D(J,K))/(DZ2D(J,K)+DZ2D(J,K-1)))
                        TMIN = Min (TMIN, (C2(J, K, II)*DZ2D(J,K-1)+C2(J, &
                       & K-1, II)*DZ2D(J,K))/(DZ2D(J,K)+DZ2D(J,K-1)))
                     Else
                        TMAX = Max (TMAX, (C2(J, K, II)*DZ2D(J,K-1)+C2(J, &
                       & K-1, II)*DZ2D(J,K))/(DZ2D(J,K)+DZ2D(J,K-1)), (C2(J, K, &
                       & II)*DZ2D(J,K+1)+C2(J, K+1, &
                       & II)*DZ2D(J,K))/(DZ2D(J,K)+DZ2D(J,K+1)))
                        TMIN = Min (TMIN, (C2(J, K, II)*DZ2D(J,K-1)+C2(J, &
                       & K-1, II)*DZ2D(J,K))/(DZ2D(J,K)+DZ2D(J,K-1)), (C2(J, K, &
                       & II)*DZ2D(J,K+1)+C2(J, K+1, &
                       & II)*DZ2D(J,K))/(DZ2D(J,K)+DZ2D(J,K+1)))
                     End If
                 !
                     If (TMIN-C2F(J, K, II) > 0.0_SP) C2F (J, K, II) = &
                    & TMIN
                     If (C2F(J, K, II)-TMAX > 0.0_SP) C2F (J, K, II) = &
                    & TMAX

                  End Do
              !
               Else
              !
                  Do K = 1, KBM1

                     If (IERR .Ne.-1) Then

                        C2F (J, K, II) = C2 (J, K, II) - 0.0014 * &    !!
                       & RAMP_TS * (C2(J, K, II)-WQOBC_TMP(I, K, II))
                    !
                     Else
                        C2F (J, K, II) = C2 (J, K, II)
                     End If
                 !
                  End Do
              !
               End If
            End Do
            Deallocate (TTMP)
         End Do !!OUTER LOOP OVER WQ VARIABLES
         Deallocate (WQOBC_TMP)
      End If
  !
  !
  !--SET BOUNDARY CONDITIONS-----------------------------------------------------|
  !
      C2F (0, :, :) = 0.0
  !
      Return
End Subroutine BCOND_WQM
!
!==============================================================================|
!
!==============================================================================|
Subroutine BCOND_NUT_PERTURBATION (T2D_NEXT, T2D, TTMP, I, J, J1)
  !==============================================================================|
  ! Calculate the OBC for temperature perturbation                               |
  !==============================================================================|
      Use MOD_SIZES, Only: NOBTY
      Use MOD_PREC, Only: SP
      Use MOD_LIMS, Only: KBM1, IINT
      Implicit None
  !
  !   INTEGER :: I1,I2,J,JN
      Integer :: I, J, J1, J2, K
      Real (SP) :: CC, CP, MU, CL
      Real (SP) :: PERT_NEXT, PERT, T2D_NEXT, T2D
      Real (SP) :: T2D_NEXT1, TM12D_NEXT2, TM12D_NEXT1, TM22D_NEXT1
      Real (SP) :: TTMP (NOBTY, KBM1)

      Return
End Subroutine BCOND_NUT_PERTURBATION
!========================================================================
