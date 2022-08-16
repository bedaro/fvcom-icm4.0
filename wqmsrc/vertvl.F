!vertvl.F
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
!Subroutine VERTVL()
!
!==============================================================================|
!   CALCULATE THE SIGMA COORDINATE VERTICAL VELOCITY FOR THE 3D MODE (omega)   |
!                                                  |
!   DETERMINED FROM EQUATION:                              |
!                                              !
!   d/dt(D) + d/dx(uD) + d/dy(uD) = d/sigma(omega)                             !
!==============================================================================|
!
Subroutine VERTVL
  !
  !------------------------------------------------------------------------------|
      Use MOD_PREC, Only: SP
      Use MOD_LIMS, Only: NCV, MLOC, MTLOC, KBM1, KB
  !
      Use MOD_TGE, Only: ISONB, NIEC, NTRG, DLTXE, DLTYE,NV
  !
      Use MOD_HYDROVARS, Only: ART1, DZ,DZ2D, H, DT1, EL, ET, DTFA, UU, VV, &
     & WTS
      Use MOD_WQM, Only: DLT
  !
      Implicit None
      Real (SP) :: XFLUX (MTLOC, KBM1)
      Real (SP) :: DIJ, UIJ, VIJ, UN, EXFLUX, TMP1
      Integer :: I, K, IA, IB, I1, J, JJ, J1, J2
  !------------------------------------------------------------------------------|
  !
  !----------------------INITIALIZE FLUX-----------------------------------------!
  !
      XFLUX = 0.0_SP
  !
  !----------------------ACCUMULATE FLUX-----------------------------------------!
  !
      Do I = 1, NCV
         I1 = NTRG (I)
         IA = NIEC (I, 1)
         IB = NIEC (I, 2)

         Do K = 1, KBM1
            DIJ = DT1 (I1) * ( DZ2D(NV(I1,1),K)+DZ2D (NV(I1,2),K)+DZ2D (NV(I1,3),K))/3.0_SP
						!to be as in vertvl_edge.F (FVCOM)
            UIJ = UU (I1, K)
            VIJ = VV (I1, K)
            EXFLUX = DIJ * (-UIJ*DLTYE(I)+VIJ*DLTXE(I))
            XFLUX (IA, K) = XFLUX (IA, K) - EXFLUX
            XFLUX (IB, K) = XFLUX (IB, K) + EXFLUX
         End Do
      End Do
  !
  !-----------------------NULLIFY BOUNDARY FLUX----------------------------------!
  !
      Do I = 1, MLOC
         If (ISONB(I) == 2) XFLUX (I, 1:KBM1) = 0.0_SP
      End Do
  !
  !---IF NO FRESH WATER INFLOW, OMEGA IS ZERO AT FREE SURFACE AND BOTTOM---------!
  !

      WTS (1:MLOC, 1) = 0.0_SP
  !
  !--------------------------CALCULATE OMEGA-------------------------------------!
  !
      Do I = 1, MLOC
         Do K = 1, KBM1
            WTS (I, K+1) = WTS (I, K) + DZ2D (I,K) * (XFLUX(I, &
           & K)/ART1(I)+(EL(I)-ET(I))/DLT)
         End Do
      End Do
  !
  !--------------------------ADJUST OMEGA----------------------------------------!
  ! IMPROVES MASS CONSERVATION
  !
      Do I = 1, MLOC
         If (Abs(WTS(I, KB)) > 1.0E-8_SP) Then
            If (ISONB(I) /= 2) Then
           !
               TMP1 = EL (I) * FLOAT (KBM1) - WTS (I, KB) * DLT / DZ2D(I, &
              & 1)!change of layer thickness due to vertical velocity
           !
               TMP1 = TMP1 / FLOAT (KBM1)
               DTFA (I) = TMP1 + H (I)
               Do K = 2, KB
                  WTS (I, K) = WTS (I, K) - FLOAT (K-1) / FLOAT (KBM1) &
                 & * WTS (I, KB)
               End Do
            End If
         End If
      End Do
  !
      Return
End Subroutine VERTVL
!==============================================================================|
